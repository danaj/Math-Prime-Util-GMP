#!/usr/bin/env python3
"""Real factoring/replay, frame integrity and interrupted-write regression tests."""
import argparse
import fcntl
import math
from pathlib import Path
import struct
import subprocess
import tempfile
import time
import zlib

import sys
sys.dont_write_bytecode = True
from _support import ROOT, RESULTS, EXE
PREAMBLE = b'MPUSIQJ\0' + struct.pack('<II', 1, 0x01020304)
FRAME = struct.Struct('<IIQIIQQQII')
SMALL = 1000000007 * 1000000009
ONE = ((1 << 89) - 1) * ((1 << 107) - 1)
TWO = ((1 << 127) - 1) * ((1 << 128) - 159)


def frames(data, partial=False):
    end = data.find(b'\n--\n')
    if end < 0 or len(data) < end+4+len(PREAMBLE):
        if partial:
            return []
        raise AssertionError('short ASCII or binary prefix')
    lines = data[:end+4].decode('ascii').splitlines()
    assert len(lines) >= 5 and lines[:3] == ['MPUSIQJ', 'Version: 1', 'Byte-order: little-endian']
    assert lines[3].startswith('N: ') and lines[3][3:].isdigit() and lines[-1] == '--'
    number = int(lines[3][3:])
    binary = end+4
    assert data[binary:binary+len(PREAMBLE)] == PREAMBLE
    result, at = [], binary+len(PREAMBLE)
    while at < len(data):
        if len(data) - at < FRAME.size:
            if partial:
                break
            raise AssertionError('short frame header')
        magic, kind, seq, length, count, run_id, utc_us, elapsed_us, header_crc, crc = FRAME.unpack_from(data, at)
        end = at + FRAME.size + length
        if end > len(data) and partial:
            break
        assert end <= len(data) and magic == 0x4b4c424a
        assert header_crc == zlib.crc32(data[at:at+48])
        assert seq == len(result)
        assert crc == zlib.crc32(data[at:at+52] + data[at+FRAME.size:end])
        if not result:
            assert kind == 1 and count == 1
            size = struct.unpack_from('<I', data, at+FRAME.size)[0]
            assert length == 4+size
            assert int.from_bytes(data[at+FRAME.size+4:end], 'big') == number
        result.append((at, end, kind, count))
        at = end
    return result


def compact_relations(data, parsed):
    """Independent structural decoder, also locating fields for bad-input checks."""
    result=[]
    for start,end,kind,event_count in parsed:
        at=start+FRAME.size
        if kind==2:
            nbytes=struct.unpack_from('<I',data,at)[0]
            at+=4+nbytes+12
            namebytes=struct.unpack_from('<I',data,at)[0]
            at+=4+namebytes
            params=struct.unpack_from('<23Q',data,at)
            bound=params[-2]
            width=max(1,(bound.bit_length()+7)//8)
            fb_size=params[4]
        if kind!=3:
            continue
        at+=192
        def varint():
            nonlocal at
            value=0;shift=0
            while True:
                assert at<end and shift<=28
                byte=data[at];at+=1
                value|=(byte&127)<<shift
                if byte<128:
                    assert shift==0 or byte!=0
                    assert value<=0xffffffff
                    return value
                shift+=7
        for _ in range(event_count):
            kind,labels=struct.unpack_from('BB',data,at);at+=2
            if kind==2:
                q=labels;assert q>0;at+=4*q
                continue
            assert kind==1 and labels in (0,1,2)
            lp_at=at
            lp1=int.from_bytes(data[at:at+width],'little') if labels==2 else 1
            if labels==2:at+=width
            lp2=int.from_bytes(data[at:at+width],'little') if labels else 1
            if labels:at+=width
            assert 1<=lp1<=lp2<=bound
            assert labels==0 or lp2>1
            assert labels!=2 or lp1>1
            ybytes=struct.unpack_from('<I',data,at)[0];at+=4+ybytes
            count=struct.unpack_from('<I',data,at)[0];at+=4
            factor_at=at;previous=0
            for i in range(count):
                code=varint();delta=code>>1
                exponent=varint() if code&1 else 1
                assert 1<=exponent<=2047 and (not (code&1) or exponent!=1)
                assert i==0 or delta>0
                previous+=delta
                assert previous<=fb_size
            result.append(dict(frame_start=start,frame_end=end,labels=labels,
                               lp_at=lp_at,width=width,factor_at=factor_at,
                               fb_size=fb_size,count=count))
        assert at==end
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--executable', type=Path, default=EXE)
    ap.add_argument('--serial-executable', type=Path)
    args = ap.parse_args()
    executable = args.executable.resolve()
    passed = 0

    def check(label):
        nonlocal passed
        passed += 1
        print('PASS', label, flush=True)

    with tempfile.TemporaryDirectory(prefix='journal-tests-', dir=RESULTS) as tmp:
        tmp = Path(tmp)
        def run(n, path=None, threads=1, ok=True, exe=executable, stdin=False, inferred=False):
            command = [str(exe), '-v']
            if exe == executable:
                command += ['-t', str(threads)]
            if path is not None:
                command += ['--journal', str(path)]
            if stdin:
                command += ['-']
            elif not inferred:
                command += [str(n)]
            p = subprocess.run(command, input=f'{n}\n' if stdin else None,
                               text=True, capture_output=True, timeout=120)
            if ok:
                assert p.returncode == 0, p.stdout + p.stderr
                lines = [line for line in p.stdout.splitlines() if line.startswith(f'{n}:')]
                assert len(lines) == 1, p.stdout + p.stderr
                factors = [int(x) for x in lines[0].split(':', 1)[1].split()]
                assert math.prod(factors) == n
                if n != 1:
                    assert all(1 < x < n for x in factors) or factors == [n]
            else:
                assert p.returncode != 0, p.stdout + p.stderr
            return p

        before = set(tmp.iterdir())
        run(SMALL)
        assert set(tmp.iterdir()) == before
        check('default path creates no journal')

        journals = {}
        for name, n, threads in [('q1', 10007*10009, 1), ('q2', 1000003*1000033, 1),
                                 ('smooth', SMALL, 1), ('one-lp', ONE, 1), ('two-lp', TWO, 4)]:
            path = tmp/(name+'.journal')
            run(n, path, threads)
            data = path.read_bytes()
            parsed = frames(data)
            assert {1, 2, 3, 4, 5} <= {f[2] for f in parsed}
            run(n, path, 1 if threads > 1 else 4)
            assert path.read_bytes() == data
            journals[name] = (n, data, parsed)
            check(name+' complete replay, cross-thread count, byte-identical journal')

        source = tmp/'two-lp.journal'
        old = source.read_bytes()
        p = run(TWO+2, source, ok=False)
        assert 'input does not match' in p.stdout+p.stderr
        assert source.read_bytes() == old
        check('wrong input rejected before modifying existing file')

        n, data, parsed = journals['two-lp']
        ascii_end = data.index(b'\n--\n')+4
        assert data[:ascii_end] == f'MPUSIQJ\nVersion: 1\nByte-order: little-endian\nN: {n}\n--\n'.encode('ascii')
        check('five-line readable ASCII prefix and matching binary INPUT')

        path = tmp/'header-notes.journal'
        notes = b'# free-form comment\n\nCollected on: a test host\n---\n-- not the delimiter\n' + b'x'*8192+b'\n'
        decorated = data[:ascii_end-3]+notes+data[ascii_end-3:]
        frames(decorated)
        path.write_bytes(decorated)
        run(n, path, inferred=True)
        assert path.read_bytes()==decorated
        check('optional header comments, metadata and long lines before exact -- delimiter')

        path = tmp/'missing-header-delimiter.journal'
        missing = data[:ascii_end-3]+notes
        path.write_bytes(missing)
        run(n, path, ok=False, inferred=True)
        assert path.read_bytes()==missing
        check('optional header notes still require a complete -- delimiter line')

        labels=set()
        for _,stored,stored_frames in journals.values():
            labels.update(r['labels'] for r in compact_relations(stored,stored_frames))
        assert labels=={0,1,2}
        check('compact smooth/1LP/2LP labels and canonical delta/exponent vectors')

        relation=next(r for r in compact_relations(data,parsed) if r['count']>=2)
        def uleb(value):
            out=bytearray()
            while value>=128:
                out.append((value&127)|128);value>>=7
            out.append(value)
            return bytes(out)
        for label,replacement in [
                ('overflow',b'\xff\xff\xff\xff\x10'),
                ('overlong',b'\x80\x00'),
                ('zero-exponent',b'\x01\x00'),
                ('explicit-one',b'\x01\x01'),
                ('wide-exponent',b'\x01'+uleb(2048)),
                ('row-overflow',uleb(2*(relation['fb_size']+1))),
                ('repeated-row',b'\x02\x00')]:
            start,end=relation['frame_start'],relation['frame_end']
            payload=bytearray(data[start+FRAME.size:end])
            offset=relation['factor_at']-start-FRAME.size
            payload[offset:offset]=replacement
            fields=list(FRAME.unpack_from(data,start))
            fields[3]=len(payload);fields[8]=fields[9]=0
            header=bytearray(FRAME.pack(*fields))
            struct.pack_into('<I',header,48,zlib.crc32(header[:48]))
            struct.pack_into('<I',header,52,zlib.crc32(header[:52]+payload))
            changed=data[:start]+header+payload+data[end:]
            path=tmp/f'bad-factors-{label}.journal'
            path.write_bytes(changed)
            run(n,path,ok=False)
            assert path.read_bytes()==changed
        check('checksummed malformed varints, exponents and row deltas rejected')

        path = tmp/'ascii-number.journal'
        changed = data.replace(f'N: {n}\n'.encode(), f'N: {n+2}\n'.encode(), 1)
        path.write_bytes(changed)
        p = run(n+2, path, ok=False)
        assert 'ASCII and binary input headers do not match' in p.stdout+p.stderr
        assert path.read_bytes() == changed
        check('ASCII input cannot disagree with the binary INPUT block')

        # Also verify the reverse: a valid, checksummed binary INPUT cannot
        # override the readable input line.
        path = tmp/'binary-number.journal'
        changed = bytearray(data)
        start,end,_,_ = parsed[0]
        changed[end-1] ^= 1
        struct.pack_into('<I', changed, start+52, zlib.crc32(changed[start:start+52]+changed[start+FRAME.size:end]))
        path.write_bytes(changed)
        p = run(n, path, ok=False)
        assert 'ASCII and binary input headers do not match' in p.stdout+p.stderr
        assert path.read_bytes() == changed
        check('binary INPUT cannot disagree with the ASCII input')

        for label, ascii_old, ascii_new, binary_offset in (
                ('version', b'Version: 1\n', b'Version: 2\n', 8),
                ('byte-order', b'Byte-order: little-endian\n', b'Byte-order: big-endian\n', 12)):
            for side in ('ascii', 'binary'):
                path = tmp/f'bad-{side}-{label}.journal'
                if side == 'ascii':
                    changed = data.replace(ascii_old, ascii_new, 1)
                else:
                    changed = bytearray(data)
                    changed[ascii_end+binary_offset] ^= 1
                path.write_bytes(changed)
                run(n, path, ok=False)
                assert path.read_bytes() == changed
            check('ASCII and binary '+label+' must agree with the supported format')

        for suffix, changed in [('truncated',data[:ascii_end-2]),
                                ('invalid-number',data.replace(f'N: {n}\n'.encode(),b'N: 12x34\n',1))]:
            path = tmp/f'bad-ascii-{suffix}.journal'
            path.write_bytes(changed)
            run(n, path, ok=False)
            assert path.read_bytes() == changed
        check('truncated and malformed ASCII headers rejected without modification')

        batch_indices = [i for i,f in enumerate(parsed) if f[2] == 3]
        first = parsed[batch_indices[2]]
        for label, cut in [('valid-prefix',first[1]),('partial-header',first[1]+11),
                           ('partial-payload',parsed[batch_indices[3]][0]+FRAME.size+19)]:
            path = tmp/(label+'.journal')
            prefix = data[:cut]
            path.write_bytes(prefix)
            run(n, path, 4)
            new = path.read_bytes()
            frames(new)
            assert new.startswith(data[:first[1]])
            assert new != prefix
            run(n, path, 1)
            check(label+' rebuilt graph and resumed with append')

        # A genuinely interrupted serial collector, then a parallel restart.
        path = tmp/'live-interruption.journal'
        with (tmp/'interrupted.log').open('w') as log:
            p = subprocess.Popen([str(executable), '-v', '-t', '1', '--journal', str(path), str(TWO)],
                                 stdout=log, stderr=log)
            deadline = time.monotonic()+30
            while time.monotonic() < deadline:
                if path.exists() and len(frames(path.read_bytes(), partial=True)) >= 5:
                    break
                assert p.poll() is None, 'collector finished before interruption'
                time.sleep(.005)
            else:
                p.kill(); p.wait(); raise AssertionError('no committed relation batch')
            p.terminate(); p.wait(timeout=10)
        prefix = path.read_bytes()
        committed = frames(prefix, partial=True)[-1][1]
        run(TWO, path, 4, inferred=True)
        assert path.read_bytes().startswith(prefix[:committed])
        frames(path.read_bytes())
        check('SIGTERM interruption; serial-to-parallel continuation')

        # Change a complete relation payload without updating its checksum.
        path = tmp/'bad-crc.journal'
        damaged = bytearray(data)
        damaged[first[0]+FRAME.size+8] ^= 1
        path.write_bytes(damaged)
        p = run(n, path, ok=False)
        assert 'checksum mismatch' in p.stdout+p.stderr
        assert path.read_bytes() == damaged
        check('complete bad checksum rejected, file unchanged')

        # A corrupted payload length must not masquerade as a short tail.
        path = tmp/'bad-length.journal'
        damaged = bytearray(data)
        struct.pack_into('<I', damaged, first[0]+16, 32*1024*1024)
        path.write_bytes(damaged)
        p = run(n, path, ok=False)
        assert 'header checksum mismatch' in p.stdout+p.stderr
        assert path.read_bytes() == damaged
        check('corrupt frame length rejected instead of discarding a valid suffix')

        # A well-checksummed context for the wrong resolved parameters must fail.
        path = tmp/'bad-context.journal'
        damaged = bytearray(data)
        start,end,_,_ = next(f for f in parsed if f[2] == 2)
        length = struct.unpack_from('<I',damaged,start+FRAME.size)[0]
        multiplier_at = start+FRAME.size+4+length
        damaged[multiplier_at] ^= 1
        struct.pack_into('<I',damaged,start+52,zlib.crc32(damaged[start:start+52]+damaged[start+FRAME.size:end]))
        path.write_bytes(damaged)
        p = run(n,path,ok=False)
        assert 'resolved parameters, factor base, or attempt differs' in p.stdout+p.stderr
        assert path.read_bytes() == damaged
        check('compatible input with different context rejected')

        path = tmp/'cofactor-calls.journal'
        run(SMALL**2, path)
        assert sum(f[2] == 5 for f in frames(path.read_bytes())) == 3
        old = path.read_bytes(); run(SMALL**2, path,4)
        assert path.read_bytes() == old
        check('perfect-power partition and repeated smaller cofactor calls')

        for n in (1, 2, 101, 3*SMALL):
            path = tmp/f'pretest-{n}.journal'
            run(n,path); old=path.read_bytes(); run(n,path,4)
            assert path.read_bytes() == old
        check('one, primes, and trial-division partitions')

        path=tmp/'stdin.journal'; run(SMALL,path,stdin=True)
        p=subprocess.run([str(executable),'--journal',str(tmp/'multiple.journal'),str(SMALL),str(ONE)],
                         capture_output=True,text=True)
        assert p.returncode and not (tmp/'multiple.journal').exists()
        p=subprocess.run([str(executable),'--journal',str(tmp/'multiple.journal')],
                         input=f'{SMALL}\n{ONE}\n',capture_output=True,text=True)
        assert p.returncode and not (tmp/'multiple.journal').exists()
        p=subprocess.run([str(executable),'--journal',str(tmp/'multiple.journal'),'-'],
                         input=f'{SMALL}\n{ONE}\n',capture_output=True,text=True)
        assert p.returncode and not (tmp/'multiple.journal').exists()
        check('single stdin input; reject multiple inputs before file creation')

        source=tmp/'smooth.journal'; old=source.read_bytes()
        run(SMALL,source,4,inferred=True)
        assert source.read_bytes()==old
        for number in (1,101):
            run(number,tmp/f'pretest-{number}.journal',inferred=True)
        check('omitted N resumes completed journals, including one and primes')

        # Keep the pipe open: inferred journal input must never wait for stdin.
        p=subprocess.Popen([str(executable),'--journal',str(source)],
                           stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        try:
            assert p.wait(timeout=5)==0
            assert p.stdout.read().startswith(f'{SMALL}:')
        finally:
            if p.poll() is None:p.kill();p.wait()
            p.stdin.close();p.stdout.close();p.stderr.close()
        check('omitted N does not read stdin or wait for EOF')

        missing=tmp/'missing-input.journal'
        p=subprocess.run([str(executable),'--journal',str(missing)],
                         input=f'{SMALL}\n',capture_output=True,text=True,timeout=5)
        assert p.returncode and 'journal does not exist' in p.stdout+p.stderr
        assert not missing.exists()
        check('omitted N with missing journal fails without creating a file')

        path=tmp/'inferred-mismatch.journal'
        changed=old.replace(f'N: {SMALL}\n'.encode(),f'N: {SMALL+2}\n'.encode(),1)
        path.write_bytes(changed)
        p=run(SMALL+2,path,ok=False,inferred=True)
        assert 'ASCII and binary input headers do not match' in p.stdout+p.stderr
        assert path.read_bytes()==changed
        check('inferred N still checks the binary INPUT block')

        path=tmp/'inferred-zero.journal'
        changed=old.replace(f'N: {SMALL}\n'.encode(),b'N: 0\n',1)
        path.write_bytes(changed)
        run(0,path,ok=False,inferred=True)
        assert path.read_bytes()==changed
        check('inferred input must be positive')

        path=tmp/'inferred-invalid-tail.journal';path.write_bytes(old+b'X'*FRAME.size)
        p=run(SMALL,path,ok=False,inferred=True)
        assert not any(line.startswith(f'{SMALL}:') for line in p.stdout.splitlines())
        assert path.read_bytes()==old+b'X'*FRAME.size
        check('inferred resume validates the journal tail before printing factors')

        path=tmp/'stdin-wrong-input.journal';path.write_bytes(old)
        run(SMALL+2,path,ok=False,stdin=True)
        assert path.read_bytes()==old
        check('explicit stdin input is checked against an existing journal')

        numbers=[SMALL,1000003*1000033]
        for stdin in (False,True):
            command=[str(executable)]+([] if stdin else list(map(str,numbers)))
            p=subprocess.run(command,input='\n'.join(map(str,numbers))+'\n' if stdin else None,
                             capture_output=True,text=True,timeout=10)
            assert p.returncode==0 and len(p.stdout.splitlines())==2
            for number,line in zip(numbers,p.stdout.splitlines()):
                assert line.startswith(f'{number}:')
                assert math.prod(map(int,line.split(':',1)[1].split()))==number
        check('ordinary multiple-input argv and implicit stdin behavior unchanged')

        path=tmp/'stdin.journal'
        with path.open('r+b') as locked:
            fcntl.lockf(locked, fcntl.LOCK_EX|fcntl.LOCK_NB)
            p=run(SMALL,path,ok=False)
            assert 'exclusive journal lock' in p.stdout+p.stderr
        check('second process cannot use an already locked journal')

        if args.serial_executable:
            serial=args.serial_executable.resolve()
            path=tmp/'serial-build.journal'; run(SMALL,path,exe=serial)
            run(SMALL,path,exe=serial,inferred=True)
            old=path.read_bytes(); run(SMALL,path,4)
            assert path.read_bytes() == old
            path=tmp/'serial-build-resume.journal'; path.write_bytes(data[:first[1]])
            run(TWO,path,exe=serial)
            check('serial-only build creation and parallel-to-serial continuation')
    print(f'All {passed}/{passed} journal tests passed.', flush=True)


if __name__ == '__main__':
    main()
