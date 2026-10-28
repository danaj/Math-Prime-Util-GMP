#!/usr/bin/env python3
"""First-frame timing checks using synthetic progress, without factoring."""
import importlib.util
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import zlib

sys.dont_write_bytecode=True
from _support import RESULTS, TOOL


def module(name, path):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result)
    return result


info=module('journal_info',TOOL)
reference=module('journal_reference',Path(__file__).with_name('test-journal.py'))
fixture=Path(sys.argv[1]) if len(sys.argv)>1 else RESULTS/'journal-attempts'


def inspect(path):
    reader=info.Journal(path,True)
    try:reader.scan();return reader.summary(300)
    finally:reader.close()


def magnitude(n):
    data=n.to_bytes((n.bit_length()+7)//8,'big')
    return struct.pack('<I',len(data))+data


def close(a, b):
    assert abs(a-b)<1e-9,(a,b)


with tempfile.TemporaryDirectory(prefix='journal-timing-',dir=RESULTS) as directory:
    tmp=Path(directory);clock=tmp/'clock.journal'
    subprocess.run([str(fixture),'--clock',str(clock)],check=True,timeout=10)
    data=clock.read_bytes()
    headers=[info.FRAME.unpack_from(data,start) for start,_,_,_ in reference.frames(data)]
    assert [h[5] for h in headers]==[0,1,1]
    assert headers[0][7]==0 and 0<headers[1][7]<=headers[2][7]
    assert [run['id'] for run in inspect(clock)['runs']]==[0,1]
    print('PASS first resumed frame includes work before serialization',flush=True)

    source=tmp/'source.journal'
    subprocess.run([str(fixture),'--events',str(source)],check=True,
                   stdout=subprocess.DEVNULL,timeout=10)
    data=source.read_bytes();parsed=reference.frames(data)
    prefix=data[:parsed[0][0]]
    def payload(item):return data[item[0]+info.FRAME.size:item[1]]
    input_payload=payload(parsed[0]);context=next(f for f in parsed if f[2]==2)
    context_payload=payload(context)
    batch=next(f for f in parsed if f[2]==3)
    records=reference.compact_relations(data,parsed)
    starts=[r['lp_at']-2 for r in records]
    family=data[batch[0]+info.FRAME.size+192:starts[0]]
    relations=[data[start:end] for start,end in zip(starts,starts[1:]+[batch[1]])]
    assert len(relations)==3

    def progress(index):
        snapshot=[0]*24;snapshot[0]=1;snapshot[2]=index+1
        snapshot[10:13]=[1,2,3]
        snapshot[13]=(1,1,2)[index];snapshot[14]=(0,1,1)[index]
        snapshot[15]=100;snapshot[16]=snapshot[13]
        return struct.pack('<24Q',*snapshot)

    def record(index, run, elapsed):
        events=(family if index==0 else b'')+relations[index]
        return 3,2 if index==0 else 1,progress(index)+events,run,elapsed

    def stream(entries):
        result=bytearray(prefix)
        for sequence,(kind,count,body,run,elapsed) in enumerate(entries):
            fields=[0x4b4c424a,kind,sequence,len(body),count,run,
                    1700000000000000+sequence*1000000,elapsed,0,0]
            header=bytearray(info.FRAME.pack(*fields))
            struct.pack_into('<I',header,48,zlib.crc32(header[:48]))
            struct.pack_into('<I',header,52,zlib.crc32(body,zlib.crc32(header[:52])))
            result.extend(header);result.extend(body)
        return bytes(result)

    initial=[(1,1,input_payload,0,0),(2,1,context_payload,0,1000000)]
    entries=initial+[record(0,0,5000000)]
    path=tmp/'progress.journal';path.write_bytes(stream(entries))
    rate=inspect(path)['contexts'][0]['rate']
    close(rate['window_seconds'],4);close(rate['raw_per_second'],.25)
    close(rate['full_per_second'],.25);close(rate['collection_eta_seconds'],396)
    print('PASS first collection batch uses the context timestamp and zero counts',flush=True)

    entries.append(record(1,0,9000000))
    previous=stream(entries);path.write_bytes(previous)
    reader=info.Journal(path,True)
    try:
        reader.scan()
        entries.append(record(2,1,10000000));resumed=stream(entries)
        with path.open('ab') as writer:writer.write(resumed[len(previous):])
        reader.scan();rate=reader.summary(300)['contexts'][0]['rate']
    finally:reader.close()
    assert rate['run']==1
    close(rate['window_seconds'],10);close(rate['raw_per_second'],.1)
    close(rate['full_per_second'],.1)
    print('PASS first resumed batch subtracts replayed counts, including incremental scans',flush=True)

    legacy=entries[:-1]+[record(2,1,0)]
    path.write_bytes(stream(legacy))
    assert inspect(path)['contexts'][0]['rate'] is None
    legacy.append((3,0,progress(2),1,10000000));path.write_bytes(stream(legacy))
    rate=inspect(path)['contexts'][0]['rate']
    close(rate['raw_per_second'],0);close(rate['full_per_second'],0)
    print('PASS zero-timestamp first batches leave their unmeasurable interval excluded',flush=True)

    n=int(inspect(source)['n'])
    ended=struct.pack('<I',0)+magnitude(n)+struct.pack('<I',1)+magnitude(n)
    later=entries+[(4,1,ended,1,11000000),(2,1,context_payload,1,20000000),
                   record(0,1,30000000)]
    path.write_bytes(stream(later));contexts=inspect(path)['contexts']
    assert len(contexts)==2
    rate=contexts[1]['rate']
    close(rate['window_seconds'],10);close(rate['raw_per_second'],.1)
    close(rate['full_per_second'],.1)
    print('PASS later contexts exclude earlier work in the same run',flush=True)

    backwards=entries+[(3,0,progress(2),1,9000000)]
    path.write_bytes(stream(backwards));result=inspect(path)
    assert result['contexts'][0]['rate'] is None
    assert any('elapsed timestamp moved backwards' in message for message in result['warnings'])
    print('PASS clock anomalies discard the affected rate interval',flush=True)

print('7/7 journal timing groups passed.',flush=True)
