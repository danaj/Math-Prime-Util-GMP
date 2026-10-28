#!/usr/bin/env python3
"""Focused two-byte event checks; synthetic streams avoid a full factor job."""
import contextlib
import importlib.util
import io
import json
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

def frame(kind, sequence, count, payload):
    fields=[0x4b4c424a,kind,sequence,len(payload),count,0,
            1700000000000000+sequence*1000000,sequence*1000000,0,0]
    header=bytearray(info.FRAME.pack(*fields))
    struct.pack_into('<I',header,48,zlib.crc32(header[:48]))
    struct.pack_into('<I',header,52,zlib.crc32(payload,zlib.crc32(header[:52])))
    return header+payload

def magnitude(n):
    data=n.to_bytes((n.bit_length()+7)//8,'big')
    return struct.pack('<I',len(data))+data

def snapshot(families=0, partials=0, target=1):
    values=[0]*24;values[0]=families;values[10:13]=[1,2,3]
    values[2]=values[3]=values[5]=partials
    values[14]=partials;values[15]=target
    return struct.pack('<24Q',*values)

with tempfile.TemporaryDirectory(prefix='event-headers-',dir=RESULTS) as directory:
    tmp=Path(directory);path=tmp/'events.journal'
    run=subprocess.run([str(fixture),'--events',str(path)],capture_output=True,text=True,timeout=10)
    assert run.returncode==0,run.stderr
    data=path.read_bytes();parsed=reference.frames(data)
    records=reference.compact_relations(data,parsed);summary=inspect(path)
    assert {v['labels'] for v in records}=={0,1,2}
    assert summary['relation_congruences_verified']==len(records)==3
    assert summary['byte_breakdown']['relation_tags_and_counts']==6*len(records)
    ctx=summary['contexts'][0]
    assert summary['byte_breakdown']['family_records']==ctx['families']*(2+4*ctx['parameters']['q_count'])
    replay=subprocess.run([str(fixture),'--events',str(path)],capture_output=True,text=True,timeout=10)
    assert replay.returncode==0 and replay.stdout==run.stdout and path.read_bytes()==data,replay.stderr
    assert json.loads(run.stdout)['families']==1
    print('PASS C smooth/1LP/2LP encoding, byte accounting, congruences and replay',flush=True)

    batch=next(f for f in parsed if f[2]==3);start,end,_,count=batch
    event=start+info.FRAME.size+192
    assert data[event]==2
    relation=next(v for v in records if start<=v['lp_at']<end)
    corruptions=[('unknown-kind',relation['lp_at']-2,255),
                 ('invalid-LP-count',relation['lp_at']-1,5),
                 ('unsupported-3LP',relation['lp_at']-1,3),
                 ('unsupported-4LP',relation['lp_at']-1,4),
                 ('zero-q',event+1,0),('wide-q',event+1,255)]
    for name,at,value in corruptions:
        payload=bytearray(data[start+info.FRAME.size:end]);payload[at-start-info.FRAME.size]=value
        changed=data[:start]+frame(3,info.FRAME.unpack_from(data,start)[2],count,payload)+data[end:]
        bad=tmp/(name+'.journal');bad.write_bytes(changed)
        run=subprocess.run([str(fixture),'--events',str(bad)],capture_output=True,text=True,timeout=10)
        assert run.returncode!=0 and bad.read_bytes()==changed,(name,run.stderr)
        try:inspect(bad)
        except info.InvalidJournal:pass
        else:raise AssertionError(name+' accepted')
    print('PASS malformed event kinds/counts/q and unsupported replay modes rejected',flush=True)

    # Synthetic, mathematically valid congruences modulo 77, not factor jobs:
    # 1^2 = 3*3*137 (mod 77); 2^2 = 3*3*3*3 (mod 77).
    prefix=b'MPUSIQJ\nVersion: 1\nByte-order: little-endian\nN: 77\n--\n'+info.PREAMBLE
    unsigned=dict.fromkeys(info.UNSIGNED,0)
    unsigned.update(bits=7,policy_first_bits=7,policy_last_bits=7,fb_size=1,
                    half_interval=1,q_count=1,poly_d=1,target_relations=1,
                    large_prime_bound=1000,smooth_bound=1000000000)
    for labels,primes,y in [(3,[3,3,137],1),(4,[3,3,3,3],2)]:
        unsigned['max_large_primes']=labels;policy=b'codec'
        context=magnitude(77)+struct.pack('<QII',1,21,len(policy))+policy
        context+=struct.pack('<23Q',*(unsigned[k] for k in info.UNSIGNED))
        context+=struct.pack('<10d',*([1.0]*10))+struct.pack('<IIII',1,2,1,1)+struct.pack('<3Q',1,2,3)
        record=bytes((1,labels))+b''.join(v.to_bytes(2,'little') for v in primes)+magnitude(y)+struct.pack('<I',0)
        path=tmp/('synthetic-'+str(labels)+'.journal')
        path.write_bytes(prefix+frame(1,0,1,magnitude(77))+frame(2,1,1,context)+frame(3,2,1,snapshot(partials=1)+record))
        report=inspect(path);assert report['relation_congruences_verified']==1
        assert report['relations'][info.RELATION_NAMES[labels]]==1
        output=io.StringIO()
        with contextlib.redirect_stdout(output):info.display(report,False)
        assert '1 '+str(labels)+'LP' in output.getvalue()
    print('PASS inspector decodes and verifies synthetic 3LP/4LP records',flush=True)

    ctx_frame=next(f for f in parsed if f[2]==2);a,b,_,_=ctx_frame
    context=bytearray(data[a+info.FRAME.size:b]);at=4+struct.unpack_from('<I',context)[0]+12
    at+=4+struct.unpack_from('<I',context,at)[0]
    struct.pack_into('<Q',context,at+6*8,16)
    target=summary['contexts'][0]['parameters']['target_relations']
    family=bytes((2,16))+struct.pack('<16I',*range(1,17))
    path=tmp/'q16.journal'
    path.write_bytes(data[:a]+frame(2,1,1,context)+frame(3,2,1,snapshot(families=1,target=target)+family))
    report=inspect(path)
    assert report['contexts'][0]['families']==1 and report['byte_breakdown']['family_records']==66
    print('PASS q=16 uses the info byte without an extra q field',flush=True)
print('All 4/4 event-header checks passed.',flush=True)
