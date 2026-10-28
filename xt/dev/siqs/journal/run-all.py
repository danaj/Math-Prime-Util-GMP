#!/usr/bin/env python3
"""Build both modes and run journal, recovery, and existing SIQS checks."""
from pathlib import Path
import json
import struct
import subprocess
import sys
import os
sys.dont_write_bytecode = True

from _support import ROOT, BUILD, RESULTS
FIXTURES = Path(__file__).resolve().parent
EXE = RESULTS/'msiqs'
os.environ['SIQS_JOURNAL_TEST_DIR'] = str(RESULTS)
os.environ['SIQS_JOURNAL_TEST_EXECUTABLE'] = str(EXE)
SOURCES = ['mpu-siqs.c','siqs.c','lanczos.c','pr32_iter.c','../squfof126.c']


def step(name, command, cwd=ROOT):
    print('START', name, flush=True)
    with (RESULTS/(name+'.log')).open('w') as log:
        p = subprocess.run(command,cwd=cwd,stdout=log,stderr=subprocess.STDOUT)
    if p.returncode:
        print('FAILED',name,'see',RESULTS/(name+'.log'),flush=True)
        raise SystemExit(p.returncode)
    print('PASS',name,flush=True)


RESULTS.mkdir(exist_ok=True)
print('Journal test artifacts:', RESULTS, flush=True)
step('eta-model',[sys.executable,'-u',str(FIXTURES/'test-eta.py')])
step('build',['perl','tools/build-siqs.pl','--output',str(EXE)]+SOURCES,BUILD)
step('build-serial',['perl','tools/build-siqs.pl','--serial','--output',str(RESULTS/'msiqs-serial'),
                     '--cflags','-O2 -falign-loops=32']+SOURCES,BUILD)
step('journal-tests',[sys.executable,'-u',str(FIXTURES/'test-journal.py'),
                      '--executable',str(EXE),'--serial-executable',str(RESULTS/'msiqs-serial')])
step('build-seeds',['cc','-O1','-DSTANDALONE','-o',str(RESULTS/'journal-seeds'),
                    str(FIXTURES/'journal-seeds.c'),str(BUILD/'lanczos.c'),
                    str(BUILD/'pr32_iter.c'),str(ROOT/'squfof126.c'),'-lgmp','-lm'])
path=RESULTS/'seeds.journal'
if path.exists():path.unlink()
step('seed-record',[str(RESULTS/'journal-seeds'),str(path)])
data=path.read_bytes()
step('seed-replay',[str(RESULTS/'journal-seeds'),str(path)])
assert path.read_bytes()==data
assert (RESULTS/'seed-record.log').read_bytes()==(RESULTS/'seed-replay.log').read_bytes()
step('seed-info',[sys.executable,str(BUILD/'tools/journal-info.py'),'--json',str(path)])
assert json.loads((RESULTS/'seed-info.log').read_text())['contexts'][0]['initial_rng']==dict(
    lanczos=0x1a6b50f924d38ce7,cofactor=0x8e721c06a9b5d43f,polynomial=0xc3f9821d4076ab5e)
step('build-attempts',['cc','-O1','-DSTANDALONE','-o',str(RESULTS/'journal-attempts'),
                       str(FIXTURES/'journal-attempts.c'),str(BUILD/'lanczos.c'),
                       str(BUILD/'pr32_iter.c'),str(ROOT/'squfof126.c'),'-lgmp','-lm'])
step('event-headers',[sys.executable,'-u',str(FIXTURES/'test-event-headers.py'),str(RESULTS/'journal-attempts')])
step('solve-flush',[sys.executable,'-u',str(FIXTURES/'test-solve-flush.py'),str(RESULTS/'journal-attempts')])
step('run-clock',[sys.executable,'-u',str(FIXTURES/'test-journal-timing.py'),str(RESULTS/'journal-attempts')])
step('codecs',[str(RESULTS/'journal-attempts'),'--codecs'])
path=RESULTS/'attempts.journal'
if path.exists():path.unlink()
step('attempt-record',[str(RESULTS/'journal-attempts'),str(path)])
data=path.read_bytes()
step('attempt-replay',[str(RESULTS/'journal-attempts'),str(path)])
assert path.read_bytes()==data
step('attempt-info',[sys.executable,str(BUILD/'tools/journal-info.py'),'--json',str(path)])
step('build-graph',['cc','-O2','-DSTANDALONE','-DPSIQS','-pthread',
                    '-o',str(RESULTS/'journal-graph'),str(FIXTURES/'journal-graph.c'),
                    str(BUILD/'lanczos.c'),str(BUILD/'pr32_iter.c'),
                    str(ROOT/'squfof126.c'),'-lgmp','-lm'])
path=RESULTS/'graph.journal'
if path.exists():path.unlink()
proof=[]
for mode in ('record','replay'):
    step('graph-'+mode,[str(RESULTS/'journal-graph'),str(path)])
    lines=(RESULTS/('graph-'+mode+'.log')).read_text().splitlines()
    proof.append(json.loads(next(line for line in lines if line.startswith('{'))))
assert proof[0]==proof[1]
step('graph-info',[sys.executable,str(BUILD/'tools/journal-info.py'),'--json','--verify-relations',str(path)])
data=path.read_bytes(); at=data.index(b'\n--\n')+4+16
while at<len(data):
    _,kind,_,length,_,_,_,_,_,_=struct.unpack_from('<IIQIIQQQII',data,at)
    if kind==4:break
    at+=56+length
prefix=RESULTS/'graph-without-result.journal';prefix.write_bytes(data[:at])
n='57896044618658097711785492504343953899242261795684735710927136105699223797919'
step('graph-solve',[str(EXE),'-v','-t','1','--journal',str(prefix),n])
assert '170141183460469231731687303715884105727 340282366920938463463374607431768211297' in (RESULTS/'graph-solve.log').read_text()
(RESULTS/'graph-proof.json').write_text(json.dumps({'live':proof[0],'replayed':proof[1],
    'solved_without_saved_partition':True},indent=2)+'\n')
path=RESULTS/'wide-header.journal'
if path.exists():path.unlink()
step('wide-header-record',[str(RESULTS/'journal-attempts'),str(path),'--wide'])
data=path.read_bytes()
step('wide-header-replay',[str(RESULTS/'journal-attempts'),str(path),'--wide'])
assert path.read_bytes()==data
step('wide-header-info',[sys.executable,str(BUILD/'tools/journal-info.py'),'--json',str(path)])
step('existing-check-siqs',['make','check-siqs'],BUILD)
step('no-journal',[sys.executable,'-u',str(FIXTURES/'test-no-journal.py')])
step('journal-info',[sys.executable,'-u',str(FIXTURES/'test-journal-info.py')])
print('All journal, recovery and existing SIQS checks passed.',flush=True)
