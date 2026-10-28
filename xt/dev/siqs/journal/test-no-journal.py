#!/usr/bin/env python3
"""Check explicit standalone and automatic Perl journal exclusion."""
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile

import sys
sys.dont_write_bytecode=True
from _support import ROOT, BUILD, RESULTS
SOURCES=['mpu-siqs.c','siqs.c','lanczos.c','pr32_iter.c','../squfof126.c']
N=((1<<89)-1)*((1<<107)-1)
FACTORS=[(1<<89)-1,(1<<107)-1]


def build_variant(name, flags):
    with (RESULTS/('build-'+name+'.log')).open('w') as log:
        subprocess.run(['perl','tools/build-siqs.pl','--cppflags','-DSIQS_NO_JOURNAL',
                        '--output',str(RESULTS/name)]+flags+SOURCES,cwd=BUILD,
                       stdout=log,stderr=subprocess.STDOUT,check=True)


def check_binary(binary, threads):
    for args,stdin in [([str(N)],None),([],str(N)+'\n')]:
        p=subprocess.run([str(binary)]+(['-t',str(threads)] if threads else [])+args,input=stdin,
                         capture_output=True,text=True,timeout=60)
        assert p.returncode==0,(binary,p.stdout,p.stderr)
        assert p.stdout.strip()==str(N)+': '+' '.join(map(str,FACTORS))
    with tempfile.TemporaryDirectory(prefix='no-journal-',dir=RESULTS) as tmp:
        path=Path(tmp)/'test.journal'
        for existing in (False,True):
            if existing:path.write_bytes(b'leave this file unchanged\n')
            old=path.read_bytes() if existing else None
            p=subprocess.run([str(binary),'--journal',str(path)]+([] if existing else [str(N)]),
                             capture_output=True,text=True,timeout=5)
            assert p.returncode!=0 and 'journal support disabled' in p.stderr+p.stdout
            assert path.read_bytes()==old if existing else not path.exists()


build_variant('msiqs-no-journal',[])
threaded='with pthreads' in (RESULTS/'build-msiqs-no-journal.log').read_text()
check_binary(RESULTS/'msiqs-no-journal',4 if threaded else 0)
print('PASS standalone no-journal factoring and disabled-option diagnostics',flush=True)
build_variant('msiqs-no-journal-serial',['--serial'])
check_binary(RESULTS/'msiqs-no-journal-serial',0)
print('PASS explicit serial no-journal build',flush=True)
build_variant('msiqs-no-journal-sanitize',[
    '--serial','--cflags','-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer',
    '--ldflags','-fsanitize=address,undefined'])
check_binary(RESULTS/'msiqs-no-journal-sanitize',0)
print('PASS ASan/UBSan no-journal factoring and rejection paths',flush=True)

perl_flags=shlex.split(subprocess.check_output(
    ['perl','-MExtUtils::Embed','-e','ccopts'],text=True))
for name,flags in [('standalone',['-DSTANDALONE','-DSIQS_NO_JOURNAL']),('perl',perl_flags)]:
    target=RESULTS/('siqs-no-journal-'+name+'.o')
    subprocess.run(['cc','-O2','-c']+flags+['siqs.c','-o',str(target)],cwd=BUILD,check=True)
    macros=subprocess.check_output(['cc','-dM','-E']+flags+['siqs.c'],cwd=BUILD,text=True)
    assert '#define SIQS_NO_JOURNAL' in macros
    text=subprocess.check_output(['cc','-E']+flags+['siqs.c'],cwd=BUILD,text=True)
    assert 'SIQS journal: frame header checksum mismatch' not in text
    assert 'SIQS journal support disabled in this build' in text
    if shutil.which('nm'):
        symbols=subprocess.check_output(['nm','-u',str(target)],text=True)
        names={line.split()[-1].lstrip('_') for line in symbols.splitlines() if line.split()}
        forbidden={'fopen','fread','fseeko','ftruncate','fcntl',
                   'clock_gettime','gettimeofday','time'}
        assert not names&forbidden,(name,names&forbidden)
    print('PASS '+name+' preprocessing/object excludes journal I/O and clocks',flush=True)
print('All no-journal checks passed.',flush=True)
