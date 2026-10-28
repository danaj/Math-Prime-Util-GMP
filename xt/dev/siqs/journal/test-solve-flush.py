#!/usr/bin/env python3
"""A tiny real matrix solve must seal its pending journal batch."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile

sys.dont_write_bytecode=True
from _support import RESULTS, TOOL

spec=importlib.util.spec_from_file_location('journal_info',TOOL)
info=importlib.util.module_from_spec(spec);spec.loader.exec_module(info)
fixture=Path(sys.argv[1]) if len(sys.argv)>1 else RESULTS/'journal-attempts'
with tempfile.TemporaryDirectory(prefix='solve-flush-',dir=RESULTS) as directory:
    path=Path(directory)/'solve.journal'
    command=[str(fixture),'--solve-events',str(path)]
    recorded=subprocess.run(command,capture_output=True,text=True,timeout=10)
    assert recorded.returncode==0,recorded.stdout+recorded.stderr
    counts=json.loads(recorded.stdout);data=path.read_bytes()
    assert counts['full']==64
    reader=info.Journal(path,True)
    try:
        reader.scan();summary=reader.summary(300)
        assert reader.last_type==3
    finally:reader.close()
    assert summary['relation_congruences_verified']==65
    assert summary['contexts'][0]['snapshot']['full_relations']==64
    replayed=subprocess.run(command,capture_output=True,text=True,timeout=10)
    assert replayed.returncode==0,replayed.stdout+replayed.stderr
    assert replayed.stdout==recorded.stdout and path.read_bytes()==data
print('PASS solver flushes a real pending batch; all 65 congruences verify and replay preserves bytes')
