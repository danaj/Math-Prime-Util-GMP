#!/usr/bin/env python3
"""Inspector checks using real journals, deliberate corruption and live writers."""
import importlib.util,json,signal,struct,subprocess,sys,tempfile,time,zlib
from pathlib import Path
sys.dont_write_bytecode=True
from _support import ROOT, RESULTS, EXE, TOOL
spec=importlib.util.spec_from_file_location('journal_info',TOOL)
info=importlib.util.module_from_spec(spec);spec.loader.exec_module(info)
ONE=((1<<89)-1)*((1<<107)-1);TWO=((1<<127)-1)*((1<<128)-159);passed=0

def check(name):
 global passed
 passed+=1;print('PASS '+name,flush=True)

def inspect(path,verify=False):
 reader=info.Journal(path,verify)
 try:reader.scan();return reader.summary(300)
 finally:reader.close()

def cli(path,*args,ok=True):
 p=subprocess.run([sys.executable,str(TOOL),*args,str(path)],capture_output=True,text=True,timeout=60)
 assert (p.returncode==0)==ok,(args,p.stdout,p.stderr)
 return p

def frames(data):
 at=data.index(b'\n--\n')+20;result=[]
 while at<len(data):
  fields=list(info.FRAME.unpack_from(data,at));end=at+info.FRAME.size+fields[3]
  result.append((at,end,fields));at=end
 return result

def rewrite(data,frame,changes=None,payload=None):
 start,end,fields=frame;fields=fields.copy()
 for key,value in (changes or {}).items():fields[key]=value
 if payload is None:payload=data[start+info.FRAME.size:end]
 fields[3]=len(payload);fields[-2:]=[0,0]
 header=bytearray(info.FRAME.pack(*fields))
 struct.pack_into('<I',header,48,zlib.crc32(header[:48]))
 struct.pack_into('<I',header,52,zlib.crc32(payload,zlib.crc32(header[:52])))
 return data[:start]+header+payload+data[end:]

with tempfile.TemporaryDirectory(prefix='journal-info-',dir=RESULTS) as directory:
 tmp=Path(directory);path=tmp/'complete.journal'
 subprocess.run([str(EXE),'-t','1','--journal',str(path),str(ONE)],stdout=subprocess.DEVNULL,check=True)
 data=path.read_bytes();ff=frames(data);batches=[v for v in ff if v[2][1]==3]
 result=inspect(path,True)
 assert result['relation_count']==result['relation_congruences_verified']>0
 assert set(result['contexts'][0]['initial_rng'])=={'lanczos','cofactor','polynomial'}
 assert len(result['runs'])==1 and result['runs'][0]['id']==0
 assert sum(result['byte_breakdown'].values())==len(data) and path.read_bytes()==data
 check('inspection, exact byte accounting and full mathematical verification')
 assert json.loads(cli(path,'--json','--verify-relations').stdout)['relation_count']==result['relation_count']
 text=cli(path,'--eta').stdout
 assert 'N: '+str(ONE) in text and ' full, ' in text and 'Current target:' in text
 assert 'ETA:' not in text and text.rstrip().endswith('Collection completed.')
 assert result['job_end']['complete'] and result['job_end']['partition']['n']==str(ONE)
 check('JSON/text parameters, compact final counters and completion')
 assert cli(path,'--follow','--eta','--interval','.01').stdout.rstrip().endswith('Collection completed.')
 assert json.loads(cli(path,'--follow','--json','--interval','.01').stdout)['job_end']['complete']
 check('text and JSON followers exit on the final job record')
 stream=tmp/'stream.journal';stream.write_bytes(data[:batches[0][1]])
 with (tmp/'progress.log').open('w') as output:
  follower=subprocess.Popen([sys.executable,str(TOOL),'--follow','--eta','--interval','.01',str(stream)],stdout=output,stderr=output)
  try:
   def await_text(fragment):
    deadline=time.monotonic()+10
    while time.monotonic()<deadline:
     content=(tmp/'progress.log').read_text()
     if fragment(content):return content
     assert follower.poll() is None,content
     time.sleep(.01)
    raise AssertionError('missing progress output')
   await_text(lambda content:'ETA:' in content)
   with stream.open('ab') as out:out.write(data[batches[0][1]:batches[1][1]])
   await_text(lambda content:content.count('ETA:')>=2)
   final=ff[-1];assert final[2][1]==6
   with stream.open('ab') as out:out.write(data[batches[1][1]:final[0]])
   await_text(lambda content:f"{result['contexts'][-1]['snapshot']['full_relations']:,} full" in content)
   assert follower.poll() is None
   with stream.open('ab') as out:out.write(data[final[0]:])
   assert follower.wait(timeout=10)==0
  finally:
   if follower.poll() is None:follower.terminate();follower.wait(timeout=10)
 progress=(tmp/'progress.log').read_text()
 assert progress.count('Journal: ')==1 and progress.count('N: ')==1
 assert progress.count('ETA:')>=2 and progress.rstrip().endswith('Collection completed.')
 final_lines=progress.rstrip().split('\n\n')[-1].splitlines()
 assert len(final_lines)==3 and final_lines[1].startswith('recent rates:')
 check('live compact updates print the header once and wait for job completion')
 for label,changed,success in [('failed',rewrite(data,ff[-1],payload=b'\0\0\0\0'+data[ff[-1][0]+info.FRAME.size+4:ff[-1][1]]),True),('bad-status',rewrite(data,ff[-1],payload=b'\x02\0\0\0'+data[ff[-1][0]+info.FRAME.size+4:ff[-1][1]]),False)]:
  altered=tmp/(label+'.journal');altered.write_bytes(changed)
  report=cli(altered,'--follow','--eta','--interval','.01',ok=success)
  if success:assert 'Job ended without a complete factorization.' in report.stdout
 extra=rewrite(data,ff[-1],{2:len(ff)})[ff[-1][0]:]
 after=tmp/'after-end.journal';after.write_bytes(data+extra)
 assert 'frame after final job record' in cli(after,ok=False).stderr
 check('failed jobs finish watching; malformed or nonterminal job markers rejected')
 context_frame=next(v for v in ff if v[2][1]==2)
 for frame,offset in ((context_frame,-8),(batches[0],96)):
  payload=bytearray(data[frame[0]+info.FRAME.size:frame[1]])
  struct.pack_into('<Q',payload,len(payload)+offset if offset<0 else offset,0)
  bad=tmp/'zero-rng.journal';bad.write_bytes(rewrite(data,frame,payload=payload))
  assert 'invalid zero RNG state' in cli(bad,ok=False).stderr
  failed=subprocess.run([str(EXE),'--journal',str(bad)],capture_output=True,text=True,timeout=10)
  assert failed.returncode!=0 and 'invalid zero RNG state' in failed.stderr
 check('zero RNG states rejected in starting and continuation checkpoints')
 for value in (1,101):
  small=tmp/('small-'+str(value)+'.journal')
  subprocess.run([str(EXE),'--journal',str(small),str(value)],stdout=subprocess.DEVNULL,check=True)
  item=inspect(small,True)
  assert item['relation_count']==0 and not item['contexts'] and item['returned_partitions']
 check('zero-relation journals for one and primes validate without contexts')
 restart=tmp/'restart.journal';restart.write_bytes(data[:batches[2][1]])
 subprocess.run([str(EXE),'-t','4','--journal',str(restart)],stdout=subprocess.DEVNULL,check=True)
 resumed=restart.read_bytes();rf=frames(resumed)
 assert [v['id'] for v in inspect(restart)['runs']]==[0,1]
 cutoff=next(v[1] for v in rf if v[2][5]==1 and v[2][1]==3)
 restart.write_bytes(resumed[:cutoff])
 subprocess.run([str(EXE),'-t','1','--journal',str(restart)],stdout=subprocess.DEVNULL,check=True)
 final=restart.read_bytes()
 assert [v['id'] for v in inspect(restart,True)['runs']]==[0,1,2]
 subprocess.run([str(EXE),'--journal',str(restart)],stdout=subprocess.DEVNULL,check=True)
 assert restart.read_bytes()==final
 check('resumed writes increment run; completed replay changes no bytes or run')
 for label,new_run in [('skip',2),('backwards',0)]:
  chosen=[v for v in frames(resumed) if v[2][5]==1][label=='backwards']
  altered=rewrite(resumed,chosen,{5:new_run});bad=tmp/(label+'.journal');bad.write_bytes(altered)
  cli(bad,ok=False)
  p=subprocess.run([str(EXE),'--journal',str(bad)],capture_output=True,text=True)
  assert p.returncode!=0 and 'run counter' in p.stderr and bad.read_bytes()==altered
 check('C/Python reject skipped or decreasing runs without changing file')
 altered=rewrite(data,batches[3],{6:0,7:0});clock=tmp/'clock.journal';clock.write_bytes(altered)
 summary=json.loads(cli(clock,'--json','--verify-relations').stdout)
 assert any('UTC' in v for v in summary['warnings']) and any('elapsed' in v for v in summary['warnings'])
 clock_prefix=tmp/'clock-prefix.journal';clock_prefix.write_bytes(altered[:batches[4][1]])
 assert inspect(clock_prefix)['contexts'][0]['rate'] is None
 subprocess.run([str(EXE),'--journal',str(clock)],stdout=subprocess.DEVNULL,check=True)
 assert clock.read_bytes()==altered
 check('non-monotonic clocks warn; C replay and mathematical checks remain valid')
 partial=tmp/'partial.journal';partial.write_bytes(data[:batches[2][1]+11]);summary=inspect(partial)
 assert summary['tail'] and summary['validated_bytes']==batches[2][1]
 cli(partial,'--strict',ok=False)
 check('incomplete tail preserves prefix; strict mode rejects it')
 changed=bytearray(data);changed[batches[0][0]+info.FRAME.size+8]^=1
 bad=tmp/'crc.journal';bad.write_bytes(changed);cli(bad,ok=False)
 payload=bytearray(data[batches[0][0]+info.FRAME.size:batches[0][1]])
 struct.pack_into('<Q',payload,16*8,struct.unpack_from('<Q',payload,16*8)[0]+1)
 bad=tmp/'histogram.journal';bad.write_bytes(rewrite(data,batches[0],payload=payload));cli(bad,ok=False)
 check('bad CRC and inconsistent histogram rejected')
 payload=bytearray(data[batches[0][0]+info.FRAME.size:batches[0][1]])
 width=(result['contexts'][0]['parameters']['large_prime_bound'].bit_length()+7)//8;at=192
 while True:
  kind,labels=struct.unpack_from('BB',payload,at);at+=2
  if kind==2:at+=4*labels
  else:break
 assert kind==1
 at+=labels*width;size=struct.unpack_from('<I',payload,at)[0];assert size>0
 payload[at+4+size-1]^=1
 bad=tmp/'congruence.journal';bad.write_bytes(rewrite(data,batches[0],payload=payload));cli(bad)
 assert 'congruence' in cli(bad,'--verify-relations',ok=False).stderr
 check('optional verification catches checksummed mathematical corruption')
 eta_data=data[:batches[8][1]]
 for frame in frames(eta_data):
  seq=frame[2][2];eta_data=rewrite(eta_data,frame,{6:1700000000000000+seq*1000000,7:seq*1000000})
 eta=tmp/'eta.journal';eta.write_bytes(eta_data);ctx=inspect(eta)['contexts'][0];rate=ctx['rate']
 assert rate and rate['run']==0 and rate['full_per_second']>0
 remaining=max(0,ctx['snapshot']['collection_target']-ctx['snapshot']['full_relations'])
 assert abs(rate['collection_eta_seconds']-remaining/rate['full_per_second'])<1e-9
 check('ETA uses same-run/context monotonic progress')
 live=tmp/'live.journal'
 with (tmp/'collector.log').open('w') as collector_log,(tmp/'follow.log').open('w') as follow_log:
  collector=subprocess.Popen([str(EXE),'-t','1','--journal',str(live),str(TWO)],stdout=collector_log,stderr=collector_log)
  follower=None
  try:
   deadline=time.monotonic()+30
   while time.monotonic()<deadline:
    try:
     if inspect(live)['frames']>=3:break
    except (OSError,info.InvalidJournal,ValueError):pass
    assert collector.poll() is None;time.sleep(.02)
   else:raise AssertionError('no live prefix')
   follower=subprocess.Popen([sys.executable,str(TOOL),'--follow','--json','--eta','--interval','.02',str(live)],stdout=follow_log,stderr=follow_log)
   while time.monotonic()<deadline:
    lines=(tmp/'follow.log').read_text().splitlines()
    if len(lines)>=2:
     reports=[json.loads(v) for v in lines]
     assert reports[-1]['frames']>reports[0]['frames'] and collector.poll() is None
     break
    assert follower.poll() is None;time.sleep(.02)
   else:raise AssertionError('no incremental output')
  finally:
   if collector.poll() is None:collector.terminate()
   collector.wait(timeout=10)
   if follower is not None and follower.poll() is None:
    follower.send_signal(signal.SIGINT);assert follower.wait(timeout=10)==0
 check('live follow reads actively locked writer and reports incremental progress')
print(f'All {passed}/{passed} journal-info checks passed.',flush=True)
