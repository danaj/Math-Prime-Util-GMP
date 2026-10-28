#!/usr/bin/env python3
"""Read-only SIQS journal validation, inspection and live collection estimates."""
import argparse
from collections import Counter, deque
from datetime import datetime, timedelta, timezone
import json
import math
import os
from pathlib import Path
import struct
import sys
import time
import zlib

FRAME = struct.Struct('<IIQIIQQQII')
PREAMBLE = b'MPUSIQJ\0'+struct.pack('<II',1,0x01020304)
MAX_FRAME = 64*1024*1024
RELATION_NAMES = ('smooth','one_lp','two_lp','three_lp','four_lp')


def recent_samples(samples, window):
    samples=list(samples)
    if len(samples)<2:return samples
    minimum=samples[-1][0]-window*1e6
    first=next((i for i,v in enumerate(samples[:-1]) if v[0]>=minimum),len(samples)-2)
    return samples[first:]


def collection_eta(samples, remaining, window=300, model='quadratic', horizon=None):
    """Forecast from full counts; measured display rates remain unchanged."""
    samples=recent_samples(samples,window)
    result=dict(seconds=None,linear_seconds=None,model='linear',fit_points=0)
    if remaining<=0:
        result.update(seconds=0.0,linear_seconds=0.0)
        return result
    if len(samples)<2:return result
    first,last=samples[0],samples[-1]
    span=(last[0]-first[0])/1e6
    if span<=0:return result
    average=(last[2]-first[2])/span
    if average<=0:return result
    result.update(seconds=remaining/average,linear_seconds=remaining/average)
    if model=='linear' or len(samples)<6 or span<min(60.0,window):return result

    # Equal-duration bins prevent dense frame bursts dominating OLS. Average
    # x, x*x and y separately so binning preserves an exact quadratic.
    bins={}
    for timestamp,_,full in samples:
        x=(timestamp-last[0])/(last[0]-first[0])
        index=min(29,max(0,int((x+1)*30)))
        sums=bins.setdefault(index,[0.0,0.0,0.0,0])
        sums[0]+=x;sums[1]+=x*x;sums[2]+=full-last[2];sums[3]+=1
    points=[(x/n,q/n,y/n) for x,q,y,n in bins.values()]
    result['fit_points']=len(points)
    if len(points)<6:return result
    mx,mq,my=(math.fsum(p[i] for p in points)/len(points) for i in range(3))
    xx=math.fsum((x-mx)**2 for x,_,_ in points)
    if xx<=1e-12:return result
    xy=math.fsum((x-mx)*(y-my) for x,_,y in points)
    xq=math.fsum((x-mx)*(q-mq) for x,q,_ in points)
    q_slope=xq/xx
    # Orthogonalize x*x against 1 and x instead of inverting normal equations.
    orth=[(q-mq-q_slope*(x-mx),y-my) for x,q,y in points]
    qq=math.fsum(q*q for q,_ in orth)
    if qq<=1e-12:return result
    a=math.fsum(q*y for q,y in orth)/qq
    b=xy/xx-a*q_slope
    if not all(map(math.isfinite,(a,b))) or b<=0 or a<=1e-10*max(1.0,abs(b)):
        return result
    velocity=b/span;coefficient=a/(span*span)
    limit=min(span,window) if horizon is None else min(span,window,horizon)
    if limit<=0:return result
    covered=velocity*limit+coefficient*limit*limit
    if remaining<=covered:
        seconds=2*remaining/(velocity+math.sqrt(velocity*velocity+4*coefficient*remaining))
    else:
        # Acceleration is local: extend it by at most one observed window,
        # then hold the resulting rate constant for the remaining forecast.
        seconds=limit+(remaining-covered)/(velocity+2*coefficient*limit)
    if math.isfinite(seconds) and seconds>=0:
        result.update(seconds=seconds,model='quadratic',current_rate=velocity,
                      quadratic_coefficient=coefficient,acceleration_horizon_seconds=limit,
                      fit_span_seconds=span)
    return result
UNSIGNED = '''bits policy_first_bits policy_last_bits max_large_primes fb_size
half_interval q_count poly_d target_relations sieve_start sieve_free_units
sieve_start_prime_floor fb_floor relation_extra multiplier_refine_divisor
fixed_half_interval sieve_initial stage1_bias sieve_raw_scores sieve_saturating
sieve_score_bound large_prime_bound smooth_bound'''.split()
DOUBLES = '''lp_multiplier residual_multiplier sieve_start_index_exponent
fb_coefficient interval_scale sieve_hit_exponent sieve_hit_bound_nominal
sieve_hit_residual_multiplier sieve_hit_bound log_scale'''.split()
PROGRESS = '''families polynomials candidates split_attempts split_square split_siqs
split_squfof split_rho split_fail split_rejected la_rng cofactor_rng poly_rng
full_relations retained_partials collection_target'''.split()
PROGRESS += ['length_'+str(i) for i in range(1,8)]+['length_8plus']


class InvalidJournal(Exception):
    pass


class FrameChecksumError(InvalidJournal):
    """A live reader may retry this frame once before reporting corruption."""


def require(condition, message):
    if not condition:
        raise InvalidJournal(message)


class Payload:
    def __init__(self, data):
        self.data=data
        self.at=0

    def take(self, count):
        require(0<=count<=len(self.data)-self.at,'truncated payload field')
        result=self.data[self.at:self.at+count]
        self.at+=count
        return result

    def u32(self):
        return struct.unpack('<I',self.take(4))[0]

    def u8(self):
        return self.take(1)[0]

    def u64(self):
        return struct.unpack('<Q',self.take(8))[0]

    def magnitude(self, maximum=MAX_FRAME):
        size=self.u32()
        require(size<=maximum,'oversized magnitude')
        raw=self.take(size)
        require(not raw or raw[0]!=0,'noncanonical magnitude')
        return int.from_bytes(raw,'big')

    def varint(self):
        value=0
        for shift in range(0,35,7):
            byte=self.take(1)[0]
            require(shift!=28 or byte<=15,'ULEB128 overflow')
            value|=(byte&127)<<shift
            if byte<128:
                require(shift==0 or byte!=0,'noncanonical ULEB128')
                return value
        raise InvalidJournal('unterminated ULEB128')

    def done(self):
        require(self.at==len(self.data),'unexpected trailing payload')


def utc_text(value):
    try:
        date=datetime.fromtimestamp(value//1000000,timezone.utc)+timedelta(microseconds=value%1000000)
        return date.isoformat(timespec='microseconds')
    except (ValueError,OverflowError,OSError):
        return 'outside supported date range'


class Journal:
    def __init__(self, path, verify=False):
        self.path=Path(path)
        self.file=self.path.open('rb')
        self.verify=verify
        self.sequence=0
        self.contexts=[]
        self.active=None
        self.primes=[]
        self.runs=[]
        self.records=Counter()
        self.verified=0
        self.factor_entries=0
        self.warnings=[]
        self.warning_count=0
        self.tail=None
        self.previous_clock=None
        self.last_type=None
        self.returned=[]
        self.job_end=None
        self.bytes=Counter()
        self.read_prefix()
        self.committed=self.file.tell()
        self.bytes['prefix']=self.committed

    def close(self):
        self.file.close()

    def warn(self, message):
        self.warning_count+=1
        if message not in self.warnings and len(self.warnings)<100:
            self.warnings.append(message)

    def read_prefix(self):
        def line():
            result=self.file.readline(MAX_FRAME+1)
            require(len(result)<=MAX_FRAME and result.endswith(b'\n'),'incomplete or oversized ASCII header')
            return result
        require(line()==b'MPUSIQJ\n','invalid ASCII magic')
        require(line()==b'Version: 1\n','unsupported ASCII version')
        require(line()==b'Byte-order: little-endian\n','unsupported ASCII byte order')
        number=line()
        require(number.startswith(b'N: ') and number[3:-1].isdigit(),'invalid ASCII input')
        self.n=int(number[3:-1])
        require(self.n>0,'input must be positive')
        size=0
        while True:
            note=line();size+=len(note)
            require(size<=MAX_FRAME,'optional ASCII header exceeds size limit')
            if note==b'--\n':break
        require(self.file.read(16)==PREAMBLE,'invalid or incomplete binary preamble')

    def context(self, r, count):
        require(count==1,'invalid context item count')
        require(self.active is None or self.active['ended'],'new context before previous attempt ended')
        n=r.magnitude()
        require(n>0 and self.n%n==0,'context cofactor does not divide original input')
        multiplier=r.u64()
        require(1<=multiplier<=2047,'invalid multiplier')
        require(r.u32()==21,'unsupported packed row width')
        name=r.take(r.u32()).decode('ascii')
        params={name:r.u64() for name in UNSIGNED}
        for key in DOUBLES:
            value=struct.unpack('<d',r.take(8))[0]
            require(math.isfinite(value),'non-finite parameter')
            params[key]=value
        fb_count=r.u32()
        require(fb_count==params['fb_size'] and 0<fb_count<(1<<21),'invalid factor-base size')
        require(1<=params['q_count']<=16 and params['poly_d'] in (1,2),'invalid polynomial geometry')
        require(params['max_large_primes']<=4,'unsupported large-prime mode')
        primes=[]
        previous=1
        for _ in range(fb_count):
            prime,root,log=r.u32(),r.u32(),r.u32()
            require(prime>previous and root<prime,'invalid ordered factor-base entry')
            primes.append(prime);previous=prime
        initial_rng={name:r.u64() for name in ('lanczos','cofactor','polynomial')}
        require(all(initial_rng.values()),'invalid zero RNG state')
        r.done()
        self.bytes['factor_base']+=12*fb_count
        self.bytes['context_metadata']+=len(r.data)-12*fb_count
        self.primes=primes
        self.active=dict(id=len(self.contexts)+1,n=str(n),bits=n.bit_length(),
                         multiplier=multiplier,policy=name,parameters=params,
                         pmax=primes[-1],initial_rng=initial_rng,relations=Counter(),families=0,
                         snapshot=None,ended=False,samples=deque(maxlen=5000),sample_run=None)
        self.contexts.append(self.active)

    def batch(self, r, count, run, elapsed, bad_clock):
        ctx=self.active
        require(ctx is not None and not ctx['ended'],'batch outside an active context')
        snapshot={name:r.u64() for name in PROGRESS}
        require(all(snapshot[name] for name in ('la_rng','cofactor_rng','poly_rng')),
                'invalid zero RNG state')
        require(count<=(len(r.data)-r.at)//2,'invalid batch event count')
        self.bytes['batch_snapshots']+=8*len(PROGRESS)
        params=ctx['parameters'];width=max(1,(params['large_prime_bound'].bit_length()+7)//8)
        n=int(ctx['n']);magnitude_limit=(n.bit_length()+7)//8
        before=sum(ctx['relations'].values())
        for _ in range(count):
            kind,info=r.u8(),r.u8()
            if kind==2:
                q=info
                require(q==params['q_count'],'family factor count differs from context')
                indices=[r.u32() for _ in range(q)]
                require(len(set(indices))==q and all(v<len(self.primes) for v in indices),'invalid A-family indices')
                ctx['families']+=1
                self.bytes['family_records']+=2+4*q
                continue
            require(kind==1,'unknown batch event')
            labels=info
            require(labels<=4,'invalid large-prime count')
            require(labels<=params['max_large_primes'],"relation exceeds the context's large-prime count")
            large_primes=[int.from_bytes(r.take(width),'little') for _ in range(labels)]
            require(all(1<p<=params['large_prime_bound'] for p in large_primes) and
                    large_primes==sorted(large_primes),'invalid large-prime labels')
            start=r.at;y=r.magnitude(magnitude_limit)
            require(y<n,'relation y is not reduced')
            self.bytes['y_magnitude_and_length']+=r.at-start
            factors=r.u32()
            require(factors<=params['fb_size']+1 and factors<=len(r.data)-r.at,'invalid factor count')
            start=r.at;row=0;rhs=math.prod(large_primes)%n
            for i in range(factors):
                code=r.varint();delta=code>>1
                exponent=r.varint() if code&1 else 1
                require(1<=exponent<=2047 and (not (code&1) or exponent!=1),'invalid exponent')
                require(i==0 or delta>0,'rows are not strictly increasing')
                row+=delta
                require(row<=params['fb_size'],'factor row outside factor base')
                if self.verify:
                    prime=-1 if row==0 else self.primes[row-1]
                    rhs=rhs*pow(prime,exponent,n)%n
            if self.verify:
                require(y*y%n==rhs,'relation congruence failed')
                self.verified+=1
            self.bytes['factor_varints']+=r.at-start
            self.factor_entries+=factors
            self.bytes['relation_tags_and_counts']+=6
            self.bytes['large_prime_labels']+=labels*width
            kind=RELATION_NAMES[labels]
            ctx['relations'][kind]+=1;self.records[kind]+=1
        r.done()
        require(snapshot['families']==ctx['families'],'family progress differs from events')
        full=snapshot['full_relations'];hist=[snapshot[key] for key in PROGRESS[-8:]]
        require(sum(hist)==full,'cycle histogram does not sum to completed relations')
        require(full<=sum(ctx['relations'].values()),'more completed relations than raw records')
        require(snapshot['retained_partials']<=sum(ctx['relations'].values()),'invalid retained-partial count')
        require(snapshot['collection_target']<=0xffffffff,'collection target out of range')
        old=ctx['snapshot']
        if old:
            for key in ('families','polynomials','candidates','full_relations'):
                require(snapshot[key]>=old[key],'decreasing '+key+' within context')
            require(all(snapshot[k]>=old[k] for k in PROGRESS[-8:]),'decreasing cycle histogram')
        ctx['snapshot']=snapshot
        if ctx['sample_run']!=run:
            ctx['samples'].clear();ctx['sample_run']=run
            # A resumed run starts after replay, with these saved counters.
            # Older journals may start at zero, leaving no first interval.
            if not bad_clock and elapsed>0:
                ctx['samples'].append((0,before,old['full_relations'] if old else 0))
        elif bad_clock:
            ctx['samples'].clear()
        if not bad_clock:
            ctx['samples'].append((elapsed,sum(ctx['relations'].values()),full))
        self.runs[-1]['relations']+=sum(ctx['relations'].values())-before

    def partition(self, r):
        n=r.magnitude()
        require(n>0 and self.n%n==0,'returned input does not divide original input')
        count=r.u32()
        require(count<=n.bit_length() and count<=(len(r.data)-r.at)//5,'invalid partition length')
        values=[r.magnitude() for _ in range(count)]
        require(all(v>=2 or (n==1 and v==1) for v in values) and math.prod(values)==n,'invalid multiplicative partition')
        return dict(n=str(n),factors=[str(v) for v in values])

    def consume(self, kind, count, payload, run, elapsed, bad_clock):
        require(self.job_end is None,'frame after final job record')
        r=Payload(payload)
        if kind==1:
            require(self.sequence==0 and count==1,'misplaced INPUT')
            require(r.magnitude()==self.n,'ASCII and binary inputs differ');r.done()
            self.bytes['input_and_returns']+=len(payload)
        elif kind==2:
            self.context(r,count)
            # A later context shares the run clock but starts its own counts.
            self.active['sample_run']=run
            if not bad_clock:self.active['samples'].append((elapsed,0,0))
        elif kind==3:
            self.batch(r,count,run,elapsed,bad_clock)
        elif kind==4:
            require(self.active is not None and not self.active['ended'],'attempt end outside context')
            require(count==1 and r.u32() in (0,1),'invalid attempt-end flag')
            self.partition(r);r.done();self.active['ended']=True
            self.bytes['input_and_returns']+=len(payload)
        elif kind==5:
            require(count==1,'invalid RETURN item count')
            require(self.active is None or self.active['ended'],'RETURN before attempt ended')
            part=self.partition(r);r.done();self.returned.append(part)
            self.bytes['input_and_returns']+=len(payload)
        elif kind==6:
            require(count==1,'invalid job-end item count')
            require(self.active is None or self.active['ended'],'job end before attempt ended')
            complete=r.u32()
            require(complete in (0,1),'invalid job-end status')
            part=self.partition(r);r.done()
            require(part['n']==str(self.n),'job-end input differs from original input')
            self.job_end=dict(complete=bool(complete),partition=part)
            self.bytes['input_and_returns']+=len(payload)
        else:
            raise InvalidJournal('unknown frame type')

    def scan(self):
        self.tail=None
        self.file.seek(self.committed)
        require(os.fstat(self.file.fileno()).st_size>=self.committed,'file was truncated inside validated prefix')
        added=0
        while True:
            at=self.file.tell();header=self.file.read(FRAME.size)
            if not header:break
            if len(header)!=FRAME.size:
                self.tail='incomplete final frame header';self.file.seek(at);break
            magic,kind,seq,length,count,run,utc,elapsed,hcrc,crc=FRAME.unpack(header)
            if zlib.crc32(header[:48])!=hcrc:
                self.file.seek(at);raise FrameChecksumError('frame header checksum mismatch at byte '+str(at))
            require(magic==0x4b4c424a and seq==self.sequence and 1<=kind<=6,'invalid frame magic/type/sequence')
            require(self.sequence!=0 or kind==1,'first frame must be INPUT')
            require(length<=MAX_FRAME,'frame exceeds size limit')
            payload=self.file.read(length)
            if len(payload)!=length:
                self.tail='incomplete final frame payload';self.file.seek(at);break
            if zlib.crc32(payload,zlib.crc32(header[:52]))!=crc:
                self.file.seek(at);raise FrameChecksumError('frame checksum mismatch at byte '+str(at))
            if not self.runs:
                require(run==0,'first run counter must be zero')
            else:
                require(run in (self.runs[-1]['id'],self.runs[-1]['id']+1),'run counter decreased or skipped')
            if not self.runs or run!=self.runs[-1]['id']:
                self.runs.append(dict(id=run,frames=0,relations=0,first_utc_us=utc,last_utc_us=utc,elapsed_us=elapsed))
            bad_clock=False
            if self.previous_clock:
                prev_run,prev_utc,prev_elapsed=self.previous_clock
                if utc<prev_utc:self.warn('UTC timestamp moved backwards (run '+str(run)+')')
                if run==prev_run and elapsed<prev_elapsed:
                    self.warn('elapsed timestamp moved backwards (run '+str(run)+'); affected rate interval excluded')
                    bad_clock=True
                    if self.active is not None:self.active['samples'].clear()
            self.consume(kind,count,payload,run,elapsed,bad_clock)
            self.runs[-1]['frames']+=1;self.runs[-1]['last_utc_us']=utc;self.runs[-1]['elapsed_us']=elapsed
            self.previous_clock=(run,utc,elapsed)
            self.bytes['frame_headers']+=FRAME.size
            self.last_type=kind;self.sequence+=1;added+=1;self.committed=self.file.tell()
        require(self.sequence>0,'missing complete INPUT frame')
        return added

    def summary(self, window, eta_model='quadratic', eta_horizon=None):
        contexts=[]
        for ctx in self.contexts:
            item={key:value for key,value in ctx.items() if key not in ('samples','sample_run')}
            item['relations']=dict(ctx['relations'])
            rate=None
            samples=recent_samples(ctx['samples'],window)
            if len(samples)>=2:
                first,last=samples[0],samples[-1]
                seconds=(last[0]-first[0])/1e6
                if seconds>0:
                    raw=(last[1]-first[1])/seconds;full=(last[2]-first[2])/seconds
                    target=ctx['snapshot']['collection_target'] or ctx['parameters']['target_relations']
                    remaining=max(0,target-last[2])
                    prediction=collection_eta(samples,remaining,window,eta_model,eta_horizon)
                    rate=dict(run=ctx['sample_run'],window_seconds=seconds,raw_per_second=raw,
                              full_per_second=full,eta=prediction,
                              collection_eta_seconds=prediction['seconds'] if not ctx['ended'] else None)
            item['rate']=rate
            contexts.append(item)
        runs=[dict(v,first_utc=utc_text(v['first_utc_us']),last_utc=utc_text(v['last_utc_us'])) for v in self.runs]
        size=os.fstat(self.file.fileno()).st_size;total=sum(self.records.values())
        require(sum(self.bytes.values())==self.committed,'internal byte accounting mismatch')
        return dict(path=str(self.path),n=str(self.n),bits=self.n.bit_length(),file_bytes=size,
                    validated_bytes=self.committed,frames=self.sequence,tail=self.tail,
                    relations=dict(self.records),relation_count=total,
                    factor_entries=self.factor_entries,
                    bytes_per_relation=self.committed/total if total else None,
                    byte_breakdown=dict(self.bytes),contexts=contexts,runs=runs,
                    relation_congruences_verified=self.verified,warnings=self.warnings,
                    warning_count=self.warning_count,last_frame_type=self.last_type,
                    returned_partitions=self.returned,job_end=self.job_end)


def duration(seconds, rough=False):
    if seconds is None:return 'insufficient progress'
    if rough and seconds<3600:
        scale,unit=(1,'second') if seconds<60 else (60,'minute')
        value=max(1,round(seconds/scale))
        return f'{value} {unit}'+('s' if value!=1 else '')
    if seconds<60:return f'{seconds:.1f} seconds'
    if seconds<3600:return f'{seconds/60:.1f} minutes'
    return f'{seconds/3600:.1f} hours'


def collection_eta_text(ctx, window=300):
    rate=ctx['rate'] if ctx else None
    seconds=rate['collection_eta_seconds'] if rate else None
    if seconds is None:return 'waiting for more data'
    if seconds<=0:return duration(seconds)
    rough=rate['window_seconds']<min(60.0,window)
    snapshot=ctx['snapshot']
    if ctx['parameters']['max_large_primes']>0 and snapshot:
        # A display cue for the early graph, not a confidence bound or model input.
        full=snapshot['full_relations']
        rough=rough or 2*(full-snapshot['length_1'])<full
    return ('roughly ' if rough else '')+duration(seconds,rough)


def display(info, eta, window=300):
    print('Journal: '+info['path'])
    print('N: '+info['n'])
    print(f"{info['bits']} bits; {info['frames']:,} valid frames; writing runs: {len(info['runs'])}; {info['validated_bytes']:,}/{info['file_bytes']:,} validated bytes")
    counts=info['relations']
    breakdown=', '.join(f"{counts.get(name,0):,} {label}" for i,(name,label) in
                       enumerate(zip(RELATION_NAMES,('smooth','1LP','2LP','3LP','4LP')))
                       if i<=2 or counts.get(name,0))
    print(f"Stored raw relations: {info['relation_count']:,} ({breakdown})")
    if info['bytes_per_relation'] is not None:print(f"Whole-prefix bytes/relation: {info['bytes_per_relation']:.2f}")
    print('Checksums and structure: OK for the complete prefix')
    if info['relation_congruences_verified']:print(f"Verified congruences: {info['relation_congruences_verified']:,} OK")
    for ctx in info['contexts']:
        p=ctx['parameters'];s=ctx['snapshot']
        print(f"Context {ctx['id']}: {ctx['bits']} bits, {ctx['policy']}, k={ctx['multiplier']}, FB={p['fb_size']:,} (pmax={ctx['pmax']}), M={p['half_interval']}, q={p['q_count']}, d={p['poly_d']}")
        print(f"  LP={p['large_prime_bound']}, residual={p['smooth_bound']}, K_L={p['lp_multiplier']:g}, K_R={p['residual_multiplier']:g}, initial target={p['target_relations']:,}" + ('; attempt ended' if ctx['ended'] else ''))
        if s:
            print(f"  Full={s['full_relations']:,}, retained partials={s['retained_partials']:,}, current target={s['collection_target']:,}, families={s['families']:,}, polys={s['polynomials']:,}")
            print('  Lengths (raw records/full), 1–7,8+: '+', '.join(str(s[k]) for k in PROGRESS[-8:]))
        rate=ctx['rate']
        if rate:
            print(f"  Run {rate['run']} recent rates: {rate['raw_per_second']:.3f} raw/s, {rate['full_per_second']:.3f} full/s over {rate['window_seconds']:.1f}s")
        if eta and not ctx['ended']:
            print('  Collection ETA: '+collection_eta_text(ctx,window))
    for run in info['runs']:
        print(f"Run {run['id']}: {run['frames']:,} frames, {run['relations']:,} raw relations; {run['first_utc']} to {run['last_utc']}; elapsed {run['elapsed_us']/1e6:.3f}s")
    if info['tail']:print('Tail: '+info['tail']+'; validated prefix remains usable')
    for warning in info['warnings']:print('Warning: '+warning)


def display_eta(info, previous, window=300):
    ctx=info['contexts'][-1] if info['contexts'] else None
    snapshot=ctx['snapshot'] if ctx else None
    target=(snapshot['collection_target'] if snapshot else ctx['parameters']['target_relations']) if ctx else None
    context=ctx['id'] if ctx else None
    if previous is None:
        print('Journal: '+info['path'])
        print('N: '+info['n'])
        print(f"{info['bits']} bits; {info['frames']:,} valid frames; {info['validated_bytes']:,}/{info['file_bytes']:,} validated bytes; writing runs: {len(info['runs'])}")
        print(f"Stored raw relations: {info['relation_count']:,}; checksums and structure OK")
    if ctx and (previous is None or previous['context']!=context):
        p=ctx['parameters']
        print(f"Context {context}: {ctx['policy']}, k={ctx['multiplier']}, FB={p['fb_size']:,}, M={p['half_interval']:,}, q={p['q_count']}, d={p['poly_d']}")
        print(f"LP={p['large_prime_bound']:,}, residual={p['smooth_bound']:,}, K_L={p['lp_multiplier']:g}, K_R={p['residual_multiplier']:g}")
    if previous is None or previous['context']!=context or previous['target']!=target:
        print(f"Current target: {target:,} relations" if target is not None else 'Current target: no collection context')
    if previous is None or previous['warnings']!=info['warning_count']:
        for warning in info['warnings']:print('Warning: '+warning)
    if info['tail'] and (previous is None or previous['tail']!=info['tail']):
        print('Tail: '+info['tail']+'; validated prefix remains usable')
    print()
    s=snapshot or dict(full_relations=0,retained_partials=0,families=0,polynomials=0)
    print(f"{s['full_relations']:,} full, {s['retained_partials']:,} partials    {s['families']:,} families, {s['polynomials']:,} polys")
    rate=ctx['rate'] if ctx else None
    if rate:
        print(f"recent rates: {rate['raw_per_second']:.3f} raw/s, {rate['full_per_second']:.3f} full/s over {rate['window_seconds']:.1f}s")
    else:print('recent rates: waiting for more data')
    if info['job_end'] is None and not (ctx and ctx['ended']):
        seconds=rate['collection_eta_seconds'] if rate else None
        print('ETA: '+collection_eta_text(ctx,window)+' until collection complete' if seconds is not None else 'ETA: waiting for more data')
    return dict(context=context,target=target,warnings=info['warning_count'],tail=info['tail'])


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('journal',type=Path)
    ap.add_argument('--verify-relations',action='store_true',help='also check every relation congruence')
    ap.add_argument('--json',action='store_true',help='emit structured results; JSON lines in follow mode')
    ap.add_argument('--eta',action='store_true',help='show compact progress with a collection ETA')
    ap.add_argument('--window',type=float,default=300,help='recent rate window in seconds (default: 300)')
    ap.add_argument('--eta-model',choices=('quadratic','linear'),default='quadratic',
                    help='ETA forecast model (default: quadratic with linear fallback)')
    ap.add_argument('--eta-horizon',type=float,
                    help='maximum seconds of future acceleration (default: one observed window)')
    ap.add_argument('--follow',action='store_true',help='follow appended frames until a final job record, without locking the writer')
    ap.add_argument('--interval',type=float,default=2,help='follow polling interval in seconds (default: 2)')
    ap.add_argument('--strict',action='store_true',help='reject an incomplete final frame')
    args=ap.parse_args()
    if not math.isfinite(args.window) or args.window<=0 or not math.isfinite(args.interval) or args.interval<=0:
        ap.error('window and interval must be positive finite values')
    if args.eta_horizon is not None and (not math.isfinite(args.eta_horizon) or args.eta_horizon<=0):
        ap.error('ETA horizon must be a positive finite value')
    journal=None;retry=None;first=True;previous=None
    try:
        journal=Journal(args.journal,args.verify_relations)
        while True:
            try:
                added=journal.scan();retry=None
            except FrameChecksumError as error:
                if not args.follow or str(error)==retry:raise
                retry=str(error);time.sleep(args.interval);continue
            if args.strict and journal.tail:raise InvalidJournal(journal.tail)
            if first or added:
                info=journal.summary(args.window,args.eta_model,args.eta_horizon)
                if args.json:print(json.dumps(info,sort_keys=True),flush=True)
                else:
                    if args.eta:previous=display_eta(info,previous,args.window)
                    else:display(info,False,args.window)
                    if info['job_end'] is not None:
                        print('Collection completed.' if info['job_end']['complete'] else 'Job ended without a complete factorization.')
                    sys.stdout.flush()
                first=False
            if not args.follow or journal.job_end is not None:break
            time.sleep(args.interval)
        return 0
    except KeyboardInterrupt:
        return 0
    except (InvalidJournal,OSError,UnicodeError,ValueError,struct.error) as error:
        print('journal-info: '+str(error),file=sys.stderr)
        return 1
    finally:
        if journal is not None:journal.close()


if __name__=='__main__':
    raise SystemExit(main())
