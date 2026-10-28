#!/usr/bin/env python3
"""Fast deterministic ETA checks; no factor jobs or journal writes."""
import importlib.util
from pathlib import Path
import sys

sys.dont_write_bytecode=True
tool=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parents[4]/'siqs/tools/journal-info.py'
spec=importlib.util.spec_from_file_location('journal_info',tool)
info=importlib.util.module_from_spec(spec);spec.loader.exec_module(info)


def samples(function, stop=300, step=10):
    return [(t*1000000,2*int(function(t)),int(function(t))) for t in range(0,stop+1,step)]


def close(actual, expected):
    assert abs(actual-expected)<1e-7*max(1,abs(expected)),(actual,expected)


linear=samples(lambda t:3*t)
result=info.collection_eta(linear,60)
assert result['model']=='linear';close(result['seconds'],20)
print('PASS constant velocity retains the measured linear ETA')

accelerating=samples(lambda t:t+t*t/20)
result=info.collection_eta(accelerating,2040)
assert result['model']=='quadratic';close(result['seconds'],60)
close(result['current_rate'],31);close(result['quadratic_coefficient'],.05)
close(result['linear_seconds'],127.5)
print('PASS local quadratic recovers positive acceleration and its future root')

covered=31*300+.05*300*300
remaining=covered+(31+2*.05*300)*600
result=info.collection_eta(accelerating,remaining)
close(result['seconds'],900);close(result['acceleration_horizon_seconds'],300)
short=info.collection_eta(accelerating,remaining,horizon=60)
assert short['seconds']>result['seconds'];close(short['acceleration_horizon_seconds'],60)
print('PASS acceleration ends at the horizon and shorter horizons are conservative')

slowing=samples(lambda t:60*t-t*t/20)
result=info.collection_eta(slowing,450)
assert result['model']=='linear';close(result['seconds'],10)
print('PASS negative curvature falls back to the measured linear rate')

assert info.collection_eta(accelerating[:5],100)['model']=='linear'
assert info.collection_eta(samples(lambda t:0),100)['seconds'] is None
assert info.collection_eta([],0)['seconds']==0
assert info.collection_eta([(0,1,1),(0,2,2)],100)['seconds'] is None
assert info.collection_eta(accelerating,2040,model='linear')['model']=='linear'
print('PASS short, flat, duplicate-time and explicit-linear cases are handled')

offset=1<<60
translated=[(t+offset,raw+offset,full+offset) for t,raw,full in accelerating]
close(info.collection_eta(translated,2040)['seconds'],60)
print('PASS integer centering preserves accuracy with large clocks and counters')

bursty=sorted(accelerating+[accelerating[10]]*1000)
close(info.collection_eta(bursty,2040)['seconds'],60)
print('PASS bursts of frames preserve the quadratic after time binning')

noisy=accelerating[:]
t,raw,full=noisy[-1];noisy[-1]=(t,raw+100,full+50)
result=info.collection_eta(noisy,500)
assert result['model']=='quadratic'
h=result['seconds']
close(result['current_rate']*h+result['quadratic_coefficient']*h*h,500)
print('PASS forecasts anchor to the actual current count rather than the fitted intercept')

short_window=info.collection_eta(accelerating,2040,window=120)
assert short_window['model']=='quadratic'
close(short_window['seconds'],60);close(short_window['fit_span_seconds'],120)
print('PASS the selected window bounds both fit history and acceleration horizon')

def context(seconds=1482, span=120, full=200, single=100, large_primes=2):
    return dict(rate=dict(collection_eta_seconds=seconds,window_seconds=span),
                parameters=dict(max_large_primes=large_primes),
                snapshot=dict(full_relations=full,length_1=single))


assert info.collection_eta_text(None)=='waiting for more data'
assert info.collection_eta_text(dict(rate=None))=='waiting for more data'
assert info.collection_eta_text(context(seconds=None))=='waiting for more data'
print('PASS unavailable estimates still wait for more data')

assert info.collection_eta_text(context(span=59))=='roughly 25 minutes'
assert info.collection_eta_text(context(span=60))=='24.7 minutes'
print('PASS short history marks the display rough without relying on model choice')

assert info.collection_eta_text(context(single=101))=='roughly 25 minutes'
assert info.collection_eta_text(context(single=100))=='24.7 minutes'
assert info.collection_eta_text(context(single=99))=='24.7 minutes'
print('PASS half the completed relations must come from combined partials')

assert info.collection_eta_text(context(single=200,large_primes=0))=='24.7 minutes'
assert info.collection_eta_text(context(span=20,single=200,large_primes=0))=='roughly 25 minutes'
print('PASS smooth-only operation skips the graph cue but still needs history')

assert info.collection_eta_text(context(span=30),window=30)=='24.7 minutes'
assert info.collection_eta_text(context(span=29),window=30)=='roughly 25 minutes'
print('PASS a shorter selected window also shortens the display history requirement')

assert info.collection_eta_text(context(seconds=24.7,span=10))=='roughly 25 seconds'
assert info.collection_eta_text(context(seconds=.2,span=10))=='roughly 1 second'
assert info.collection_eta_text(context(seconds=60,span=10))=='roughly 1 minute'
assert info.collection_eta_text(context(seconds=0,span=10))=='0.0 seconds'
print('PASS rough seconds/minutes use whole units and zero remains exact')

import copy
ctx=context(single=150)
before=copy.deepcopy(ctx)
assert info.collection_eta_text(ctx)=='roughly 25 minutes'
assert ctx==before
print('PASS display formatting leaves numerical estimates and progress unchanged')

print('16/16 ETA checks passed.')
