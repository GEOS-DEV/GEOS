#!/usr/bin/env python3
"""Compare two GEOS MPM thermal runs (MPI/CPU/GPU parity or convergence).

Uses per-tracer particle histories rather than rank-local nodal fields. Outputs
absolute/relative mismatch and fails if a tracer or thermal variable is absent.
"""
import argparse
import csv
import json
import math
from pathlib import Path

from thermal_post import read_tracers


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--candidate',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,default=Path('comparison'))
    p.add_argument('--rtol',type=float,default=1e-4)
    p.add_argument('--atol',type=float,default=0.005)
    args=p.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    try:
        a=read_tracers(args.reference); b=read_tracers(args.candidate)
        if len(a)!=len(b):raise ValueError('Tracer count mismatch {} != {}'.format(len(a),len(b)))
        diffs=[]
        for i,((na,ra),(nb,rb)) in enumerate(zip(a,b)):
            if not ra or not rb:raise ValueError('Empty tracer')
            # Compare last output time; no sample-time interpolation because
            # identical run configurations must reach the same end time.
            x,y=ra[-1],rb[-1]
            if abs(x['t']-y['t'])>1e-8*max(1.,abs(x['t'])):
                raise ValueError('Different end times in tracer '+str(i))
            if abs(x['x']-y['x'])>1e-6:
                raise ValueError('Different tracer coordinates in tracer '+str(i))
            d=y['temperature']-x['temperature']
            diffs.append(dict(tracer=i,x=x['x'],time=x['t'],reference=x['temperature'],
                              candidate=y['temperature'],difference=d,
                              tolerance=args.atol+args.rtol*abs(x['temperature'])))
        l2=math.sqrt(sum(d['difference']**2 for d in diffs)/len(diffs))
        maxerr=max(abs(d['difference']) for d in diffs)
        passed=all(abs(d['difference'])<=d['tolerance'] for d in diffs)
        summary=dict(status='PASS' if passed else 'FAIL',n=len(diffs),rms_error_K=l2,
                     max_error_K=maxerr,reference=str(args.reference),candidate=str(args.candidate),
                     atol=args.atol,rtol=args.rtol)
        with (args.output_dir/'comparison.csv').open('w',newline='') as h:
            w=csv.DictWriter(h,fieldnames=list(diffs[0]));w.writeheader();w.writerows(diffs)
    except Exception as e:
        summary=dict(status='FAIL',error=str(e))
    (args.output_dir/'comparison_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    return 0 if summary['status']=='PASS' else 2

if __name__=='__main__':raise SystemExit(main())
