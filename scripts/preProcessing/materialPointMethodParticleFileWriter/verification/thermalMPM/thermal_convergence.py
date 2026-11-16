#!/usr/bin/env python3
"""Require decreasing Fourier-mode error under 18/34/66 node refinement."""
import argparse
import json
import math
from pathlib import Path
from thermal_cases import CASES, temperature,AMP
from thermal_post import read_tracers

def err(folder,case_id):
    traces=read_tracers(folder)
    case=CASES[case_id]
    ds=[]
    for _, rows in traces:
        r=rows[-1]
        ds.append(r['temperature']-temperature('cosine',r['x'],r['t'],case=case))
    return math.sqrt(sum(x*x for x in ds)/len(ds))/AMP

def main():
    p=argparse.ArgumentParser()
    for c in ('cosine_coarse','cosine','cosine_fine'):
        p.add_argument('--'+c,required=True,type=Path)
    p.add_argument('--output',type=Path,default=Path('thermal_convergence.json'))
    a=p.parse_args()
    try:
        errors={k:err(getattr(a,k),k) for k in ('cosine_coarse','cosine','cosine_fine')}
        status='PASS' if errors['cosine_fine']<errors['cosine']<errors['cosine_coarse'] else 'FAIL'
        result=dict(status=status,errors=errors,criterion='strict monotonically decreasing normalized RMSE')
    except Exception as exc:result=dict(status='FAIL',error=str(exc))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    return 0 if result['status']=='PASS' else 2

if __name__=='__main__':raise SystemExit(main())
