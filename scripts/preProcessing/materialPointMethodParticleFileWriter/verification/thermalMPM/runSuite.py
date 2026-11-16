#!/usr/bin/env python3
"""Dispatch individual GEOS verification case folders via their runTest."""
import argparse
import subprocess
import sys
from pathlib import Path
from thermal_cases import CASES

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--case',action='append',choices=sorted(CASES))
    parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--force',action='store_true')
    args=parser.parse_args()
    root=Path(__file__).resolve().parent
    names=args.case or list(CASES)
    failures=[]
    for name in names:
        cmd=[sys.executable,str(root/name/'runTest')]
        if args.prepare_only:cmd.append('--prepare-only')
        if args.force:cmd.append('--force')
        print('Running:', ' '.join(cmd),flush=True)
        rc=subprocess.call(cmd,cwd=root/name)
        if rc!=0:failures.append((name,rc))
    if failures:print('Non-passing cases:',failures)
    return 0 if not failures else 1

if __name__=='__main__':raise SystemExit(main())
