#!/usr/bin/env python3
"""Reproducible search/state checks. No GEOS build or mechanics is implied."""
from pathlib import Path
import argparse,json,os,platform,shutil,subprocess,sys
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--cxx',default='g++')
p.add_argument('--sanitize',action='store_true')
p.add_argument('--benchmark',action='store_true')
p.add_argument('--mpi',action='store_true')
a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
here=Path(__file__).resolve().parent
flags=['-std=c++14','-O3','-Wall','-Wextra','-Werror']
if a.sanitize: flags=['-std=c++14','-O1','-g','-fsanitize=address,undefined','-fno-omit-frame-pointer','-fno-pie','-no-pie']
report={'scope':'standalone/extracted-block checks, not full GEOS', 'platform':platform.platform(),'sanitizers':a.sanitize,'tests':{},'commands':[]}
def run(cmd):
 report['commands'].append(' '.join(map(str,cmd)))
 r=subprocess.run(list(map(str,cmd)),check=True,text=True,capture_output=True,timeout=120)
 return r.stdout+r.stderr
report['compiler']=run([a.cxx,'--version']).splitlines()[0]
for name in ['testIndexedSearch','testStateValidation']:
 run([a.cxx,*flags,here/(name+'.cpp'),'-o',out/name])
 text=run([out/name]);(out/(name+'.log')).write_text(text);report['tests'][name]=text.strip()
cmd=[sys.executable,here/'testProductionSearchAdapter.py','--output',out/'adapter','--cxx',a.cxx]
if a.sanitize:cmd+=['--sanitize']
text=run(cmd);(out/'adapter.log').write_text(text);report['tests']['productionSearchFitBlock']=text.strip()
text=run([sys.executable,here/'testOctantSearchOption.py'])
(out/'octant_search_option.log').write_text(text);report['tests']['octantSearchMode']=text.strip()
if a.benchmark:
 report['benchmarks']=[]
 for repeat in range(4):
  text=run([out/'testIndexedSearch','--benchmark']);value=json.loads(text);value['repetition']=repeat;value['warmup']=repeat==0;report['benchmarks'].append(value)
if a.mpi:
 cxx,launcher=shutil.which('mpicxx'),shutil.which('mpiexec')
 if cxx and launcher:
  run([cxx,*flags,here/'testMPIState.cpp','-o',out/'testMPIState'])
  report['mpi']=[{'ranks':n,'log':run([launcher,'-n',n,out/'testMPIState'])}for n in [1,2,4,8]]
 else:report['mpi']='NOT RUN: mpicxx/mpiexec unavailable; emulated ownership is not MPI coverage'
else:report['mpi']='NOT REQUESTED'
(out/'results.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items()if k!='commands'},indent=2))
