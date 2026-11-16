#!/usr/bin/env python3
"""Regenerate per-case pfw_input / postProcess / runTest wrappers from CASES.

Run from any directory. Edits to a specific generated wrapper will be replaced;
put shared logic in thermal_pfw.py and thermal_post.py instead.
"""
from pathlib import Path
from thermal_cases import CASES

ROOT=Path(__file__).resolve().parent
for name,case in CASES.items():
    d=ROOT/name
    d.mkdir(parents=True,exist_ok=True)
    (d/('pfw_input_'+name+'.py')).write_text('''# -*- coding: utf-8 -*-
"""{desc}.

Units: meters, kg, seconds, Kelvin, watts. See ../README.md.
Requires thermal PFW integration from geos_mpm_thermal_extension.patch.
"""
#[pfw_dependency] input:../thermal_pfw.py
#[pfw_dependency] input:../thermal_cases.py
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from thermal_pfw import build
pfw = build({name!r})
'''.format(desc=case['description'],name=name))
    post=d/('postProcess_'+name+'.py')
    post.write_text('''#!/usr/bin/env python3
"""Post-process {name}: analytical/conservation comparison and plots."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from thermal_post import run_case
if __name__=='__main__':raise SystemExit(run_case({name!r}))
'''.format(name=name))
    run=d/'runTest'
    run.write_text('''#!/usr/bin/env python3
"""GEOS MPM verification runner for {name}."""
from pathlib import Path
import sys
SOURCE_DIR = Path(__file__).resolve().parent
PFW_ROOT = SOURCE_DIR.parent.parent.parent
sys.path.insert(0,str(PFW_ROOT))
from mpm_vv_folder_runner import run_folder
VARIANTS=[{{"name":"default","label":{label!r},"case_name":"thermalVV_{name}","env":{{}}}}]
if __name__=='__main__':
    raise SystemExit(run_folder(source_dir=SOURCE_DIR,
        input_file='pfw_input_{name}.py',post_script='postProcess_{name}.py',
        variants=VARIANTS,default_case_id='thermal_{name}',
        description={label!r}))
'''.format(name=name,label=case['description']))
    testtex=d/'test.tex'
    escaped=name.replace('_',r'\_')
    testtex.write_text(
        r'\subsection{Thermal MPM: '+escaped+'}'+'\n'+
        case['description']+'. The reference solution and acceptance metric are in '
        +r'\texttt{thermal\_post.py}. See the generated '
        +r'\texttt{'+escaped+r'\_results.tex} for the status.'+'\n')
    run.chmod(0o755);post.chmod(0o755)
print('Generated',len(CASES),'case folders.')
