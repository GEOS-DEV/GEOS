"""Emit and check four axisymmetric contact/DFG cases with FLIP and FMPM.

Uses emulated collectives to check deterministic PFW output. Does not run
GEOS or test distributed contact dynamics. Use a fresh --output directory.
"""
import argparse
import getpass
import os
from pathlib import Path
import subprocess
import sys

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--ranks', type=int, nargs='+', default=[1, 3])
args = parser.parse_args()
base = Path(__file__).resolve().parents[1]
case = base/'examples/axisymmetricDFG'

for mode, update in ((m, u) for m in ('intact', 'damage', 'contact', 'combined') for u in ('FLIP', 'FMPM')):
    reference = None
    for size in args.ranks:
        folder = args.output.resolve()/f'{mode}_{update}_r{size}'
        folder.mkdir(parents=True, exist_ok=False)
        env = dict(os.environ, AXIS_DFG_UPDATE=update, AXIS_DFG_CASE=mode,
                   OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        (folder/f'userDefs_{getpass.getuser()}.py').write_text(
            'geosPath="/bin/true"\npfwPath='+repr(str(base))+'\ndefaultBank="test"\n')
        subprocess.run([sys.executable, str(case/'caseSettings.py'), '--write', 'case.json'],
                       cwd=folder, env=env, check=True)
        logs = []
        jobs = []
        try:
            for rank in range(size):
                log = (folder/f'rank{rank}.log').open('w')
                logs.append(log)
                jobs.append(subprocess.Popen([
                    sys.executable, str(base/'tests/_collective_test_worker.py'),
                    str(base/'particleFileWriter.py'), str(case/'pfw_input_axisymmetricDFG.py'),
                    str(rank), str(size), str(folder/'collectives')],
                    cwd=folder, env=env, stdout=log, stderr=subprocess.STDOUT))
            for job in jobs:
                if job.wait(timeout=120) != 0:
                    raise RuntimeError('Writer failed; inspect '+str(folder))
        finally:
            for job in jobs:
                if job.poll() is None:
                    job.kill()
                    job.wait()
            for log in logs:
                log.close()
        subprocess.run([sys.executable, str(case/'checkParticles.py')],
                       cwd=folder, env=env, check=True)
        content = (folder/'mpmParticleFile_axisymmetricDFG').read_bytes()
        if reference is None:
            reference = content
        assert content == reference, 'Particle file depends on writer rank count'
        print(mode, update, size, 'emulated ranks passed', flush=True)
