"""Generate/check both jacketed-sample variants with emulated collectives.

Run with --output NEW_DIRECTORY; this is not a real MPI or GEOS dynamics test.
"""
import argparse
import getpass
import os
from pathlib import Path
import subprocess
import sys

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--ranks', nargs='+', type=int, default=[1, 3])
args = parser.parse_args()
base = Path(__file__).resolve().parents[1]
case = base/'examples/axisymmetricJacketedSample'
for variant in ('jacketed', 'homogeneous'):
    reference = None
    for size in args.ranks:
        folder = args.output.resolve()/f'{variant}_r{size}'
        folder.mkdir(parents=True, exist_ok=False)
        env = dict(os.environ, JACKET_CASE=variant, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        (folder/f'userDefs_{getpass.getuser()}.py').write_text(
            'geosPath="/bin/true"\npfwPath='+repr(str(base))+'\ndefaultBank="test"\n')
        subprocess.run([sys.executable, str(case/'caseSettings.py'), '--write', 'case.json'], cwd=folder, env=env, check=True)
        logs = []; jobs = []
        try:
            for rank in range(size):
                log = (folder/f'rank{rank}.log').open('w'); logs.append(log)
                jobs.append(subprocess.Popen([sys.executable, str(base/'tests/_collective_test_worker.py'),
                    str(base/'particleFileWriter.py'), str(case/'pfw_input_axisymmetricJacketedSample.py'),
                    str(rank), str(size), str(folder/'collectives')],
                    cwd=folder, env=env, stdout=log, stderr=subprocess.STDOUT))
            for job in jobs:
                if job.wait(timeout=120) != 0:
                    raise RuntimeError('Writer failed; inspect '+str(folder))
        finally:
            for job in jobs:
                if job.poll() is None: job.kill(); job.wait()
            for log in logs: log.close()
        subprocess.run([sys.executable, str(case/'checkParticles.py')], cwd=folder, env=env, check=True)
        content = (folder/'mpmParticleFile_axisymmetricJacketedSample').read_bytes()
        if reference is None: reference = content
        assert content == reference, 'Particle output depends on writer rank count'
        print(variant, size, 'emulated ranks passed', flush=True)
