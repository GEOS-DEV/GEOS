"""Check contact inputs and full writer output with 1/3 emulated ranks.

Usage: python tests/run_axisymmetric_contact_writer_tests.py --output NEW_DIR
Emulated collectives check deterministic file emission, not GEOS MPI dynamics.
"""
import argparse
import getpass
import os
from pathlib import Path
import subprocess
import sys
import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--ranks', type=int, nargs='+', default=[1, 3])
args = parser.parse_args()
base = Path(__file__).resolve().parents[1]
case = base/'examples/axisymmetricContact'
sys.path.insert(0, str(base))
import pfw_axisymmetric as axis

# Input validation admits ordinary contact but preserves unsupported-feature
# errors. These direct calls also exercise fractional/negative group inputs.
data = dict(particleType=2, group=0, surfaceFlag=2, czTag=0, velocity=[0, 1, 0])
for group in (0, 1, 7, 99):
    axis.validate_particle([2, 0, 0], dict(data, group=group), (3, np.eye(2)))
for change in (dict(group=-1), dict(group=.5), dict(group=100),
               dict(czTag=1), dict(surfaceFlag=3), dict(surfaceFlag=4),
               dict(velocity=[0, 0, 1])):
    try:
        axis.validate_particle([2, 0, 0], dict(data, **change), (3, np.eye(2)))
    except ValueError:
        pass
    else:
        raise AssertionError('Unsupported particle accepted: '+repr(change))
for change in (dict(damageFieldPartitioning=1), dict(overlapCorrection='normalForce'),
               dict(computeParticleSurfaceNormalsAndPositions=1), dict(computeNodalArea=1)):
    try:
        axis.validate_options(dict(xmin=0, xmax=2, **change))
    except ValueError:
        pass
    else:
        raise AssertionError('Unsupported option accepted: '+repr(change))
axis.validate_options(dict(xmin=0, xmax=2, enableContact=1))
print('Contact-group and unsupported-feature validation passed.', flush=True)

for update in ('FLIP', 'FMPM'):
    reference = None
    for size in args.ranks:
        folder = args.output.resolve()/f'{update}_r{size}'
        folder.mkdir(parents=True, exist_ok=False)
        env = dict(os.environ, AXIS_CONTACT_UPDATE=update,
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
                    str(base/'particleFileWriter.py'), str(case/'pfw_input_axisymmetricContact.py'),
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
        content = (folder/'mpmParticleFile_axisymmetricContact').read_bytes()
        if reference is None:
            reference = content
        assert content == reference, 'Particle file depends on writer rank count'
        print(update, size, 'emulated ranks passed', flush=True)
