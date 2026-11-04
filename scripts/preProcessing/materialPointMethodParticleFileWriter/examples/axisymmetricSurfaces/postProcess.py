"""Check the initial full-ring particle volume and report GEOS area diagnostics.

Run in the output folder. Choose the same AXIS_SURFACE_CASE as generation.
No dynamics acceptance is inferred from successful input generation.
"""
from pathlib import Path
import os
import re
import numpy as np
case=os.environ.get('AXIS_SURFACE_CASE','cylinder')
p=Path('mpmParticleFile_case');names=p.read_text().splitlines()[0].split()
data=np.loadtxt(p,skiprows=1,ndmin=2)
volume=data[:,names.index('ParticleVolume')].sum()
expected=np.pi if case=='cylinder' else 4*np.pi/3
area=2*np.pi if case=='cylinder' else 4*np.pi
print(f'Initial ring volume {volume:.9g}; analytical {expected:.9g}')
assert np.isclose(volume,expected,rtol=.02 if case=='sphere' else 1e-10)
print(f'Initial surface area reference {area:.9g}; Laplace pressure {(.01 if case=="cylinder" else .02):.9g} when tension is enabled')
log=Path('geos.log')
if log.exists():
    values=[float(v) for v in re.findall(r'Total surface area for field 0\s*=\s*([-+\d.eE]+)',log.read_text(errors='replace'))]
    if values:
        print(f'GEOS area: first={values[0]:.9g}, last={values[-1]:.9g}; initial relative error={(values[0]/area-1):.3%}')
        assert np.all(np.isfinite(values)) and min(values)>0
    else: raise SystemExit('No nodal-area output found; verify computeNodalArea=1 and a completed timestep')
