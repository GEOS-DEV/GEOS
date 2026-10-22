"""Check generated cylinder measures; not a substitute for a dynamics benchmark."""
import sys
import numpy as np
path=sys.argv[1]
with open(path) as f: names=f.readline().split()
a=np.loadtxt(path,skiprows=1,ndmin=2)
radius=a[:,names.index('PositionX')]
volume=a[:,names.index('ParticleVolume')]
shape=a[:,names.index('GeometryShape')].astype(int)
g=a[:,[names.index('Geometry'+ij) for ij in ('XX','XY','YX','YY')]].reshape(-1,2,2)
expected=2*np.pi*radius*np.linalg.det(g)*np.where(shape==4,0.5,1.0)
assert np.allclose(volume,expected,rtol=1e-12,atol=0)
assert set(shape)=={3,4}
assert np.isclose(volume.sum(),np.pi,rtol=1e-12,atol=0)
print('Full-ring cylinder volume:',volume.sum(),'mass at rho=1000:',1000*volume.sum())
