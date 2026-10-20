"""Small PFW-only fixture; never submits or runs GEOS.

PFW_TEST_DIM=2|3, PFW_TEST_MODE=cartesian|simplex,
PFW_TEST_GEOMETRY=plane|separated, PFW_TEST_LEVEL=1|2|...
The make_objects factory deliberately filters by the rank interval.
"""
import os
import numpy as np
import pfw_geometryObjects as geom

DIM = int(os.environ.get('PFW_TEST_DIM', '3'))
MODE = os.environ.get('PFW_TEST_MODE', 'simplex')
SCENE = os.environ.get('PFW_TEST_GEOMETRY', 'plane')
LEVEL = int(os.environ.get('PFW_TEST_LEVEL', '2'))

pfw = dict(
    mSubmitJobs=False, mBatch=True, mBank='test', mWallTime='00:01:00',
    planeStrain=int(DIM == 2), periodic=[False, False, False],
    xmin=0.0, xmax=1.0, ymin=0.0, ymax=1.0, zmin=0.0, zmax=1.0,
    nI=7, nJ=6, nK=3 if DIM == 2 else 5, ppc=1,
    xpar=1, ypar=1, zpar=1, outputType='silo',
    endTime=0.01, plotInterval=0.01, restartInterval=1.0,
    surfaceRefinementMethod=MODE, surfaceRefinementLevel=LEVEL,
    particleReferenceGeometry=os.environ.get('PFW_TEST_REFERENCE_GEOMETRY', 'equivalentBox'),
    sortObjects=True,
    materials=['mat0', 'mat1'],
    particleFileFields=['Velocity', 'MaterialType', 'ContactGroup',
                        'SurfaceFlag', 'RVector', 'SurfaceNormal', 'SurfacePosition'],
    materialPropertyString='\n'.join(
        '<ElasticIsotropic name="mat%d" defaultDensity="1" '
        'defaultYoungModulus="%d" defaultPoissonRatio="0.25"/>' % (i, i + 1)
        for i in (0, 1)),
)


class PlanarSide(geom.Geometry):
    def __init__(self, sign, material):
        super().__init__('side%d' % material, mat=material, group=material, particleType=1)
        self.sign = sign
        self.normal = np.array([1.0, 0.31, 0.19 if DIM == 3 else 0.0])
        self.normal /= np.linalg.norm(self.normal)

    def isInterior(self, point, depth):
        if np.any(np.asarray(point) < 0) or np.any(np.asarray(point) >= 1):
            return -1
        distance = self.sign * (np.dot(self.normal, point) - 0.49)
        return -1 if distance > 0 else (2 if abs(distance) < depth else 0)

    def getSurfaceNormal(self, point):
        return self.sign * self.normal

    def getSurfacePosition(self, point):
        return -(np.dot(self.normal, point) - 0.49) * self.normal

    def xMin(self):
        return 0.0

    def xMax(self):
        return 1.0


def make_objects(rankxmin, rankxmax):
    if SCENE == 'plane':
        return [PlanarSide(1, 0), PlanarSide(-1, 1)]
    descriptions = [([0.0, 0.08, 0.08], [0.32, 0.93, 0.93], 0),
                    ([0.58, 0.08, 0.08], [1.0, 0.93, 0.93], 1)]
    result = []
    for lower, upper, material in descriptions:
        if lower[0] <= rankxmax and upper[0] >= rankxmin:
            result.append(geom.box('box%d' % material, lower[:DIM], upper[:DIM],
                                   mat=material, group=material, particleType=1, dim=DIM))
    return result
