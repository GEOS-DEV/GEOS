"""Unified geometry example: box-like bulk particles and cut surface tets.

Use the normal PFW userDefs configuration, then run particleFileWriter.py with
this input. It writes files and a launch script; mSubmitJobs=False prevents job
submission. All outputs use canonical geometry schema 2.
"""
import pfw_geometryObjects as geom

pfw = dict(
    mSubmitJobs=False, mBatch=True, mBank='test', mWallTime='00:05:00',
    planeStrain=0, periodic=[False, False, False],
    xmin=0.0, xmax=1.0, ymin=0.0, ymax=1.0, zmin=0.0, zmax=1.0,
    nI=10, nJ=10, nK=10, ppc=1, xpar=1, ypar=1, zpar=1,
    endTime=0.001, plotInterval=0.001, restartInterval=1.0, outputType='silo',
    materials=['solid'],
    materialPropertyString='<ElasticIsotropic name="solid" defaultDensity="1000" '
                           'defaultYoungModulus="1e6" defaultPoissonRatio="0.25"/>',
    particleFileFields=['Velocity', 'MaterialType', 'ContactGroup', 'SurfaceFlag',
                        'SurfaceNormal', 'SurfacePosition'],
    # Retain cut geometry while using point B-spline interpolation.
    particleGeometry='fromRefinement', geometryReset='FOnly',
    surfaceRefinementMethod='simplex', surfaceRefinementLevel=2,
    subdivideParticles=0, oversizedParticleTreatment='none',
)

# The interpolation algorithm remains SinglePointBSpline, including cut tets.
# Runtime tet/wedge splitting is available but disabled in this example; CPTI
# remains future work. See ../geometrySplitting/pfw_input_geometry_splitting.py.
pfw['objects'] = [geom.sphere('ball', [0.5, 0.5, 0.5], 0.3,
                             mat=0, group=0, particleType=1)]
