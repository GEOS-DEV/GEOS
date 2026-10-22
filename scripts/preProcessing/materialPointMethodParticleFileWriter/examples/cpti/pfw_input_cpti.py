"""All-CPTI sphere: six tets per Cartesian piece; native cut tets at the surface.

Run with the normal PFW userDefs configuration. This writes inputs/launch files
but does not submit a job. CPTI requires 3D; ppc counts the seed lattice.
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
    particleGeometry='fromRefinement', geometryReset='FOnly',
    surfaceRefinementMethod='simplex', surfaceRefinementLevel=2,
    cpdiDomainScaling=1, cpdiDomainScalingType='homel',
    # Initial discretization only; runtime splitting can be enabled separately.
    subdivideParticles=0, oversizedParticleTreatment='none',
)
pfw['objects'] = [geom.sphere('ball', [0.5, 0.5, 0.5], 0.3,
                             mat=0, group=0, particleType=3)]
