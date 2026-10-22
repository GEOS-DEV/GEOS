# Hex bulk with tetrahedral surface refinement, using CPDI and CPTI.
# Writes particle/XML/launch files; does not submit a job.
# See HEX_BULK_TET_SURFACE_PLAN.md for configuration and diagnostics.
import pfw_geometryObjects as geom

pfw = dict(
    # Normal PFW execution controls; supply your usual userDefs configuration.
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

    # hexTet selects native geometry and CPTI for actual cuts. Uncut/fallback
    # pieces keep the object's CPDI type. Level is a factor per direction;
    # level 1 still cuts surfaces, while level 2 first divides candidate parents.
    surfaceRefinementMethod='hexTet', surfaceRefinementLevel=2,
    surfaceRefinementFallback='cartesian', surfaceRefinementPlaneTolerance=0.05,

    # Physical geometry evolves independently of a constitutive F-only reset.
    # Mapping may temporarily scale CPDI/CPTI domains using their own algorithms.
    geometryReset='FOnly', cpdiDomainScaling=1, cpdiDomainScalingType='homel',

    # This example isolates initial refinement. Enable runtime splitting
    # separately if wanted; both shapes already have native split support.
    subdivideParticles=0, oversizedParticleTreatment='none',
)

# Type 2 gives CPDI hex bulk. Do not set this object to type 3: that would
# request an all-tetrahedron object, which the hexTet preset rejects.
pfw['objects'] = [geom.sphere('ball', [0.5, 0.5, 0.5], 0.3,
                             mat=0, group=0, particleType=2)]
# Inspect mpmParticleFile_<case>_refinement.json / particle_layout after writing.
# Its shape/type counts and volumes describe final rows, including fallback.
