"""Axisymmetric rubber-jacketed specimen with steel end plates.

x=r, y=axial position, z=0. The full height is modeled. All interfaces share
contact group 0, so they remain bonded. A virtual pressure bath loads the closed
outer boundary, including the jacket caps. No contact/CZ, Weibull, fracture,
explicit fluid, pore-pressure evolution or prescribed platen motion is used.

Requires the axisymmetric stage-1 and uniform-confinement example patches.
"""
#[pfw_dependency] input:caseSettings.py
import pfw_geometryObjects as geom
import pfw_tracerPoints as tracers
from caseSettings import settings

case = settings()
r = case['radius']
a = case['sample_half_height']
tp = case['plate_thickness']
tj = case['jacket_thickness']
refine = case['refine']

# ---------------------------------------------------------------------------
# Mesh and solver: genuine r-y mesh, full-ring particle volumes, no thickness.
# nI/nJ include the writer's two outer buffer cells. The jacket is two particles
# thick at the default resolution; JACKET_REFINE=2 doubles cells per direction.
# PIC provides numerical damping for this pressure-ramp/settling smoke test.
# ---------------------------------------------------------------------------
pfw = dict(
    mSubmitJobs=True, mBatch=True, mWallTime='00:20:00',
    formulation='Axisymmetric', xmin=0., xmax=0.014, ymin=-0.015, ymax=0.015,
    nI=14*refine+2, nJ=30*refine+2, ppc=2, xpar=1, ypar=1, zpar=1,
    periodic=[False, False], boundaryConditionTypes=[1, 0, 0, 0],
    enableContact=0, damageFieldPartitioning=0,
    timeIntegrationOption='ExplicitDynamic', updateMethod='PIC',
    cflFactor=0.25, initialDt=1.e-10,
    endTime=case['end_time'], plotInterval=case['end_time']/40.,
    restartInterval=case['end_time']/2., outputType='silo',
    particleGeometry='fromRefinement', geometryReset='FOnly',
    cpdiDomainScaling=1, cpdiDomainScalingType='homel',
    particleFileFields=['Velocity', 'MaterialType', 'ContactGroup', 'ParticleColor',
                        'SurfaceFlag', 'SurfaceNormal', 'SurfacePosition'],
    boxAverageHistory=1, boxAverageWriteInterval=case['end_time']/200.,
    boxAverageMin='{ 0, -0.009, -0.001 }', boxAverageMax='{ 0.009, 0.009, 0.001 }',
)

pfw['runDebug'] = True

# ---------------------------------------------------------------------------
# Three elastic materials. Rubber uses a nearly incompressible *elastic*
# surrogate at small strains. Steel is deformable with its actual stiffness
# scale; it therefore controls the explicit timestep. No mass scaling is used.
# The homogeneous variant replaces all three parameter sets by sample values.
# ---------------------------------------------------------------------------
pfw['materials'] = [m['name'] for m in case['materials']]
pfw['materialPropertyString'] = '\n'.join(
    '<ElasticIsotropic name="{name}" defaultDensity="{density}" '
    'defaultYoungModulus="{E}" defaultPoissonRatio="{nu}"/>'.format(**m)
    for m in case['materials'])

# ---------------------------------------------------------------------------
# Nonoverlapping meridional rectangles, revolved implicitly about r=0.
# CPDI quads fill sample and plates; CPTI triangles fill the thin jacket/caps.
# Surface flags mark only the outside of the jacket. Internal seams and the
# symmetry axis are not physical free surfaces. Box flag order: r-, y-, r+, y+.
# ---------------------------------------------------------------------------
def piece(name, lo, hi, material, particle_type, flags):
    obj = geom.box(name, lo, hi, dim=2, mat=material, group=0,
                   particleType=particle_type, flaggedSurfaces=flags)
    obj.particleColor = material
    return obj

pfw['objects'] = [
    piece('sample', [0., -a], [r, a], 0, 2, [False]*4),
    piece('lower_plate', [0., -a-tp], [r, -a], 1, 2, [False]*4),
    piece('upper_plate', [0., a], [r, a+tp], 1, 2, [False]*4),
    piece('jacket_side', [r, -a-tp], [r+tj, a+tp], 2, 3, [False, False, True, False]),
    piece('jacket_lower_cap', [0., -a-tp-tj], [r+tj, -a-tp], 2, 3, [False, True, True, False]),
    piece('jacket_upper_cap', [0., a+tp], [r+tj, a+tp+tj], 2, 3, [False, False, True, True]),
]

# ---------------------------------------------------------------------------
# Fluid pressure without fluid particles: a uniform background stress -p I.
# The cylindrical weak form includes the hoop pressure contribution. A cosine
# ramp avoids an initial step load. Pressure is held after the event finishes.
# Do not supply box/sphere exclusions or a second explicit surface traction.
# ---------------------------------------------------------------------------
pfw['useEvents'] = 1
pfw['mpmEventsString'] = f'''
<ConfiningPressure name="pressureRamp" startTime="0" endTime="{case['ramp_time']}"
                   startPressure="0" endPressure="{case['pressure']}"
                   interpType="1" uniformPressure="1"/>
'''

# ---------------------------------------------------------------------------
# Histories: sample-interior box averages and five Lagrangian tracers.
# stressXX=radial, stressYY=axial, stressZZ=hoop. All volumes are full-ring values.
# ---------------------------------------------------------------------------
tracers.set_tracers(
    pfw,
    points=[[0.00525, 0.00025, 0.], [0.00525, 0.00975, 0.],
            [0.00525, 0.01125, 0.], [0.01025, 0.00025, 0.], [0.00525, 0.01275, 0.]],
    variables=['particleID', 'materialType', 'volume', 'density',
               'velocityX', 'velocityY', 'stressXX', 'stressYY', 'stressZZ', 'stressXY'],
    write_interval=case['end_time']/200., output_prefix='jacketTracer')
