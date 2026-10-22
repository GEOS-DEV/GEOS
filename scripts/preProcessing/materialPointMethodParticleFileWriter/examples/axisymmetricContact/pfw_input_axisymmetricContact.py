"""Two elastic cylinders approach across a gap: ordinary axisymmetric contact.

The lower and upper cylinders have different contact groups. Both use CPDI
cores and CPTI outer layers. There are no CZ, DFG, pressure loads or moving
boundaries. The material's zero Poisson ratio isolates axial impact initially.
"""
#[pfw_dependency] input:caseSettings.py
import pfw_geometryObjects as geom
import pfw_tracerPoints as tracers
from caseSettings import settings

case = settings()

# ---------------------------------------------------------------------------
# Native r-y grid. nI/nJ include two buffer cells; spacing is 1 mm, with two
# particles per cell direction. Only the radial axis is a symmetry boundary.
# ---------------------------------------------------------------------------
pfw = dict(
    mSubmitJobs=True, mBatch=True, mWallTime='00:10:00',
    formulation='Axisymmetric', xmin=0., xmax=0.012, ymin=-0.014, ymax=0.014,
    nI=14, nJ=30, ppc=2, xpar=1, ypar=1, zpar=1,
    periodic=[False, False], boundaryConditionTypes=[1, 0, 0, 0],
    timeIntegrationOption='ExplicitDynamic', updateMethod=case['update'], updateOrder=2,
    cflFactor=0.25, initialDt=1.e-9, endTime=case['end_time'],
    plotInterval=case['end_time']/60., restartInterval=case['end_time']/2.,
    outputType='silo', particleGeometry='fromRefinement', geometryReset='FOnly',
    cpdiDomainScaling=1, cpdiDomainScalingType='homel',
    particleFileFields=['Velocity', 'MaterialType', 'ContactGroup', 'ParticleColor',
                        'SurfaceFlag', 'SurfaceNormal', 'SurfacePosition'],
    boxAverageHistory=1, boxAverageWriteInterval=case['end_time']/200.,
    plotGridFields=1, gridFieldNames=['gridMass', 'gridVelocity', 'gridSurfaceNormal',
                                    'gridContactForce'],
)

pfw['runDebug'] = True

# ---------------------------------------------------------------------------
# Ordinary body contact. Gap gating uses the explicitly supplied interface
# positions. Normal directions use the corrected volume-gradient estimator;
# no explicit-normal bias or radius-dependent friction parameter is needed.
# ---------------------------------------------------------------------------
pfw.update(enableContact=1, damageFieldPartitioning=0, frictionCoefficient=0.,
           contactNormalType='MassWeighted', contactGapCorrection='Implicit',
           useSurfacePositionForContact=1, explicitSurfaceNormalInfluence=0.,
           computeParticleSurfaceNormalsAndPositions=0, overlapCorrection='Off',
           disableSurfaceNormalsAndPositionsOnCPDIScaling=0)
pfw['materials'] = ['solid']
pfw['materialPropertyString'] = (
    '<ElasticIsotropic name="solid" defaultDensity="{density}" '
    'defaultYoungModulus="{E}" defaultPoissonRatio="{nu}"/>'.format(**case))

# ---------------------------------------------------------------------------
# ContactGroup identifies a body; ParticleColor is only for visualization.
# Flag only the opposing axial faces. The CPDI/CPTI seam within each cylinder
# shares its contact group and therefore does not become a contact interface.
# ---------------------------------------------------------------------------
class MovingPiece(geom.box):
    def getVelocity(self, point):
        return [0., case['speed'] if self.group == 0 else -case['speed'], 0.]


pfw['objects'] = []
for group, lo, hi, flags in [
        (0, -0.009, -0.001, [False, False, False, True]),
        (1,  0.001,  0.009, [False, True, False, False])]:
    for part, rlo, rhi, particle_type in [('core', 0., 0.008, 2), ('skin', 0.008, 0.010, 3)]:
        obj = MovingPiece(f'body{group}_{part}', [rlo, lo], [rhi, hi], dim=2,
                          mat=0, group=group, particleType=particle_type, flaggedSurfaces=flags)
        obj.particleColor = group
        pfw['objects'].append(obj)

# ---------------------------------------------------------------------------
# Lagrangian histories resolve arrival of compression and subsequent rebound.
# The initial geometric gap closes at 1 ms in free flight. Tracer positions
# are particle centres, so their separation is not the physical surface gap.
# ---------------------------------------------------------------------------
tracers.set_tracers(
    pfw, points=[[0.00425, -0.00125, 0.], [0.00425, 0.00125, 0.],
                 [0.00425, -0.00875, 0.], [0.00425, 0.00875, 0.]],
    variables=['particleID', 'velocityX', 'velocityY', 'velocityZ', 'stressYY', 'stressZZ'],
    write_interval=case['end_time']/200., output_prefix='contactTracer')
