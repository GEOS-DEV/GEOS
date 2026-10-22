"""Two touching cylinders initially pull apart against a prescribed CZ.

Mixed CPDI cores and CPTI skins share one contact group per cylinder. The
elastic variant remains tied; the debond variant uses a low opening limit.
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
    formulation='Axisymmetric', xmin=0., xmax=0.012, ymin=-0.012, ymax=0.012,
    nI=14, nJ=26, ppc=2, xpar=1, ypar=1, zpar=1,
    periodic=[False, False], boundaryConditionTypes=[1, 0, 0, 0],
    timeIntegrationOption='ExplicitDynamic', updateMethod=case['update'], updateOrder=2,
    cflFactor=0.25, initialDt=1.e-9, endTime=case['end_time'],
    plotInterval=case['end_time']/60., restartInterval=case['end_time']/2.,
    outputType='silo', particleGeometry='fromRefinement', geometryReset='FOnly',
    cpdiDomainScaling=1, cpdiDomainScalingType='homel',
    particleFileFields=['Velocity', 'MaterialType', 'ContactGroup', 'ParticleColor',
                        'SurfaceFlag', 'SurfaceNormal', 'SurfacePosition', 'CZTag'],
    boxAverageHistory=1, boxAverageWriteInterval=case['end_time']/200.,
    plotGridFields=1, gridFieldNames=['gridMass', 'gridVelocity', 'gridSurfaceNormal',
                                    'gridContactForce', 'gridCohesiveForce', 'gridCohesiveArea', 'gridCohesiveNode',
                                    'gridDisplacement'],
)

pfw['runDebug'] = True

# ---------------------------------------------------------------------------
# Contact supplies compression across an intact CZ and ordinary contact after
# complete debonding. Friction is zero; cohesive shear remains in the CZ law.
# ---------------------------------------------------------------------------
pfw.update(preventCZInterpenetration=1, enableContact=1, damageFieldPartitioning=0, frictionCoefficient=0.,
           contactNormalType='MassWeighted', contactGapCorrection='Implicit',
           useSurfacePositionForContact=1, explicitSurfaceNormalInfluence=0.,
           computeParticleSurfaceNormalsAndPositions=0, overlapCorrection='Off',
           disableSurfaceNormalsAndPositionsOnCPDIScaling=0)
pfw['materials'] = ['solid']
pfw['materialPropertyString'] = (
    '<ElasticIsotropic name="solid" defaultDensity="{density}" '
    'defaultYoungModulus="{E}" defaultPoissonRatio="{nu}"/>'.format(**case))

# ---------------------------------------------------------------------------
# Traction-separation stiffness has units Pa/m. Reference area is the full
# surface of revolution. No plane-strain thickness or extra 2*pi*r is supplied.
# CZTag selects the interface, not a material or a new contact group.
# ---------------------------------------------------------------------------
pfw['materialPropertyString'] += (
    '<UncoupledCohesiveZone name="interfaceLaw" normalForceConstant="{stiffness}" '
    'shearForceConstant="{stiffness}" maxNormalDisplacement="{failure_opening}" '
    'maxTangentialDisplacement="1.0"/>'.format(**case))
pfw['cohesiveZoneRegions'] = """
<CohesiveZoneRegion name="bond" constitutiveModel="interfaceLaw" tag="1"
                    czVolumeNormalization="1"/>"""
# Defining the activation event does not enable the solver's event dispatcher.
# SurfaceFlag/CZTag identify candidates; the event creates their active bond.
pfw['useEvents'] = 1
pfw['mpmEventsString'] = (
    '<ReferenceCohesiveZones name="activateBond" startTime="0" regionNames="{bond}" '
    'computeNormalsAndPositions="0" czSurfaceDisplacementUpdate="' + case['displacement'] + '"/>')

# ---------------------------------------------------------------------------
# ContactGroup identifies a body; ParticleColor is only for visualization.
# Flag only the opposing axial faces. The CPDI/CPTI seam within each cylinder
# shares its contact group and therefore does not become a contact interface.
# ---------------------------------------------------------------------------
class MovingPiece(geom.box):
    def getVelocity(self, point):
        return [0., -case['speed'] if self.group == 0 else case['speed'], 0.]


pfw['objects'] = []
for group, lo, hi, flags in [
        (0, -0.008, 0.0, [False, False, False, True]),
        (1,  0.0,  0.008, [False, True, False, False])]:
    for part, rlo, rhi, particle_type in [('core', 0., 0.008, 2), ('skin', 0.008, 0.010, 3)]:
        obj = MovingPiece(f'body{group}_{part}', [rlo, lo], [rhi, hi], dim=2,
                          mat=0, group=group, particleType=particle_type, flaggedSurfaces=flags)
        obj.particleColor = group
        obj.czTag = 1
        pfw['objects'].append(geom.surfaceFlagWrapper(obj.name+'_cz', obj, 3))

# ---------------------------------------------------------------------------
# Lagrangian histories show the tensile pulse and subsequent return or
# separation. Tracer-centre spacing is not the CZ opening: it includes bulk
# particle deformation and must not be used directly in the cohesive law.
# ---------------------------------------------------------------------------
tracers.set_tracers(
    pfw, points=[[0.00425, -0.00025, 0.], [0.00425, 0.00025, 0.],
                 [0.00425, -0.00775, 0.], [0.00425, 0.00775, 0.]],
    variables=['particleID', 'velocityX', 'velocityY', 'velocityZ', 'stressYY', 'stressZZ'],
    write_interval=case['end_time']/200., output_prefix='czTracer')
