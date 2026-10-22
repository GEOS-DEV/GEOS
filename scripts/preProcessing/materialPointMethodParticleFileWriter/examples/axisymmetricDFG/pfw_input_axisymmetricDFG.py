"""Three-body contact and same-body DFG opening, with no cohesive zones.

intact/damage: a continuous cylinder initially moves away from y=0. Only
damage has a pre-damaged band; neither prescribes an internal surface normal.
contact: three cylinders approach with independent contact groups.
combined: the three-body case also has a damaged band in its central body.
"""
#[pfw_dependency] input:caseSettings.py
import pfw_geometryObjects as geom
import pfw_tracerPoints as tracers
from caseSettings import settings

case = settings()
multi = case['groups'] == 3

# A true r-y grid with 1 mm cells, CPDI bulk and CPTI outer skin. The radial
# axis is symmetric; all other mesh boundaries are remote outflow boundaries.
pfw = dict(
    mSubmitJobs=True, mBatch=True, mWallTime='00:10:00',
    formulation='Axisymmetric', xmin=0., xmax=.009, ymin=-.012, ymax=.012,
    nI=11, nJ=26, ppc=2, xpar=1, ypar=1, zpar=1,
    periodic=[False, False], boundaryConditionTypes=[1, 0, 0, 0],
    timeIntegrationOption='ExplicitDynamic', updateMethod=case['update'], updateOrder=2,
    cflFactor=.2, initialDt=1.e-9, endTime=case['end_time'],
    plotInterval=case['end_time']/50., restartInterval=case['end_time']/2.,
    outputType='silo', particleGeometry='fromRefinement', geometryReset='FOnly',
    cpdiDomainScaling=1, cpdiDomainScalingType='homel',
    particleFileFields=['Velocity', 'MaterialType', 'ContactGroup', 'ParticleColor',
                        'Damage', 'SurfaceFlag', 'SurfaceNormal', 'SurfacePosition'],
    boxAverageHistory=1, boxAverageWriteInterval=case['end_time']/200.,
    plotGridFields=1, gridFieldNames=['gridMass', 'gridVelocity', 'gridDamage',
                                    'gridMaxDamage', 'gridDamageGradient',
                                    'gridSurfaceNormal', 'gridContactForce'],
)
pfw['runDebug'] = True
# DFG selects side*numberOfGroups+group. Its existing damage and quality
# criteria decide whether same-group fields tie or can separate/contact.
# No explicit surface positions are available for the internal damage band.
pfw.update(enableContact=1, damageFieldPartitioning=case['dfg'], neighborRadius=.0015,
           separabilityMinDamage=.5, surfaceQualityThreshold=.5,
           maxSingleFieldStateFractionForSeparability=-1., frictionCoefficient=.2,
           contactNormalType='MassWeighted', contactGapCorrection='Simple',
           useSurfacePositionForContact=0, explicitSurfaceNormalInfluence=0.,
           computeParticleSurfaceNormalsAndPositions=0, overlapCorrection='Off',
           resetDefGradForFullyDamagedParticles=0)
pfw['materials'] = ['solid']
pfw['materialPropertyString'] = (
    '<ElasticIsotropic name="solid" defaultDensity="{density}" '
    'defaultYoungModulus="{E}" defaultPoissonRatio="0"/>'.format(**case))


class Piece(geom.box):
    def getVelocity(self, point):
        # Outer bodies approach the centre; the damaged central body's two
        # halves initially separate slowly and are subsequently compressed.
        if multi and self.group != 1:
            vy = .5 if self.group == 0 else -.5
        else:
            vy = case['speed'] if point[1] > 0 else -case['speed']
            if case['mode'] == 'contact':
                vy = 0.
        return [0., vy, 0.]

    def getDamage(self, point):
        # This is a prescribed initial damage field, not a new failure law.
        # No gap, contact-group split or explicit crack normal is inserted.
        return float(case['mode'] in ('damage', 'combined') and abs(point[1]) < .00051)


pfw['objects'] = []
pieces = [(0, -.008, -.003), (1, -.002, .002), (2, .003, .008)] if multi else [(0, -.006, .006)]
for group, lo, hi in pieces:
    for part, rlo, rhi, particle_type in [('core', 0., .005, 2), ('skin', .005, .006, 3)]:
        obj = Piece(f'body{group}_{part}', [rlo, lo], [rhi, hi], dim=2,
                    mat=0, group=group, particleType=particle_type,
                    flaggedSurfaces=[False, False, False, False])
        obj.particleColor = group
        pfw['objects'].append(obj)

# Tracers compare the interior tensile response with outer-body motion.
end = .00775 if multi else .00575
tracers.set_tracers(
    pfw, points=[[.00325, -.00075, 0.], [.00325, .00075, 0.],
                 [.00325, -end, 0.], [.00325, end, 0.]],
    variables=['particleID', 'velocityX', 'velocityY', 'velocityZ', 'stressYY', 'stressZZ'],
    write_interval=case['end_time']/200., output_prefix='dfgTracer')
