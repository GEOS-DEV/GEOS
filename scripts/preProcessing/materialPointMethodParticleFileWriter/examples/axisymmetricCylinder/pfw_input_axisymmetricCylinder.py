"""Full solid cylinder: CPDI core, CPTI outer annulus, radial contraction pulse.

Coordinates are (r, axial position); physical radius=1, length=1, density=1000.
No contact, external tractions, stochastic strength or moving-grid controls.
"""
import pfw_geometryObjects as geom

pfw=dict(
    mSubmitJobs=False, mBatch=True, mBank='test', mWallTime='00:05:00',
    formulation='Axisymmetric', xmin=0.0, xmax=1.2, ymin=0.0, ymax=1.2,
    nI=14, nJ=14, ppc=2, xpar=1, ypar=1, zpar=1, periodic=[False,False],
    boundaryConditionTypes=[1,0,1,0], enableContact=0, damageFieldPartitioning=0,
    endTime=0.002, plotInterval=0.0002, restartInterval=0.001, outputType='silo',
    materials=['solid'],
    materialPropertyString='<ElasticIsotropic name="solid" defaultDensity="1000" '
                           'defaultYoungModulus="1e6" defaultPoissonRatio="0.25"/>',
    particleFileFields=['Velocity','MaterialType','ContactGroup','ParticleColor',
                        'SurfaceFlag','SurfaceNormal','SurfacePosition'],
    particleGeometry='fromRefinement', geometryReset='FOnly',
    cpdiDomainScaling=1, cpdiDomainScalingType='homel',
    timeIntegrationOption='ExplicitDynamic', updateMethod='FLIP', boxAverageHistory=1,
)
class CylinderPiece(geom.box):
    def getVelocity(self, point):
        return [-0.1*point[0],0.0,0.0]
    def getParticleColor(self, point):
        return self.particleType

pfw['objects']=[CylinderPiece('core',[0.,0.],[0.6,1.],dim=2,mat=0,group=0,particleType=2),
                CylinderPiece('skin',[0.6,0.],[1.,1.],dim=2,mat=0,group=0,particleType=3)]
