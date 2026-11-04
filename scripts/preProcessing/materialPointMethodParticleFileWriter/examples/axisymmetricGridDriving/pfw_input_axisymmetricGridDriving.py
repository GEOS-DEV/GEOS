"""Cylinder in a moving cylindrical wall and axial platens.

Select AXIS_GRID_MODE=stretch (default), stress, or mixed. Optional
AXIS_GRID_WALL=contact uses unilateral walls (FLIP); Moving is bilateral.
The plane is (r,y), radius=height=1; quantities represent complete rings.
"""
import os
import pfw_geometryObjects as geom

mode = os.environ.get('AXIS_GRID_MODE', 'stretch')
wall = os.environ.get('AXIS_GRID_WALL', 'moving')
update = os.environ.get('AXIS_GRID_UPDATE', 'FLIP')
if mode not in ('stretch', 'stress', 'mixed') or wall not in ('moving', 'contact'):
    raise ValueError('AXIS_GRID_MODE: stretch/stress/mixed; AXIS_GRID_WALL: moving/contact')
if update not in ('FLIP', 'PIC', 'XPIC', 'FMPM'):
    raise ValueError('AXIS_GRID_UPDATE: FLIP/PIC/XPIC/FMPM')
face_type = 2 if wall == 'moving' else 3
pfw = dict(
    mSubmitJobs=False, mBatch=True, mBank='test', mWallTime='00:10:00',
    formulation='Axisymmetric', xmin=0., xmax=1., ymin=0., ymax=1.,
    nI=22, nJ=22, ppc=2, xpar=1, ypar=1, zpar=1, periodic=[False,False],
    boundaryConditionTypes=[1,face_type,1,2], enableContact=0,
    boundaryFaceFrictionCoefficients=[0,0,0,0,0,0],
    prescribedBoundaryFTable=1,
    # Grid stretches: time, Frr, Fyy, inactive coordinate (always 1).
    fTable=[[0.,1.,1.,1.], [.2,.99,1.,1.], [.3,.99,1.,1.]],
    fTableInterpType='Smoothstep',
    endTime=.3, plotInterval=.01, restartInterval=.1, outputType='silo',
    materials=['solid'], materialPropertyString='<ElasticIsotropic name="solid" defaultDensity="1000" '
      'defaultYoungModulus="1e6" defaultPoissonRatio="0.25"/>',
    particleFileFields=['Velocity','MaterialType','ContactGroup','ParticleColor',
                        'SurfaceFlag','SurfaceNormal','SurfacePosition'],
    particleGeometry='fromRefinement', cpdiDomainScaling=1, cpdiDomainScalingType='homel',
    timeIntegrationOption='ExplicitDynamic', updateMethod=update, updateOrder=2,
    damageFieldPartitioning=0, reactionHistory=1, reactionWriteInterval=.001,
    boxAverageHistory=1, boxAverageWriteInterval=.001, boxAverageResizeWithDomain=1,
    boxAverageMin=[0,0,0], boxAverageMax=[1,1,0],
)
if mode in ('stress','mixed'):
    pfw.update(stressControl=[1,0,0], stressControlKp=.01,
               stressTable=[[0.,0.,0.,0.], [.2,-1e4,0.,0.], [.3,-1e4,0.,0.]],
               stressTableInterpType='Smoothstep')
    # The controller replaces the radial F-table column. Axial F is either
    # held at 1 (stress) or independently ramped to .995 (mixed).
    axial = .995 if mode == 'mixed' else 1.
    pfw['fTable'] = [[0.,1.,1.,1.], [.2,1.,axial,1.], [.3,1.,axial,1.]]
pfw['objects'] = [geom.box('cylinder',[0.,0.],[1.,1.],dim=2,mat=0,group=0,particleType=2)]
