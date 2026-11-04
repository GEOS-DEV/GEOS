"""Reconstructed surfaces, full-ring area and capillary-pressure smoke cases.

AXIS_SURFACE_CASE=cylinder (default) or sphere; AXIS_SURFACE_INTERPOLATION=
 mixed (CPDI core, CPTI skin), cpdi, cpti or bspline. x=r; y=axis; zz=hoop.
Use AXIS_SURFACE_TENSION=0 to isolate geometric reconstruction and area.
"""
import os
import pfw_geometryObjects as geom

case=os.environ.get('AXIS_SURFACE_CASE','cylinder')
interpolation=os.environ.get('AXIS_SURFACE_INTERPOLATION','mixed')
if case not in ('cylinder','sphere'): raise ValueError('Choose cylinder or sphere')
if interpolation not in ('mixed','cpdi','cpti','bspline'): raise ValueError('Unknown interpolation')
tension=int(os.environ.get('AXIS_SURFACE_TENSION','1'))
if tension not in (0,1): raise ValueError('AXIS_SURFACE_TENSION must be 0 or 1')
resolution=int(os.environ.get('AXIS_SURFACE_RESOLUTION','24'))
if resolution<8: raise ValueError('Use at least eight radial cells per unit length')
R=1.;H=1.;gamma=.01;rho=1.;bulk=10.
pfw=dict(
    mSubmitJobs=False,mBatch=True,mBank='test',mWallTime='00:10:00',
    formulation='Axisymmetric',xmin=0.,xmax=1.25,ymin=0.,ymax=H,
    nI=int(1.25*resolution)+2,nJ=resolution+2,ppc=2,
    xpar=1,ypar=1,zpar=1,periodic=[False,False],boundaryConditionTypes=[1,0,1,1],
    materials=['liquid'],materialPropertyString=f'<Liquid name="liquid" defaultDensity="{rho}" '
        f'defaultBulkModulus="{bulk}" defaultViscosity="0"/>',
    particleFileFields=['Velocity','MaterialType','ContactGroup','ParticleColor',
                        'SurfaceFlag','SurfaceNormal','SurfacePosition'],
    particleGeometry='fromRefinement',geometryReset='FOnly',cpdiDomainScaling=0,
    timeIntegrationOption='ExplicitDynamic',updateMethod='FLIP',
    enableContact=0,damageFieldPartitioning=0,overlapCorrection='Off',
    computeParticleSurfaceNormalsAndPositions=1,overwriteExistingNormalsAndPositions=1,
    normalAndPositionMethod='DFGAndVolumeIntegration',computeNodalArea=1,
    enableSurfaceTension=tension,surfaceTensionPairs=[[0,-1,gamma]],
    endTime=.05,plotInterval=.005,restartInterval=.025,outputType='silo',
    plotGridFields=1,gridFieldNames=['gridSurfaceArea','gridBasedSurfaceNormal',
      'gridBasedSurfacePosition','gridSurfaceMaterialVolume','gridSurfaceTensionForce'],
    boxAverageHistory=1,boxAverageWriteInterval=.001,
    boxAverageMin=[0,-2,0],boxAverageMax=[2,2,0],
    useEvents=1,
)
# Constant initial Laplace pressure provides a refinement diagnostic. The
# diffuse discretization is not exactly well-balanced; residual waves remain.
pressure=tension*gamma/R*(2 if case=='sphere' else 1)
pfw['mpmEventsString']=f'<InitializeStress name="laplacePressure" pressure="{pressure}" targetRegion="all"/>'

class UnpaintedBox(geom.box):
    def getSurfaceNormal(self,point): return [0.,0.,0.]
    def getSurfacePosition(self,point): return [0.,0.,0.]

if case=='cylinder':
    # Only r=R is a physical free surface. The CPDI/CPTI seam is the same body.
    pieces=[(0.,.75,2,[False]*4),(.75,R,3,[False,False,True,False])] if interpolation=='mixed' else [
        (0.,R,{'cpdi':2,'cpti':3,'bspline':1}[interpolation],[False,False,True,False])]
    pfw['objects']=[UnpaintedBox(f'piece{k}',[lo,0.],[hi,H],dim=2,mat=0,group=0,
        particleType=ptype,flaggedSurfaces=flags) for k,(lo,hi,ptype,flags) in enumerate(pieces)]
else:
    pfw.update(ymin=-1.25,ymax=1.25,nJ=int(2.5*resolution)+2,boundaryConditionTypes=[1,0,0,0])
    # The curved-surface cutter needs the object's analytical normals while
    # generating triangles. GEOS overwrites those normals in its first pass.
    ptype={'mixed':2,'cpdi':2,'cpti':3,'bspline':1}[interpolation]
    pfw['objects']=[geom.disk('sphere',[0.,0.],R,mat=0,group=0,particleType=ptype)]
    if interpolation=='mixed':
        pfw.update(surfaceRefinementMethod='simplex',surfaceRefinementLevel=1,simplexParticleType=3,
                   surfaceRefinementFallback='cartesian',surfaceRefinementPlaneTolerance=.05)
