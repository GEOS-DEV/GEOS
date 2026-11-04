"""Spatial pressure loads: annular cylinder (default) or solid sphere.

AXIS_PRESSURE_CASE=borehole: independent inner and outer pressure ramps.
AXIS_PRESSURE_CASE=sphere: external pressure selected by r*r+y*y > Rcut*Rcut.
All pressure values are Pa, positive in compression; x=r, y=axial, zz=hoop.
"""
import os
import pfw_geometryObjects as geom

case = os.environ.get('AXIS_PRESSURE_CASE', 'borehole')
if case not in ('borehole', 'sphere'):
    raise ValueError('AXIS_PRESSURE_CASE must be borehole or sphere')

pfw = dict(
    mSubmitJobs=True, mBatch=True, mWallTime='00:10:00',
    formulation='Axisymmetric',
    xmin=.005, xmax=.025, ymin=-.02, ymax=.02, nI=22, nJ=42,
    ppc=2, xpar=1, ypar=1, zpar=1, periodic=[False,False],
    # Unconstrained radial motion; axial symmetry planes impose plane strain.
    boundaryConditionTypes=[0,0,1,1], enableContact=0,
    materials=['solid'], materialPropertyString='<ElasticIsotropic name="solid" defaultDensity="2000" '
      'defaultYoungModulus="1e8" defaultPoissonRatio="0.25"/>',
    particleFileFields=['Velocity','MaterialType','ContactGroup','ParticleColor'],
    particleGeometry='fromRefinement', cpdiDomainScaling=1, cpdiDomainScalingType='homel',
    timeIntegrationOption='ExplicitDynamic', updateMethod='FLIP', damageFieldPartitioning=0,
    endTime=.007, plotInterval=.00025, restartInterval=.0035, outputType='silo',
    boxAverageHistory=1, boxAverageWriteInterval=.00005,
    boxAverageMin=[0,-1,0], boxAverageMax=[1,1,0],
    useEvents=1,
)
pfw['runDebug'] = True

if case == 'borehole':
    # Material a=.01 < r < b=.02. The selected nodes at r=a and r=b load
    # the inner/outer boundaries. Keep a < r0 < r1 < b; cutoffs are fixed in
    # physical space and must still cover the boundary support as it moves.
    pfw['objects'] = [geom.box('annulus',[.01,-.02],[.02,.02],dim=2,mat=0,group=0,particleType=2)]
    pfw['mpmEventsString'] = '''
    <BoreholePressure name="innerPressure" startTime="0" endTime="0.005"
      boreholeRadius="0.0105" startPressure="0" endPressure="100000" interpType="2"/>
    <ConfiningPressure name="outerPressure" startTime="0" endTime="0.005"
      confiningPressureBoxMin="{-1,-1,0}" confiningPressureBoxMax="{0.0195,1,0}"
      startPressure="0" endPressure="20000" interpType="2"/>
    '''
else:
    # Positive-r half of a meridional disk becomes a full sphere. Only the
    # revolution axis is constrained. The curved boundary is particle-resolved.
    pfw.update(xmin=0., xmax=.025, ymin=-.025, ymax=.025, nI=27, nJ=52,
               boundaryConditionTypes=[1,0,0,0])
    pfw['objects'] = [geom.disk('sphere',[0.,0.],.02,mat=0,group=0,particleType=2)]
    pfw['mpmEventsString'] = '''
    <ConfiningPressure name="spherePressure" startTime="0" endTime="0.005"
      confiningPressureSphereRadius="0.0195"
      startPressure="0" endPressure="20000" interpType="2"/>
    '''
