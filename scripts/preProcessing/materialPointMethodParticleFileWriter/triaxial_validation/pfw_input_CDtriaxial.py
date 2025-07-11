import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree   


pfw = {}
pfw["runDebug"] = False
stopTime = 1000.0
confiningPressure = 0.0032  # confining pressure GPa (make sure this p < -p0/3)
maxCompressiveStrain = 0.035

# PID Control Parameters
kp = 0.1
ki = 0.00
kd = 0.01

# DOMAIN ---------------------------------------------------------------------------------

sampleWidth = 1.0  # mm
sampleHeight = 1.0 # mm
sampleLength = 1.0 # mm

domainWidth = sampleWidth  # This would be increased for unconfined compression.
domainHeight = sampleHeight
domainLength = sampleLength

pfw["xmin"] = 0.0             # mm
pfw["xmax"] = domainWidth    # mm
pfw["ymin"] = 0.0 # mm
pfw["ymax"] = domainHeight # mm
pfw["zmin"] = 0.0 # mm
pfw["zmax"] = domainLength # mm

refine=1  # partitions in each direction
cpp=3     # cells per partition in each direction

pfw["xpar"]=refine
pfw["ypar"]=refine
pfw["zpar"]=refine

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cpp  # grid cells in the z-direction
pfw["ppc"]=2               # particles per cell in each direction

# BATCH PARAMETERS  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"]="00:05:00"
pfw["mSubmitJobs"]=False

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/100
pfw["restartInterval"]=5*stopTime
pfw['lastRestartBufferInSeconds'] = 0.

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.025
pfw["initialDt"]=1e-16
pfw["reactionHistory"]=1
pfw["reactionWriteInterval"] = stopTime/1000
pfw["boxAverageHistory"] = 1
pfw["boxAverageWriteInterval"] = stopTime/1000

pfw["solverProfiling"]=0         
pfw["frictionCoefficient"]=0.25  

pfw["updateMethod"]="PIC" # best for single element.
pfw["updateOrder"]=2      # This only needs to be specified with XPIC or FMPM.
pfw["useAPIC"]=0          # This option should only be used for single-element simulations, and will improve accuracy of reactions to be more consistent with box data.

#pfw["singleElementDriver"] = 0
#pfw["useConstantTimeStep"] = 1
#pfw["constantTimeStepValue"] = stopTime/10000

pfw["useEvents"] = 1


# MATERIAL PROPERTIES For concrete in the ceramic damage model---------------------------------------

density = 2.75
bulk = 275./( 3.*( 1. - 2.*0.3) )
shear = 275./( 2.*( 1. + 0.3) )
tensileStrength = 0.005
compressiveStrength = 8.*tensileStrength
maximumStrength = 5.0*compressiveStrength
crackSpeed = 1.8
thirdInvariantDependence=1


pfw["materials"] = [ "concrete" ]
pfw["materialPropertyString"]="""
 <CeramicDamage

    name="concrete"
    defaultDensity="""+'"'+str(density)+'"'+"""
    defaultBulkModulus="""+'"'+str(bulk)+'"'+"""
    defaultShearModulus="""+'"'+str(shear)+'"'+"""
    tensileStrength="""+'"'+str(tensileStrength)+'"'+"""
    compressiveStrength="""+'"'+str(compressiveStrength)+'"'+"""
    maximumStrength="""+'"'+str(maximumStrength)+'"'+"""
    crackSpeed="""+'"'+str(crackSpeed)+'"'+"""
    thirdInvariantDependence=""" + '"' + str(thirdInvariantDependence) + '"' + """/>

"""


# GEOMETRY OBJECTS -------------------------------------------------------
# single block filling domain for single-element test.

block = geom.box('block',[pfw["xmin"],pfw["ymin"],pfw["zmin"]],[pfw["xmax"],pfw["ymax"],pfw["zmax"]],v=[0.0,0.0,0.0],mat=0,group=0)
pfw["objects"]=[block]

# DEFORMATION -----------------------------------------------------------------------------

# We initialize the domain with p=confiningPressure, if p > -p0/3 this will fail.
# we then apply a stress boundary condition on x,y faces to be in equilibrium, and do a 
# strain-controlled motion in z.

# Ftable only controls z-direction
pfw["boundaryConditionTypes"]=[ 2, 2, 2, 2, 2, 2 ]
pfw["fTableInterpType"]='Cosine'
pfw["prescribedBoundaryFTable"] = 1
pfw["fTable"]=[
    [0.,	       1.0, 1.0, 1.0],
    [.6*stopTime,  1.0, 1.0, np.exp(-maxCompressiveStrain) ],
    [stopTime,	   1.0, 1.0, np.exp(-.77*maxCompressiveStrain) ]
    ]

# stress table controls x- and y- directions.
pfw["stressControl"]=[ 1, 1, 0]
pfw["stressTableInterpType"] = 'Cosine'
pfw["stressControlKp"] = kp
pfw["stressControlKi"] = ki
pfw["stressControlKd"] = kd
pfw["stressTable"]=[[0.0,      	   -confiningPressure, -confiningPressure, -confiningPressure],
					[stopTime,     -confiningPressure, -confiningPressure, -confiningPressure]]

# This should initialize stress in equilibrium with boundary conditions:
pfw["mpmEventsString"]="""
<MPMEvents>
    <InitializeStress 
        time="0.0"
        interval="0.1"
        targetRegion="all"
        pressure=""" + '"' + str(confiningPressure) + '"' + """
        />
</MPMEvents>
"""