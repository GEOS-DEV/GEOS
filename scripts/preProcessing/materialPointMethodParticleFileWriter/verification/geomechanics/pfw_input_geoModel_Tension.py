# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

# This is currently just a smoke test to see if the geomechanics model is implemented
# successfully and runs. 
#
# Initializes material with a confining stress (should be p < -p0/3 ) and then ramps
# load in 1 direction with that confinement maintained through stress control on
# lateral boundaries.
# TODO: add some actual verification, so make sure the response is correct.


pfw = {}
pfw["runDebug"] = False
stopTime = 1000.0

physicalStopTime = 360*24*3600.e6 # actual end of 1-year test time, micro-seconds.
timeScale = physicalStopTime / stopTime  # strain rate multiplier to match 1-year creep test in stopTime

#confiningPressure = 0.005  # confining pressure GPa
#maxCompressiveStrain = 0.05 # stress control end point (GPa)
confiningPressure = 0.00  
maxCompressiveStrain = -0.1

fractureEnergyReleaseRate = 1.e-10
fractureStress = 0.005

kp = 0.05
ki = 0.00
kd = 0.005
pfw["cflFactor"]=0.005



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
pfw["mSubmitJobs"]=True

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/100
pfw["restartInterval"]=5*stopTime
pfw['lastRestartBufferInSeconds'] = 0.

pfw["timeIntegrationOption"]="ExplicitDynamic"

pfw["initialDt"]=1e-16
pfw["reactionHistory"]=1
pfw["reactionWriteInterval"] = stopTime/1000
pfw["boxAverageHistory"] = 1
pfw["boxAverageWriteInterval"] = stopTime/1000

pfw["solverProfiling"]=0         
pfw["frictionCoefficient"]=0.25  

pfw["updateMethod"]="PIC"
pfw["updateOrder"]=2

#pfw["singleElementDriver"] = 0
#pfw["useConstantTimeStep"] = 1
#pfw["constantTimeStepValue"] = stopTime/10000

pfw["useEvents"] = 1

# MATERIAL PROPERTIES --------------------------------------------------------------------

# Material model parameters.  Will be read from input file for plotting.
MPa = 0.001 # GPa
density = 2.3

# nonlinear bulk modulus
# b0 = 2.003 
# b1 = 36.33
# b2 = 650.0*MPa
# b3 = 1.0
# b4 = 4.0e-3 
b0 = 1.67 #low pressure bulk modulus
b1 = 30.0 #high pressure bulk modulus
b2 = 0.3 #shape parameter
b3 = 1.42 
b4 = 0.015

# nonlinear shear modulus defined with a pressure-dependent poisson's ratio.
# Two options, 
#   (i) constant shear modulus ( g0>0, g1=g2=0)
#   (ii) pressure-depenendent poissons ratio,  ( 0 < g1+g2 < 0.5 )
#        low pressure poisson's ratio: g1+g2
#        high-pressure poisson's ratio: g1
#
#        G = g0 for plastically dilated states, for I1<0 and evp <=0,
#        G = 1.5*kp*(1-2*nu)/(1+nu), where kp=b0+b1*exp(b2/I1), and nu = g1+g2*exp(g3/I1) 
#
#        in this mode, g0 will be the initial or dilated shear modules, for consistency
#        with poisson control, set:
#        nu = m_g1 + m_g2;
#        g0 = 1.5 * b0 * ( 1.0 - 2.0 * nu ) / ( 1.0 + nu );
#
# g0 = 0.40128  # dilated shear modulus (g0>0), units of GPa
# g1 = 0.25    # low-pressure poisson's ratio is nu = g1+g2 (should be 0<nu<0.5)
# g2 = -0.15   # high-pressure poisson's ratio is nu = g1
# g3 = 0.0     # unused for now.
# g4 = 0.      # unused for now.

g1 = 0.3    # low-pressure poisson's ratio is nu = g1+g2 (should be 0<nu<0.5)
g2 = -0.15   # high-pressure poisson's ratio is nu = g1
nu = g1 + g2
g0 = 1.5 * b0 * ( 1.0 - 2.0 * nu ) / ( 1.0 + nu )  # dilated shear modulus (g0>0), units of GPa
g3 = 0.3     #
g4 = 0.      # unused for now.

# crush curve
BETA_nonassociativity = 1.16 # 0.5<beta<2.
# p0 = -0.059
# p1 = 4.38  # GPa^-1
# p2 = 0.0
# p3 = 0.255
# p4 = 0.0
p0 = -0.03
p1 = 0.65  # GPa^-1
p2 = 0.0
p3 = 0.255
p4 = 0.0
CR = 0.9  # ZZZZ 0.2<CR<0.8  ####For TP1, kep CR close to 0.2-.4. 0.8 is too high (8/9)
Kf = 0.0   # fluid bulk modulus (2.2 GPa)
pfi = 0.0  # initial fluid pressure (GPa)

# yield surface
#These following values are from Figure 5 of Kikibas et al
PEAKI1 = 0.0161
FSLOPE = 0.613
STREN = 0.0135
YSLOPE = 0.0001
BETA_nonassociativity = 1.16

# fluid pressure (not used)
Kf = 0.0   # fluid bulk modulus (2.2 GPa)
pfi = 0.0  # initial fluid pressure (GPa)

# overstress model (not used)
t1RateDependence=0.0
t2RateDependence=0.0

# If the input value for STREN above is treated as an initial value STREN_i
# the strength will be a function of deviatoric plastic strain:
# Note STREN goes to STREN_i + K, as n*gamma_p is large:

# set K=0 to disable hardening:
strainHardeningK = (0.)*0.0075 # (GPa), STREN = STREN_i + K*(1.-exp(-n*gamma_p))
strainHardeningn = 350.0 # STREN = STREN_i + K*(1.-exp(-n*gamma_p))

# set creep flag to 0 to disable all creep options.
enableCreep = 0

# Porous compaction creep:
A = 0.238    # creep: phi_e = A*exp(-p/B), note that A is 0-pressure equilbrium porosity and phi_i = 1 -exp(-p3)
B = 0.08     # creep: phi_e = A*exp(-p/B), B is GPa
C = timeScale*5.7e-09   # creep rate: dphi/dt = -p*C*(phi-phi_e)
creepD = 1.
creepE = 0.
creepF = 1.

# deviatoric creep:
creepc0 = timeScale*2.0e-11  # dev creep rate term: creep strain sincrement = c0 * sigma_vm * std::pow( elasticShearStrain - equilibriumShearStrain, c1 ) * dt;
creepc1 = 1.0      # dev creep rate term: creep strain sincrement = c0 * sigma_vm * std::pow( elasticShearStrain - equilibriumShearStrain, c1 ) * dt;

# bulk and shear are used for initial time-step calculations and 
# some normalizations and are defined from the geo input parameters
bulk = b0+b1
shear = g0

# These parameters are for the fracture propagation.  We don't have any data
# to constrain them so these are best guesses.
KIc = 1.591*(.001)*np.sqrt(1000.)       # fracture toughness GPa mm^1/2
fractureRadius = 0.5*(KIc**2)/(np.pi*PEAKI1**2)  # radius of fracture process zone
constrainedModulus = bulk + 4./3.*shear
Gf = KIc/constrainedModulus

pfw["materials"] = [ "ghareb" ]
pfw["materialPropertyString"]="""
<Geomechanics
   name="ghareb"
   defaultDensity="""+'"'+str(density)+'"'+"""
   b0="""+'"'+str(b0)+'"'+"""
   b1="""+'"'+str(b1)+'"'+"""
   b2="""+'"'+str(b2)+'"'+"""
   b3="""+'"'+str(b3)+'"'+"""
   b4="""+'"'+str(b4)+'"'+"""
   g0="""+'"'+str(g0)+'"'+"""
   g1="""+'"'+str(g1)+'"'+"""
   g2="""+'"'+str(g2)+'"'+"""
   g3="""+'"'+str(g3)+'"'+"""
   g4="""+'"'+str(g4)+'"'+"""
   p0="""+'"'+str(p0)+'"'+"""
   p1="""+'"'+str(p1)+'"'+"""
   p2="""+'"'+str(p2)+'"'+"""
   p3="""+'"'+str(p3)+'"'+"""
   p4="""+'"'+str(p4)+'"'+"""
   cr="""+'"'+str(CR)+'"'+"""
   fluidBulkModulus="""+'"'+str(Kf)+'"'+"""
   fluidInitialPressure="""+'"'+str(pfi)+'"'+"""
   t1RateDependence="""+'"'+str(t1RateDependence)+'"'+"""
   t2RateDependence="""+'"'+str(t2RateDependence)+'"'+"""
   peakI1="""+'"'+str(PEAKI1)+'"'+"""
   fSlope="""+'"'+str(FSLOPE)+'"'+"""
   stren="""+'"'+str(STREN)+'"'+"""
   ySlope="""+'"'+str(YSLOPE)+'"'+"""
   beta="""+'"'+str(BETA_nonassociativity)+'"'+"""
   enableCreep="""+'"'+str(enableCreep)+'"'+"""
   creepC0="""+'"'+str(creepc0)+'"'+"""
   creepC1="""+'"'+str(creepc1)+'"'+"""
   creepA="""+'"'+str(A)+'"'+"""
   creepB="""+'"'+str(B)+'"'+"""
   creepC="""+'"'+str(C)+'"'+"""
   creepD="""+'"'+str(creepD)+'"'+"""
   creepE="""+'"'+str(creepE)+'"'+"""
   creepF="""+'"'+str(creepF)+'"'+"""
   strainHardeningN="""+'"'+str(strainHardeningn)+'"'+"""
   strainHardeningK="""+'"'+str(strainHardeningK)+'"'+"""
   fractureEnergyReleaseRate="""+'"'+str(fractureEnergyReleaseRate)+'"'+"""
   fractureStress="""+'"'+str(fractureStress)+'"'+"""
   
   />
"""

# GEOMETRY OBJECTS -------------------------------------------------------

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
    [0,	          1.000,	    1.,	1.],
    [stopTime,	    1.0,1.0, np.exp(-maxCompressiveStrain) ]
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