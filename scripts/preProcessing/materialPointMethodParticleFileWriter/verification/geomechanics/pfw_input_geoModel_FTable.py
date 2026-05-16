# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

# This is currently just a smoke test to see if the geomechanics model is implemented
# successfully and runs. 
#
# TODO: add some actual verification, so make sure the response is correct.


pfw = {}
pfw["runDebug"] = False
stopTime = 100.0

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
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mNodes"]=int(np.ceil(float(pfw["mCores"])/36.)) 
pfw["mSubmitJobs"]=True

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/100
pfw["restartInterval"]=stopTime
pfw['lastRestartBufferInSeconds'] = 0.

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25  
pfw["initialDt"]=1e-16
pfw["reactionHistory"]=1
pfw["reactionWriteInterval"]=stopTime/2000
pfw["boxAverageHistory"]=1
pfw["boxAverageWriteInterval"]=stopTime/2000

pfw["solverProfiling"]=0         
pfw["frictionCoefficient"]=0.25  

pfw["updateMethod"]="XPIC"
pfw["updateOrder"]=2

# MATERIAL PROPERTIES --------------------------------------------------------------------

# Will be read from input file for plotting.
density = 2.648

# Material model parameters.  Will be read from input file for plotting.
MPa = 0.001 # GPa
density = 2.3

# nonlinear bulk modulus
b0 = 2.003 
b1 = 36.33
b2 = 650.0*MPa
b3 = 1.0
b4 = 4.0e-3 

# nonlinear shear modulus
g0 = 0.40128  # (0.924<shear modulus<1.502)
g1 = 0.088
g2 = 0.0001
g3 = 0.0
g4 = 0.

# crush curve
BETA_nonassociativity = 1.16 # 0.5<beta<2.
p0 = -0.059
p1 = 4.38  # GPa^-1
p2 = 0.0
p3 = 0.255
p4 = 0.0

CR = 0.5  # ZZZZ 0.2<CR<0.8  ####For TP1, kep CR close to 0.2-.4. 0.8 is too high (8/9)
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
strainHardeningK = 0.0075 # (GPa), STREN = STREN_i + K*(1.-exp(-n*gamma_p))
strainHardeningn = 350.0 # STREN = STREN_i + K*(1.-exp(-n*gamma_p))

# set creep flag to 0 to disable all creep options.
enableCreep = 0

# Porous compaction creep:
A = 0.238    # creep: phi_e = A*exp(-p/B), note that A is 0-pressure equilbrium porosity and phi_i = 1 -exp(-p3)
B = 0.08    # creep: phi_e = A*exp(-p/B), B is GPa
C = 5.7e-09   # creep rate: dphi/dt = -p*C*(phi-phi_e)

# deviatoric creep:
creepc0 = 2.0e-11  # dev creep rate term: creep strain sincrement = c0 * sigma_vm * std::pow( elasticShearStrain - equilibriumShearStrain, c1 ) * dt;
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
   strainHardeningN="""+'"'+str(strainHardeningn)+'"'+"""
   strainHardeningK="""+'"'+str(strainHardeningK)+'"'+"""
   />
"""

# GEOMETRY OBJECTS -------------------------------------------------------

block = geom.box('block',[pfw["xmin"],pfw["ymin"],pfw["zmin"]],[pfw["xmax"],pfw["ymax"],pfw["zmax"]],v=[0.0,0.0,0.0],mat=0,group=0)
pfw["objects"]=[block]

# DEFORMATION -----------------------------------------------------------------------------

pfw["fTableInterpType"]='Smoothstep'
pfw["prescribedBoundaryFTable"]=1
# pfw["fTable"]=[[0,	 1,	    1,	1],
#                [.4*stopTime,	 0.998,	0.998,	0.998],
#                [.8*stopTime,	 0.996,	0.998,	0.998],
#                [stopTime,	 1,	    1,	1]
#                ]
pfw["fTable"]=[[0,	          1.000,	    1,	1],
[.1*stopTime,	 0.995,0.995,0.995]
[.2*stopTime,	 0.998,0.998,0.998]
[.3*stopTime,	 0.980,0.980,0.980]
[.4*stopTime,	 0.985,0.985,0.985]
[.5*stopTime,	 0.960,0.960,0.960]
[.6*stopTime,	 0.965,0.965,0.965]
[.7*stopTime,	 0.940,0.940,0.940]
[.8*stopTime,	 0.945,0.945,0.945]
[.9*stopTime,	 0.920,0.920,0.920]
[stopTime,	     0.925,0.925,0.925]
]

# prescribed deformation (moving pistons) at all faces
pfw["boundaryConditionTypes"]=[2, 2, 2, 2, 2, 2]