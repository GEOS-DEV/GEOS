# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree
import pfw_materials as matdb

# This is currently just a smoke test to see if the geomechanics model is implemented
# successfully and runs. 
#
# Initializes material with a confining stress (should be p < -p0/3 ) and then ramps
# load in 1 direction with that confinement maintained through stress control on
# lateral boundaries.
# TODO: add some actual verification, so make sure the response is correct.


#V3: Mike added damage into the code. 
#V$: Got better values for CR and related yield envelope variables
#V6; has variables for the fracture damage updated. like kappa_i, damageEvolutionCriteron(s)

pfw = {}
pfw["runDebug"] = False
stopTime = 1000.0



confiningPressure = 0.00203   # confining pressure GPa
maxCompressiveStrain = 0.0704 # F-table (strain) control end point (GPa)

# PID Control Parameters
kp = 0.05
ki = 0.00
kd = 0.005




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

pfw["updateMethod"]="PIC" # best for single element or quasistatic loading where dissipation isn't an issue.
pfw["updateOrder"]=2      # This only needs to be specified with XPIC or FMPM.

pfw["useEvents"] = 1

# MATERIAL PROPERTIES --------------------------------------------------------------------

# Material model parameters.  Will be read from input file for plotting.
MPa = 0.001     # GPa
density = 1.57  # from excel sheet (only affects dynamics)

# nonlinear bulk modulus (tune to unload slopes in HYD compression)
b0 = 1.67
b1 = 30.0 # high pressure bulk modulus ia b0+b1 (GPa)
b2 = 0.3  # shape parameter for pressure dependence
b3 = 1.42 # elastic-plastic coupling, reduction in bulk modulus with plastic vol. strain
b4 = 0.015 # elastic-plastic coupling shape parameter.

# nonlinear shear modulus defined with a pressure-dependent poisson's ratio.
# Two options, 
#   (i) constant shear modulus ( g0>0, g1=g2=0)
#   (ii) pressure-depenndent poissons ratio,  ( 0 < g1+g2 < 0.5 )
#        low pressure poisson's ratio: g1+g2
#        high-pressure poisson's ratio: g1
#
#        G = g0 for plastically dilated states, for I1<0 and evp <=0,
#        G = 1.5*kp*(1-2*nu)/(1+nu), where kp=b0+b1*exp(b2/I1), and nu = g1+g2*exp(g3/I1), I1<0
#
#        in this mode, g0 will be the initial or dilated shear modules, for consistency
#        with poisson control, set:
#        nu = m_g1 + m_g2;
#        g0 = 1.5 * b0 * ( 1.0 - 2.0 * nu ) / ( 1.0 + nu );
g1 = 0.1
g2=-0.0001# high-pressure poisson's ratio is nu = g1+g2 (should be 0<nu<0.5)
nu = g1 + g2
g0 = 1.5 * b0 * ( 1.0 - 2.0 * nu ) / ( 1.0 + nu )  # dilated shear modulus (g0>0), units of GPa
g3 = 0.3     #
g4 = 0.      # unused for now.

# crush curve
p0 = min(-.03,-3.0001*confiningPressure) # initial hydrostatic compressive strength, force it to be larger than confining pressure.  VAlue is I1=-3p, negative in compression
p1 = 0.75  # Compaction rate parameter: GPa^-1
p2 = 0.0   # unused
p3 = 0.255 # maximum achievable compressive volumetric plastic strain, relaties to initial porosity:  poro0 = 1-exp(-p3)
p4 = 0.0   # unused

# Yield surface parameters, tuned to initial yield stress (inflection point) in TXC data
CR = 0.2
STREN = 0.50
FSLOPE = 0.18
FSLOPEFAILED = FSLOPE-(0.05)
PEAKI1 = 0.0322/2.
YSLOPE = 0.004/2. #yslope should be less than fslope

# dilation control
BETA_nonassociativity = 1.4 # 0.5<beta<2.

# fluid pressure (not used), set Kf=pfi=0 to disable the legacy fluid effects models, which probably need to be updated.
Kf = 0.0   # fluid bulk modulus (2.2 GPa)
pfi = 0.0  # initial fluid pressure (GPa)

# overstress model (not used).  This was for rate-dependent strengthening, but hasn't been tested since the model
# was ported from Uintah. 
t1RateDependence=0.0
t2RateDependence=0.0

# If the input value for STREN above is treated as an initial value STREN_i
# the strength will be a function of deviatoric plastic strain:
# STREN = STREN_i + K*(1.-exp(-n*gamma_p))
# Note STREN goes to STREN_i + K, as n*gamma_p is large:
# set K=0 to disable hardening:
strainHardeningK = 0.017
strainHardeningn = 250.0 


# set creep flag to 0 to disable all creep options.
enableCreep = 0

# if enableCreep = 1, you need to specify the physical time represented by the test time.
physicalStopTime = 360*24*3600.e6        # actual end of 1-year test time, micro-seconds.
timeScale = physicalStopTime / stopTime  # strain rate multiplier to match 1-year creep test in stopTime

# Porous compaction creep (Tuned to HYDcreep data)
creepA = 0.238    # creep: phi_e = A*exp(-p/B), note that A is 0-pressure equilbrium porosity and phi_i = 1 -exp(-p3)
creepB = 0.08     # creep: phi_e = A*exp(-p/B), B is GPa
creepC = timeScale*5.7e-09   # creep rate: dphi/dt = -p*C*(phi-phi_e)
creepD = 1.
creepE = 0.
creepF = 1.

# Deviatoric creep: (tuned to TXC creep data)
creepc0 = timeScale*2.0e-11  # dev creep rate term: creep strain sincrement = c0 * sigma_vm * std::pow( elasticShearStrain - equilibriumShearStrain, c1 ) * dt;
creepc1 = 1.0      # dev creep rate term: creep strain sincrement = c0 * sigma_vm * std::pow( elasticShearStrain - equilibriumShearStrain, c1 ) * dt;
creepc2 = 0

# Parameters for the geomechanics model damage evolution (tuned to TXCco data)
fractureEnergyReleaseRate =1.5e-8
fractureStress= 0.0184   #GPa
fractureSofteningExponent = 0.75   # controls how quickly the strength drops to resiscancel -u malendadual value as damage -> 1
damageEvolutionCriterion = 1       # user-specied brittle-ductile transition pressure
brittleDuctileTransition = 0.020


pfw["materials"] = [matdb.ghareb["name"]]
pfw["materialPropertyString"] = matdb.ghareb["materialString"]

# GEOMETRY OBJECTS -------------------------------------------------------
# single block filling domain for single-element test.

block = geom.box('block',[pfw["xmin"],pfw["ymin"],pfw["zmin"]],[pfw["xmax"],pfw["ymax"],pfw["zmax"]],vel=[0.0,0.0,0.0],mat=0,group=0)
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
<InitializeStress 
		time="0.0"
		interval="0.1"
		targetRegion="all"
		pressure=
""" + '"' + str(confiningPressure) + '"' + """
		/>
</MPMEvents>
"""