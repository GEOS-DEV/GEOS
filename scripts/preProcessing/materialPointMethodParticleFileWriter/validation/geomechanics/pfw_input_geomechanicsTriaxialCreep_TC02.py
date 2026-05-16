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


#V3: Mike added damage into the code. 
#V$: Got better values for CR and related yield envelope variables
#V6; has variables for the fracture damage updated. like kappa_i, damageEvolutionCriteron(s)

pfw = {}
pfw["runDebug"] = False
stopTime = 1000.0
#physicalStopTime = 360*24*3600.e6        # actual end of 1-year test time, micro-seconds.
physicalStopTime = 7851*1.e6     # the TP02 test lasted 7851 seconds,  converting this to microseconds.

confiningPressure = 1.070101939000000016e-02 # confining pressure GPa
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
pfw["mWallTime"]="00:08:00"
pfw["mSubmitJobs"]=True

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/100
pfw["restartInterval"]=5*stopTime
pfw['lastRestartBufferInSeconds'] = 0.

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.005
pfw["initialDt"]=1e-16
pfw["reactionHistory"]=1
pfw["reactionWriteInterval"] = stopTime/1000
pfw["boxAverageHistory"] = 1
pfw["boxAverageWriteInterval"] = stopTime/1000

pfw["solverProfiling"]=0         
pfw["frictionCoefficient"]=0.25  

pfw["updateMethod"]="PIC" # best for single element.
pfw["updateOrder"]=2      # This only needs to be specified with XPIC or FMPM.
pfw["useAPIC"]=1          # This option should only be used for single-element simulations, and will improve accuracy of reactions to be more consistent with box data.

#pfw["singleElementDriver"] = 0
#pfw["useConstantTimeStep"] = 1
#pfw["constantTimeStepValue"] = stopTime/10000

pfw["useEvents"] = 0
# use this if you are not ramping up the hydrostatic stress.

# MATERIAL PROPERTIES --------------------------------------------------------------------

# Material model parameters.  Will be read from input file for plotting.
MPa = 0.001 # GPa
density = 2.3

# nonlinear bulk modulus (tune to unload slopes in HYD compression - HL01
b0 = (0.667)*2.003
b1 = 36.33  # high pressure bulk modulus ia b0+b1 (GPa)
b2 = 0.65   # shape parameter for pressure dependence
b3 = 1.00   # elastic-plastic coupling, reduction in bulk modulus with plastic vol. strain
b4 = 0.004  # elastic-plastic coupling shape parameter.

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
#
# g1 = 0.1  #0.3       # low-pressure poisson's ratio is nu = g1+g2 (should be 0<nu<0.5)
# g2 = -0.0001 #-0.15     # high-pressure poisson's ratio is nu = g1
# nu = g1 + g2
# g0 = 1.5 * b0 * ( 1.0 - 2.0 * nu ) / ( 1.0 + nu )  # dilated shear modulus (g0>0), units of GPa
# g3 = 0.3     #
# g4 = 0.      # unused for now.
g0=1.2
g1=0
g2=0
g3=0
g4=0

# crush curve
p0 = min(-.057,-3.0001*confiningPressure) # initial hydrostatic compressive strength, force it to be larger than confining pressure.  VAlue is I1=-3p, negative in compression
p1 = 4.38  # Compaction rate parameter: GPa^-1
p2 = 0.0   # unused
p3 = 0.255 # maximum achievable compressive volumetric plastic strain, relaties to initial porosity:  poro0 = 1-exp(-p3)
p4 = 0.0   # unused


BETA_nonassociativity = 1.6 #1.16 # 0.5<beta<2.
CR = 0.233   # 0.233 #0.9   # Cap branch point parameter:  ZZZZ 0.2<CR<0.8  ####For TP1, kep CR close to 0.2-.4. 0.8 is too high (8/9)

# yield surface
PEAKI1 = 0.0015       #0.0161 #0.008
FSLOPE = 0.3166       #0.18 #0.613
FSLOPEFAILED = 0.2666 #0.13 #0.305
STREN = 0.0745        #0.50 #0.0135
YSLOPE = 0.0002       #0.002 #0.0001

# If the input value for STREN above is treated as an initial value STREN_i
# the strength will be a function of deviatoric plastic strain:
# Note STREN goes to STREN_i + K, as n*gamma_p is large:
# set K=0 to disable hardening:
strainHardeningK = 0.015 #0.017# 0.0075 # (GPa), STREN = STREN_i + K*(1.-exp(-n*gamma_p))
strainHardeningn = 250.   #250.0 #350.0 # STREN = STREN_i + K*(1.-exp(-n*gamma_p))

# set creep flag to 0 to disable all creep options.
# if enableCreep = 1, you need to specify the physical time represented by the test time.
enableCreep = 1
timeScale = physicalStopTime / stopTime  # strain rate multiplier to match 1-year creep test in stopTime

# Porous compaction creep:
creepA = 0.2205    # creep: phi_e = A*exp(-p/B), note that A is 0-pressure equilbrium porosity and phi_i = 1 -exp(-p3) = 0.225
creepB = 0.000628     # creep: phi_e = A*exp(-p/B), B is GPa
creepC = timeScale*3.08e-11   # creep rate: dphi/dt = -p*C*(phi-phi_e)
creepD = 2.307  # equilibriumPorosityPressureExponent
creepE = 0
creepF = 0.564  # compactionRateExponent

# Deviatoric creep: (tuned to TXC creep data)
#   equilibriumShearStrainConstant = c0;
#   equilibriumShearStrainExponent = c1;
#   shearStrainRateConstant = c2;
#   real64 equilibriumVMShearStrain = equilibriumShearStrainConstant * std :: pow(vonMisesStress_old , equilibriumShearStrainExponent);
#   real64 creepVMStrainIncrement =  Dt * shearStrainRateConstant * std::max(equilibriumVMShearStrain - plasticVMshearStrain , 0.0);
creepc0 = 10.0               # equilibriumShearStrainConstant (1/GPa)
creepc1 = 1.229              # equilibriumShearStrainExponent (dimensionless)
creepc2 = timeScale*5.0e-11 # shearStrainRateConstant (1/microsecond), 7.45e-8 was from the mathematica notebook fit.

# Parameters for the geomechanics model damage evolution.
fractureEnergyReleaseRate =1.5e-8
fractureStress= 0.015   #GPa
fractureSofteningExponent = 0.75   # controls how quickly the strength drops to resiscancel -u malendadual value as damage -> 1
damageEvolutionCriterion = 1       # user-specied brittle-ductile transition pressure
brittleDuctileTransition = 0.020

# fluid pressure (not used), set Kf=pfi=0 to disable the legacy fluid effects models, which probably need to be updated.
Kf = 0.0   # fluid bulk modulus (2.2 GPa)
pfi = 0.0  # initial fluid pressure (GPa)

# overstress model (not used).  This was for rate-dependent strengthening, but hasn't been tested since the model
# was ported from Uintah.
t1RateDependence=0.0
t2RateDependence=0.0
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
   fSlopeFailed="""+'"'+str(FSLOPEFAILED)+'"'+"""
   stren="""+'"'+str(STREN)+'"'+"""
   ySlope="""+'"'+str(YSLOPE)+'"'+"""
   beta="""+'"'+str(BETA_nonassociativity)+'"'+"""
   enableCreep="""+'"'+str(enableCreep)+'"'+"""
   creepC0="""+'"'+str(creepc0)+'"'+"""
   creepC1="""+'"'+str(creepc1)+'"'+"""
   creepC2="""+'"'+str(creepc2)+'"'+"""
   creepA="""+'"'+str(creepA)+'"'+"""
   creepB="""+'"'+str(creepB)+'"'+"""
   creepC="""+'"'+str(creepC)+'"'+"""
   creepD="""+'"'+str(creepD)+'"'+"""
   creepE="""+'"'+str(creepE)+'"'+"""
   creepF="""+'"'+str(creepF)+'"'+"""
   strainHardeningN="""+'"'+str(strainHardeningn)+'"'+"""
   strainHardeningK="""+'"'+str(strainHardeningK)+'"'+"""
   fractureEnergyReleaseRate="""+'"'+str(fractureEnergyReleaseRate)+'"'+"""
   fractureSofteningExponent="""+'"'+str(fractureSofteningExponent)+'"'+"""
   fractureStress="""+'"'+str(fractureStress)+'"'+"""
   damageEvolutionCriterion="""+'"'+str(damageEvolutionCriterion)+'"'+"""
   />
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
pfw["prescribedBoundaryFTable"] = 0
pfw["fTable"]=[
	[0,	          1.000,	    1.,	1.],
   [stopTime,	          1.000,	    1.,	1.]
	]

# stress table controls x- and y- directions.
pfw["stressControl"]=[ 1, 1, 1]
pfw["stressTableInterpType"] = 'Cosine'
pfw["stressControlKp"] = kp
pfw["stressControlKi"] = ki
pfw["stressControlKd"] = kd
pfw["stressTable"]=[[2.487418596875372521e-05*stopTime, 0, 0, 0],
[1.012563303257513030e-02*stopTime, -2.108803608999999975e-03, -2.108803608999999975e-03, -2.108803608999999975e-03],
[2.022639187918150866e-02*stopTime, -4.593658936999999755e-03, -4.593658936999999755e-03, -4.593658936999999755e-03],
[3.032715072578788529e-02*stopTime, -7.453454796999999833e-03, -7.453454796999999833e-03, -7.453454796999999833e-03],
[5.052866841900063855e-02*stopTime, -7.630466507000000821e-03, -7.630466507000000821e-03, -7.630466507000000821e-03],
[6.062942726560701517e-02*stopTime, -confiningPressure, - confiningPressure, -1.070101939000000016e-02],
[6.062942726560701517e-02*stopTime, -confiningPressure, - confiningPressure, -1.070101939000000016e-02],
[7.073018611221337792e-02*stopTime, -confiningPressure, - confiningPressure, -1.509734128000000085e-02],
[9.093170380542614506e-02*stopTime, -confiningPressure, - confiningPressure, -1.547444327999999980e-02],
[1.010324626520325147e-01*stopTime, -confiningPressure, - confiningPressure, -1.795881333000000010e-02],
[1.010324626520325147e-01*stopTime, -confiningPressure, - confiningPressure, -1.795881333000000010e-02],
[1.111332214986388844e-01*stopTime, -confiningPressure, - confiningPressure, -2.113893686000000230e-02],
[1.313347391918516516e-01*stopTime, -confiningPressure, - confiningPressure, -2.104603639999999901e-02],
[1.5e-01*stopTime, -confiningPressure, - confiningPressure, -2.68e-02],
[1.5e-01*stopTime, -confiningPressure, - confiningPressure, -2.68e-02],
[stopTime, -confiningPressure, - confiningPressure, -2.615328664999999955e-02]]

# This should initialize stress in equilibrium with boundary conditions:
# pfw["mpmEventsString"]="""
# <MPMEvents>
# 	<InitializeStress 
# 		time="0.0"
# 		interval="0.1"
# 		targetRegion="all"
# 		pressure=""" + '"' + str(confiningPressure) + '"' + """
# 		/>
# </MPMEvents>
# """