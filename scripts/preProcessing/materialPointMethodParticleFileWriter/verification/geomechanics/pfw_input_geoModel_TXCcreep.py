# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

# This is currently just a smoke test to see if the geomechanics model is implemented
# successfully and runs. 
#
# TODO: add some actual verification, so make sure the response is correct.


pfw = {}
pfw["runDebug"] = True
stopTime = 1000.0

physicalStopTime = 360*24*3600.e6 # actual end of 1-year test time, micro-seconds.
#physicalStopTime = 24*3600.e6 # actual end of 1-day test time, micro-seconds.
timeScale = physicalStopTime / stopTime  # strain rate multiplier to match 1-year creep test in stopTime

confiningPressure = 0.005 # GPa
stressDifference = 0.010 # GPa


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
pfw["ppc"]=1               # particles per cell in each direction

# BATCH PARAMETERS  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"]="00:05:00"
pfw["mSubmitJobs"]=True

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/100
pfw["restartInterval"]=stopTime
pfw['lastRestartBufferInSeconds'] = 0.

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.005  # this needs to be really small, or the integraion error in m_domainL to get m_domainF is large since there is noisy in stress controlled m_domainL
pfw["initialDt"]=1e-6
pfw["reactionHistory"]=1
pfw["reactionWriteInterval"]=stopTime/1000
pfw["boxAverageHistory"]=1
pfw["boxAverageWriteInterval"]=stopTime/1000

pfw["solverProfiling"]=0         
pfw["frictionCoefficient"]=0.25  

pfw["updateMethod"]="PIC"
pfw["updateOrder"]=2

pfw["useAPIC"]=1
pfw["useInteralForceAsFaceReaction"]=1

# MATERIAL PROPERTIES --------------------------------------------------------------------

# Material model parameters.  Will be read from input file for plotting.
MPa = 0.001     # GPa
density = 1.57  # from excel sheet (only affects dynamics)

# nonlinear bulk modulus (tune to unload slopes in HYD compression)
b0 = 2.8
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
# g1 = 0.1
# g2=-0.0001# high-pressure poisson's ratio is nu = g1+g2 (should be 0<nu<0.5)
# nu = g1 + g2
# g0 = 1.5 * b0 * ( 1.0 - 2.0 * nu ) / ( 1.0 + nu )  # dilated shear modulus (g0>0), units of GPa
# g3 = 0.3     #
# g4 = 0.      # unused for now.

g0=1.0
g1=0
g2=0
g3=0
g4=0

# crush curve
p0 = min(-.059,-3.0001*confiningPressure) # initial hydrostatic compressive strength, force it to be larger than confining pressure.  VAlue is I1=-3p, negative in compression
p1 = 0.65  # Compaction rate parameter: GPa^-1
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
enableCreep = 1

# Porous compaction creep (Tuned to HYDcreep data)
# real64 A = m_creepA;  // volumetric creep rate parameter
# real64 B = m_creepB;  // volumetric creep rate parameter
# real64 C = m_creepC;  // volumetric creep rate parameter
# real64 equilibriumPorosityPressureExponent = m_creepD;  // volumetric creep rate parameter
# real64 equilibriumPorosityOffset = m_creepE;  // volumetric creep rate parameter
# real64 compactionRatePressureExponent = m_creepF;  // volumetric creep rate parameter
# phi_e = std::max(1.e-10 , A * exp( -std::pow(p,equilibriumPorosityPressureExponent) / B ) + equilibriumPorosityOffset );    
# dphidt = -1.0*std::pow(p,compactionRatePressureExponent)*C*( phi_p - phi_e );
creepA = 0.2205              # 0-pressure equilbrium porosity and phi_i = 1 -exp(-p3)
creepB = 0.000628            # pressure term
creepC = timeScale*0.308e-10  # creep rate multiplier: dphi/dt = -p*C*(phi-phi_e)
creepD = 2.307              # equilibrium porosity exponent
creepE = 0.                 # equilibriumPorosityOffset
creepF = 0.564              # compactionRatePressureExponent

# Deviatoric creep: (tuned to TXC creep data)
#   equilibriumShearStrainConstant = c0;
#   equilibriumShearStrainExponent = c1;
#   shearStrainRateConstant = c2;
#   real64 equilibriumVMShearStrain = equilibriumShearStrainConstant * std :: pow(vonMisesStress_old , equilibriumShearStrainExponent);
#   real64 creepVMStrainIncrement =  Dt * shearStrainRateConstant * std::max(equilibriumVMShearStrain - plasticVMshearStrain , 0.0);
creepc0 = 10.0               # equilibriumShearStrainConstant (1/GPa)
creepc1 = 1.229              # equilibriumShearStrainExponent (dimensionless)
creepc2 = timeScale*7.45e-08 # shearStrainRateConstant (1/microsecond)

# Parameters for the geomechanics model damage evolution (tuned to TXCco data)
fractureEnergyReleaseRate = 1.5e-8  # fracture energy release per unit area, normalizes damage evolution
fractureStress= 0.0184             # GPa
fractureSofteningExponent = 0.75   # controls how quickly the strength drops to residual value as damage -> 1
damageEvolutionCriterion = 1       # user-specied brittle-ductile transition pressure
brittleDuctileTransition = 0.020   # pressure (GPa) below which damage accumulates with damageEvolutionCriterion


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
   brittleDuctileTransition="""+'"'+str(brittleDuctileTransition)+'"'+"""
   />
"""

# GEOMETRY OBJECTS -------------------------------------------------------

block = geom.box('block',[pfw["xmin"],pfw["ymin"],pfw["zmin"]],[pfw["xmax"],pfw["ymax"],pfw["zmax"]],v=[0.0,0.0,0.0],mat=0,group=0)
pfw["objects"]=[block]

# DEFORMATION -----------------------------------------------------------------------------

# Ftable only controls z-direction
pfw["boundaryConditionTypes"]=[ 2, 2, 2, 2, 2, 2 ]
pfw["fTableInterpType"]='Cosine'
pfw["prescribedBoundaryFTable"] = 1

pfw["fTable"]=[
    [0,	          1., 1., 1.],
    [0.05*stopTime,	  0.995, 0.995, 0.995],
    [0.25*stopTime,	  1.0025, 1.0025, 0.98],
    [stopTime,	  1.0025, 1.0025, 0.98]
]

# # stress table controls x- and y- directions.
# pfw["stressControl"]=[ 1, 1, 0]
# pfw["stressTableInterpType"] = 'Cosine'
# pfw["stressControlKp"] = 0.05
# pfw["stressControlKi"] = 0.0
# pfw["stressControlKd"] = 0.0
# pfw["stressTable"]=[[0.,      	   -confiningPressure, -confiningPressure, -confiningPressure],
# 					[stopTime,     -confiningPressure, -confiningPressure, -confiningPressure]
#      ]

# # This should initialize stress in equilibrium with boundary conditions:
# pfw["useEvents"] = 1
# pfw["mpmEventsString"]="""
# <MPMEvents>
#     <InitializeStress 
#         time="0.0"
#         interval="0.1"
#         targetRegion="all"
#         pressure=""" + '"' + str(confiningPressure) + '"' + """
#         />
# </MPMEvents>
# """