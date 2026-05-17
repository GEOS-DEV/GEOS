# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

# This is currently just a smoke test to see if the geomechanics model is implemented
# successfully and runs using the buckling capability with non-monotonic cascading
# compaction.
#

pfw = {}
pfw["runDebug"] = False
stopTime = 100.0

maxCompressiveStrain = 0.40

# material direction:
matDir = np.array([1.,1., 0])
matDir = matDir / np.linalg.norm(matDir)

# DOMAIN ---------------------------------------------------------------------------------

sampleX = 0.8  # mm
sampleY = 1.0 # mm

domainX = sampleY  # This would be increased for unconfined compression.
domainY = sampleY

pfw["xmin"] = -0.5*domainX             # mm
pfw["xmax"] = 0.5*domainX    # mm
pfw["ymin"] = 0.0 # mm
pfw["ymax"] = domainY # mm

refine=1  # partitions in each direction
cpp=20    # cells per partition in each direction

pfw["xpar"]=refine
pfw["ypar"]=refine

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["ppc"]=2               # particles per cell in each direction


# Set domain Thickness to be 1 grid cell.
DX = (pfw["xmax"] - pfw["xmin"])/pfw["nI"]
pfw["planeStrain"] = 1
pfw["zmin"] = -0.5*DX # mm
pfw["zmax"] = 0.5*DX
pfw["zpar"] = 1
pfw["nK"] = 3  # grid cells in the z-direction

# These fields are needed
pfw["particleFileFields"] = ["Velocity",
                             "MaterialType",
                             "ContactGroup",
                             "StrengthScale",       # needed for weibull
                             "SurfaceFlag",         # needed for CZ and contact
                             "MaterialDirection",   # needed for graphite model
                             "RVector"             # needed for cpdi
]     # needed for CZ and enhanced contact.

# BATCH PARAMETERS  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"]="12:00:00"
pfw["mSubmitJobs"]=True

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/100
pfw["restartInterval"]=5*stopTime
pfw['lastRestartBufferInSeconds'] = 0.

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25
pfw["initialDt"]=1e-16
pfw["reactionHistory"]=1
pfw["reactionWriteInterval"] = stopTime/1000
pfw["boxAverageHistory"] = 1
pfw["boxAverageWriteInterval"] = stopTime/1000

pfw["solverProfiling"]=0         
pfw["frictionCoefficient"]=0.25  

pfw["updateMethod"]="XPIC"
pfw["updateOrder"]=2

# MATERIAL PROPERTIES --------------------------------------------------------------------

# Material model parameters.  Will be read from input file for plotting.
MPa = 0.001     # GPa
density = 1.57  # from excel sheet (only affects dynamics)

# nonlinear bulk modulus (tune to unload slopes in HYD compression - ISR01
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
p0 = -0.03 # initial hydrostatic compressive strength, force it to be larger than confining pressure.  VAlue is I1=-3p, negative in compression
p1 = 5.0  # Compaction rate parameter: GPa^-1
p2 = 0.0   # unused
p3 = 0.65 # maximum achievable compressive volumetric plastic strain, relaties to initial porosity:  poro0 = 1-exp(-p3)
p4 = 0.0   # unused

# Yield surface parameters, tuned to initial yield stress (inflection point) in TXC data
BETA_nonassociativity = 1.4 # 0.5<beta<2. # dilation control
CR = 0.2

STREN = 0.50
FSLOPE = 0.18
FSLOPEFAILED = FSLOPE-(0.05)
PEAKI1 = 0.0322/2.
YSLOPE = 0.002               #yslope should be less than fslopefailed

# If the input value for STREN above is treated as an initial value STREN_i
# the strength will be a function of deviatoric plastic strain:
# Note STREN goes to STREN_i + K, as n*gamma_p is large:
# set K=0 to disable hardening:
strainHardeningK = 0.017 #0.017# 0.0075 # (GPa), STREN = STREN_i + K*(1.-exp(-n*gamma_p))
strainHardeningn = 250.   #250.0 #350.0 # STREN = STREN_i + K*(1.-exp(-n*gamma_p))

# set creep flag to 0 to disable all creep options.
enableCreep = 0
timeScale = 1.    # ratio of physical time to simulation time.  

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
fractureSofteningExponent = 0.75   # controls how quickly the strength drops to resiscancel -u malendadual value as damage -> 1
damageEvolutionCriterion = 1       # user-specied brittle-ductile transition pressure
brittleDuctileTransition = 0.020

# fluid pressure (not used), set Kf=pfi=0 to disable the legacy fluid effects models, which probably need to be updated.
Kf = 0.0   # fluid bulk modulus (2.2 GPa) - set to 0
pfi = 0.0  # initial fluid pressure (GPa) - set to 0
# overstress model (not used).  This was for rate-dependent strengthening, but hasn't been tested since the model
# was ported from Uintah. 
t1RateDependence=0.0 # - set to 0
t2RateDependence=0.0 # - set to 0

# Buckling parameters to allow cascading compaction (non-monotonic)
enableBuckling = 2                   # 0: no buckling, 1: isotropic, 2: directional (requires material direction in particleFields)
bucklingLength = 0.2*sampleY         # physical length, if smaller than element size, multiple collapse will occur
bucklingAmplitude = 0.8              # reduction in compaction strength during collapse, 1: no reduction, 2: complete loss of strength


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
   enableBuckling="""+'"'+str(enableBuckling)+'"'+"""
   bucklingLength="""+'"'+str(bucklingLength)+'"'+"""
   bucklingAmplitude="""+'"'+str(bucklingAmplitude)+'"'+"""
   />
"""

# GEOMETRY OBJECTS -------------------------------------------------------

block = geom.box('block',[-0.5*sampleX,pfw["ymin"],pfw["zmin"]],[0.5*sampleX,pfw["ymax"],pfw["zmax"]],vel=[0.0,0.0,0.0],mat=0,group=0)

blockWDir = geom.materialDirectionWrapper(name='blockWDir',subObject=block,matDir=matDir)

pfw["objects"]=[blockWDir]

# DEFORMATION -----------------------------------------------------------------------------

# We initialize the domain with p=confiningPressure, if p > -p0/3 this will fail.
# we then apply a stress boundary condition on x,y faces to be in equilibrium, and do a 
# strain-controlled motion in z.

# Ftable only controls z-direction
pfw["boundaryConditionTypes"]=[ 0, 0, 2, 2, 1, 1 ]
pfw["fTableInterpType"]='Cosine'
pfw["prescribedBoundaryFTable"] = 1
pfw["fTable"]=[
    [0,	       1.0,	1.0, 1.0],
    [stopTime, 1.0, np.exp(-maxCompressiveStrain), 1.0 ]
    ]
