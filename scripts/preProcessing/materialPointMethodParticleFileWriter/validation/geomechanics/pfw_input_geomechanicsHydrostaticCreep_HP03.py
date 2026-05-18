# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree
import re
import pfw_materials as matdb

# Hydrostatic compression test for comparison to HPO3 creep test on Ghareb chalk.

pfw = {}
pfw["runDebug"] = True

# Prescribed deformation and BC tables:
stopTime = 1000.  #simulation end time
physicalStopTime= 90.78*24*60*60*1000000 # for micro seconds
timeScale = physicalStopTime/stopTime  # this will be used to scale creep-rate term.

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
pfw["mSubmitJobs"]=False

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

pfw["useAPIC"]=0
pfw["useInteralForceAsFaceReaction"]=1

# MATERIAL PROPERTIES --------------------------------------------------------------------

# Material model parameters.  Will be read from input file for plotting.
MPa = 0.001 # GPa
density = 2.3

# nonlinear bulk modulus (tune to unload slopes in HYD compression - HL01
b0 = 2.003  # initial bulk modulus.  this had been reduced by 2/3 when matching TXC creep data.
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
g1 = 0.1  #0.3       # low-pressure poisson's ratio is nu = g1+g2 (should be 0<nu<0.5)
g2 = -0.0001 #-0.15     # high-pressure poisson's ratio is nu = g1
nu = g1 + g2
g0 = 1.5 * b0 * ( 1.0 - 2.0 * nu ) / ( 1.0 + nu )  # dilated shear modulus (g0>0), units of GPa
g3 = 0.3     #
g4 = 0.      # unused for now.

# crush curve
p0 = -0.057 # initial hydrostatic compressive strength, force it to be larger than confining pressure.  VAlue is I1=-3p, negative in compression
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
pfw["materials"] = [matdb.ghareb["name"]]
pfw["materialPropertyString"] = matdb.ghareb["materialString"]


# # MATERIAL PROPERTIES --------------------------------------------------------------------

# # Will be read from input file for plotting.
# density = 2.648

# # Material model parameters.  Will be read from input file for plotting.
# MPa = 0.001 # GPa
# density = 2.3

# # nonlinear bulk modulus
# b0 = 2.003 
# b1 = 36.33
# b2 = 650.0*MPa
# b3 = 1.0
# b4 = 4.0e-3 

# # nonlinear shear modulus
# g0 = 0.40128  # (0.924<shear modulus<1.502)
# g1 = 0.088
# g2 = 0.0001
# g3 = 0.0
# g4 = 0.

# # crush curve
# BETA_nonassociativity = 1.16 # 0.5<beta<2.
# p0 = -0.059  #-0.059
# p1 = 4.38  # GPa^-1
# p2 = 0.0
# p3 = 0.255
# p4 = 0.0

# CR = 0.5  # ZZZZ 0.2<CR<0.8  ####For TP1, kep CR close to 0.2-.4. 0.8 is too high (8/9)
# Kf = 0.0   # fluid bulk modulus (2.2 GPa)
# pfi = 0.0  # initial fluid pressure (GPa)

# # yield surface
# #These following values are from Figure 5 of Kikibas et al
# PEAKI1 = 0.0161
# FSLOPE = 0.613
# STREN = 0.0135
# YSLOPE = 0.0001
# BETA_nonassociativity = 1.16

# # fluid pressure (not used)
# Kf = 0.0   # fluid bulk modulus (2.2 GPa)
# pfi = 0.0  # initial fluid pressure (GPa)

# # overstress model (not used)
# t1RateDependence=0.0
# t2RateDependence=0.0

# # If the input value for STREN above is treated as an initial value STREN_i
# # the strength will be a function of deviatoric plastic strain:
# # Note STREN goes to STREN_i + K, as n*gamma_p is large:

# # set K=0 to disable hardening:
# strainHardeningK = (0.)*0.0075 # (GPa), STREN = STREN_i + K*(1.-exp(-n*gamma_p))
# strainHardeningn = 350.0 # STREN = STREN_i + K*(1.-exp(-n*gamma_p))

# # deviatoric creep:
# creepc0 = timeScale*2.0e-11  # dev creep rate term: creep strain sincrement = c0 * sigma_vm * std::pow( elasticShearStrain - equilibriumShearStrain, c1 ) * dt;
# creepc1 = 1.0      # dev creep rate term: creep strain sincrement = c0 * sigma_vm * std::pow( elasticShearStrain - equilibriumShearStrain, c1 ) * dt;

# # bulk and shear are used for initial time-step calculations and 
# # some normalizations and are defined from the geo input parameters
# bulk = b0+b1
# shear = g0

# # These parameters are for the fracture propagation.  We don't have any data
# # to constrain them so these are best guesses.
# KIc = 1.591*(.001)*np.sqrt(1000.)       # fracture toughness GPa mm^1/2
# fractureRadius = 0.5*(KIc**2)/(np.pi*PEAKI1**2)  # radius of fracture process zone
# constrainedModulus = bulk + 4./3.*shear
# Gf = KIc/constrainedModulus

# pfw["materials"] = [matdb.ghareb["name"]]
# pfw["materialPropertyString"]="""
# <Geomechanics
#    name="ghareb"
#    defaultDensity="""+'"'+str(density)+'"'+"""
#    b0="""+'"'+str(b0)+'"'+"""
#    b1="""+'"'+str(b1)+'"'+"""
#    b2="""+'"'+str(b2)+'"'+"""
#    b3="""+'"'+str(b3)+'"'+"""
#    b4="""+'"'+str(b4)+'"'+"""
#    g0="""+'"'+str(g0)+'"'+"""
#    g1="""+'"'+str(g1)+'"'+"""
#    g2="""+'"'+str(g2)+'"'+"""
#    g3="""+'"'+str(g3)+'"'+"""
#    g4="""+'"'+str(g4)+'"'+"""
#    p0="""+'"'+str(p0)+'"'+"""
#    p1="""+'"'+str(p1)+'"'+"""
#    p2="""+'"'+str(p2)+'"'+"""
#    p3="""+'"'+str(p3)+'"'+"""
#    p4="""+'"'+str(p4)+'"'+"""
#    cr="""+'"'+str(CR)+'"'+"""
#    fluidBulkModulus="""+'"'+str(Kf)+'"'+"""
#    fluidInitialPressure="""+'"'+str(pfi)+'"'+"""
#    t1RateDependence="""+'"'+str(t1RateDependence)+'"'+"""
#    t2RateDependence="""+'"'+str(t2RateDependence)+'"'+"""
#    peakI1="""+'"'+str(PEAKI1)+'"'+"""
#    fSlope="""+'"'+str(FSLOPE)+'"'+"""
#    stren="""+'"'+str(STREN)+'"'+"""
#    ySlope="""+'"'+str(YSLOPE)+'"'+"""
#    beta="""+'"'+str(BETA_nonassociativity)+'"'+"""
#    enableCreep="""+'"'+str(enableCreep)+'"'+"""
#    creepC0="""+'"'+str(creepc0)+'"'+"""
#    creepC1="""+'"'+str(creepc1)+'"'+"""
#    creepA="""+'"'+str(A)+'"'+"""
#    creepB="""+'"'+str(B)+'"'+"""
#    creepC="""+'"'+str(C)+'"'+"""
#    creepD="""+'"'+str(equilibriumPorosityPressureExponent)+'"'+"""
#    creepE="""+'"'+str(0.)+'"'+"""
#    creepF="""+'"'+str(compactionRateExponent)+'"'+"""
#    strainHardeningN="""+'"'+str(strainHardeningn)+'"'+"""
#    strainHardeningK="""+'"'+str(strainHardeningK)+'"'+"""
#    />
# """

   # plasticStrainTolerance="1e-10"
   # stressReturnTolerance="1e-6"
   # maxAllowedSubcycles="256"
   # failedStepResponse="3"

# GEOMETRY OBJECTS -------------------------------------------------------

block = geom.box('block',[pfw["xmin"],pfw["ymin"],pfw["zmin"]],[pfw["xmax"],pfw["ymax"],pfw["zmax"]],vel=[0.0,0.0,0.0],mat=0,group=0)
pfw["objects"]=[block]

# DEFORMATION -----------------------------------------------------------------------------

# Ftable only controls z-direction
pfw["boundaryConditionTypes"]=[ 2, 2, 2, 2, 2, 2 ]
pfw["fTableInterpType"]='Cosine'
pfw["prescribedBoundaryFTable"] = 1
pfw["fTable"] = [[0,1,1,1],[stopTime,1,1,1]]

#stress table controls x- and y- directions.
pfw["stressControl"]=[ 1, 1, 1]
pfw["stressTableInterpType"] = 'Cosine'
pfw["stressControlKp"] = 0.05
pfw["stressControlKi"] = 0.0
pfw["stressControlKd"] = 0.0

pfw["stressTable"]=[[0.0,  0.0, 0.0, 0.0],
[(stopTime*0.000000000000000000e+00),8.9631840000000004100e-06,8.9631840000000004100e-06,8.9631840000000004100e-06],
[(stopTime*1.010101010101010187e-02),-2.078079760000000036e-03,-2.078079760000000036e-03,-2.078079760000000036e-03],
[(stopTime*1.212121212121212155e-01),-2.078079760000000036e-03,-2.078079760000000036e-03,-2.078079760000000036e-03],
[(stopTime*1.313131313131313260e-01),-4.174085887999999175e-03,-4.174085887999999175e-03,-4.174085887999999175e-03],
[(stopTime*2.121212121212121271e-01),-4.154091092999999575e-03,-4.154091092999999575e-03,-4.154091092999999575e-03],
[(stopTime*2.222222222222222376e-01),-6.946467678000000333e-03,-6.946467678000000333e-03,-6.946467678000000333e-03],
[(stopTime*3.030303030303030387e-01),-6.954741386000000587e-03,-6.954741386000000587e-03,-6.954741386000000587e-03],
[(stopTime*3.131313131313131493e-01),-1.037591980999999844e-02,-1.037591980999999844e-02,-1.037591980999999844e-02],
[(stopTime*3.939393939393939781e-01),-1.038074613999999916e-02,-1.038074613999999916e-02,-1.038074613999999916e-02],
[(stopTime*4.040404040404040886e-01),-1.380950880000000013e-02,-1.380950880000000013e-02,-1.380950880000000013e-02],
[(stopTime*4.444444444444444753e-01),-1.378399819000000054e-02,-1.378399819000000054e-02,-1.378399819000000054e-02],
[(stopTime*4.545454545454545858e-01),-2.076838703999999911e-02,-2.076838703999999911e-02,-2.076838703999999911e-02],
[(stopTime*6.363636363636364646e-01),-2.073184481999999773e-02,-2.073184481999999773e-02,-2.073184481999999773e-02],
[(stopTime*6.464646464646465196e-01),-1.723275565000000065e-02,-1.723275565000000065e-02,-1.723275565000000065e-02],
[(stopTime*7.070707070707071829e-01),-1.728308737000000192e-02,-1.728308737000000192e-02,-1.728308737000000192e-02],
[(stopTime*7.171717171717172379e-01),-1.380399298999999969e-02,-1.380399298999999969e-02,-1.380399298999999969e-02],
[(stopTime*7.676767676767677351e-01),-1.380812983999999896e-02,-1.380812983999999896e-02,-1.380812983999999896e-02],
[(stopTime*7.777777777777777901e-01),-6.920957076999999409e-03,-6.920957076999999409e-03,-6.920957076999999409e-03],
[(stopTime*8.686868686868687295e-01),-6.903030707999999877e-03,-6.903030707999999877e-03,-6.903030707999999877e-03],
[(stopTime*8.787878787878788955e-01),-3.443241646000000118e-03,-3.443241646000000118e-03,-3.443241646000000118e-03],
[(stopTime*9.494949494949496138e-01),-3.476336478999999830e-03,-3.476336478999999830e-03,-3.476336478999999830e-03],
[(stopTime*9.595959595959596689e-01),-2.408338620000000010e-03,-2.408338620000000010e-03,-2.408338620000000010e-03],
[(stopTime*9.797979797979798899e-01),-2.422128133999999853e-03,-2.422128133999999853e-03,-2.422128133999999853e-03],
[(stopTime*9.898989898989899450e-01),-1.367230313000000145e-03,-1.367230313000000145e-03,-1.367230313000000145e-03],
[(stopTime*1.0),-1.381019826999999988e-03,-1.381019826999999988e-03,-1.381019826999999988e-03]]
