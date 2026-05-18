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
stopTime = 1.0
# DOMAIN ---------------------------------------------------------------------------------

sampleX = 1.0  # mm
sampleY = 1.0 # mm
sampleZ = 1.0 # mm

domainX = 1.2*sampleX  # This would be increased for unconfined compression.
domainY = 1.2*sampleY
domainZ = sampleZ

pfw["xmin"] = 0.0             # mm
pfw["xmax"] = domainX    # mm
pfw["ymin"] = 0.0 # mm
pfw["ymax"] = domainY # mm
pfw["zmin"] = 0.0 # mm
pfw["zmax"] = domainZ # mm

refine=1  # partitions in each direction
cpp=4    # cells per partition in each direction

pfw["xpar"]=refine
pfw["ypar"]=refine
pfw["zpar"]=refine

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cpp  # grid cells in the z-direction
pfw["ppc"]=2               # particles per cell in each direction

# BATCH PARAMETERS  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"]="00:2:00"
pfw["mSubmitJobs"]=False

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------
pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/100
pfw["restartInterval"]=5*stopTime
pfw['lastRestartBufferInSeconds'] = 0.
pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25
pfw["initialDt"]=1e-16
pfw["updateMethod"]="PIC"

# MATERIAL PROPERTIES --------------------------------------------------------------------
# The Ghareb model has temperature as a wrapper, so we should see the temp
# table profile in both the particle and material model field in paraview

# Material model parameters.  Will be read from input file for plotting.
MPa = 0.001     # GPa
density = 1.57  # from excel sheet (only affects dynamics)

# nonlinear bulk modulus (tune to unload slopes in HYD compression - ISR01
b0 = 1.67
b1 = 30.0 # high pressure bulk modulus ia b0+b1 (GPa)
b2 = 0.3  # shape parameter for pressure dependence
b3 = 1.42 # elastic-plastic coupling, reduction in bulk modulus with plastic vol. strain
b4 = 0.015 # elastic-plastic coupling shape parameter.
g1 = 0.1
g2=-0.0001# high-pressure poisson's ratio is nu = g1+g2 (should be 0<nu<0.5)
nu = g1 + g2
g0 = 1.5 * b0 * ( 1.0 - 2.0 * nu ) / ( 1.0 + nu )  # dilated shear modulus (g0>0), units of GPa
g3 = 0.3     #
g4 = 0.      # unused for now.
p0 = -0.03 # initial hydrostatic compressive strength, force it to be larger than confining pressure.  VAlue is I1=-3p, negative in compression
p1 = 5.0  # Compaction rate parameter: GPa^-1
p2 = 0.0   # unused
p3 = 0.65 # maximum achievable compressive volumetric plastic strain, relaties to initial porosity:  poro0 = 1-exp(-p3)
p4 = 0.0   # unused
BETA_nonassociativity = 1.4 # 0.5<beta<2. # dilation control
CR = 0.2
STREN = 0.50
FSLOPE = 0.18
FSLOPEFAILED = FSLOPE-(0.05)
PEAKI1 = 0.0322/2.
YSLOPE = 0.002               #yslope should be less than fslopefailed
strainHardeningK = 0.017 #0.017# 0.0075 # (GPa), STREN = STREN_i + K*(1.-exp(-n*gamma_p))
strainHardeningn = 250.   #250.0 #350.0 # STREN = STREN_i + K*(1.-exp(-n*gamma_p))
enableCreep = 0
timeScale = 1.    # ratio of physical time to simulation time.  
creepA = 0.2205              # 0-pressure equilbrium porosity and phi_i = 1 -exp(-p3)
creepB = 0.000628            # pressure term
creepC = timeScale*0.308e-10  # creep rate multiplier: dphi/dt = -p*C*(phi-phi_e)
creepD = 2.307              # equilibrium porosity exponent
creepE = 0.                 # equilibriumPorosityOffset
creepF = 0.564              # compactionRatePressureExponent
creepc0 = 10.0               # equilibriumShearStrainConstant (1/GPa)
creepc1 = 1.229              # equilibriumShearStrainExponent (dimensionless)
creepc2 = timeScale*7.45e-08 # shearStrainRateConstant (1/microsecond)
fractureEnergyReleaseRate = 1.5e-8  # fracture energy release per unit area, normalizes damage evolution
fractureStress= 0.0184             # GPa
fractureSofteningExponent = 0.75   # controls how quickly the strength drops to resiscancel -u malendadual value as damage -> 1
damageEvolutionCriterion = 1       # user-specied brittle-ductile transition pressure
brittleDuctileTransition = 0.020
Kf = 0.0   # fluid bulk modulus (2.2 GPa) - set to 0
pfi = 0.0  # initial fluid pressure (GPa) - set to 0
t1RateDependence=0.0 # - set to 0
t2RateDependence=0.0 # - set to 0
enableBuckling = 1                   # 0: no buckling, 1: isotropic, 2: directional (requires material direction in particleFields)
bucklingLength = 0.2*sampleZ         # physical length, if smaller than element size, multiple collapse will occur
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

block = geom.box('block',[-0.5*sampleX,-0.5*sampleY,pfw["zmin"]],[0.5*sampleX,0.5*sampleY,pfw["zmax"]],vel=[0.0,0.0,0.0],mat=0,group=0)
pfw["objects"]=[block]

# DEFORMATION -----------------------------------------------------------------------------
# no deformation just a prescribed temp we can view in paraview.

# MPM EVENTS -------------------------------------------------------------------------------
pfw["temperatureTable"]=[
    [0,	       300.],
    [0.5*stopTime, 200.],
    [stopTime, 400.]
    ]
pfw["temperatureTableInterpType"] = 'Linear'

pfw["useEvents"]=1
pfw["mpmEventsString"]="""
<TemperatureProfile
        time=
""" + '"' + str(0.0) + '"' + """
        interval=""" + '"' + str(stopTime) + '"' + """
        />
</MPMEvents>
"""