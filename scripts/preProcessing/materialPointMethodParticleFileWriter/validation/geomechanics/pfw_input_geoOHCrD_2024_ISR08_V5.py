# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree
import re
import pfw_materials as matdb
############################################################
#       this cript is associated with HP2/IS2 data      
############################################################

# Test of the prescibed F-table and boundary condition table capability.
# An elastic block is loaded in tensile uniaxial strain, using prescibed F loading.
# The lateral boundaries are then freed using outflow BCs so that subsequent loading is uniaxial stress.
# Material is then unloaded to zero strain.
# Reactions are plotted against the expected result for an elastic material
# Initial slope is constrained modulus, 2nd leg is young's modulus, 3rd leg unloads with young's modulus.
# The Mathematica notebook 'rotateParticles.nb' can be used to rotate particles in-place to see the effect of grid-alignment.


#V5 has new temp table capabilites that Mike put in

pfw = {}
pfw["runDebug"] = False
stopTime = 2000.  # simulated test time in microseconds
#stopTime = 1000.
#the stopTime on the order of 100us 
#strain versus time will be off

#this will need to be varied when we run the creep simulation and the creep + collapse cannot be separated


#The ISR08 test was 18 days long, or (90.78*24) 2187.898 hours
physicalStopTime = 18.*24*60*60*1000000# physical stop time in microsectonds (e.g. for a 30 h test)

timeScale = physicalStopTime/stopTime  # ratio of actual test timet o simulated time.
#timeScale = 0.018356

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
#***in the 2023 code, refine is set to 1***#
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
pfw["mWallTime"]="00:20:00"
#in the 2023 code, mWall time is 00:10:00
pfw["mSubmitJobs"]=False

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------
#this section updated on 10/15/24

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

pfw["useEvents"]=1
#ISR08 temperature was 20.3 Celsius
initialTemperature = 20.3+273.15
pfw["initialTemperature"]=initialTemperature

pfw["useAPIC"]=0
pfw["useInteralForceAsFaceReaction"]=1


# MATERIAL PROPERTIES --------------------------------------------------------------------
#this section updated on 10/15/24

density = 2.3
b0 = 2.003 #use the same values as HL01
b1 = 36.33 # high pressure bulk modulus ia b0+b1 (GPa); use the same values as HL01
b2 = 0.65  # shape parameter for pressure dependence; use the same values as HL01
b3 = 1.0 # elastic-plastic coupling, reduction in bulk modulus with plastic vol. strain; use the same values as HL01
b4 = 0.004 # elastic-plastic coupling shape parameter; use the same values as HL01
C1a = 0.
C2a = 0.
C3 = 0.
C4 = 0.
C5 = 0.
g1 = 0.1
g2=-0.0001# high-pressure poisson's ratio is nu = g1+g2 (should be 0<nu<0.5)
nu = g1 + g2
g0 = 1.5 * b0 * ( 1.0 - 2.0 * nu ) / ( 1.0 + nu )  # dilated shear modulus (g0>0), units of GPa
g3 = 0.3     #
g4 = 0.      # unused for now.
p0 = -0.057 # initial hydrostatic compressive strength, force it to be larger than confining pressure.  VAlue is I1=-3p, negative in compression; use the same values as HL01
p1 = 4.38  # Compaction rate parameter: GPa^-1 ; use the same values as HL01
p2 = 0.0   # unused
p3 = 0.31 # maximum achievable compressive volumetric plastic strain, relaties to initial porosity:  poro0 = 1-exp(-p3) ;use the same values as HL01


p4 = 0.0   # unused
BETA_nonassociativity = 1.4 # 0.5<beta<2. # dilation control
CR = 0.2
STREN = 0.0745 #using the TP1-TP3 values because these are also outcrop
FSLOPE = 0.3166#using the TP1-TP3 values because these are also outcrop
FSLOPEFAILED = FSLOPE-(0.05)
PEAKI1 = 0.0015 #using the TP1-TP3 values because these are also outcrop
YSLOPE = 0.00002               #yslope should be less than fslopefailed; #using the TP1-TP3 values because these are also outcrop
strainHardeningK = 0.015 #0.017# 0.0075 # (GPa), STREN = STREN_i + K*(1.-exp(-n*gamma_p))
strainHardeningn = 250.   #250.0 #350.0 # STREN = STREN_i + K*(1.-exp(-n*gamma_p))
enableCreep = 1
#timeScale = 1.    # ratio of physical time to simulation time.  


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ hydrostatic creep ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A = 0.24
B = 0.08
C = timeScale*7.0e-13
#D 
equilibriumPorosityPressureExponent= 1.13
#E
equilibriumPorosityOffset = 0.0
#F
compactionRatePressureExponent = -1.75
p3 = 0.26



creepA = A             # 0-pressure equilbrium porosity and phi_i = 1 -exp(-p3)
creepB = B            # pressure term #from HC01 best fit
#CreepC is meant to = creepC_old * (b_0)^creepF_new    creepC_new = (7e-12*tstop)*(2.003^-1.5) = 0.0301    = 1.5e-5 * tstop
creepC = C # creep rate multiplier: dphi/dt = - M * p/b0 *e^(C*(phi-phi_e))
creepD = equilibriumPorosityPressureExponent      # equilibrium porosity exponent
creepE = equilibriumPorosityOffset              # equilibriumPorosityOffset
#creepFnew is meant to be = (F_old*ln(p))/(ln(p/K)) = (-1.5 * ln (0.008))/(ln (0.008/2.003)) = -1.31
creepF = compactionRatePressureExponent          # compactionRatePressureExponent

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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
bucklingLength = 0.2*sampleLength       # physical length, if smaller than element size, multiple collapse will occur
bucklingAmplitude = 0.8     
stressControlKp = 0.01
stressControlKi = 0.0
stressControlKd = 0.001

initialTemperature = 20.3 + 273.15 #in celsius put to kelvin
pfw["initialTemperature"]=initialTemperature

Q = 2200. #activation energy in J/mol

pfw["materials"] = [matdb.ghareb["name"]]
pfw["materialPropertyString"]="""
<ElasticIsotropic
	name="aluminum"
	defaultDensity="2.7"
	defaultBulkModulus="70.0"
	defaultShearModulus="24.0"/>

<Geomechanics
   name="ghareb"
   defaultDensity="""+'"'+str(density)+'"'+"""
   b0="""+'"'+str(b0)+'"'+"""
   b1="""+'"'+str(b1)+'"'+"""
   b2="""+'"'+str(b2)+'"'+"""
   b3="""+'"'+str(b3)+'"'+"""
   b4="""+'"'+str(b4)+'"'+"""
   C1a="""+'"'+str(C1a)+'"'+"""
   C2a="""+'"'+str(C2a)+'"'+"""
   C3="""+'"'+str(C3)+'"'+"""
   C4="""+'"'+str(C4)+'"'+"""
   C5="""+'"'+str(C5)+'"'+"""
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
   peakI1="""+'"'+str(PEAKI1)+'"'+"""
   fSlope="""+'"'+str(FSLOPE)+'"'+"""
   stren="""+'"'+str(STREN)+'"'+"""
   ySlope="""+'"'+str(YSLOPE)+'"'+"""
   beta="""+'"'+str(BETA_nonassociativity)+'"'+"""
   enableCreep="""+'"'+str(enableCreep)+'"'+"""
   C1a="""+'"'+str(C1a)+'"'+"""
   C2a="""+'"'+str(C2a)+'"'+"""
   C3="""+'"'+str(C3)+'"'+"""
   C4="""+'"'+str(C4)+'"'+"""
   C5="""+'"'+str(C5)+'"'+"""
   creepC0="""+'"'+str(creepc0)+'"'+"""
   creepC1="""+'"'+str(creepc1)+'"'+"""
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
   initialTemperature="""+'"'+str(initialTemperature)+'"'+"""
   damageEvolutionCriterion="""+'"'+str(damageEvolutionCriterion)+'"'+"""
   brittleDuctileTransition="""+'"'+str(brittleDuctileTransition)+'"'+"""
   enableBuckling="""+'"'+str(enableBuckling)+'"'+"""
   bucklingLength="""+'"'+str(bucklingLength)+'"'+"""
   bucklingAmplitude="""+'"'+str(bucklingAmplitude)+'"'+"""
   />
"""


# DEFORMATION -----------------------------------------------------------------------------
#We don't use F tables for creep tests. Only stress tables
pfw["fTableInterpType"]='Smoothstep'
pfw["prescribedBoundaryFTable"]=1
#probably don't need the following line
pfw["fTable"]=[[0,1,1,1,],[stopTime,1,1,1]]

# stress table controls x- and y- directions.
# stress table controls x- and y- directions.
pfw["stressControl"]=[ 1, 1, 1]
pfw["stressTableInterpType"] = 'Smoothstep'
pfw["stressControlKp"] = stressControlKp
pfw["stressControlKi"] = stressControlKi
pfw["stressControlKd"] = stressControlKd

#stress table values can be found in: /Users/malenda1/Desktop/creep-compaction/organized by test type/OHCrW/IS2_CombinedTests_0to3000psi_MM_calculations.xlsx
#in the 'pressure vs. time stress table tab

pfw["stressTable"]=[[0.0,  0.0, 0.0, 0.0],
[(stopTime*0.000000000000000000e+00),-0.000000000000000000e+00,-0.000000000000000000e+00, -0.000000000000000000e+00],
[(stopTime*1.351351351351351426e-02),-8.004315999999999168e-03,-8.004315999999999168e-03, -8.004315999999999168e-03],
[(stopTime*9.864864864864865135e-01),-7.993767000000000580e-03,-7.993767000000000580e-03, -7.993767000000000580e-03],
[(stopTime*1.000000000000000000e+00),1.86157999999999997e-06,1.86157999999999997e-06, 1.86157999999999997e-06]]


# prescribed deformation (moving pistons) at all faces
pfw["boundaryConditionTypes"]=[2, 2, 2, 2, 2, 2]


# GEOMETRY OBJECTS -------------------------------------------------------

block = geom.box('block',[pfw["xmin"],pfw["ymin"],pfw["zmin"]],[pfw["xmax"],pfw["ymax"],pfw["zmax"]],vel=[0.0,0.0,0.0],mat=0,group=0)
pfw["objects"]=[block]

# GEOS MPM Events -----------------------------------------------------------------------



pfw["temperatureTable"]=[[(stopTime*0.00e+00), (2.03e+01)+273.15],
[(stopTime*1.51e-01), (2.03e+01)+273.15],
[(stopTime*1.64e-01), (4.08e+01)+273.15],
[(stopTime*3.79e-01), (4.10e+01)+273.15],
[(stopTime*3.92e-01), (6.37e+01)+273.15],
[(stopTime*5.44e-01), (6.40e+01)+273.15],
[(stopTime*5.56e-01), (8.07e+01)+273.15],
[(stopTime*7.72e-01), (8.05e+01)+273.15],
[(stopTime*7.84e-01), (1.03e+02)+273.15],
[(stopTime*9.36e-01), (1.00e+02)+273.15],
[(stopTime*9.49e-01), (4.65e+01)+273.15],
[(stopTime*9.62e-01), (2.63e+01)+273.15],
[(stopTime*9.87e-01), (2.10e+01)+273.15]]

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