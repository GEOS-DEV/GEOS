# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree
import math

pfw = {} 
pfw['particleFileFields'] = ["Velocity",  # +3
                            "MaterialType", # +1
                            "ContactGroup", # +1
                            "SurfaceFlag", # +1
                            "RVector",
                            "StrengthScale",
                            "SurfaceNormal",
                            "SurfacePosition",]

# Geometry
pdcDiameter = 10 # mm
pdcBackRakeAngle = (15*3.14159/180.) # rad
pdcLength = 10 # mm
pdcDepth = 2 # mm
pdcGap = 1 # mm gap between domain boundary.
pdcVelocity = 0.01 # mm/us = km/s

sampleX = 40 # mm, cutting translation direction
sampleY = 20 # mm, sample surface normal direction
sampleZ = 10 # mm, half depth, with sym at z=0

# sample will cut domain width 2x, with 2x depth of cut per pass.
stopTime = 2.0*sampleX/pdcVelocity

# MATERIAL PROPERTIES --------------------------------------------------------------------

density = 2.648
bulk = 36.3
shear = 26.0
tensileStrength = 0.449
compressiveStrength = 2.27
maximumStrength = 5.0
crackSpeed = 1.8
damagedMaterialFrictionSlope = 0.5773502691896258

# Domain ---------------------------------------------------------------------------------

refine=2
cpp=12
pfw["xpar"]=4*refine  # grid partitions
pfw["ypar"]=3*refine
pfw["zpar"]=1*refine

pfw["nI"]=pfw["xpar"]*cpp   	# grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp 	# grid cells in the y-direction
pfw["nK"]=pfw["zpar"]*cpp  			# grid cells in the z-direction
pfw["ppc"]=2   		# particles per cell in each direction

domainX = sampleX
domainY = sampleY + pdcDiameter*np.cos(pdcBackRakeAngle)
domainZ = sampleZ
                
# Define all the geometric objects -------------------------------------------------------
pfw["xmin"] = 0.0
pfw["xmax"] = domainX
pfw["ymin"] = 0.0	# m
pfw["ymax"] = pfw["ymin"] + domainY 
pfw["zmin"] = 0.0
pfw["zmax"] = domainZ

# Batch parameters for GEOS runs.  --------------------------------------------------------
pfw["mBatch"]=True
pfw["mWallTime"]="12:00:00"
pfw["mCores"]=pfw["xpar"]*pfw["ypar"]*pfw["zpar"]
pfw["mSubmitJobs"]=True
pfw["autoRestart"]=False

# GEOS MPM i/o parameters ---------------------------------------------------------------

pfw["outputType"]="silo"

# GEOSX MPM PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime            
pfw["plotInterval"]=stopTime/100
pfw["restartInterval"]=stopTime/5 # Don't need restarts for now

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25 
pfw["initialDt"]=1e-16
pfw["cpdiDomainScaling"]=1
pfw["damageFieldPartitioning"]=1

pfw["frictionCoefficient"]=0.25

# Exact contact surface options.
pfw["contactGapCorrection"]="Implicit"
pfw["explicitSurfaceNormalInfluence"]= 1000.  # 1000
pfw["useSurfacePositionForContact"]= 1  # 1


# END GEOSX MPM PARAMETERS ---------------------------------------------------------------

# Deformation ---------------------------------------------------------------------------------
pfw["periodic"] = [ True, False, False]
pfw["boundaryConditionTypes"]=[ 0, 0, 2, 2, 1, 1 ]

# "enablePrescribedBoundaryTransverseVelocities must be of length 6. "
# "The 6 entries correspond to transverse velocity BCs on the x-, x+, y-, y+, z- and z+ faces." 
pfw['enablePrescribedBoundaryTransverseVelocities'] = [0, 0, 1, 1, 0, 0]
pfw['prescribedBoundaryTransverseVelocities'] = [ [0., 0.], [0.,0.], [0.0, -pdcVelocity], [0.0, 0.0], [0., 0.], [0.,0.] ]
pfw["fTableInterpType"]="Smoothstep"
pfw["prescribedBoundaryFTable"]=1
pfw["fTable"]=[[0,        1.00, 1.00, 1.00],
               [stopTime, 1.00, (domainY - 2.0*pdcDepth)/domainY, 1.00]
               ]

# Define all the geometric objects -------------------------------------------------------

pdc_x1 = np.array( [ 0.5*domainX + 0.5*pdcLength*np.cos(pdcBackRakeAngle), pfw["ymin"] + sampleY + 0.5*pdcDiameter*np.cos(pdcBackRakeAngle) , 0.0 ] )
pdc_x2 = pdc_x1 + pdcLength*np.array([-np.cos(pdcBackRakeAngle), + np.sin(pdcBackRakeAngle), 0.0 ] )
                      
cutter = geom.cylinder('pdc',
                      x1 = pdc_x1,
                      x2 = pdc_x2,
                      r = pdcDiameter*0.5,
                      vel = [0.0,0,0],
                      mat = 0,
                      group = 0)

substrate = geom.box('substrate',
    [pfw["xmin"], pfw["ymin"], pfw["zmin"]],
    [pfw["xmax"], pfw["ymin"] + sampleY, pfw["zmax"]],
    vel=[-pdcVelocity, 0.0, 0.0], mat=1, group=0, dim=3, flaggedSurfaces=[False, False, False, True, False, False])

DX = domainX/pfw["nI"]
weibullFlawSize = 6.0*DX
weibullVolume = 1000 # mm^3
grainWeibullModulus = 6.0

weibullSubstrate = geom.voronoiWeibullBoxWrapper('weibullSubstrate',                                                     
    subObject=substrate,       
    x0 = np.array( [ pfw["xmin"] - weibullFlawSize, pfw["ymin"] - weibullFlawSize, pfw["zmin"] - weibullFlawSize ] ),
    x1 = np.array( [ pfw["xmax"] + weibullFlawSize, pfw["ymin"] + sampleY + weibullFlawSize, pfw["zmax"] + weibullFlawSize ] ),
    flawSize=weibullFlawSize,
    weibullVolume=weibullVolume,
    weibullModulus=grainWeibullModulus,
    weibullSeed=1,
    vMin=(DX)**3.,
    vpts=None,
    dim=3,
    randomMatDir = False
    )

pfw["objects"]=[cutter, weibullSubstrate]

pfw["materials"] = [ "diamond", "sand" ]
pfw["materialPropertyString"]="""
<CeramicDamage
	name="sand"
	defaultDensity="""+'"'+str(density)+'"'+"""
	defaultBulkModulus="""+'"'+str(bulk)+'"'+"""
	defaultShearModulus="""+'"'+str(shear)+'"'+"""
	tensileStrength="""+'"'+str(tensileStrength)+'"'+"""
	compressiveStrength="""+'"'+str(compressiveStrength)+'"'+"""
	maximumStrength="""+'"'+str(maximumStrength)+'"'+"""
	crackSpeed="""+'"'+str(crackSpeed)+'"'+"""
	damagedMaterialFrictionSlope="""+'"' + str(damagedMaterialFrictionSlope) + '"' + """
	/>
 <VonMisesJ
  name="diamond"
  defaultDensity="3.3"
  defaultBulkModulus="54.0"
  defaultShearModulus="44.0"
  defaultYieldStrength="35.0"
  />
"""