# -*- coding: utf-8 -*-
import os
import math
import numpy as np                   # math stuff
import pfw_geometryObjects as geom
import importlib

#[pfw_dependency] pfw:pfw_materials.py
#[pfw_dependency] /usr/workspace/crook5/pfwx/validation_materials.py

import pfw_materials as material_db
import validation_materials as valMats

pfw = {}

# INPUTS --------------------------------------------------------------------------------

temp = 273.15+23
maxCompressiveStrain = 0.2
radialBias = 0.0

case_name = os.environ.get("DISK_CASE_NAME", "disk")

g2pMethod = os.environ.get("DISK_METHOD", "PIC")
usePrecomputed = os.environ.get("DISK_PRECOMPUTE", "Precompute")
seed = 1436359

variant_label = os.environ.get("DISK_VARIANT_LABEL", f"{g2pMethod}_{usePrecomputed}")

# BATCH PARAMETERS _----------------------------------------------------------------------

pfw["runDebug"] = True
pfw["mBatch"]=True
pfw["mBank"]="mahem"
pfw["mWallTime"]="00:05:00"
pfw["mSubmitJobs"]=True
pfw["autoRestart"]=False

# DOMAIN ---------------------------------------------------------------------------------

pfw["planeStrain"]=1
pfw["periodic"]=[True, False, False]

refine=3
pfw["xpar"] = refine
pfw["ypar"] = refine
pfw["zpar"] = 1

cpp=20
pfw["nI"]=pfw["xpar"]*cpp          # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp          # grid cells in the y-direction
pfw["nK"]=3                        # grid cells in the z-direction
pfw["ppc"]=2

diameter = 25.4/2
radius = diameter/2

domainHeight = diameter
domainWidth = domainHeight

dx = domainHeight/(pfw["nI"]-2)
dy = domainWidth/(pfw["nJ"]-2)

domainLength = 3*dx

dz = domainLength/(pfw["nK"]-2)

pfw["xmin"] = -0.5*domainWidth   # mm
pfw["xmax"] =  0.5*domainWidth   # mm
pfw["ymin"] = -0.5*domainHeight  # mm
pfw["ymax"] =  0.5*domainHeight  # mm
pfw["zmin"] = -0.5*domainLength  # mm
pfw["zmax"] =  0.5*domainLength  # mm

# GEOSX MPM SOLVER PARAMETERS -------------------------------------------------------------------

stopTime = 20.0*domainHeight/math.sqrt(15./1.92)*(maxCompressiveStrain/0.25) # 200.0*domainHeight/math.sqrt(15./1.92)*(maxCompressiveStrain/0.25)
pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/200.0
pfw["restartInterval"]=stopTime
pfw['lastRestartBufferInSeconds']=600

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25
pfw["initialDt"]=1e-16

pfw["plotUnscaledParticles"]=1
pfw["cpdiDomainScaling"]=1
pfw["damageFieldPartitioning"]=0

pfw["reactionHistory"]=1
pfw["reactionWriteInterval"]=stopTime/2000.0
pfw["boxAverageHistory"]=1
pfw["boxAverageWriteInterval"]=stopTime/2000.0

pfw["maxParticleVelocity"]=10.0
pfw["minParticleJacobian"]=0.01
pfw["maxParticleJacobian"]=10.0
pfw["FSubcycles"]=10

pfw["gridToParticleMapping"]="Precomputed"
pfw["updateMethod"]=g2pMethod
pfw["updateOrder"]="2"

pfw["contactGapCorrection"]="Implicit"
pfw["explicitSurfaceNormalInfluence"]=1
pfw["useSurfacePositionForContact"]= 1 # Turned off for now because of issues with large deformation

pfw["plottableFields"]=[
    "particleRank",
    "particleMass",
    "particleVelocity",
    "particleSurfaceFlag",
    "particleStrengthScale",
    "particleDamage",
    "particleDamageGradient",
    "particleCohesiveForce",
    "particleStress",
    "particleSurfaceNormal",
    "particleSurfacePosition",
    "particlePlasticStrain"
]

pfw['particleFileFields'] = [
    "Velocity",
    "MaterialType",
    "ContactGroup",
    "Damage",
    "StrengthScale",
    "SurfaceFlag",
    "SurfaceNormal",
    "SurfacePosition",
    "RVector"
]	

# MATERIAL PROPERTIES ----------------------------------------------------------------
matFile = importlib.import_module("pfw_materials")
# quartz = matFile.elasticAluminumSI
quartz = matFile.quartz

pfw["materials"] = [quartz["name"]]
pfw["materialPropertyString"] = quartz["materialString"]

# GEOMETRY OBJECTS -------------------------------------------------------

# Weibull variability is applied per each grain
grainDiameter = diameter*0.3
def make_objects():
    disk1 = geom.cylinder(
        "disk1",
        [0.0, 0.0, pfw["zmin"]],
        [0.0, 0.0, pfw["zmax"]],
        domainHeight / 2.0,
        vel = [0.0, 0.0, 0.0],
        mat = 0,
        group = 0,
    )

    weibullFlawSize = 6.0 * dx

    weibullSample = geom.voronoiWeibullBoxWrapper(
        "weibullSubstrate",
        subObject = disk1,
        x0 = np.array([
            pfw["xmin"] - weibullFlawSize,
            pfw["ymin"] - weibullFlawSize,
            pfw["zmin"] - weibullFlawSize,
        ]),
        x1 = np.array([
            pfw["xmax"] + weibullFlawSize,
            pfw["ymax"] + weibullFlawSize,
            pfw["zmax"] + weibullFlawSize,
        ]),
        flawSize = weibullFlawSize,
        weibullVolume = quartz["weibullReferenceVolume"],
        weibullModulus = quartz["weibullModulus"], 
        # weibullVolume = 1.0,
        # weibullModulus = 5.0, 
        weibullSeed = 1,
        vMin = dx**3,
        vpts = None,
        dim = 3,
        randomMatDir = False,
    )
    # damageWrapper = geom.damageBoxWrapper("dmg",
    #            weibullSample,
    #            [-domainHeight / 8.0, -domainHeight / 8.0, -100.0],
    #            [domainHeight / 8.0, domainHeight / 8.0, 100.0],
    #            1.0)

    return [weibullSample]


# DEFORMATION ---------------------------------------------------------------------------------

pfw["boundaryConditionTypes"]= [ 0, 0, 2, 2, 1, 1 ]

pfw["fTableInterpType"]='Cosine'
pfw["prescribedBoundaryFTable"]=0
pfw["fTable"]=[[0.0,      1.0, 1.0,                      1.0],
               [stopTime, 1.0, 1.0-maxCompressiveStrain, 1.0]]

pfw_expected = {
    "variant_label": variant_label,
    "case_name": case_name
}