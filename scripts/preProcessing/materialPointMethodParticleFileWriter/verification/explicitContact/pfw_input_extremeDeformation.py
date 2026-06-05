# -*- coding: utf-8 -*-
import pfw_geometryObjects as geom   # this contains all the geometry object functions for pfw
import numpy as np                   # math stuff
from sklearn.neighbors import KDTree          # nearest neighbor search with KDTree

pfw = {}
pfw["runDebug"] = False
stopTime = 30

# Domain ---------------------------------------------------------------------------------
domainHeight = 1.0
domainWidth = 1.0

pfw["xmin"] = -0.5*domainWidth # mm
pfw["xmax"] = 0.5*domainWidth # mm
pfw["ymin"] = -0.5*domainHeight # mm
pfw["ymax"] = 0.5*domainHeight # mm

pfw["planeStrain"] = 1

pfw["periodic"] = [True, False, False]

refine = 5
cpp = 20 # cells per partition in each direction

pfw["xpar"]=refine  # grid partitions
pfw["ypar"]=refine
pfw["zpar"]=1

pfw["nI"]=pfw["xpar"]*cpp  # grid cells in the x-direction
pfw["nJ"]=pfw["ypar"]*cpp  # grid cells in the y-direction
pfw["nK"]=3          # grid cells in the z-direction
pfw["ppcx"]=2         # particles per cell in each direction
pfw["ppcy"]=5         # particles per cell in each direction

domainLength = domainHeight/refine/cpp

pfw["zmin"] =-0.5*domainLength # mm
pfw["zmax"] = 0.5*domainLength # mm

# Batch parameters for GEOS runs.  --------------------------------------------------------

pfw["mBatch"]=True
pfw["mWallTime"] = "01:00:00"
pfw["mSubmitJobs"]=True
pfw["autoRestart"]=False
pfw["lastRestartBufferInSeconds"] = 600

# GEOSX MPM PARAMETERS -------------------------------------------------------------------

pfw["endTime"]=stopTime
pfw["plotInterval"]=stopTime/200
pfw["restartInterval"]=stopTime/20

pfw["timeIntegrationOption"]="ExplicitDynamic"
pfw["cflFactor"]=0.25
pfw["initialDt"]=1e-16
pfw["cpdiDomainScaling"]=1
pfw["damageFieldPartitioning"]=1

pfw["solverProfiling"]=0
pfw["needsNeighborList"]=0
pfw["reactionHistory"]=1
pfw["boxAverageHistory"]=1
pfw["useEvents"]=0 #1
pfw["frictionCoefficient"]=0.0

pfw["plotUnscaledParticles"]=1
# pfw["overlapCorrection"]=2
# pfw["overlapThreshold1"]=1.05
# pfw["overlapThreshold2"]=1.10

pfw["maxParticleVelocity"]=10.0
pfw["minParticleJacobian"]=0.01
pfw["maxParticleJacobian"]=10.0

pfw["updateMethod"]="FMPM"
pfw["updateOrder"]=2

pfw["contactGapCorrection"] = "Implicit"
pfw["useSurfacePositionForContact"]=1
pfw["explicitSurfaceNormalInfluence"]=1000
pfw["disableSurfaceNormalsAndPositionsOnCPDIScaling"]=1

pfw["particleFileFields"] = ["Velocity",
                             "MaterialType",
                             "ContactGroup",
							 "Damage",
							 "StrengthScale",
                             "SurfaceFlag",
                             "RVector",
							 "SurfaceNormal",
							 "SurfacePosition"]

# END GEOSX MPM PARAMETERS ---------------------------------------------------------------

# Deformation ----------------------------------------------------------------------------

pfw["prescribedBcTable"]=0
pfw["boundaryConditionTypes"]=[ 0, 0, 2, 2, 1, 1 ]

pfw["prescribedBoundaryFTable"]=1
pfw["fTableInterpType"] = "Cosine"
pfw["fTable"]=[[0,            1.0,  1.0,    1.0],
               [stopTime,     1.0,  0.4,    1.0]]

# Define all the geometric objects -------------------------------------------------------

# MATERIAL PROPERTIES --------------------------------------------------------------------

density = 1.00 # mass density mg/mm^3
E = 10.0   # stiffness (GPa)
nu = 0.22   # poisson's ratio (-)
strength = 0.5

pfw["materials"] = ["elasticPlastic"]
pfw["materialPropertyString"]="""
<VonMisesJ
    name="elasticPlastic"
    defaultDensity=""" + '"' + str(density) + '"' + """
    defaultYoungModulus=""" + '"' + str(E) + '"' + """
    defaultPoissonRatio=""" + '"' + str(nu) + '"' + """
    defaultYieldStrength=""" + '"' + str(strength) + '"' + """/>"""

# GEOMETRY -------------------------------------------------------------------------------

def make_objects():
    # Prepopulate domain using poisson disc sampling
    s_radius = 0.1
    s_spacing = 1.05*s_radius
    seed = 12345

    x0 = [pfw["xmin"], pfw["ymin"]]
    dx = [pfw["xmax"]-pfw["xmin"], pfw["ymax"]-pfw["ymin"]]

    circles = geom.poisson(2*s_spacing, x0=x0, dx=dx, seed=12345, numTrials=100, dim=2, periodic=[True, False])
    out = geom.add_pore_images(circles, x0=x0, dx=dx, dim=2, periodic=[True, False])
    circles = out[0]

    objects= []
    for n, c in enumerate(circles):
        circle = geom.cylinder("circle_" + str(n),
                               [c[0], c[1], 10*pfw["zmin"]],
                               [c[0], c[1], 10*pfw["zmax"]],
                               s_radius, vel=[0.0, 0.0, 0.0], mat=0, group=0)
        objects.append(circle)

    return objects
