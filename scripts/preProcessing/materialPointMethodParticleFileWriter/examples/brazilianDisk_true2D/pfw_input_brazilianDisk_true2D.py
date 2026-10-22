# -*- coding: utf-8 -*-
# ---- GEOS-MPM example input metadata ----
# Purpose: Native-2D Brazilian disk, CPDI bulk and CPTI surface triangles.
# Solver/PFW features: FMPM(2), quartz damage, 2D Weibull strength, moving y
# boundaries, periodic x translation, native surface refinement and histories.
# Requires: GEOS-native-2D.patch and its updated particle file writer.
# ---- end example input metadata ----
"""Native-2D mixed CPDI/CPTI version of the Brazilian-disk FMPM example.

A unit-diameter quartz disk translates in periodic x while moving y boundaries
compress the specimen by 20% using the existing cosine F-table loading. Interior
particles are CPDI parallelograms. Accepted surface cuts are CPTI triangles;
uncut or Cartesian-fallback pieces retain CPDI interpolation.

The mesh is genuinely XY: no z bounds, z cells, or extruded object end caps.
One physical thickness supplies volume, mass, interface area and Weibull sample
volume. Constitutive tensors retain the plane-strain out-of-plane stress.

This input generates files and a batch script without submitting a job. Change
mSubmitJobs to True when ready to run with your usual userDefs configuration.
"""

import numpy as np
import pfw_geometryObjects as geom
# [pfw_dependency] pfw:pfw_materials.py
import pfw_materials as matdb

pfw = {}

# -----------------------------------------------------------------------------
# 1. Physical specimen and native XY background grid
# -----------------------------------------------------------------------------
# The default grid and 3x3 solver partitioning match brazilianDisk_FMPM.
# refine/cpp change resolution and partitioning; physical thickness is fixed.
refine = 3
cpp = 24
ppc = 2

domainHeight = 1.0
domainWidth = 1.5 * domainHeight
diskRadius = 0.5 * domainHeight
physicalThickness = domainHeight / 24.0
# This matches the old example's default represented thickness, but unlike an
# artificial z extent it does not change when refine/cpp are changed.

pfw['formulation'] = 'PlaneStrain'
pfw['planeStrainThickness'] = physicalThickness
pfw['xmin'], pfw['xmax'] = -0.5 * domainWidth, 0.5 * domainWidth
pfw['ymin'], pfw['ymax'] = 0.0, domainHeight
pfw['xpar'], pfw['ypar'], pfw['zpar'] = refine, refine, 1
pfw['nI'] = round(1.5 * refine * cpp)
pfw['nJ'] = refine * cpp
pfw['ppc'] = ppc
pfw['periodic'] = [True, False]

# nI includes no x padding because x is periodic. Nonperiodic nJ includes two
# padding cells, as in the existing PFW convention. DX sets the Weibull spacing.
DX = domainWidth / pfw['nI']

# -----------------------------------------------------------------------------
# 2. Quartz material and physical 2D Weibull strength field
# -----------------------------------------------------------------------------
# Keep the existing Brazilian-disk material preset and seed. Native dim=2
# Voronoi cells have statistical volume area*physicalThickness. The realization
# consequently differs from the old 3D Weibull construction.
quartz = matdb.quartz
pfw['materials'] = [quartz['name']]
pfw['materialPropertyString'] = quartz['materialString']

# A native disk needs only an XY center and radius. Its normals and surface
# offsets cannot select an unwanted z end cap. Type 2 means CPDI bulk.
disk = geom.disk(
    'disk', center=[0.0, 0.5 * domainHeight], r=diskRadius,
    vel=[1.0, 0.0, 0.0], mat=0, group=0, particleType=2,
)
weibullFlawSize = 6.0 * DX
weibullDisk = geom.voronoiWeibullBoxWrapper(
    'weibullDisk', subObject=disk,
    x0=np.array([pfw['xmin'] - weibullFlawSize,
                 pfw['ymin'] - weibullFlawSize]),
    x1=np.array([pfw['xmax'] + weibullFlawSize,
                 pfw['ymax'] + weibullFlawSize]),
    flawSize=weibullFlawSize,
    weibullVolume=quartz['weibullReferenceVolume'],
    weibullModulus=quartz['weibullModulus'], weibullSeed=1,
    vMin=DX**2 * physicalThickness,
    dim=2, randomMatDir=False,
)
pfw['objects'] = [weibullDisk]

# -----------------------------------------------------------------------------
# 3. CPDI interior with a refined CPTI surface band
# -----------------------------------------------------------------------------
# The coarse lattice has ppc^2 particles per occupied XY cell. Candidate surface
# domains are subdivided by surfaceRefinementLevel in each in-plane direction.
# Actual cut pieces retain triangle geometry and use type 3 (CPTI). Uncut pieces
# in the band are smaller CPDI quads; fallback pieces also remain CPDI quads.
# The hexTet preset is 3D only, so use this explicit native-2D combination.
pfw['particleGeometry'] = 'fromRefinement'
pfw['surfaceRefinementMethod'] = 'simplex'
pfw['surfaceRefinementLevel'] = 2  # Increase to 4 for a finer surface band.
pfw['simplexParticleType'] = 3
pfw['surfaceRefinementFallback'] = 'cartesian'
pfw['surfaceRefinementPlaneTolerance'] = 0.05

# Surface normal/position fields are included explicitly for geometry/contact
# diagnostics. Unified G/G0 geometry is always written; no RVector is needed.
pfw['particleFileFields'] = [
    'Velocity', 'MaterialType', 'ContactGroup', 'StrengthScale',
    'SurfaceFlag', 'SurfaceNormal', 'SurfacePosition',
]

# -----------------------------------------------------------------------------
# 4. Time integration, mapping and geometry controls
# -----------------------------------------------------------------------------
stopTime = 1.0
pfw['endTime'] = stopTime
pfw['timeIntegrationOption'] = 'ExplicitDynamic'
pfw['cflFactor'] = 0.25
pfw['initialDt'] = 1.0e-16
pfw['updateMethod'] = 'FMPM'
pfw['updateOrder'] = 2
pfw['damageFieldPartitioning'] = 1
pfw['frictionCoefficient'] = 0.25
pfw['needsNeighborList'] = 1
pfw['solverProfiling'] = 1

# CPDI uses its Homel domain projection; CPTI uses uniform diameter clipping.
# These affect temporary integration geometry, not the displayed material G.
pfw['cpdiDomainScaling'] = 1
pfw['cpdiDomainScalingType'] = 'homel'
pfw['geometryReset'] = 'FOnly'

# Disable runtime splitting for an initial comparison of the mixed layout.
# To exercise native runtime splitting too, set oversizedParticleTreatment to
# 'split', defGradResetMaxParticleDomainToGridCellRatio to e.g. 0.9, and
# oversizedParticleResetRankParticleCountThreshold to -1 (no count fallback).
pfw['subdivideParticles'] = 0
pfw['oversizedParticleTreatment'] = 'none'

# -----------------------------------------------------------------------------
# 5. Periodic translation and Brazilian-disk compression
# -----------------------------------------------------------------------------
# Face order: x-, x+, y-, y+. Outflow=0, Symmetry=1, Moving=2, Contact=3.
# Periodicity handles the x pair. The moving y pair follows the cosine F table.
# This preserves the loading of the existing example; it has no separate platens.
pfw['boundaryConditionTypes'] = [0, 0, 2, 2]
pfw['prescribedBcTable'] = 0
pfw['prescribedBoundaryFTable'] = 1
pfw['fTableInterpType'] = 'Cosine'
# Rows remain [time, Fxx, Fyy, Fzz]; plane strain requires Fzz=1.
pfw['fTable'] = [
    [0.0,      1.0, 1.0, 1.0],
    [stopTime, 1.0, 0.8, 1.0],
]

# -----------------------------------------------------------------------------
# 6. Native particle output, histories and job preparation
# -----------------------------------------------------------------------------
# Silo contains native quads/triangles. Use outputType='vtk' for ParaView.
# Plot 200 intervals through loading; keep reaction histories at finer cadence.
pfw['outputType'] = 'silo'
pfw['plotInterval'] = stopTime / 200.0
pfw['restartInterval'] = 2.0 * stopTime
pfw['reactionHistory'] = 1
pfw['reactionWriteInterval'] = stopTime / 2000.0
pfw['boxAverageHistory'] = 1
pfw['boxAverageWriteInterval'] = stopTime / 2000.0
pfw['runDebug'] = True
pfw['mBatch'] = True
pfw['mWallTime'] = '00:05:00'
pfw['mSubmitJobs'] = True
