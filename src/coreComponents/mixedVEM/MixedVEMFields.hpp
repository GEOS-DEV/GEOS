/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MixedVEMFields.hpp
 *
 * Mesh fields of the mixed VEM.
 *
 * The primary unknowns do not live where a visualization writer expects them: the stress
 * is a face quantity and the displacement is a six dimensional rigid motion per cell, so
 * neither is directly a cell array. The cell fields below are the writable projections,
 * produced by MixedVEMCellOutput.hpp:
 *
 *   displacement  u_h(x_E), the translation part of u_h|_E, comparable to a nodal
 *                 displacement sampled at the cell center,
 *   rotation      the infinitesimal rotation omega of u_h|_E,
 *   stress        Pi_E sigma_h, the constant symmetric stress of the cell, comparable to
 *                 a cell averaged finite element stress.
 */

#ifndef GEOS_MIXEDVEM_MIXEDVEMFIELDS_HPP_
#define GEOS_MIXEDVEM_MIXEDVEMFIELDS_HPP_

#include "common/DataLayouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace mixedVEM
{

DECLARE_FIELD( faceStress,
               "faceStress",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Stress degrees of freedom of a face, the six traction modes of T_h(f)" );

DECLARE_FIELD( multiplier,
               "displacementMultiplier",
               array2d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "Displacement trace moments on the interior faces, the hybridization multiplier" );

DECLARE_FIELD( displacementTrace,
               "displacementTrace",
               array2d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Prescribed displacement on a boundary face, the natural condition of the mixed form" );

DECLARE_FIELD( traction,
               "faceTraction",
               array2d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Prescribed traction on a boundary face, the essential condition of the mixed form" );

DECLARE_FIELD( boundaryType,
               "faceBoundaryType",
               array1d< integer >,
               0,
               NOPLOT,
               NO_WRITE,
               "Boundary role of a face: 0 interior, 1 on the boundary of the domain" );

DECLARE_FIELD( displacementMask,
               "displacementComponentMask",
               array1d< integer >,
               0,
               NOPLOT,
               NO_WRITE,
               "Bit k is set when the displacement component k is prescribed on the face; "
               "the complementary components carry a prescribed traction" );

DECLARE_FIELD( rigidMotion,
               "rigidMotion",
               array2d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Coefficients of u_h on the rigid body motion basis of the cell" );

DECLARE_FIELD( displacement,
               "displacement",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Displacement of the cell evaluated at its center" );

DECLARE_FIELD( rotation,
               "rotation",
               array2d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "Infinitesimal rotation of the cell" );

DECLARE_FIELD( stress,
               "stress",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Projected constant stress of the cell, in the order (xx, yy, zz, yz, xz, xy)" );

}

}

}

#endif // GEOS_MIXEDVEM_MIXEDVEMFIELDS_HPP_
