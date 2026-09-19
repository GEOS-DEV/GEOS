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
 * @file MixedMimeticFields.hpp
 */

#ifndef GEOS_MIXEDMIMETIC_MIXEDMIMETICFIELDS_HPP_
#define GEOS_MIXEDMIMETIC_MIXEDMIMETICFIELDS_HPP_

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace mixedMimetic
{

DECLARE_FIELD( faceMassFlux,
               "faceMassFlux",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Face mass flux unknown of the mixed mimetic formulation (positive in the direction of the global face normal)" );

DECLARE_FIELD( faceMassFlux_n,
               "faceMassFlux_n",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Face mass flux at the previous converged time step" );

DECLARE_FIELD( mfdFlag,
               "mfdFlag",
               array1d< integer >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Inner product selector eta of the cell: 1 = consistent (MFD) product, 0 = diagonal (TPFA) product" );

DECLARE_FIELD( prescribedMfdFlag,
               "prescribedMfdFlag",
               array1d< integer >,
               -1,
               LEVEL_0,
               WRITE_AND_READ,
               "Constraint on eta read from the mesh, not the flag itself (mfdFlag is the answer): -1 = free, the layers "
               "decide; 0 = diagonal (TPFA) product; 1 = consistent (MFD) product; a prescribed value is final; "
               "any other value is rejected" );

DECLARE_FIELD( consistencyIndicator,
               "consistencyIndicator",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Residual-based consistency indicator per cell" );

DECLARE_FIELD( degeneracyIndicator,
               "degeneracyIndicator",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Cell volume as a percentage of the total volume of its node star (admissibility of the consistent (MFD) product)" );

DECLARE_FIELD( faceOrientationCell,
               "faceOrientationCell",
               array1d< globalIndex >,
               -1,
               NOPLOT,
               NO_WRITE,
               "Global index of the cell E_min whose outward normal orients the face flux: sigma_{E,f} = +1 if E = E_min, -1 otherwise" );

DECLARE_FIELD( faceDofScale,
               "faceDofScale",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Characteristic scale of the face flux unknown, s_f = p_scale / |M_ff|, used by the residual-norm weights" );

DECLARE_FIELD( faceStencilLabel,
               "faceStencilLabel",
               array1d< integer >,
               0,
               NOPLOT,
               NO_WRITE,
               "Face classification for the adaptive solver (0 = all adjacent cells TPFA-compatible: the flux is condensed "
               "into a two-point expression; 1 = adjacent to at least one MFD-compatible cell: the flux stays a live unknown)" );

DECLARE_FIELD( faceResidual,
               "faceResidual",
               array1d< real64 >,
               0,
               LEVEL_1,
               NO_WRITE,
               "Face-assembled normalized residual of the consistency layer" );

}

}

}

#endif // GEOS_MIXEDMIMETIC_MIXEDMIMETICFIELDS_HPP_
