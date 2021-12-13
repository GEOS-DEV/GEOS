/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PermeabilityExtrinsicData.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_PERMEABILITYEXTRINSICDATA_HPP_
#define GEOSX_CONSTITUTIVE_PERMEABILITYEXTRINSICDATA_HPP_

#include "common/DataLayouts.hpp"
#include "mesh/ExtrinsicMeshData.hpp"

namespace geosx
{

namespace extrinsicMeshData
{

namespace permeability
{

EXTRINSIC_MESH_DATA_TRAIT( permeability,
                           "permeability",
                           array3d< real64 >,
                           -1,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Rock permeability" );

EXTRINSIC_MESH_DATA_TRAIT( dPerm_dPressure,
                           "dPerm_dPressure",
                           array3d< real64 >,
                           0,
                           LEVEL_3,
                           WRITE_AND_READ,
                           "Derivative of rock permeability with respect to pressure" );

EXTRINSIC_MESH_DATA_TRAIT( dPerm_dAperture,
                           "dPerm_dAperture",
                           array3d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Derivative of rock permeability with respect to aperture" );

EXTRINSIC_MESH_DATA_TRAIT( permeabilityMultiplier,
                           "permeabilituMultiplier",
                           array3d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Rock permeability multiplier" );

}

}

}

#endif // GEOSX_CONSTITUTIVE_PERMEABILITYEXTRINSICDATA_HPP_
