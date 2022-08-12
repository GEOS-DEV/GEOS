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
 * @file MultiFluidExtrinsicData.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDEXTRINSICDATA_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDEXTRINSICDATA_HPP_

#include "constitutive/fluid/layouts.hpp"
#include "mesh/ExtrinsicMeshData.hpp"

namespace geosx
{

namespace extrinsicMeshData
{

namespace reactivefluid
{

using array2dLayoutComp = array2d< real64, compflow::LAYOUT_COMP >;

EXTRINSIC_MESH_DATA_TRAIT( primarySpeciesConcentration,
                           "primarySpeciesConcentration",
                           array2dLayoutComp,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "primarySpeciesConcentration" );

EXTRINSIC_MESH_DATA_TRAIT( primarySpeciesTotalConcentration,
                           "primarySpeciesTotalConcentration",
                           array2dLayoutComp,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "primarySpeciesTotalConcentration" );

EXTRINSIC_MESH_DATA_TRAIT( secondarySpeciesConcentration,
                           "secondarySpeciesConcentration",
                           array2dLayoutComp,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "secondarySpeciesConcentration" );

EXTRINSIC_MESH_DATA_TRAIT( kineticReactionRates,
                           "kineticReactionRates",
                           array2dLayoutComp,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "kineticReactionRates" );
}

}

}

#endif // GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDEXTRINSICDATA_HPP_
