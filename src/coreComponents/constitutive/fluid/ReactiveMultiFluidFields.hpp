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

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUIDFIELDS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUIDFIELDS_HPP_

#include "constitutive/fluid/layouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace reactivefluid
{

using array2dLayoutComp = array2d< real64, compflow::LAYOUT_COMP >;

DECLARE_FIELD( primarySpeciesConcentration,
               "primarySpeciesConcentration",
               array2dLayoutComp,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "primarySpeciesConcentration" );

DECLARE_FIELD( primarySpeciesTotalConcentration,
               "primarySpeciesTotalConcentration",
               array2dLayoutComp,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "primarySpeciesTotalConcentration" );

DECLARE_FIELD( secondarySpeciesConcentration,
               "secondarySpeciesConcentration",
               array2dLayoutComp,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "secondarySpeciesConcentration" );

DECLARE_FIELD( kineticReactionRates,
               "kineticReactionRates",
               array2dLayoutComp,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "kineticReactionRates" );
}

}

}

#endif // GEOS_CONSTITUTIVE_FLUID_MULTIFLUIDFIELDS_HPP_
