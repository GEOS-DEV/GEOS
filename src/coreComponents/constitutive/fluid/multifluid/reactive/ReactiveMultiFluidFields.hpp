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
 * @file ReactiveMultiFluidFields.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_REACTIVE_MULTIFLUIDFIELDS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_REACTIVE_MULTIFLUIDFIELDS_HPP_

#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace reactivefluid
{

using array2dLayoutComp = array2d< real64, compflow::LAYOUT_COMP >;
using array3dLayoutComp_dC = array3d< real64, compflow::LAYOUT_COMP_DC >;

DECLARE_FIELD( primarySpeciesConcentration,
               "primarySpeciesConcentration",
               array2dLayoutComp,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "primarySpeciesConcentration" );

DECLARE_FIELD( primarySpeciesAggregateConcentration,
               "primarySpeciesAggregateConcentration",
               array2dLayoutComp,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "primarySpeciesAggregateConcentration" );

DECLARE_FIELD( primarySpeciesAggregateConcentration_n,
               "primarySpeciesAggregateConcentration_n",
               array2dLayoutComp,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "primarySpeciesAggregateConcentration at the previous timestep" );

DECLARE_FIELD( dPrimarySpeciesAggregateConcentration_dLogPrimaryConc,
               "dPrimarySpeciesAggregateConcentration_dLogPrimaryConc",
               array3dLayoutComp_dC,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Deivatives of primarySpeciesAggregateConcentration w.r.t log primary species concentration" );

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
