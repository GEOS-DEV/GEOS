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
 * @file ReactiveFluidFields.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_REACTIVEFLUID_REACTIVEFLUIDFIELDS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_REACTIVEFLUID_REACTIVEFLUIDFIELDS_HPP_

#include "constitutive/fluid/reactivefluid/ReactiveFluidLayouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace reactivefluid
{

using array3dLayoutSpecies = array3d< real64, constitutive::reactivefluid::LAYOUT_SPECIES >;
using array4dLayoutSpecies_dC = array4d< real64, constitutive::reactivefluid::LAYOUT_SPECIES_DC >;

DECLARE_FIELD( initialPrimarySpeciesConcentration,
               "initialPrimarySpeciesConcentration",
               array3dLayoutSpecies,
               1e-16,
               LEVEL_0,
               WRITE_AND_READ,
               "initialPrimarySpeciesConcentration" );

DECLARE_FIELD( primarySpeciesAggregateConcentration,
               "primarySpeciesAggregateConcentration",
               array3dLayoutSpecies,
               1e-16,
               LEVEL_0,
               WRITE_AND_READ,
               "primarySpeciesAggregateConcentration" );

DECLARE_FIELD( primarySpeciesAggregateConcentration_n,
               "primarySpeciesAggregateConcentration_n",
               array3dLayoutSpecies,
               1e-16,
               LEVEL_0,
               WRITE_AND_READ,
               "primarySpeciesAggregateConcentration at the previous timestep" );

DECLARE_FIELD( primarySpeciesMobileAggregateConcentration,
               "primarySpeciesMobileAggregateConcentration",
               array3dLayoutSpecies,
               1e-16,
               LEVEL_0,
               WRITE_AND_READ,
               "primarySpeciesMobileAggregateConcentration" );

DECLARE_FIELD( dPrimarySpeciesAggregateConcentration_dLogPrimarySpeciesConcentrations,
               "dPrimarySpeciesAggregateConcentration_dLogPrimarySpeciesConcentrations",
               array4dLayoutSpecies_dC,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Deivatives of primarySpeciesAggregateConcentration w.r.t log primary species concentration" );

DECLARE_FIELD( dPrimarySpeciesMobileAggregateConcentration_dLogPrimarySpeciesConcentrations,
               "dPrimarySpeciesMobileAggregateConcentration_dLogPrimarySpeciesConcentrations",
               array4dLayoutSpecies_dC,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Deivatives of primarySpeciesMobileAggregateConcentration w.r.t log primary species concentration" );

DECLARE_FIELD( secondarySpeciesConcentration,
               "secondarySpeciesConcentration",
               array3dLayoutSpecies,
               1e-16,
               NOPLOT,
               WRITE_AND_READ,
               "secondarySpeciesConcentration" );

DECLARE_FIELD( kineticReactionRates,
               "kineticReactionRates",
               array3dLayoutSpecies,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "kineticReactionRates" );

DECLARE_FIELD( aggregateSpeciesRates,
               "aggregateSpeciesRates",
               array3dLayoutSpecies,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "aggregateSpeciesRates" );

DECLARE_FIELD( dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations,
               "dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations",
               array4dLayoutSpecies_dC,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Deivatives of aggregate concentration rates w.r.t log primary species concentration" );
}

}

}

#endif // GEOS_CONSTITUTIVE_FLUID_REACTIVEFLUID_REACTIVEFLUIDFIELDS_HPP_
