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
 * @file SinglePhaseReactiveTransportFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVETRANSPORTFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVETRANSPORTFIELDS_HPP_

#include "mesh/MeshFields.hpp"

namespace geos
{
/**
 * A scope for field traits.
 */
namespace fields
{

namespace flow
{
using array2dLayoutComp = array2d< real64, compflow::LAYOUT_COMP >;
using array2dLayoutFluid_dC = array2d< real64, compflow::LAYOUT_FLUID_DC >;

DECLARE_FIELD( logPrimarySpeciesConcentration,
               "logPrimarySpeciesConcentration",
               array2dLayoutComp,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Natural log of primary species concentration (molarity)" );

DECLARE_FIELD( logPrimarySpeciesConcentration_n,
               "logPrimarySpeciesConcentration_n",
               array2dLayoutComp,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Natural log of primary species concentration (molarity) at the previous converged time step" );

DECLARE_FIELD( bcLogPrimarySpeciesConcentration,
               "bcLogPrimarySpeciesConcentration",
               array2dLayoutComp,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Boundary condition for natural log of primary species concentration (molarity)" );

DECLARE_FIELD( totalPrimarySpeciesAmount,
               "totalPrimarySpeciesAmount",
               array2dLayoutComp,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Total amount of species (both primary and secondary) that contain the ion in the primary species (mole)" );

DECLARE_FIELD( totalPrimarySpeciesAmount_n,
               "totalPrimarySpeciesAmount_n",
               array2dLayoutComp,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Total amount of species (both primary and secondary) that contain the ion in the primary species (mole) at the previous converged time step" );

DECLARE_FIELD( dMobility_dLogPrimaryConc,
               "dMobility_dLogPrimaryConc",
               array2dLayoutFluid_dC,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of fluid mobility with respect to log of primary species concentration" );

}

}

}

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVETRANSPORTFIELDS_HPP_
