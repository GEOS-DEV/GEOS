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
 * @file ProppantTransportFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_PROPPANT_PROPPANTTRANSPORTFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_PROPPANT_PROPPANTTRANSPORTFIELDS_HPP_

#include "mesh/MeshFields.hpp"

namespace geos
{
/**
 * A scope for field traits.
 */
namespace fields
{

namespace proppant
{

DECLARE_FIELD_WITH_NAMESPACE( proppantConcentration,
                              "proppantConcentration",
                              array1d< real64 >,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Proppant concentration",
                              "proppant" );

DECLARE_FIELD_WITH_NAMESPACE( proppantConcentration_n,
                              "proppantConcentration_n",
                              array1d< real64 >,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Proppant concentration at the previous converged time step",
                              "proppant" );

DECLARE_FIELD_WITH_NAMESPACE( componentConcentration,
                              "componentConcentration",
                              array2d< real64 >,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Component concentration",
                              "proppant" );

DECLARE_FIELD_WITH_NAMESPACE( componentConcentration_n,
                              "componentConcentration_n",
                              array2d< real64 >,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Component concentration at the previous converged time step",
                              "proppant" );

DECLARE_FIELD_WITH_NAMESPACE( bcComponentConcentration,
                              "bcComponentConcentration",
                              array2d< real64 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Boundary component concentration",
                              "proppant" );

DECLARE_FIELD_WITH_NAMESPACE( cellBasedFlux,
                              "cellBasedFlux",
                              array2d< real64 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Cell-based flux",
                              "proppant" );

DECLARE_FIELD_WITH_NAMESPACE( isProppantBoundary,
                              "isProppantBoundary",
                              array1d< integer >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Flag denoting the proppant boundary",
                              "proppant" );

DECLARE_FIELD_WITH_NAMESPACE( isProppantMobile,
                              "isProppantMobile",
                              array1d< integer >,
                              1.0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Flag indicating whether proppant is mobile",
                              "proppant" );

DECLARE_FIELD_WITH_NAMESPACE( componentDensity_n,
                              "componentDensity_n",
                              array2d< real64 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Component density at the previous converged time step",
                              "proppant" );

DECLARE_FIELD_WITH_NAMESPACE( proppantPackVolumeFraction,
                              "proppantPackVolumeFraction",
                              array1d< real64 >,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Proppant pack volume fraction",
                              "proppant" );

DECLARE_FIELD_WITH_NAMESPACE( proppantExcessPackVolume,
                              "proppantExcessPackVolume",
                              array1d< real64 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Proppant pack volume fraction",
                              "proppant" );

DECLARE_FIELD_WITH_NAMESPACE( proppantLiftFlux,
                              "proppantLiftFlux",
                              array1d< real64 >,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Proppant lift flux",
                              "proppant" );

}

}

}

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_PROPPANT_PROPPANTTRANSPORTFIELDS_HPP_
