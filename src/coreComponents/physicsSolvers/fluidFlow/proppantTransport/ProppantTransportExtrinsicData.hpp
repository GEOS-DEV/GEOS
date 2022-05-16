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
 * @file ProppantTransportExtrinsicData.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_PROPPANT_PROPPANTTRANSPORTEXTRINSICDATA_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_PROPPANT_PROPPANTTRANSPORTEXTRINSICDATA_HPP_

#include "mesh/ExtrinsicMeshData.hpp"

namespace geosx
{
/**
 * A scope for extrinsic mesh data traits.
 */
namespace extrinsicMeshData
{

namespace proppant
{

EXTRINSIC_MESH_DATA_TRAIT( proppantConcentration,
                           "proppantConcentration",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Proppant concentration" );

EXTRINSIC_MESH_DATA_TRAIT( proppantConcentration_n,
                           "proppantConcentration_n",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Proppant concentration at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( componentConcentration,
                           "componentConcentration",
                           array2d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Component concentration" );

EXTRINSIC_MESH_DATA_TRAIT( componentConcentration_n,
                           "componentConcentration_n",
                           array2d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Component concentration at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( bcComponentConcentration,
                           "bcComponentConcentration",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Boundary component concentration" );

EXTRINSIC_MESH_DATA_TRAIT( cellBasedFlux,
                           "cellBasedFlux",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Cell-based flux" );

EXTRINSIC_MESH_DATA_TRAIT( isProppantBoundary,
                           "isProppantBoundary",
                           array1d< integer >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Flag denoting the proppant boundary" );

EXTRINSIC_MESH_DATA_TRAIT( isProppantMobile,
                           "isProppantMobile",
                           array1d< integer >,
                           1.0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Flag indicating whether proppant is mobile" );

EXTRINSIC_MESH_DATA_TRAIT( componentDensity_n,
                           "componentDensity_n",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Component density at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( proppantPackVolumeFraction,
                           "proppantPackVolumeFraction",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Proppant pack volume fraction" );

EXTRINSIC_MESH_DATA_TRAIT( proppantExcessPackVolume,
                           "proppantExcessPackVolume",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Proppant pack volume fraction" );

EXTRINSIC_MESH_DATA_TRAIT( proppantLiftFlux,
                           "proppantLiftFlux",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Proppant lift flux" );


}

}

}

#endif // GEOSX_PHYSICSSOLVERS_FLUIDFLOW_PROPPANT_PROPPANTTRANSPORTEXTRINSICDATA_HPP_
