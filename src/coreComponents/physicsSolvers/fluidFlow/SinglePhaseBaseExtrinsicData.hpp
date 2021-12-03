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
 * @file SinglePhaseBaseExtrinsicData.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASEEXTRINSICDATA_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASEEXTRINSICDATA_HPP_

#include "mesh/ExtrinsicMeshData.hpp"

namespace geosx
{
/**
 * A scope for extrinsic mesh data traits.
 */
namespace extrinsicMeshData
{

EXTRINSIC_MESH_DATA_TRAIT( mobility,
                           "mobility",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Mobility" );

EXTRINSIC_MESH_DATA_TRAIT( dMobility_dPressure,
                           "dMobility_dPressure",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of mobility with respect to pressure" );

EXTRINSIC_MESH_DATA_TRAIT( densityOld,
                           "densityOld",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Density at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( facePressure,
                           "facePressure",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Face pressure" );

}

}

#endif // GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASEEXTRINSICDATA_HPP_
