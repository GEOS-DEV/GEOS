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
 * @file WellSolverBaseExtrinsicData.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLSOLVERBASEEXTRINSICDATA_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLSOLVERBASEEXTRINSICDATA_HPP_

#include "mesh/ExtrinsicMeshData.hpp"

namespace geosx
{
/**
 * A scope for extrinsic mesh data traits.
 */
namespace extrinsicMeshData
{

namespace well
{

EXTRINSIC_MESH_DATA_TRAIT( pressure,
                           "pressure",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Pressure" );

EXTRINSIC_MESH_DATA_TRAIT( pressureOld,
                           "pressureOld",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Pressure at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( gravityCoefficient,
                           "gravityCoefficient",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Gravity coefficient (dot product of gravity acceleration by gravity vector)" );

}

}

}

#endif // GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLSOLVERBASEEXTRINSICDATA_HPP_
