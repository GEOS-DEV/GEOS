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
 * @file SinglePhaseWellExtrinsicData.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELLEXTRINSICDATA_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELLEXTRINSICDATA_HPP_

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

EXTRINSIC_MESH_DATA_TRAIT( connectionRate,
                           "connectionRate",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Connection rate" );

EXTRINSIC_MESH_DATA_TRAIT( deltaConnectionRate,
                           "deltaConnectionRate",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Accumulated connection rate updates" );

EXTRINSIC_MESH_DATA_TRAIT( densityOld,
                           "densityOld",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Density at the previous converged time step" );

EXTRINSIC_MESH_DATA_TRAIT( perforationRate,
                           "perforationRate",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Perforation rate" );

EXTRINSIC_MESH_DATA_TRAIT( dPerforationRate_dPres,
                           "dPerforationRate_dPres",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Derivative of perforation rate with respect to pressure" );

}

}

}

#endif // GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELLEXTRINSICDATA_HPP_
