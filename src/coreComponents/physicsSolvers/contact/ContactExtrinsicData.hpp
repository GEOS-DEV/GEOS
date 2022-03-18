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
 * @file ContactExtrinsicData.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_CONTACT_CONTACTEXTRINSICDATA_HPP_
#define GEOSX_PHYSICSSOLVERS_CONTACT_CONTACTEXTRINSICDATA_HPP_

#include "mesh/ExtrinsicMeshData.hpp"

namespace geosx
{
/**
 * A scope for extrinsic mesh data traits.
 */
namespace extrinsicMeshData
{

namespace contact
{

EXTRINSIC_MESH_DATA_TRAIT( dispJump,
                           "displacementJump",
                           array2d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Displacement jump vector" );

EXTRINSIC_MESH_DATA_TRAIT( deltaDispJump,
                           "deltaDisplacementJump",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Delta displacement jump vector" );

EXTRINSIC_MESH_DATA_TRAIT( oldDispJump,
                           "oldDisplacementJump",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Displacement jump vector at the previous time-step" );

EXTRINSIC_MESH_DATA_TRAIT( traction,
                           "traction",
                           array2d< real64 >,
                           0,
                           LEVEL_0,
                           NO_WRITE,
                           "Fracture traction vector" );

EXTRINSIC_MESH_DATA_TRAIT( deltaTraction,
                           "traction",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "An array that holds the traction increments on the fracture." );


EXTRINSIC_MESH_DATA_TRAIT( dTraction_dJump,
                           "dTraction_dJump",
                           array3d< real64 >,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of the traction w.r.t. the displacement jump." );

EXTRINSIC_MESH_DATA_TRAIT( dTraction_dPressure,
                           "dTraction_dPressure",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           NO_WRITE,
                           "Derivative of the traction w.r.t. to the fluid pressure." );

}

}

}

#endif // GEOSX_PHYSICSSOLVERS_CONTACT_CONTACTEXTRINSICDATA_HPP_
