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
 * @file SolidMechanicsExtrinsicData.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSEXTRINSICDATA_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSEXTRINSICDATA_HPP_

#include "common/DataLayouts.hpp"
#include "mesh/ExtrinsicMeshData.hpp"

namespace geosx
{
/**
 * A scope for extrinsic mesh data traits.
 */
namespace extrinsicMeshData
{

namespace solidMechanics
{

using array2dLayoutTotalDisplacement = array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM >;
using arrayView2dLayoutTotalDisplacement = arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD >;
using arrayViewConst2dLayoutTotalDisplacement = arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD >;

using array2dLayoutIncrDisplacement = array2d< real64, nodes::INCR_DISPLACEMENT_PERM >;
using arrayView2dLayoutIncrDisplacement = arrayView2d< real64, nodes::INCR_DISPLACEMENT_USD >;
using arrayViewConst2dLayoutIncrDisplacement = arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD >;

using array2dLayoutVelocity = array2d< real64, nodes::VELOCITY_PERM >;
using arrayView2dLayoutVelocity = arrayView2d< real64, nodes::VELOCITY_USD >;
using arrayViewConst2dLayoutVelocity = arrayView2d< real64 const, nodes::VELOCITY_USD >;

using array2dLayoutAcceleration = array2d< real64, nodes::ACCELERATION_PERM >;
using arrayView2dLayoutAcceleration = arrayView2d< real64, nodes::ACCELERATION_USD >;
using arrayViewConst2dLayoutAcceleration = arrayView2d< real64 const, nodes::ACCELERATION_USD >;


EXTRINSIC_MESH_DATA_TRAIT( totalDisplacement,
                           "totalDisplacement",
                           array2dLayoutTotalDisplacement,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Total displacements at the nodes" );

EXTRINSIC_MESH_DATA_TRAIT( incrementalDisplacement,
                           "incrementalDisplacement",
                           array2dLayoutIncrDisplacement,
                           0,
                           LEVEL_3,
                           WRITE_AND_READ,
                           "Incremental displacements for the current time step on the nodes" );

EXTRINSIC_MESH_DATA_TRAIT( velocity,
                           "velocity",
                           array2dLayoutVelocity,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Current velocity on the nodes" );

EXTRINSIC_MESH_DATA_TRAIT( acceleration,
                           "acceleration",
                           array2dLayoutAcceleration,
                           0,
                           LEVEL_1,
                           WRITE_AND_READ,
                           "Current acceleration on the nodes. This array also is used "
                           "to hold the summation of nodal forces resulting from the governing equations" );

EXTRINSIC_MESH_DATA_TRAIT( externalForce,
                           "externalForce",
                           array2d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "External forces on the nodes. This includes any boundary"
                           " conditions as well as coupling forces such as hydraulic forces" );

EXTRINSIC_MESH_DATA_TRAIT( mass,
                           "mass",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Mass on the nodes" );

EXTRINSIC_MESH_DATA_TRAIT( velocityTilde,
                           "velocityTilde",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Velocity predictors on the nodes" );

EXTRINSIC_MESH_DATA_TRAIT( uhatTilde,
                           "uhatTilde",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Incremental displacement predictors on the nodes" );

EXTRINSIC_MESH_DATA_TRAIT( contactForce,
                           "contactForce",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Contact force" );

}

}

}

#endif // GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSEXTRINSICDATA_HPP_
