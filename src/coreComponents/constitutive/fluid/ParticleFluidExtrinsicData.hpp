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
 * @file ParticleFluidExtrinsicData.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PARTICLEFLUIDEXTRINSICDATA_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PARTICLEFLUIDEXTRINSICDATA_HPP_

#include "mesh/ExtrinsicMeshData.hpp"

namespace geosx
{

namespace extrinsicMeshData
{

namespace particlefluid
{

EXTRINSIC_MESH_DATA_TRAIT( settlingFactor,
                           "settlingFactor",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Settling factor" );

EXTRINSIC_MESH_DATA_TRAIT( dSettlingFactor_dPressure,
                           "dSettlingFactor_dPressure",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Derivative of settling factor with respect to pressure" );

EXTRINSIC_MESH_DATA_TRAIT( dSettlingFactor_dProppantConcentration,
                           "dSettlingFactor_dProppantConcentration",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Derivative of settling factor with respect to proppant concentration" );

EXTRINSIC_MESH_DATA_TRAIT( dSettlingFactor_dComponentConcentration,
                           "dSettlingFactor_dComponentConcentration",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Derivative of settling factor with respect to component concentration" );

EXTRINSIC_MESH_DATA_TRAIT( collisionFactor,
                           "collisionFactor",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Collision factor" );

EXTRINSIC_MESH_DATA_TRAIT( dCollisionFactor_dProppantConcentration,
                           "dCollisionFactor_dProppantConcentration",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Derivative of collision factor with respect to proppant concentration" );

EXTRINSIC_MESH_DATA_TRAIT( proppantPackPermeability,
                           "proppantPackPermeability",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Proppant pack permeability" );


}

}

}

#endif // GEOSX_CONSTITUTIVE_FLUID_PARTICLEFLUIDEXTRINSICDATA_HPP_
