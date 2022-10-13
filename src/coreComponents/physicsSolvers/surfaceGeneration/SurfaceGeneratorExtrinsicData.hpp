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
 * @file SurfaceGeneratorExtrinsicData.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SURFACEGENERATION_SURFACEGENERATOREXTRINSICDATA_HPP_
#define GEOSX_PHYSICSSOLVERS_SURFACEGENERATION_SURFACEGENERATOREXTRINSICDATA_HPP_

#include "common/DataLayouts.hpp"
#include "mesh/ExtrinsicMeshData.hpp"

namespace geosx
{
/**
 * A scope for extrinsic mesh data traits.
 */
namespace extrinsicMeshData
{

namespace surfaceGeneration
{

EXTRINSIC_MESH_DATA_TRAIT( degreeFromCrack,
                           "degreeFromCrack",
                           array1d< integer >,
                           -1,
                           LEVEL_1,
                           WRITE_AND_READ,
                           "Distance to the crack in terms of topological distance. "
                           "(i.e. how many nodes are along the path to the closest "
                           "node that is on the crack surface." );

EXTRINSIC_MESH_DATA_TRAIT( degreeFromCrackTip,
                           "degreeFromCrackTip",
                           array1d< integer >,
                           100000,
                           LEVEL_1,
                           WRITE_AND_READ,
                           "Distance to the crack tip in terms of topological distance. "
                           "(i.e. how many nodes are along the path to the closest "
                           "node that is on the crack surface." );

EXTRINSIC_MESH_DATA_TRAIT( SIFNode,
                           "SIFNode",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Calculated Stress Intensity Factor on the node." );

EXTRINSIC_MESH_DATA_TRAIT( ruptureRate,
                           "ruptureRate",
                           array1d< real64 >,
                           1.0e99,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Rate of rupture in terms of number of objects split per time." );

EXTRINSIC_MESH_DATA_TRAIT( SIF_I,
                           "SIF_I",
                           array1d< real64 >,
                           -1,
                           LEVEL_1,
                           WRITE_AND_READ,
                           "Calculated mode 1 Stress Intensity Factor on the node." );

EXTRINSIC_MESH_DATA_TRAIT( SIF_II,
                           "SIF_II",
                           array1d< real64 >,
                           -1,
                           LEVEL_1,
                           WRITE_AND_READ,
                           "Calculated mode 2 Stress Intensity Factor on the node." );

EXTRINSIC_MESH_DATA_TRAIT( SIF_III,
                           "SIF_III",
                           array1d< real64 >,
                           -1,
                           LEVEL_1,
                           WRITE_AND_READ,
                           "Calculated mode 3 Stress Intensity Factor on the node." );

EXTRINSIC_MESH_DATA_TRAIT( ruptureState,
                           "ruptureState",
                           array1d< integer >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Rupture state of the face: \n 0=not ready for rupture \n 1=ready for rupture \n 2=ruptured." );

EXTRINSIC_MESH_DATA_TRAIT( SIFonFace,
                           "SIFonFace",
                           array1d< real64 >,
                           1,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Calculated Stress Intensity Factor on the face." );


/// The template type T for registration of a container<T>.

EXTRINSIC_MESH_DATA_TRAIT( K_IC,
                           "K_IC",
                           array2d< real64 >,
                           1e99,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Critical Stress Intensity Factor :math:`K_{IC}` in the plane of the face." );

EXTRINSIC_MESH_DATA_TRAIT( K_IC_00,
                           "K_IC_00",
                           array1d< real64 >,
                           -1,
                           NOPLOT,
                           WRITE_AND_READ,
                           ":math:`K_{IC}` on 0-plane, in 0-direction." );

EXTRINSIC_MESH_DATA_TRAIT( K_IC_01,
                           "K_IC_01",
                           array1d< real64 >,
                           -1,
                           NOPLOT,
                           WRITE_AND_READ,
                           ":math:`K_{IC}` on 0-plane, in 1-direction." );

EXTRINSIC_MESH_DATA_TRAIT( K_IC_02,
                           "K_IC_02",
                           array1d< real64 >,
                           -1,
                           NOPLOT,
                           WRITE_AND_READ,
                           ":math:`K_{IC}` on 0-plane, in 2-direction." );

EXTRINSIC_MESH_DATA_TRAIT( K_IC_10,
                           "K_IC_10",
                           array1d< real64 >,
                           -1,
                           NOPLOT,
                           WRITE_AND_READ,
                           ":math:`K_{IC}` on 1-plane, in 0-direction." );

EXTRINSIC_MESH_DATA_TRAIT( K_IC_11,
                           "K_IC_11",
                           array1d< real64 >,
                           -1,
                           NOPLOT,
                           WRITE_AND_READ,
                           ":math:`K_{IC}` on 1-plane, in 1-direction." );

EXTRINSIC_MESH_DATA_TRAIT( K_IC_12,
                           "K_IC_12",
                           array1d< real64 >,
                           -1,
                           NOPLOT,
                           WRITE_AND_READ,
                           ":math:`K_{IC}` on 1-plane, in 2-direction." );

EXTRINSIC_MESH_DATA_TRAIT( K_IC_20,
                           "K_IC_20",
                           array1d< real64 >,
                           -1,
                           NOPLOT,
                           WRITE_AND_READ,
                           ":math:`K_{IC}` on 2-plane, in 0-direction." );

EXTRINSIC_MESH_DATA_TRAIT( K_IC_21,
                           "K_IC_21",
                           array1d< real64 >,
                           -1,
                           NOPLOT,
                           WRITE_AND_READ,
                           ":math:`K_{IC}` on 2-plane, in 1-direction." );

EXTRINSIC_MESH_DATA_TRAIT( K_IC_22,
                           "K_IC_22",
                           array1d< real64 >,
                           -1,
                           NOPLOT,
                           WRITE_AND_READ,
                           ":math:`K_{IC}` on 2-plane, in 2-direction." );

EXTRINSIC_MESH_DATA_TRAIT( primaryCandidateFace,
                           "primaryCandidateFace",
                           array1d< localIndex >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "??" );

EXTRINSIC_MESH_DATA_TRAIT( isFaceSeparable,
                           "isFaceSeparable",
                           array1d< integer >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "A flag to mark if the face is separable." );

}

}

}

#endif // GEOSX_PHYSICSSOLVERS_SURFACEGENERATION_SURFACEGENERATOREXTRINSICDATA_HPP_
