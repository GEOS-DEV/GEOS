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
 * @file PerforationExtrinsicData.hpp
 */

#ifndef GEOSX_MESH_PERFORATIONEXTRINSICDATA_HPP_
#define GEOSX_MESH_PERFORATIONEXTRINSICDATA_HPP_

#include "mesh/ExtrinsicMeshData.hpp"

namespace geosx
{
/**
 * A scope for extrinsic mesh data traits.
 */
namespace extrinsicMeshData
{

namespace perforation
{

///@cond DO_NOT_DOCUMENT
EXTRINSIC_MESH_DATA_TRAIT( reservoirElementRegion,
                           "reservoirElementRegion",
                           array1d< localIndex >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "For each perforation, elementRegion index of the perforated element" );

EXTRINSIC_MESH_DATA_TRAIT( reservoirElementSubRegion,
                           "reservoirElementSubregion",
                           array1d< localIndex >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "For each perforation, elementSubRegion index of the perforated element" );

EXTRINSIC_MESH_DATA_TRAIT( reservoirElementIndex,
                           "reservoirElementIndex",
                           array1d< localIndex >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "For each perforation, element index of the perforated element" );

EXTRINSIC_MESH_DATA_TRAIT( wellElementIndex,
                           "wellElementIndex",
                           array1d< localIndex >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "For each perforation, index of the well element" );

EXTRINSIC_MESH_DATA_TRAIT( wellTransmissibility,
                           "wellTransmissibility",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "For each perforation, well transmissibility" );

EXTRINSIC_MESH_DATA_TRAIT( location,
                           "location",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "For each perforation, physical location (x,y,z coordinates)" );

///@endcond DO_NOT_DOCUMENT
}

}

}

#endif // GEOSX_MESH_PERFORATIONEXTRINSICDATA_HPP_
