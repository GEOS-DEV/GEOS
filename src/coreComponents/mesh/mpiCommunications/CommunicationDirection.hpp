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
 * @file CommunicationDirection.hpp
 */

#ifndef GEOS_MESH_MPICOMMUNICATIONS_COMMUNICATIONDIRECTION_HPP_
#define GEOS_MESH_MPICOMMUNICATIONS_COMMUNICATIONDIRECTION_HPP_

namespace geos
{

/**
 * @brief Direction to use when packing and unpacking persistent ghost-sync lists.
 *
 * The default GEOS synchronization direction is OwnerToGhost:
 *   - pack ghostsToSend
 *   - unpack into ghostsToReceive
 *
 * Some algorithms, including MPM grid accumulation, need the reverse direction:
 *   - pack ghostsToReceive
 *   - unpack into ghostsToSend
 *
 * Keeping the direction explicit avoids physically swapping the stored
 * ghostsToSend/ghostsToReceive arrays.
 */
enum class CommunicationDirection
{
  OwnerToGhost,
  GhostToOwner
};

} /* namespace geos */

#endif /* GEOS_MESH_MPICOMMUNICATIONS_COMMUNICATIONDIRECTION_HPP_ */
