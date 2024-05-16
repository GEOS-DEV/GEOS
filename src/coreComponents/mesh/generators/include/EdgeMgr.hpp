/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2020-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_EDGEMGR_HPP
#define GEOS_EDGEMGR_HPP

#include "common/DataTypes.hpp"

namespace geos::generators
{

class EdgeMgr
{
public:
  /**
   * @brief Total number of edges across all the cell blocks.
   * @return The total number of edges.
   */
  [[nodiscard]] virtual localIndex numEdges() const = 0;

  /**
   * @brief Returns the edge to nodes mapping.
   * @return A 1 to 2 relationship. The result is meant to have size (numEdges, 2).
   */
  [[nodiscard]] virtual array2d< localIndex > getEdgeToNodes() const = 0;

  /**
   * @brief Returns the edge to faces mapping.
   * @return A one to many relationship.
   */
  [[nodiscard]] virtual ArrayOfArrays< localIndex > getEdgeToFaces() const = 0;

  /**
   * @brief Returns the ghost rank mapping. Index is an edge index local to the MPI rank.
   * @return A @c numEdges length array.
   */
  [[nodiscard]] virtual array1d< integer > getGhostRank() const = 0;

  // TODO Use inheritance?
  [[nodiscard]] virtual array1d< globalIndex > getLocalToGlobal() const = 0;
};

}

#endif //GEOS_EDGEMGR_HPP
