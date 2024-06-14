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

#ifndef GEOS_NODEMGR_HPP
#define GEOS_NODEMGR_HPP


#include "../CellBlockUtilities.hpp" // TODO At least part of it should become public.

#include "GhostExchange.hpp"

#include "common/DataTypes.hpp"

namespace geos::generators
{

class NodeMgr: public GhostExchange
{
public:
  /**
 * @brief Total number of nodes across all the cell blocks.
 * @return The total number of nodes.
 *
 * Nodes shared by multiple cell blocks are counted only once.
 */
  [[nodiscard]] virtual localIndex numNodes() const = 0;

  /**
 * @brief Returns the node coordinates in a (numNodes, 3) 2d array.
 * @return A const view to the array.
 */
  [[nodiscard]] virtual array2d< real64, nodes::REFERENCE_POSITION_PERM > getNodePositions() const = 0;

  /**
   * @brief Returns the node to edges mapping.
   * @return The one to many relationship.
   */
  [[nodiscard]] virtual ArrayOfArrays< localIndex > getNodeToEdges() const = 0;

  /**
   * @brief Returns the face to nodes mappings.
   * @return The one to many relationship.
   */
  [[nodiscard]] virtual ArrayOfArrays< localIndex > getNodeToFaces() const = 0;

  /**
   * @brief Returns the node to elements mapping.
   * @return A one to many relationship.
   */
  [[nodiscard]] virtual ToCellRelation< ArrayOfArrays< localIndex > > getNodeToElements() const = 0;

  /**
   * @brief Returns the node sets. Key of the map is the name of the set.
   * @return A reference to constant map.
   */
  [[nodiscard]] virtual std::map< string, SortedArray< localIndex > > const & getNodeSets() const = 0;
};

}

#endif //GEOS_NODEMGR_HPP
