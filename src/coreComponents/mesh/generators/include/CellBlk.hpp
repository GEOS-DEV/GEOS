/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_CELLBLK_HPP
#define GEOS_CELLBLK_HPP

#include "GhostExchange.hpp"

#include "mesh/ElementType.hpp"
#include "common/DataTypes.hpp"


namespace geos::generators
{

class CellBlk: public generators::GhostExchange
{
public:

  /**
   * @brief Get the type of element in this subregion.
   * @return a string specifying the type of element in this subregion
   *
   * See class FiniteElementBase for possible element type.
   */
  [[nodiscard]] virtual ElementType getElementType() const = 0;

  /**
   * @brief Get the number of nodes per element.
   * @return number of nodes per element
   */
  [[nodiscard]] virtual localIndex numNodesPerElement() const = 0;

  /**
   * @brief Get the number of edges per element.
   * @return number of edges per element
   */
  [[nodiscard]] virtual localIndex numEdgesPerElement() const = 0;

  /**
   * @brief Get the number of faces per element.
   * @return number of faces per element
   */
  [[nodiscard]] virtual localIndex numFacesPerElement() const = 0;

  /**
   * @brief Get the number of elements.
   * @return number of elements in the cell block
   */
  [[nodiscard]] virtual localIndex numElements() const = 0;

  /**
   * @brief Get the element-to-nodes map.
   * @return The mapping relationship as a 2d-array.
   */
  [[nodiscard]] virtual array2d< localIndex, cells::NODE_MAP_PERMUTATION > getElemToNodes() const = 0;

  /**
   * @brief Get the element-to-edges map.
   * @return The mapping relationship as a 2d-array.
   */
  [[nodiscard]] virtual array2d< localIndex > getElemToEdges() const = 0;

  /**
   * @brief Get the element-to-faces map.
   * @return The mapping relationship as a 2d-array.
   */
  [[nodiscard]] virtual array2d< localIndex > getElemToFaces() const = 0;
};

} // end of namespace

#endif // include guard
