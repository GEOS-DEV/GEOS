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

#ifndef GEOS_GHOSTEXCHANGE_HPP
#define GEOS_GHOSTEXCHANGE_HPP

#include "common/DataTypes.hpp"

namespace geos::generators
{

class GhostExchange
{
public:
  /**
   * @brief Get global to local map for the contained geometrical quantity (nodes, edges...).
   * @return The mapping relationship as a map (global indexing is not contiguous).
   */
  [[nodiscard]] virtual unordered_map< globalIndex, localIndex > getGlobalToLocal() const = 0;

  /**
   * @brief Get the list of geometrical quantity (nodes, edges...) that need to be sent to neighbors.
   * @return A mapping with the neighbor rank as key and the list of entities as value.
   */
  [[nodiscard]] virtual std::map< integer, array1d< localIndex > > getSend() const = 0;

  /**
   * @brief Get the list of geometrical quantity (nodes, edges...) that will be received from neighbors.
   * @return A mapping with the neighbor rank as key and the list of entities as value.
   */
  [[nodiscard]] virtual std::map< integer, array1d< localIndex > > getRecv() const = 0;
};

}

#endif //GEOS_GHOSTEXCHANGE_HPP
