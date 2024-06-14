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

#ifndef GEOS_MESHMAPPINGS_HPP
#define GEOS_MESHMAPPINGS_HPP

#include "CellMgr.hpp"
#include "EdgeMgr.hpp"
#include "FaceMgr.hpp"
#include "NodeMgr.hpp"

#include "dataRepository/Group.hpp"

namespace geos::generators
{

class MeshMappings : public dataRepository::Group
{
public:
  /**
   * @brief Constructor
   * @param name The name of this Group.
   * @param parent The parent Group.
   */
  MeshMappings( string const & name, Group * const parent ):
  Group( name, parent )
  {
    // Left blank
  }

  /**
   * @brief Extra space for node to edges mapping.
   * @return Number of extra values as an integer.
   */
  static constexpr localIndex edgeMapExtraSpacePerNode()
  { return 8; }

  /**
   * @brief Extra space for node to faces mapping.
   * @return Number of extra values as an integer.
   */
  static constexpr localIndex faceMapExtraSpacePerNode()
  { return 8; }

  /**
   * @brief Extra space for node to elements mapping.
   * @return Number of extra values as an integer.
   */
  static constexpr localIndex elemMapExtraSpacePerNode()
  { return 8; }

  /**
   * @brief Extra space for extra nodes.
   * @return Number of extra values as an integer.
   */
  static constexpr localIndex nodeMapExtraSpacePerFace()
  { return 4; }

  /**
   * @brief Extra space for extra faces.
   * @return Number of extra values as an integer.
   */
  static constexpr localIndex edgeMapExtraSpacePerFace()
  { return 4; }

  /**
   * @brief Extra space for extra edges.
   * @return Number of extra values as an integer.
   */
  static constexpr localIndex faceMapExtraSpacePerEdge()
  { return 4; }

  virtual CellMgr const & getCellMgr() const = 0 ;
  virtual EdgeMgr const & getEdgeMgr() const = 0 ;
  virtual FaceMgr const & getFaceMgr() const = 0 ;
  virtual NodeMgr const & getNodeMgr() const = 0 ;

  virtual std::set< integer > const & getNeighbors() const = 0;
};

}

#endif //GEOS_MESHMAPPINGS_HPP
