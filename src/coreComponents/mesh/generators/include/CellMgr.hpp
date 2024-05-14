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

#ifndef GEOS_CELLMGR_HPP
#define GEOS_CELLMGR_HPP

#include "LineBlockABC.hpp"

#include "dataRepository/Group.hpp"

#include "common/DataTypes.hpp"


namespace geos::generators
{

class CellMgr
{
public:
  /**
   * @brief Returns a group containing the cell blocks as @p CellBlockABC instances.
   * @return Mutable reference to the cell blocks group.
   *
   * @note It should probably be better not to expose a non-const accessor here.
   */
  virtual dataRepository::Group & getCellBlocks() = 0;

  /**
   * @brief Returns a group containing the face blocks as @p FaceBlockABC instances.
   * @return Mutable reference to the face blocks group.
   *
   * @note It should probably be better not to expose a non-const accessor here.
   */
  virtual dataRepository::Group & getFaceBlocks() = 0;

  /**
   * @brief Returns LineBlockABC corresponding to the given identifier
   * @param name the name of the required LineBlockABC
   * @return The LineBlockABC associated with the given name
   */
  virtual LineBlockABC const & getLineBlock( string name ) const = 0;

  /**
   * @brief Returns a group containing the cell blocks as CellBlockABC instances
   * @return Const reference to the Group instance.
   */
  virtual const dataRepository::Group & getCellBlocks() const = 0;

  /**
   * @brief Returns a group containing the face blocks as FaceBlockABC instances
   * @return Const reference to the Group instance.
   */
  virtual const dataRepository::Group & getFaceBlocks() const = 0;
};

}

#endif //GEOS_CELLMGR_HPP
