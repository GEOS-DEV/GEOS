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

#ifndef GEOS_FACEMGR_HPP
#define GEOS_FACEMGR_HPP

#include "../CellBlockUtilities.hpp" // TODO At least part of it should become public.

#include "common/DataTypes.hpp"

namespace geos::generators
{

class FaceMgr
{
public:
  /**
   * @brief Total number of faces across all the cell blocks.
   * @return The total number of faces.
   */
  virtual localIndex numFaces() const = 0;

  /**
   * @brief Returns the face to nodes mapping.
   * @return The one to many relationship.
   */
  virtual ArrayOfArrays< localIndex > getFaceToNodes() const = 0;

  /**
   * @brief Returns the face to edges mapping.
   * @return A one to many relationship.
   */
  virtual ArrayOfArrays< localIndex > getFaceToEdges() const = 0;

  /**
   * @brief Returns the face to elements mapping.
   * @return A 1 to 2 relationship. The result is meant to have size (numFaces, 2).
   *
   * In case the face only belongs to one single element, the second value of the table is -1.
   */
  virtual ToCellRelation< array2d< localIndex > > getFaceToElements() const = 0;
};

}

#endif //GEOS_FACEMGR_HPP
