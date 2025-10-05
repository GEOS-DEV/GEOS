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
 * @file GeometricPartitioner.hpp
 * @brief Defines the abstract base class for partitioners based on spatial decomposition.
 */

#ifndef GEOS_PARTITIONER_GEOMETRICPARTITIONER_HPP_
#define GEOS_PARTITIONER_GEOMETRICPARTITIONER_HPP_

#include "PartitionerBase.hpp"

namespace geos
{

/**
 * @class GeometricPartitioner
 * @brief An abstract class for partitioners that operate on spatial coordinates.
 *
 * This class provides an interface for any partitioner that divides the domain
 * geometrically (e.g., Cartesian grids, octrees, etc.).
 */
class GeometricPartitioner : public PartitionerBase
{
public:
  /// Inherit the constructor from PartitionerBase.
  using PartitionerBase::PartitionerBase;

  virtual ~GeometricPartitioner() = default;

  /**
   * @brief Initializes the partitioner with the global domain's bounding box.
   *
   * This method provides the overall spatial extent of the problem.
   *
   * @param globalMin The minimum coordinates (x, y, z) of the global bounding box.
   * @param globalMax The maximum coordinates (x, y, z) of the global bounding box.
   */
  virtual void initializeDomain( const R1Tensor & globalMin, const R1Tensor & globalMax ) = 0;

  /**
   * @brief Returns the destination rank for a given coordinate.
   *
   * Determine which partition an object belongs to based on its position.
   *
   * @param coords The coordinates of the object to be placed.
   * @return The rank ID of the destination partition, or -1 if the coordinate
   *         is outside the entire simulation domain.
   */
  virtual int getDestinationRank( const R1Tensor & coords ) const = 0;
};

} // namespace geos

#endif // GEOS_PARTITIONER_GEOMETRICPARTITIONER_HPP_
