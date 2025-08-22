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
 * @file ScatterDataProvider.hpp
 */

#ifndef GEOS_FILEIO_TIMEHISTORY_SCATTERDATAPROVIDER_HPP_
#define GEOS_FILEIO_TIMEHISTORY_SCATTERDATAPROVIDER_HPP_

#include "common/DataTypes.hpp"

namespace geos
{

/**
 * @brief Interface for classes that provide scattered data (coordinates and values)
 */
class ScatterDataProvider
{
public:

  /// Virtual destructor
  virtual ~ScatterDataProvider() = default;

  /**
   * @brief Get the number of scatter points
   * @return Number of points where scattered data is computed
   */
  virtual localIndex getNumScatterPoints() const = 0;

  /**
   * @brief Get the scattered data values
   * @return Array of data values at scatter points
   */
  virtual array1d< real64 > const & getScatterData() const = 0;

  /**
   * @brief Get the coordinates of scatter points
   * @return 2D array of coordinates [nPoints x nDim]
   */
  virtual array2d< real64 > const & getScatterCoordinates() const
  {
    static array2d< real64 > empty;
    return empty;
  }

  /**
   * @brief Get optional metadata for scatter points
   * @return String array with metadata for each point (e.g., station names)
   */
  virtual string_array const & getScatterMetadata() const
  {
    static string_array empty;
    return empty;
  }

protected:

  /// Protected constructor
  ScatterDataProvider() = default;

  /// Non-copyable
  ScatterDataProvider( ScatterDataProvider const & ) = delete;
  ScatterDataProvider & operator=( ScatterDataProvider const & ) = delete;

  /// Non-movable
  ScatterDataProvider( ScatterDataProvider && ) = delete;
  ScatterDataProvider & operator=( ScatterDataProvider && ) = delete;
};

} // namespace geos

#endif /* GEOS_FILEIO_TIMEHISTORY_SCATTERDATAPROVIDER_HPP_ */
