/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CommID.hpp
 */

#ifndef GEOSX_MESH_MPICOMMUNICATIONS_COMMID_HPP_
#define GEOSX_MESH_MPICOMMUNICATIONS_COMMID_HPP_

#include <set>

namespace geosx
{

/**
 * This class provides management of a set of integers which can be used to
 *
 *
 */
class CommID
{
public:
  /**
   * Constructor
   * @param freeIDs
   */
  CommID( std::set< int > & freeIDs );

  /**
   * Destructor
   */
  ~CommID();

  /// default copy constructor
  CommID( CommID const & ) = default;

  /**
   * Move constructor
   * @param src The source to move data from.
   */
  CommID( CommID && src );

  /// deleted copy assignment operator
  CommID & operator=( CommID const & ) = delete;

  /// deleted move assignment operator
  CommID & operator=( CommID && ) = delete;

  /// user defined conversion operator to int
  constexpr operator int()
  { return m_id; }

private:
  /// Reference to the set of free ID's
  std::set< int > & m_freeIDs;

  /// current ID
  int m_id = -1;
};

} /* namespace geosx */

#endif /* GEOSX_MESH_MPICOMMUNICATIONS_COMMID_HPP_ */
