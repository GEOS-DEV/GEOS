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

/**
 * @file CommID.hpp
 */

#ifndef GEOSX_MESH_MPICOMMUNICATIONS_COMMID_HPP_
#define GEOSX_MESH_MPICOMMUNICATIONS_COMMID_HPP_

#include <set>

namespace geosx
{


class CommID
{
public:
  CommID( std::set< int > & freeIDs );

  CommID( CommID && src );

  ~CommID();

  CommID( CommID const & ) = delete;
  CommID & operator=( CommID const & ) = delete;
  CommID & operator=( CommID && ) = delete;

  constexpr operator int()
  { return m_id; }

private:
  std::set< int > & m_freeIDs;
  int m_id = -1;
};

} /* namespace geosx */

#endif /* GEOSX_MESH_MPICOMMUNICATIONS_COMMID_HPP_ */
