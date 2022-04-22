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

#ifndef GEOSX_MESH_MPICOMMUNICATIONS_MPICOMMUNICATIONUTILS_HPP_
#define GEOSX_MESH_MPICOMMUNICATIONS_MPICOMMUNICATIONUTILS_HPP_

#include "common/DataTypes.hpp"

namespace geosx
{

struct SyncFieldsID
{
  FieldLocation location;
  std::vector< string > fieldNames;
  std::vector< string > regionNames;

  SyncFieldsID( FieldLocation const location_, std::vector< string > const & fields, arrayView1d< string const > const & regions ):
    location( location_ ),
    fieldNames( fields ),
    regionNames( regions.begin(), regions.end() )
  {}

  SyncFieldsID( FieldLocation const location_, std::vector< string > const & fields, std::vector< string > const & regions ):
    location( location_ ),
    fieldNames( fields ),
    regionNames( regions )
  {}

  SyncFieldsID( FieldLocation const location_, std::vector< string > const & fields ):
    SyncFieldsID( location_, fields, std::vector< string >{} )
  {}

  array1d< string > getFieldNames() const
  {
    array1d< string > tmp;
    tmp.reserve( fieldNames.size() );

    for( auto const & field : fieldNames )
    {
      tmp.emplace_back( field );
    }

    return tmp;
  }
};

} /* namespace geosx */

#endif /* GEOSX_MESH_MPICOMMUNICATIONS_MPICOMMUNICATIONUTILS_HPP_ */
