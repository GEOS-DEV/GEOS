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
  array1d< string > regionNames;
  array1d< string > fieldNames;

  SyncFieldsID( FieldLocation const location, arrayView1d< string const > const & regions, arrayView1d< string const > const & fields ):
    location( location )
  {
    fillRegions( regions );
    fillFields( fields );
  }

  SyncFieldsID( FieldLocation const location, arrayView1d< string const > const & regions, std::initializer_list< string > fields ):
    location( location )
  {
    fillRegions( regions );
    fillFields(fields);
  }

  SyncFieldsID( FieldLocation const location, std::vector< string > const & regions, std::initializer_list< string > fields ):
    location( location )
  {
    fillRegions( regions );
    fillFields( fields );  
  }

  SyncFieldsID( FieldLocation const location, string const & region, std::initializer_list< string > fields ):
    location( location )
  {
    regionNames.emplace_back(region);
    fillFields(fields);
  }

  template< typename T >
  void fillRegions( T const & regions )
  {
    regionNames.resize( regions.size());
    for( integer i = 0; i < regionNames.size(); i++ )
    {
      regionNames[i] = regions[i];
    }
  }

  void fillFields( arrayView1d<string const> const & fields )
  {
    for( auto field : fields )
    {
      fieldNames.emplace_back(field);
    }
  }

  void fillFields( std::initializer_list< string > fields )
  {
    for( auto field : fields )
    {
      fieldNames.emplace_back(field);
    }
  }
};

} /* namespace geosx */

#endif /* GEOSX_MESH_MPICOMMUNICATIONS_MPICOMMUNICATIONUTILS_HPP_ */
