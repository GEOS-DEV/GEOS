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
  void addFields(FieldLocation const location, std::vector<string> const fieldNames )
  {
    string const key = generateKey(location);
    addFields(fieldNames, key);
  }

  void addFields(FieldLocation const location, std::vector<string> const fieldNames, std::vector<string> const regionNames )
  {
    for (auto const & regionName : regionNames)
    {
      string const key = generateKey(location, regionName);
      addFields(fieldNames, key);
    }
  }

  std::map< string, std::vector<string> > const & getFields() const
  {
    return m_fields;
  }

  array1d< string > getFieldNames() const
  {
    array1d< string > tmp;

    for( auto const & field : fields.second )
    {
      tmp.emplace_back( field );
    }

    return tmp;
  }
  
  static string const getRegionName( string const & key )
  {
    string regionName(key);
    regionName.erase(0, 5);
    return regionName;
  }

  static FieldLocation const getLocation( string const & key )
  {
    FiedlLocation location;
    if (key.contains("nodes"))
    {
      location = FieldLocation::Node;
    }else if ( key.contains("edges") )
    {
      location = FieldLocation::Edge;
    }
    else if ( key.contains("faces") )
    {
      location = FieldLocation::Face;
    }
    else if ( key.contains("elems") )
    {
      location = FieldLocation::Elem;
    }
    return location;
  }

  private:
  
  std::map< string, std::vector<string> > m_fields;

  static string const generateKey( FieldLocation const location ) const
  
  {
    string key;
    switch(location)
    {
      case FieldLocation::Node:
      {
        key = "nodes"; 
        break;
      }
      case FieldLocation::Edge:
      {
       key = "edges"
        break;
      }
      case FieldLocation::Face:
      {
        key = "faces"
        break;
      }
      case FieldLocation::Elem:
      {
        GEOSX_ERROR("An element located field also requires a region name to be specified.")
        break;
      } 
      return key;
    }
  }

  static string const generateKey( FieldLocation const location, string const regionName ) const
  {
    if ( location == FieldLocation::Elem ) 
    {
      return strcat("elems/", regionName);
    }else
    {
      return generateKey(location);
    }
  }

  void addFields(std::vector<string> const fieldNames, string const key)
  {
    for ( auto const & field : filedNames )
    {
      fields[key].emplace_back(field);
    }
  }
};

} /* namespace geosx */

#endif /* GEOSX_MESH_MPICOMMUNICATIONS_MPICOMMUNICATIONUTILS_HPP_ */
