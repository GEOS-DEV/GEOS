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

class FieldIdentifiers
{
  public: 
  
  void addFields( FieldLocation const location, std::vector<string> const fieldNames )
  {
    string const key = generateKey(location);
    addFields(fieldNames, key);
  }

  void addElementFields( std::vector<string> const fieldNames, std::vector<string> const regionNames )
  {
    for (string const & regionName : regionNames)
    {
      string const key = generateKey( regionName );
      addFields(fieldNames, key);
    }
  }

  std::map< string, array1d<string> > const & getFields() const
  {
    return m_fields;
  }
  
  static string const getRegionName( string const & key )
  {
    string regionName(key);
    regionName.erase(0, locationKeys.elemsKey().length());
    return regionName;
  }

  static FieldLocation const getLocation( string const & key )
  {
    FiedlLocation location;
    if (key.contains(locationKeys.nodesKey()))
    {
      location = FieldLocation::Node;
    }else if ( key.contains(locationKeys.edgesKey()) )
    {
      location = FieldLocation::Edge;
    }
    else if ( key.contains(locationKeys.facesKey()) )
    {
      location = FieldLocation::Face;
    }
    else if ( key.contains(locationKeys.elemsKey()) )
    {
      location = FieldLocation::Elem;
    }
    return location;
  }

  private:
  
  std::map< string, array1d<string> > m_fields;

  struct keysStruct 
  {
    /// @return String key for 
    static constexpr char const * nodesKey() { return "nodes"; }
    /// @return String key for 
    static constexpr char const * edgesKey() { return "edges"; }
    /// @return String key f
    static constexpr char const * facesKey() { return "faces"; }
    /// @return String key 
    static constexpr char const * elemsKey() { return "elems/"; }
  } locationKeys;

  static string const generateKey( FieldLocation const location ) const
  
  {
    string key;
    switch(location)
    {
      case FieldLocation::Node:
      {
        key = locationKeys.nodesKey(); 
        break;
      }
      case FieldLocation::Edge:
      {
       key = locationKeys.edgesKey();
        break;
      }
      case FieldLocation::Face:
      {
        key = locationKeys.facesKey();
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

  static string const generateKey( string const regionName ) const
  {
      return strcat(locationKeys.elemsKey(), regionName);
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
