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
 * @file FieldIdentifiers.hpp
 */

#ifndef GEOSX_MESH_MPICOMMUNICATIONS_FIELDIDENTIFIERS_HPP_
#define GEOSX_MESH_MPICOMMUNICATIONS_FIELDIDENTIFIERS_HPP_

#include "common/DataTypes.hpp"
#include "codingUtilities/StringUtilities.hpp"

namespace geosx
{

enum class FieldLocation
{
  Elem,   //!< location is element (like pressure in finite volumes)
  Face,   //!< location is face (like flux in mixed finite elements)
  Edge,   //!< location is edge (like flux between fracture elements)
  Node    //!< location is node (like displacements in finite elements)
};

class FieldIdentifiers
{
  public: 
  
  void addFields( FieldLocation const location, std::vector<string> const & fieldNames )
  {
    string key;
    generateKey(location, key);
    addFields(fieldNames, key);
  }

  void addElementFields( std::vector<string> const & fieldNames, std::vector<string> const & regionNames )
  {
    for (string const & regionName : regionNames)
    {
      addFields(fieldNames, generateKey( regionName ) );
    }
  }

  void addElementFields( std::vector<string> const & fieldNames, arrayView1d< string const > const & regionNames )
  {
    std::vector<string> regions(regionNames.begin(), regionNames.end());
    addElementFields( fieldNames, regions );
  }

  std::map< string, array1d<string> > const & getFields() const
  {
    return m_fields;
  }
  
  string getRegionName( string const & key ) const
  {
    string regionName(key);
    regionName.erase(0, std::string(m_locationKeys.elemsKey()).length());
    return regionName;
  }

   void getLocation( string const & key,
                     FieldLocation & location ) const
  {
    if ( key.find(m_locationKeys.nodesKey() ) != string::npos )
    {
      location = FieldLocation::Node;
    }else if ( key.find(m_locationKeys.edgesKey()) != string::npos )
    {
      location = FieldLocation::Edge;
    }
    else if ( key.find(m_locationKeys.facesKey()) != string::npos )
    {
      location = FieldLocation::Face;
    }
    else if ( key.find(m_locationKeys.elemsKey()) != string::npos )
    {
      location = FieldLocation::Elem;
    }else
    {
      GEOSX_ERROR( GEOSX_FMT( "Invalid key, {}, was provided. Location cannot be retrieved.", key) );
    }
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
  } m_locationKeys;

  void generateKey( FieldLocation const location,
                    string & key ) const
  
  {
    switch(location)
    {
      case FieldLocation::Node:
      {
        key = m_locationKeys.nodesKey(); 
        break;
      }
      case FieldLocation::Edge:
      {
       key = m_locationKeys.edgesKey();
        break;
      }
      case FieldLocation::Face:
      {
        key = m_locationKeys.facesKey();
        break;
      }
      case FieldLocation::Elem:
      {
        GEOSX_ERROR("An element located field also requires a region name to be specified.");
        break;
      } 
    }
  }

  string generateKey( string const & regionName ) const
  {
    return stringutilities::concat("", m_locationKeys.elemsKey(), regionName);  
  }

  void addFields(std::vector<string> const fieldNames, string const key)
  {
    for ( string const & field : fieldNames )
    {
      m_fields[key].emplace_back(field);
    }
  }
};

} /* namespace geosx */

#endif /* GEOSX_MESH_MPICOMMUNICATIONS_FIELDIDENTIFIERS_HPP_ */
