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
 * @file MeshObjectPath.hpp
 */

#ifndef GEOSX_MESH_MESHOBJECTPATH_HPP_
#define GEOSX_MESH_MESHOBJECTPATH_HPP_


#include "codingUtilities/EnumStrings.hpp"
#include "MeshLevel.hpp"

namespace geosx
{
class MeshBody;

class MeshObjectPath
{
public:

  using permutationMapType = std::map< string, std::map< string, std::map< string, std::vector< string > > > >;

  enum class ObjectTypes : int
  {
    nodes,
    edges,
    faces,
    elems,
    invalid
  };


  MeshObjectPath( string const path,
                  dataRepository::Group const & meshBodies );

  virtual ~MeshObjectPath();



  void processPath( string const path,
                    dataRepository::Group const & meshBodies );

  ObjectTypes const & getObjectType() const
  {
    return m_objectType;
  }

  permutationMapType const & pathPermutations() const
  {
    return m_pathPermutations;
  }


  template< typename OBJECT_TYPE = ObjectManagerBase,
            typename FUNC = std::function< void( dataRepository::Group &) > >
  void forObjectsInPath( dataRepository::Group & meshBodies,
                         FUNC && func );

  template< typename OBJECT_TYPE = ObjectManagerBase,
            typename FUNC = std::function< void( dataRepository::Group const &) > >
  void forObjectsInPath( dataRepository::Group const & meshBodies,
                         FUNC && func );

#if defined(MESH_OBJECT_PATH_PRIVATE_FUNCTION_UNIT_TESTING)
  template< typename OBJECT_TYPE >
  void testCheckObjectTypeConsistency()
  {
    checkObjectTypeConsistency< OBJECT_TYPE >();
  }

  std::vector< string > testFillPathTokens( string const & path,
                                            dataRepository::Group const & meshBodies )
  {
    return fillPathTokens( path, meshBodies );
  }
#endif

private:
  template< typename OBJECT_TYPE >
  void checkObjectTypeConsistency();

  void printPermutations() const;

  std::vector< string > fillPathTokens( string const & path,
                                        dataRepository::Group const & meshBodies );

  void processPathTokens( std::vector< string > const & pathTokens,
                          dataRepository::Group const & meshBodies );

  ObjectTypes const m_objectType;
  permutationMapType m_pathPermutations;

};

} /* namespace geosx */

#include "MeshBody.hpp"

namespace geosx
{


template< typename OBJECT_TYPE >
void MeshObjectPath::checkObjectTypeConsistency()
{
  bool consistent = false;
  if( m_objectType == ObjectTypes::nodes )
  {
    consistent = std::is_same< NodeManager, OBJECT_TYPE >::value;
  }
  else if( m_objectType == ObjectTypes::edges )
  {
    consistent = std::is_same< EdgeManager, OBJECT_TYPE >::value;
  }
  else if( m_objectType == ObjectTypes::faces )
  {
    consistent = std::is_same< FaceManager, OBJECT_TYPE >::value;
  }
  else if( m_objectType == ObjectTypes::elems )
  {
    consistent = std::is_same< ElementRegionManager, OBJECT_TYPE >::value ||
                 std::is_same< ElementRegionBase, OBJECT_TYPE >::value ||
                 std::is_base_of< ElementRegionBase, OBJECT_TYPE >::value ||
                 std::is_same< ElementSubRegionBase, OBJECT_TYPE >::value ||
                 std::is_base_of< ElementSubRegionBase, OBJECT_TYPE >::value;
  }

  GEOSX_ERROR_IF( !consistent,
                  GEOSX_FMT( "Inconsistent type specified. Type {} is not consistent with m_objectType of {}",
                             OBJECT_TYPE::catalogName(),
                             m_objectType ) );
}

template< typename OBJECT_TYPE,
          typename FUNC >
void MeshObjectPath::forObjectsInPath( dataRepository::Group & meshBodies,
                                       FUNC && func )
{
  forObjectsInPath< OBJECT_TYPE >( const_cast< dataRepository::Group const & >(meshBodies),
                                   [&]( auto const & object )
  {
    func( const_cast< OBJECT_TYPE & >(object) );
  } );
}

template< typename OBJECT_TYPE, typename FUNC >
void MeshObjectPath::forObjectsInPath( dataRepository::Group const & meshBodies,
                                       FUNC && func )
{
  checkObjectTypeConsistency< OBJECT_TYPE >();
  for( auto const & meshBodyPair : m_pathPermutations )
  {
    MeshBody const & meshBody = meshBodies.getGroup< MeshBody >( meshBodyPair.first );
    for( auto const & meshLevelPair : meshBodyPair.second )
    {
      MeshLevel const & meshLevel = meshBody.getMeshLevel( meshLevelPair.first );

      if( m_objectType == ObjectTypes::nodes )
      {
        func( dynamic_cast< OBJECT_TYPE const & >(meshLevel.getNodeManager() ) );
      }
      else if( m_objectType == ObjectTypes::edges )
      {
        func( dynamic_cast< OBJECT_TYPE const & >(meshLevel.getEdgeManager()) );
      }
      else if( m_objectType == ObjectTypes::faces )
      {
        func( dynamic_cast< OBJECT_TYPE const & >(meshLevel.getFaceManager()) );
      }
      else if( m_objectType == ObjectTypes::elems )
      {
        ElementRegionManager const & elemRegionMan = meshLevel.getElemManager();
        for( auto const & elemRegionPair : meshLevelPair.second )
        {
          ElementRegionBase const & elemRegion = elemRegionMan.getRegion( elemRegionPair.first );
          if( std::is_same< ElementRegionBase, OBJECT_TYPE >::value ||
              std::is_base_of< ElementRegionBase, OBJECT_TYPE >::value )
          {
            func( dynamic_cast< OBJECT_TYPE const & >(elemRegion) );
          }
          else
          {
            for( auto const & elemSubRegionName : elemRegionPair.second )
            {
              ElementSubRegionBase const & subRegion = elemRegion.getSubRegion( elemSubRegionName );
              if( std::is_same< ElementSubRegionBase, OBJECT_TYPE >::value ||
                  std::is_base_of< ElementSubRegionBase, OBJECT_TYPE >::value )
              {
                func( dynamic_cast< OBJECT_TYPE const & >(subRegion) );
              }
              else
              {
                GEOSX_ERROR( "You shouldn't be here" );
              }
            }
          }
        }
      }
    }
  }
}

ENUM_STRINGS( MeshObjectPath::ObjectTypes,
              MeshLevel::groupStructKeys::nodeManagerString(),
              MeshLevel::groupStructKeys::edgeManagerString(),
              MeshLevel::groupStructKeys::faceManagerString(),
              MeshLevel::groupStructKeys::elemManagerString(),
              "invalid" );

} /* namespace geosx */

#endif /* GEOSX_MESH_MESHOBJECTPATH_HPP_ */
