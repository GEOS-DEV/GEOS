/*
 * MeshPath.hpp
 *
 *  Created on: Jun 21, 2022
 *      Author: settgast
 */
#ifndef GEOSX_MESH_MESHOBJECTPATH_HPP_
#define GEOSX_MESH_MESHOBJECTPATH_HPP_


#include "common/DataTypes.hpp"
#include "dataRepository/Group.hpp"

namespace geosx
{
class MeshBody;
class MeshLevel;
class ElementRegionBase;
class ObjectManagerBase;

class MeshObjectPath
{
public:

  using permutationMapType = std::map< string, std::map< string, std::map< string, std::vector< string > > > >;

  MeshObjectPath( string const path,
                  dataRepository::Group const & meshBodies );

  virtual ~MeshObjectPath();


  static string getObjectType( string const & path );


  void fillPathTokens( string const & path,
                       std::vector<string> & pathTokens,
                       dataRepository::Group const & meshBodies );

  void ExpandPathTokens( std::vector<string> & pathTokens,
                         dataRepository::Group const & meshBodies );


  void processPath( string const path,
                    dataRepository::Group const & meshBodies );

  void printPermutations() const;

//  std::vector<string> const & getPaths() const { return m_splitPath; }

  string const & getObjectType() const
  {
    return m_objectType;
  }

  bool isElementPath() const;

  permutationMapType const & pathPermutations() const
  {
    return m_pathPermutations;
  }

  template< typename OBJECT_TYPE >
  void checkObjectTypeConsistency();

  template< typename OBJECT_TYPE = ObjectManagerBase, 
            typename FUNC = std::function< void( dataRepository::Group &)> >
  void forObjectsInPath( dataRepository::Group & meshBodies,
                         FUNC&& func )
  {
    forObjectsInPath<OBJECT_TYPE>( const_cast<dataRepository::Group const &>(meshBodies),
                                   [&]( auto const & object )
    {
      func( const_cast<OBJECT_TYPE &>(object) );
    } );
  }

  template< typename OBJECT_TYPE = ObjectManagerBase, 
            typename FUNC = std::function< void( dataRepository::Group const &)> >
  void forObjectsInPath( dataRepository::Group const & meshBodies,
                         FUNC&& func );

private:

  string m_objectType;
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
    if( m_objectType == MeshLevel::groupStructKeys::nodeManagerString() )
    {
      consistent = std::is_same< NodeManager, OBJECT_TYPE >::value;
    }
    else if( m_objectType == MeshLevel::groupStructKeys::edgeManagerString() )
    {
      consistent = std::is_same< EdgeManager, OBJECT_TYPE >::value;
    }
    else if( m_objectType == MeshLevel::groupStructKeys::faceManagerString() )
    {
      consistent = std::is_same< FaceManager, OBJECT_TYPE >::value;
    }
    else if( m_objectType == MeshLevel::groupStructKeys::elemManagerString() )
    {
      consistent = std::is_same< ElementRegionManager, OBJECT_TYPE >::value || 
                   std::is_same< ElementRegionBase, OBJECT_TYPE >::value ||
                   std::is_base_of< ElementRegionBase, OBJECT_TYPE >::value ||
                   std::is_same< ElementSubRegionBase, OBJECT_TYPE >::value ||
                   std::is_base_of< ElementSubRegionBase, OBJECT_TYPE >::value ;
    }

    GEOSX_ERROR_IF( !consistent,
                    GEOSX_FMT( "Inconsistent type specified. Type {} is not consistent with m_objectType of {}",
                               OBJECT_TYPE::catalogName(),
                               m_objectType) );
  }

  template< typename OBJECT_TYPE, typename FUNC > 
  void MeshObjectPath::forObjectsInPath( dataRepository::Group const & meshBodies,
                                         FUNC&& func )
  {
    checkObjectTypeConsistency<OBJECT_TYPE>();
    for( auto const & meshBodyPair : m_pathPermutations )
    {
      MeshBody const & meshBody = meshBodies.getGroup<MeshBody>( meshBodyPair.first );
      for( auto const & meshLevelPair : meshBodyPair.second )
      {
        MeshLevel const & meshLevel = meshBody.getMeshLevel( meshLevelPair.first );

        if( m_objectType == MeshLevel::groupStructKeys::nodeManagerString() )
        {
          func( dynamic_cast<OBJECT_TYPE const &>(meshLevel.getNodeManager() ) );
        }
        else if( m_objectType == MeshLevel::groupStructKeys::edgeManagerString() )
        {
          func( dynamic_cast<OBJECT_TYPE const &>(meshLevel.getEdgeManager()) );
        }
        else if( m_objectType == MeshLevel::groupStructKeys::faceManagerString() )
        {
          func( dynamic_cast<OBJECT_TYPE const &>(meshLevel.getFaceManager()) );
        }
        else if( m_objectType == MeshLevel::groupStructKeys::elemManagerString() )
        {
          ElementRegionManager const & elemRegionMan = meshLevel.getElemManager();
          for( auto const & elemRegionPair : meshLevelPair.second )
          {
            ElementRegionBase const & elemRegion = elemRegionMan.getRegion( elemRegionPair.first );
            if( std::is_same< ElementRegionBase, OBJECT_TYPE >::value ||
                std::is_base_of< ElementRegionBase, OBJECT_TYPE >::value )
            {
              func( dynamic_cast<OBJECT_TYPE const &>(elemRegion) );
            }
            else
            {
              for( auto const & elemSubRegionName : elemRegionPair.second )
              {
                ElementSubRegionBase const & subRegion = elemRegion.getSubRegion( elemSubRegionName );
                if( std::is_same< ElementSubRegionBase, OBJECT_TYPE >::value ||
                    std::is_base_of< ElementSubRegionBase, OBJECT_TYPE >::value )
                {
                  func( dynamic_cast<OBJECT_TYPE const &>(subRegion) );
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
} /* namespace geosx */

#endif /* GEOSX_MESH_MESHOBJECTPATH_HPP_ */
