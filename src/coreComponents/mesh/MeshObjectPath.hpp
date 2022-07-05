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

  permutationMapType const & pathPermutations() const
  {
    return m_pathPermutations;
  }


  template< typename LAMBDA >
  void forObjectsInPath( LAMBDA lambda )
  {
    lambda();
  }


private:

  string m_objectType;
  permutationMapType m_pathPermutations;

};

} /* namespace geosx */

#endif /* GEOSX_MESH_MESHOBJECTPATH_HPP_ */
