/*
 * MeshBody.hpp
 *
 *  Created on: Sep 13, 2017
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MANAGERS_MESHBODY_HPP_
#define SRC_COMPONENTS_CORE_SRC_MANAGERS_MESHBODY_HPP_

#include "MeshLevel.hpp"


namespace geosx
{



class MeshLevel;

class MeshBody : public dataRepository::ManagedGroup
{
public:
  MeshBody( string const & name,
             ManagedGroup * const parent );
  virtual ~MeshBody();

  MeshLevel * CreateMeshLevel( int32 const newLevel );

  MeshLevel * getMeshLevel( int32 const level ) { return this->GetGroup<MeshLevel>(level); }
  MeshLevel const * getMeshLevel( int32 const level ) const { return this->GetGroup<MeshLevel>(level); }

  struct viewKeysStruct
  {


    dataRepository::ViewKey meshLevels                = { "meshLevels" };
  }viewKeys;

  struct groupStructKeys
  {
  }groupKeys;


};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_MESHLEVEL_HPP_ */
