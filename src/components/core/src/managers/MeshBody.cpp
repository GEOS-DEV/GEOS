/*
 * MeshBody.cpp
 *
 *  Created on: Sep 13, 2017
 *      Author: settgast
 */

#include "MeshBody.hpp"
#include "MeshLevel.hpp"

namespace geosx
{
using namespace dataRepository;

MeshBody::MeshBody( string const & name,
                      ManagedGroup * const parent ):
  ManagedGroup(name,parent)
{
  RegisterViewWrapper<int32>( viewKeys.meshLevels );
}

MeshBody::~MeshBody()
{
  // TODO Auto-generated destructor stub
}



MeshLevel * MeshBody::CreateMeshLevel( int32 const newLevel )
{
  this->RegisterGroup<MeshLevel>( "Level0" );
  this->RegisterGroup<MeshLevel>( "Level1" );
}

} /* namespace geosx */
