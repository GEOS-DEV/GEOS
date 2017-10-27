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
  RegisterViewWrapper<integer>( viewKeys.meshLevels );
}

MeshBody::~MeshBody()
{
  // TODO Auto-generated destructor stub
}



MeshLevel * MeshBody::CreateMeshLevel( integer const newLevel )
{
  return this->RegisterGroup<MeshLevel>( "Level0" );
}

} /* namespace geosx */
