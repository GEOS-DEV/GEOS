/*
 * MeshLevel.cpp
 *
 *  Created on: Sep 13, 2017
 *      Author: settgast
 */

#include "MeshLevel.hpp"
#include "NodeManager.hpp"
//#include "EdgeManager.hpp"
#include "FaceManager.hpp"
#include "ElementRegionManager.hpp"

namespace geosx
{
using namespace dataRepository;

MeshLevel::MeshLevel( string const & name,
                      ManagedGroup * const parent ):
  ManagedGroup(name,parent)//,
//  m_nodeManager( groupStructKeys::nodeManagerString,this)
{
//  RegisterGroup<NodeManager>( groupKeys.vertexManager );
//  RegisterGroup<EdgeManager>( groupKeys.cellManager );


  RegisterGroup<NodeManager>( groupKeys.nodeManager );
//  RegisterGroup<EdgeManager>( groupKeys.edgeManager );
  RegisterGroup<FaceManager>( groupKeys.faceManager );
  RegisterGroup<ElementRegionManager>( groupKeys.elemManager );

//  RegisterGroup( groupStructKeys::nodeManagerString, &m_nodeManager, false );
//  RegisterGroup<EdgeManager>( groupKeys.edgeManager );
//  RegisterGroup<FaceManager>( groupKeys.faceManager );
//  RegisterGroup<ElementRegionManager>( groupKeys.elemManager );


  RegisterViewWrapper<integer>( viewKeys.meshLevel );
}

MeshLevel::~MeshLevel()
{
  // TODO Auto-generated destructor stub
}

} /* namespace geosx */
