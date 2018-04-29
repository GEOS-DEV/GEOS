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
  ManagedGroup(name,parent),
  m_nodeManager( groupStructKeys::nodeManagerString,this),
  m_faceManager( groupStructKeys::faceManagerString,this),
  m_elementManager( groupStructKeys::elemManagerString,this)
{

  RegisterGroup( groupStructKeys::nodeManagerString, &m_nodeManager, false );
//  RegisterGroup<EdgeManager>( groupKeys.edgeManager );
  RegisterGroup<FaceManager>( groupStructKeys::faceManagerString, &m_faceManager, false );
  RegisterGroup<ElementRegionManager>( groupStructKeys::elemManagerString, &m_elementManager, false );


  RegisterViewWrapper<integer>( viewKeys.meshLevel );
}

MeshLevel::~MeshLevel()
{
  // TODO Auto-generated destructor stub
}

} /* namespace geosx */
