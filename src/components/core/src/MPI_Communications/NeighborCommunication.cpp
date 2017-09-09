//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @file NeighborCommunication.cpp
 * @author settgast1
 * @date Mar 9, 2011
 */

#include "NeighborCommunication.hpp"

#include "legacy/IO/BinStream.h"
//#include "../ObjectManagers/PhysicalDomainT.h"
//#include "Utilities/Utilities.h"


#if ATK_FOUND
#include <slic/slic.hpp>
#endif


static void GetModifiedNeighborIndices( const ObjectDataStructureBaseT& object,
                                        const lArray1d& existingNeighborIndices,
                                        const lSet& modifiedIndices,
                                        lArray1d& modifiedNeighborIndices );

static void GetModifiedNeighborIndices( const ObjectDataStructureBaseT& object,
                                        const lArray1d& existingNeighborIndices,
                                        const lArray1d& modifiedIndices,
                                        lArray1d& modifiedNeighborIndices );





static void GetModifiedNeighborIndices( const ObjectDataStructureBaseT& object,
                                        const lArray1d& existingNeighborIndices,
                                        const lSet& modifiedIndices,
                                        lArray1d& modifiedNeighborIndices )
{

//  if( !(modifiedIndices.empty()) )
//  {
//    lSet existingNeighborIndicesSet( existingNeighborIndices.begin(), existingNeighborIndices.end() );
//
//    for( lSet::const_iterator a=modifiedIndices.begin() ; a!=modifiedIndices.end() ; ++a )
//    {
//
//      const localIndex parentIndex = object.GetParentIndex(*a);
//      if( existingNeighborIndicesSet.count( parentIndex ) == 1 )
//      {
//        modifiedNeighborIndices.push_back(*a);
//      }
//
//    }
//  }

}


static void GetModifiedNeighborIndices( const ObjectDataStructureBaseT& object,
                                        const lArray1d& existingNeighborIndices,
                                        const lArray1d& modifiedIndices,
                                        lArray1d& modifiedNeighborIndices )
{
  lSet modifiedIndicesSet( modifiedIndices.begin(), modifiedIndices.end() );
  GetModifiedNeighborIndices( object, existingNeighborIndices, modifiedIndicesSet, modifiedNeighborIndices );

}




namespace geosx
{

NeighborCommunication::NeighborCommunication():
tempNeighborData(),
m_neighborRank(-1),
m_receiveLocalIndices(),
//m_domain(NULL),
m_rank(-1),
m_size(-1),
m_rankOfNeighborNeighbors(),
m_sendLocalIndices(),
m_elementRegionsSendLocalIndices(),
m_elementRegionsReceiveLocalIndices(),
m_sendBuffer(),
m_receiveBuffer()
{
  for( int i=0 ; i<CommRegistry::maxComm ; ++i )
  {
    mpiSendSizeRequest[i] = MPI_REQUEST_NULL;
    mpiRecvSizeRequest[i] = MPI_REQUEST_NULL;
    mpiSendBufferRequest[i] = MPI_REQUEST_NULL;
    mpiRecvBufferRequest[i] = MPI_REQUEST_NULL;
    m_sendSize[i] = 0;
    m_receiveSize[i] = 0;
  }
}

NeighborCommunication::NeighborCommunication( const NeighborCommunication& init ):
tempNeighborData(init.tempNeighborData),
m_neighborRank(init.m_neighborRank),
m_receiveLocalIndices(init.m_receiveLocalIndices),
//m_domain(init.m_domain),
m_rank(init.m_rank),
m_size(init.m_size),
m_rankOfNeighborNeighbors(init.m_rankOfNeighborNeighbors),
m_sendLocalIndices(init.m_sendLocalIndices),
m_elementRegionsSendLocalIndices(init.m_elementRegionsSendLocalIndices),
m_elementRegionsReceiveLocalIndices(init.m_elementRegionsReceiveLocalIndices),
m_sendBuffer(init.m_sendBuffer),
m_receiveBuffer(init.m_receiveBuffer)
{
 for( int i=0 ; i<CommRegistry::maxComm ; ++i )
  {
    mpiSendSizeRequest[i] = init.mpiSendSizeRequest[i];
    mpiRecvSizeRequest[i] = init.mpiRecvSizeRequest[i];
    mpiSendBufferRequest[i] = init.mpiSendBufferRequest[i];
    mpiRecvBufferRequest[i] = init.mpiRecvBufferRequest[i];
    m_sendSize[i] = init.m_sendSize[i];
    m_receiveSize[i] = init.m_receiveSize[i];
  }
}


NeighborCommunication& NeighborCommunication::operator=( const NeighborCommunication& init )
{
//  SetDomain( *(init.m_domain) );

  tempNeighborData = init.tempNeighborData;
  m_neighborRank = init.m_neighborRank;
  m_receiveLocalIndices = init.m_receiveLocalIndices;
  m_rank = init.m_rank;
  m_size = init.m_size;
  m_rankOfNeighborNeighbors = init.m_rankOfNeighborNeighbors;
  m_sendLocalIndices = init.m_sendLocalIndices;
  m_elementRegionsSendLocalIndices = init.m_elementRegionsSendLocalIndices;
  m_elementRegionsReceiveLocalIndices = init.m_elementRegionsReceiveLocalIndices;
  m_sendBuffer = init.m_sendBuffer;
  m_receiveBuffer = init.m_receiveBuffer;


 for( int i=0 ; i<CommRegistry::maxComm ; ++i )
  {
    mpiSendSizeRequest[i] = init.mpiSendSizeRequest[i];
    mpiRecvSizeRequest[i] = init.mpiRecvSizeRequest[i];
    mpiSendBufferRequest[i] = init.mpiSendBufferRequest[i];
    mpiRecvBufferRequest[i] = init.mpiRecvBufferRequest[i];
    m_sendSize[i] = init.m_sendSize[i];
    m_receiveSize[i] = init.m_receiveSize[i];
  }
 return *this;
}


NeighborCommunication::~NeighborCommunication()
{}

void NeighborCommunication::Initialize( const int neighborRank,
                                        const int rank,
                                        const int size )
{
  m_neighborRank = neighborRank;
  m_rank = rank;
  m_size = size;
}

void NeighborCommunication::Clear()
{
  tempNeighborData.clear();

  std::map<string, lArray1d>::iterator it;
  for(it = m_sendLocalIndices.begin(); it != m_sendLocalIndices.end(); ++it)
  {
    it->second.clear();
  }
  for( std::map< std::string, lArray1d>::iterator i=m_elementRegionsSendLocalIndices.begin() ;
       i!=m_elementRegionsSendLocalIndices.end() ; ++i )
  {
    i->second.clear();
  }

  for(it = m_receiveLocalIndices.begin(); it != m_receiveLocalIndices.end(); ++it)
  {
    it->second.clear();
  }
  for( std::map< std::string, lArray1d>::iterator i=m_elementRegionsReceiveLocalIndices.begin() ;
      i!=m_elementRegionsReceiveLocalIndices.end() ; ++i )
  {
    i->second.clear();
  }
}


//void NeighborCommunication::SetDomain( PhysicalDomainT& domain )
//{
//  // set the const pointer "m_domain" by casting away the the const
//  // on the address of the pointer, and modifying what the address of
//  // the pointer.
////  PhysicalDomainT** temp = const_cast<PhysicalDomainT**>(&m_domain);
////  *temp = &domain;
//
//}

void NeighborCommunication::CommunicatePackedObjectBufferSizes( )
{
  const size_t n = tempNeighborData.objectsToSend.size();
  if(tempNeighborData.objectsToReceive.size() != n) {
#if ATK_FOUND
    SLIC_ERROR("Cannot have number of object types received differ from those sent");
#endif
  }

  bufvector::size_type* sendSizes = new bufvector::size_type[n];
  bufvector::size_type* receiveSizes = new bufvector::size_type[n];
  {
    std::map<string, bufvector>::const_iterator it;
    localIndex i = 0;
    for(it = tempNeighborData.objectsToSend.begin(); it != tempNeighborData.objectsToSend.end(); ++it, ++i)
    {
      sendSizes[i] = it->second.size();
    }
  }

  MPI_Request mpiRequest[2];
  MPI_Status mpiStatus;

  SendReceive( reinterpret_cast<char*>(&receiveSizes), n*sizeof(bufvector::size_type),
               reinterpret_cast<char*>(&sendSizes), n*sizeof(bufvector::size_type),
               mpiRequest[0], mpiRequest[1] );

  MPI_Waitall( 2, mpiRequest, &mpiStatus );

  delete[] sendSizes;
  {
    std::map<string, bufvector>::iterator it;
    localIndex i = 0;
    for(it = tempNeighborData.objectsToReceive.begin(); it != tempNeighborData.objectsToReceive.end(); ++it, ++i)
    {
      it->second.resize( receiveSizes[i]);
    }
  }
  delete[] receiveSizes;
}

/**
 * @brief Fill the matchedObjectsIndices and matchedObjectsNumbers
 *
 * Use data in tempNeighborData.neighborNumbers, localObjectNumbers, and the global-to-local map in object
 * to fill the tempNeighborData.matchedObjectsIndices and tempNeighborData.matchedObjectsNumbers
 *
 * @param[in] object ObjectDataStructureBaseT manager whose globalToLocal map structure corresponds to the query
 * @param[in] name ObjectDataStructureBaseT manager name in domain
 * @param[in] localObjectNumbers Array of local numbers to query
 */
void NeighborCommunication::DetermineMatchedBoundaryObject( const ObjectDataStructureBaseT& object,
                                                            const string name,
                                                            const gArray1d& localObjectNumbers)
{
//  //get temporary reference to the (already created) numbers list for the given object type
//  std::map<string, gArray1d>::const_iterator cgit;
//  cgit = tempNeighborData.neighborNumbers.find(name);
//  if (cgit == tempNeighborData.neighborNumbers.end())
//  {
////#if ATK_FOUND
////    SLIC_ERROR(
////        "Cannot find name " + toString<int>(name)
////            + " in neighborNumbers in NeighborCommunication::DetermineMatchedBoundaryObject");
////#endif
//  }
//  const gArray1d& neighborObjectNumbers = cgit->second;
//
//  //instantiate indices and numbers lists for the given object type
//  //and assign a temporary reference
//  lArray1d& matchedObjectsIndices = tempNeighborData.matchedIndices[name];
//  gArray1d& matchedObjectsNumbers = tempNeighborData.matchedNumbers[name];
//
//  // now compare local and neighbor lists. They are sorted, so that makes it a lot easier.
//  // the looping here is a little funny. We are going to loop over all entries in the localNodes array, but the
//  // neighborNodes array is also being iterated over...just not as part of the loop counter.
//  gArray1d::const_iterator localObjectNumber = localObjectNumbers.begin();
//  gArray1d::const_iterator neighborObjectNumber = neighborObjectNumbers.begin();
//
//  // loop over local objects. exit loop if either the local or neighbor object is at end()
//  for( ; (localObjectNumber!=localObjectNumbers.end() && neighborObjectNumber!=neighborObjectNumbers.end()) ;
//      ++localObjectNumber)
//  {
//
//    // a loop to increment the neighbor objects
//    while( neighborObjectNumber!=neighborObjectNumbers.end() )
//    {
//      // if the object indices match, then
//      if( *localObjectNumber == *neighborObjectNumber )
//      {
//        localIndex localMatchedObjectIndex = stlMapLookup( object.m_globalToLocalMap, *localObjectNumber );
//
//        matchedObjectsIndices.push_back(localMatchedObjectIndex);
//        matchedObjectsNumbers.push_back( *localObjectNumber );
//
//      }
//      else if ( *localObjectNumber < *neighborObjectNumber )
//        break;
//
//      ++neighborObjectNumber;
//    }
//  }

  // TODO probably should check to see that the matchedBoundaryNodesGlobal array is the same on both processes
}

void NeighborCommunication::FindGhosts( const bool contactActive,
                                        const int depth )
{
//  //Currently, this MUST fill the objectsToSend associated with the following
//  //(see NeighborCommunication::SyncNames for the current list)
//  // 0: NodeManager -> PackNodes
//  // 1: EdgeManager -> PackEdges
//  // 2: FaceManager -> PackFaces
//  // 3: FiniteElementElementManager -> PackElements
//  // 4: DiscreteElementManager ->
//  // 5: DiscreteElementNodeManager ->
//  // 6: DiscreteElementFaceManager ->
//  // 7: EllipsoidalDiscreteElementManager ->
//
//  // so now we know which nodes are "shared" between the neighbors. what we do now is to collect all the objects
//  // attached to these nodes. These will be "ghost" objects on the neighbor.
//
//
//  lSet& allNodes = tempNeighborData.objectLocalIndicesToSend[PhysicalDomainT::FiniteElementNodeManager];
//  lSet& allEdges = tempNeighborData.objectLocalIndicesToSend[PhysicalDomainT::FiniteElementEdgeManager];
//  lSet& allFaces = tempNeighborData.objectLocalIndicesToSend[PhysicalDomainT::FiniteElementFaceManager];
//
//  // get list of elements connected to the matched nodes into "m_elementRegionsSendLocalIndices"
//  this->m_elementRegionsSendLocalIndices.clear();
//  m_domain->m_feElementManager.ConstructListOfIndexesFromMap( this->m_domain->m_feNodeManager.m_toElementsRelation,
//                                                              this->tempNeighborData.matchedIndices[PhysicalDomainT::FiniteElementNodeManager],
//                                                              this->m_elementRegionsSendLocalIndices,
//                                                              depth );
//
//  // now pack up the elements that are going over to the neighbor...also get the nodes that we are going to
//  // need to send.
//  m_domain->m_feElementManager.PackElements( this->tempNeighborData.objectsToSend[PhysicalDomainT::FiniteElementElementManager],
//                                           allNodes,
//                                           allFaces,
//                                           this->m_elementRegionsSendLocalIndices,
//                                           this->m_domain->m_feNodeManager,
//                                           this->m_domain->m_feFaceManager,
//                                           true, true, true, true);
//
//  // add the "external" faces if the contact is on
//  if( contactActive )
//  {
//    const iArray1d& isExternalFace = m_domain->m_feFaceManager.m_isExternal;
//
//    for( iArray1d::size_type a=0 ; a<m_domain->m_feFaceManager.m_numFaces ; ++a )
//    {
//      if( isExternalFace[a] == 1 )
//      {
//        bool allValidNodes = true;
//        for( lArray1d::const_iterator i=m_domain->m_feFaceManager.m_toNodesRelation[a].begin() ;
//             i!=m_domain->m_feFaceManager.m_toNodesRelation[a].end() ; ++i )
//        {
//          const globalIndex gnode = m_domain->m_feNodeManager.m_localToGlobalMap[*i];
//          const int owningRank = GlobalIndexManager::OwningRank( gnode );
//          if( !(m_rankOfNeighborNeighbors.count(owningRank)) )
//          {
//            allValidNodes = false;
//          }
//        }
//        if( allValidNodes )
//        {
//          allFaces.insert(a);
//        }
//      }
//    }
//  }
//
//
//  for( lSet::const_iterator faceIndex=allFaces.begin() ; faceIndex!=allFaces.end() ; ++faceIndex )
//  {
//    for( lArray1d::const_iterator edgeIndex=m_domain->m_feFaceManager.m_toEdgesRelation[*faceIndex].begin() ;
//         edgeIndex!=m_domain->m_feFaceManager.m_toEdgesRelation[*faceIndex].end() ; ++edgeIndex )
//    {
//      const globalIndex gi = m_domain->m_feEdgeManager.m_localToGlobalMap[*edgeIndex];
//      if( m_rankOfNeighborNeighbors.count( GlobalIndexManager::OwningRank(gi) ) )
//      {
//        allEdges.insert(*edgeIndex);
//      }
//    }
//
//    for( lArray1d::const_iterator i=m_domain->m_feFaceManager.m_toNodesRelation[*faceIndex].begin() ;
//         i!=m_domain->m_feFaceManager.m_toNodesRelation[*faceIndex].end() ; ++i )
//    {
//      allNodes.insert(*i);
//    }
//  }
//
//
//  const PhysicalDomainT::ObjectDataStructureKeys keys[3] = { PhysicalDomainT::FiniteElementNodeManager,
//                                                             PhysicalDomainT::FiniteElementEdgeManager,
//                                                             PhysicalDomainT::FiniteElementFaceManager };
//
//  for( int i=0 ; i<3 ; ++i )
//  {
//    const lSet& localSends = tempNeighborData.objectLocalIndicesToSend[keys[i]];
//    gArray1d& globalSends = tempNeighborData.objectGlobalIndicesToSend[keys[i]];
//
//    for( lSet::const_iterator a=localSends.begin() ; a!=localSends.end() ; ++a )
//    {
//      const globalIndex gIndex = m_domain->GetObjectDataStructure(keys[i]).m_localToGlobalMap[*a];
//      if( GlobalIndexManager::OwningRank( gIndex ) != m_rank && GlobalIndexManager::OwningRank( gIndex ) != m_neighborRank )
//      {
//        globalSends.push_back( gIndex );
//      }
//    }
//  }
//}
//
//void NeighborCommunication::FindPackGhosts_Step2(   )
//{
//  const lSet& allNodes = tempNeighborData.objectLocalIndicesToSend[PhysicalDomainT::FiniteElementNodeManager];
//  const lSet& allEdges = tempNeighborData.objectLocalIndicesToSend[PhysicalDomainT::FiniteElementEdgeManager];
//  const lSet& allFaces = tempNeighborData.objectLocalIndicesToSend[PhysicalDomainT::FiniteElementFaceManager];
//
//  lArray1d nodalSendList;
//  for( lSet::const_iterator a=allNodes.begin() ; a!=allNodes.end() ; ++a )
//  {
//    const globalIndex gnode = m_domain->m_feNodeManager.m_localToGlobalMap[*a];
//    if( GlobalIndexManager::OwningRank( gnode ) == m_rank )
//    {
//      nodalSendList.push_back( *a );
//    }
//  }
//  m_domain->m_feNodeManager.PackNodes( nodalSendList, m_domain->m_feFaceManager, tempNeighborData.objectsToSend[PhysicalDomainT::FiniteElementNodeManager], true, true, true, true );
//  m_sendLocalIndices[PhysicalDomainT::FiniteElementNodeManager] = nodalSendList;
//
//
//  //***** FACES *****
//  lArray1d faceSendList;
//  for( lSet::const_iterator a=allFaces.begin() ; a!=allFaces.end() ; ++a )
//  {
//    const globalIndex gface = m_domain->m_feFaceManager.m_localToGlobalMap[*a];
//    if( GlobalIndexManager::OwningRank( gface ) == m_rank )
//    {
//      faceSendList.push_back( *a );
//    }
//  }
//
//  m_domain->m_feFaceManager.PackFaces( faceSendList,
//                                     m_domain->m_feNodeManager,
//                                     &m_domain->m_feEdgeManager,
//                                     tempNeighborData.objectsToSend[PhysicalDomainT::FiniteElementFaceManager], true, true, true, true );
//  m_sendLocalIndices[PhysicalDomainT::FiniteElementFaceManager] = faceSendList;
//  m_sendLocalIndices[PhysicalDomainT::VirtualFaceManager] = faceSendList;
//
//
//
//  lArray1d edgeSendList;
//  for( lSet::const_iterator a=allEdges.begin() ; a!=allEdges.end() ; ++a )
//  {
//    const globalIndex gEdge = m_domain->m_feEdgeManager.m_localToGlobalMap[*a];
//    if( GlobalIndexManager::OwningRank( gEdge ) == m_rank )
//    {
//      edgeSendList.push_back( *a );
//    }
//  }
//  m_domain->m_feEdgeManager.PackEdges( edgeSendList,
//                                     m_domain->m_feNodeManager, m_domain->m_feFaceManager,
//                                     tempNeighborData.objectsToSend[PhysicalDomainT::FiniteElementEdgeManager], true, true, true, true );
//  m_sendLocalIndices[PhysicalDomainT::FiniteElementEdgeManager] = edgeSendList;
//
//
//

}



void NeighborCommunication::SetAllOwned(const gArray1d& localToGlobal, lArray1d& indices) const
{
//  indices.clear();
//  localIndex i = 0;
//  for( gArray1d::const_iterator a=localToGlobal.begin() ; a!=localToGlobal.end() ; ++a, ++i )
//    if( GlobalIndexManager::OwningRank( *a ) == m_rank)
//      indices.push_back( i );
}










template< typename TYPE >
void NeighborCommunication::PackTopologyModifications( const string key,
                                                       const TYPE& newObjects,
                                                       const TYPE& modifiedObjects,
                                                       const bool packConnectivityToGlobal,
                                                       const bool reverseOp )
{

//  const typename TYPE::size_type numNew = newObjects.size();
//  m_sendBuffer.Pack( numNew );
//  if( numNew > 0 )
//  {
//    lArray1d newObjectsOnNeighbor;
//
//    GetModifiedNeighborIndices( m_domain->GetObjectDataStructure(key),
//                                reverseOp ? ReceiveLocalIndices(key) : SendLocalIndices(key),
//                                newObjects,
//                                newObjectsOnNeighbor );
//
//    // Pack New Node Data
//    m_domain->Pack( key, newObjectsOnNeighbor, m_sendBuffer, packConnectivityToGlobal, true, true, true );
//
//    if( !reverseOp )
//    {
//      m_sendLocalIndices[key].insert( m_sendLocalIndices[key].end(), newObjectsOnNeighbor.begin(), newObjectsOnNeighbor.end() );
//    }
//    else
//    {
//      m_receiveLocalIndices[key].insert( m_receiveLocalIndices[key].end(), newObjectsOnNeighbor.begin(), newObjectsOnNeighbor.end() );
//    }
//  }
//
//  const typename TYPE::size_type numModified = modifiedObjects.size();
//  m_sendBuffer.Pack( numModified );
//  if( numModified > 0 )
//  {
//    lArray1d modifiedSendIndices;
//    GetModifiedNeighborIndices(  m_domain->GetObjectDataStructure(key),
//                                 reverseOp ? ReceiveLocalIndices(key) : SendLocalIndices(key),
//                                 modifiedObjects,
//                                 modifiedSendIndices );
//    m_domain->Pack( key, modifiedSendIndices, m_sendBuffer, packConnectivityToGlobal, true, true, true );
//  }
}
template void NeighborCommunication::PackTopologyModifications( const string, const lSet&, const lSet&, const bool, const bool );
template void NeighborCommunication::PackTopologyModifications( const string, const lArray1d&, const lArray1d&, const bool, const bool );
template<>
void NeighborCommunication::PackTopologyModifications( const string key,
                                                       const std::map< std::string, lSet >&,
                                                       const std::map< std::string, lSet >& modifiedObjects,
                                                       const bool packConnectivityToGlobal,
                                                       const bool reverseOp )
{
//
//  const std::map< std::string, lSet >::size_type numModified = modifiedObjects.size();
//  m_sendBuffer.Pack( numModified );
//  if( numModified > 0 )
//  {
//    if( key == PhysicalDomainT::FiniteElementElementManager )
//    {
//      std::map< std::string, lArray1d > modifiedSendIndices;
//      for( std::map< std::string, lSet >::const_iterator iter_mod=modifiedObjects.begin() ; iter_mod!=modifiedObjects.end() ; ++iter_mod )
//      {
//        const std::string& name = iter_mod->first;
//        const lSet& modifiedIndices = iter_mod->second;
//
//        GetModifiedNeighborIndices( m_domain->GetObjectDataStructure(PhysicalDomainT::FiniteElementElementRegion,name),
//                                    reverseOp ? ElementRegionsReceiveLocalIndices(name) : ElementRegionsSendLocalIndices(name),
//                                    modifiedIndices,
//                                    modifiedSendIndices[name] );
//      }
//
//      m_domain->Pack( PhysicalDomainT::FiniteElementElementManager, modifiedSendIndices, m_sendBuffer, packConnectivityToGlobal, false, true, true );
//    }
//    else
//    {
//#if ATK_FOUND
//      SLIC_ERROR("NeighborCommunication::PackTopologyModifications: inappropriate type for PhysicalDomainT::Unpack " + toString<int>(key));
//#endif
//    }
//  }
}












void NeighborCommunication::UnpackTopologyModifications( const string key ,
                                                         const char*& pbuffer,
                                                         lArray1d& newIndices,
                                                         lArray1d& modifiedIndices,
                                                         const bool reverseOp )
{
//  if( key!=PhysicalDomainT::FiniteElementElementManager )
//  {
//
//    lArray1d::size_type numNew;
//    bufvector::Unpack( pbuffer, numNew );
//
//    if( numNew > 0 )
//    {
//      m_domain->Unpack( key, pbuffer, newIndices, false, true, true, true );
//
//      if( !reverseOp )
//      {
//        m_receiveLocalIndices[key].insert( m_receiveLocalIndices[key].end(), newIndices.begin(), newIndices.end() );
//      }
//      else
//      {
//        m_sendLocalIndices[key].insert( m_sendLocalIndices[key].end(), newIndices.begin(), newIndices.end() );
//      }
//    }
//
//    lArray1d::size_type numModified;
//    bufvector::Unpack( pbuffer, numModified );
//
//    if( numModified > 0 )
//    {
//      m_domain->Unpack( key, pbuffer, modifiedIndices , false, true, true, true );
//    }
//  }
//  else
//  {
//#if ATK_FOUND
//    SLIC_ERROR("NeighborCommunication::UnpackTopologyModifications: inappropriate type for PhysicalDomainT::Unpack " + toString<int>(key));
//#endif
//  }

}

void NeighborCommunication::UnpackTopologyModifications( const string key ,
                                                         const char*& pbuffer,
                                                         std::map< std::string, lArray1d>& modifiedIndices )
{
//  if( key==PhysicalDomainT::FiniteElementElementManager )
//  {
//    lArray1d::size_type numModified;
//    bufvector::Unpack( pbuffer, numModified );
//
//    if( numModified > 0 )
//    {
//      m_domain->Unpack( key, pbuffer, modifiedIndices, false, false, true, true );
//    }
//  }
//  else
//  {
//#if ATK_FOUND
//    SLIC_ERROR("NeighborCommunication::UnpackTopologyModifications: inappropriate type for PhysicalDomainT::Unpack " + toString<int>(key));
//#endif
//  }

}





void NeighborCommunication::PackNewGlobalIndexRequests( const ModifiedObjectLists& modifiedObjects )
{
//  for( lSet::const_iterator i=modifiedObjects.newNodes.begin() ; i!=modifiedObjects.newNodes.end() ; ++i )
//  {
//    const globalIndex gIndex = m_domain->m_feNodeManager.m_localToGlobalMap[*i];
//    const globalIndex parent_gIndex = m_domain->m_feNodeManager.m_localToGlobalMap[m_domain->m_feNodeManager.GetParentIndex(*i)];
//    if( gIndex == GLOBALINDEX_MAX && GlobalIndexManager::OwningRank( parent_gIndex ) == this->m_neighborRank )
//    {
//      tempNeighborData.objectWithUnassignedGlobal[PhysicalDomainT::FiniteElementNodeManager].push_back(*i);
//    }
//  }
//
//  for( lSet::const_iterator i=modifiedObjects.newEdges.begin() ; i!=modifiedObjects.newEdges.end() ; ++i )
//  {
//    const globalIndex gIndex = m_domain->m_feEdgeManager.m_localToGlobalMap[*i];
//    const globalIndex parent_gIndex = m_domain->m_feEdgeManager.m_localToGlobalMap[m_domain->m_feEdgeManager.GetParentIndex(*i)];
//    if( gIndex == GLOBALINDEX_MAX && GlobalIndexManager::OwningRank( parent_gIndex ) == this->m_neighborRank )
//    {
//      tempNeighborData.objectWithUnassignedGlobal[PhysicalDomainT::FiniteElementEdgeManager].push_back(*i);
//    }
//  }
//
//  for( lSet::const_iterator i=modifiedObjects.newFaces.begin() ; i!=modifiedObjects.newFaces.end() ; ++i )
//  {
//    const globalIndex gIndex = m_domain->m_feFaceManager.m_localToGlobalMap[*i];
//    const globalIndex parent_gIndex = m_domain->m_feFaceManager.m_localToGlobalMap[m_domain->m_feFaceManager.GetParentIndex(*i)];
//    if( gIndex == GLOBALINDEX_MAX && GlobalIndexManager::OwningRank( parent_gIndex ) == this->m_neighborRank )
//    {
//      tempNeighborData.objectWithUnassignedGlobal[PhysicalDomainT::FiniteElementFaceManager].push_back(*i);
//    }
//  }
//
//
//  tempNeighborData.sendGlobalIndexRequests[0] = tempNeighborData.objectWithUnassignedGlobal[PhysicalDomainT::FiniteElementNodeManager].size();
//  tempNeighborData.sendGlobalIndexRequests[1] = tempNeighborData.objectWithUnassignedGlobal[PhysicalDomainT::FiniteElementEdgeManager].size();
//  tempNeighborData.sendGlobalIndexRequests[2] = tempNeighborData.objectWithUnassignedGlobal[PhysicalDomainT::FiniteElementFaceManager].size();

}





void NeighborCommunication::UnpackNewGlobalIndices(  )
{
//
//  {
//    globalIndex num = 0;
//    for( lArray1d::const_iterator i=tempNeighborData.objectWithUnassignedGlobal[PhysicalDomainT::FiniteElementNodeManager].begin() ;
//         i!=tempNeighborData.objectWithUnassignedGlobal[PhysicalDomainT::FiniteElementNodeManager].end() ; ++i, ++num )
//    {
//      m_domain->m_feNodeManager.m_localToGlobalMap[*i] = tempNeighborData.recvFirstNewGlobalIndices[0] + num;
//      m_domain->m_feNodeManager.m_globalToLocalMap[ m_domain->m_feNodeManager.m_localToGlobalMap[*i] ] = *i;
//    }
//  }
//
//  {
//    globalIndex num = 0;
//    for( lArray1d::const_iterator i=tempNeighborData.objectWithUnassignedGlobal[PhysicalDomainT::FiniteElementEdgeManager].begin() ;
//         i!=tempNeighborData.objectWithUnassignedGlobal[PhysicalDomainT::FiniteElementEdgeManager].end() ; ++i, ++num )
//    {
//      m_domain->m_feEdgeManager.m_localToGlobalMap[*i] = tempNeighborData.recvFirstNewGlobalIndices[1] + num;
//      m_domain->m_feEdgeManager.m_globalToLocalMap[ m_domain->m_feEdgeManager.m_localToGlobalMap[*i] ] = *i;
//    }
//  }
//
//  {
//    globalIndex num = 0;
//    for( lArray1d::const_iterator i=tempNeighborData.objectWithUnassignedGlobal[PhysicalDomainT::FiniteElementFaceManager].begin() ;
//         i!=tempNeighborData.objectWithUnassignedGlobal[PhysicalDomainT::FiniteElementFaceManager].end() ; ++i, ++num )
//    {
//      m_domain->m_feFaceManager.m_localToGlobalMap[*i] = tempNeighborData.recvFirstNewGlobalIndices[2] + num;
//      m_domain->m_feFaceManager.m_globalToLocalMap[ m_domain->m_feFaceManager.m_localToGlobalMap[*i] ] = *i;
//    }
//  }
//

}


void NeighborCommunication::ProcessNewGlobalIndexRequests( lSet& newNodeGlobals, lSet& newEdgeGlobals, lSet& newFaceGlobals )
{
//  const lArray1d::size_type numGlobalNodeRequests = tempNeighborData.recvGlobalIndexRequests[0];
//  const lArray1d::size_type numGlobalEdgeRequests = tempNeighborData.recvGlobalIndexRequests[1];
//  const lArray1d::size_type numGlobalFaceRequests = tempNeighborData.recvGlobalIndexRequests[2];
//
//  if( numGlobalNodeRequests > 0 )
//  {
//    for( lArray1d::size_type i=0 ; i<numGlobalNodeRequests ; ++i )
//    {
//      newNodeGlobals.insert( m_domain->m_feNodeManager.DataLengths() + i );
//    }
//    tempNeighborData.sendFirstNewGlobalIndices[0] = m_domain->m_feNodeManager.resize( m_domain->m_feNodeManager.DataLengths() + numGlobalNodeRequests, true );
//  }
//
//  if( numGlobalEdgeRequests > 0 )
//  {
//    for( lArray1d::size_type i=0 ; i<numGlobalEdgeRequests ; ++i )
//    {
//      newEdgeGlobals.insert( m_domain->m_feNodeManager.DataLengths() + i );
//    }
//    tempNeighborData.sendFirstNewGlobalIndices[1] = m_domain->m_feEdgeManager.resize( m_domain->m_feEdgeManager.DataLengths() + numGlobalEdgeRequests, true );
//  }
//
//  if( numGlobalFaceRequests > 0 )
//  {
//    for( lArray1d::size_type i=0 ; i<numGlobalFaceRequests ; ++i )
//    {
//      newFaceGlobals.insert( m_domain->m_feNodeManager.DataLengths() + i );
//    }
//    tempNeighborData.sendFirstNewGlobalIndices[2] = m_domain->m_feFaceManager.resize( m_domain->m_feFaceManager.DataLengths() + numGlobalFaceRequests, true );
//  }
}






int NeighborCommunication::UnpackGhosts(const string name)
{
  std::map<string, bufvector>::iterator it = tempNeighborData.objectsToReceive.find(name);
  if(it == tempNeighborData.objectsToReceive.end())
    return 0;
  const char* pbuffer = reinterpret_cast<const char*>(it->second.data());


//  m_domain->Unpack(name, pbuffer, m_receiveLocalIndices[name], true, true, true, true );
  return 1;
}

void NeighborCommunication::UnpackGhostElements( std::map<std::string,lArray1d>& newElementIndices )
{
//
//
//  m_domain->m_feElementManager.UnpackElements( tempNeighborData.objectsToReceive[PhysicalDomainT::FiniteElementElementManager],
//                                             m_domain->m_feNodeManager,
//                                             m_domain->m_feFaceManager,
//                                             m_elementRegionsReceiveLocalIndices,
//                                             true, true, true, true );
//
//  for( std::map< std::string, lArray1d>::iterator i=m_elementRegionsReceiveLocalIndices.begin();
//      i!=m_elementRegionsReceiveLocalIndices.end() ;
//      ++i )
//  {
//    newElementIndices[i->first].insert( newElementIndices[i->first].end(), i->second.begin(), i->second.end() );
//  }
}

int NeighborCommunication::UnpackGhosts(const string name, lArray1d& newIndices )
{
  //for edge, face, and de_face
  int i = UnpackGhosts(name);
  if(i == 0)
    return i;

//  newIndices.insert(newIndices.end(),m_receiveLocalIndices[name].begin(),m_receiveLocalIndices[name].end());
//  if(name == PhysicalDomainT::FiniteElementFaceManager)
//  {
//    m_receiveLocalIndices[PhysicalDomainT::VirtualFaceManager] = m_receiveLocalIndices[name];
//    std::sort(newIndices.begin(),newIndices.end());
//    // now remove the duplicates
//    lArray1d::iterator iend = std::unique(newIndices.begin(),newIndices.end());
//    newIndices.resize( iend - newIndices.begin() );
//  }

  return i;
}


bufvector::size_type NeighborCommunication::PackBuffer( const std::map<string, sArray1d>& fieldNames,
                                                        const CommRegistry::commID commID,
                                                        const bool doBufferPacking )
{
  //TODO though this is now somewhat checked it is still ridiculously unsafe ... we should change it
  //doBufferPacking==0 implies that the m_sendBuffer will not be written to, which is why we cannot resize
  //it here ... why this was working before is unclear, but m_sendSize _is_ uninitialized during
  //ProblemManagerT::CompleteObjectInitialization() when this is called and doBufferPacking=0
  if(doBufferPacking)
    m_sendBuffer.resize(this->m_sendSize[commID]);

  char* buffer = m_sendBuffer.data();
  bufvector::size_type bufferSize = 0;

  // pack node contributions to the buffer
  // pack element contributions to the buffer
  // pack face contributions to the buffer
  // pack external face contributions to the buffer
  // pack virtual face contributions to the buffer

//  const Array1dT<ObjectDataStructureBaseT*>& omap = m_domain->GetObjectDataStructures();
//  std::map< PhysicalDomainT::ObjectDataStructureKeys, sArray1d>::const_iterator it;
//  for(it = fieldNames.begin(); it != fieldNames.end(); ++it)
//  {
//    if(it->first == PhysicalDomainT::FiniteElementElementManager)
//    {
//      for( std::map< ElementManagerT::RegKeyType, ElementRegionT >::iterator ielemReg=m_domain->m_feElementManager.m_ElementRegions.begin() ;
//          ielemReg!=m_domain->m_feElementManager.m_ElementRegions.end() ; ++ielemReg )
//      {
//        const std::string& elementRegionName = ielemReg->first;
//        const ElementRegionT& elementRegion = ielemReg->second;
//        const lArray1d* const elementRegionSendLocalIndices = stlMapLookupPointer( m_elementRegionsSendLocalIndices,
//                                                                                   elementRegionName );
//        if( elementRegionSendLocalIndices )
//          bufferSize += elementRegion.PackFieldsIntoBuffer( buffer,
//                                                            it->second,
//                                                            *elementRegionSendLocalIndices,
//                                                            doBufferPacking );
//      }
//    }
//    else
//    {
//      bufferSize += m_domain->GetObjectDataStructure(it->first).PackFieldsIntoBuffer( buffer,
//                                                                                      it->second,
//                                                                                      m_sendLocalIndices[it->first],
//                                                                                      doBufferPacking );
//#if ATK_FOUND
//      //SLIC_ERROR("Trying for name (" + toString<int>(it->first) + ") size=" + toString<bufvector::size_type>(bufferSize));
//#endif
//    }
//  }
//
//  if( doBufferPacking && bufferSize != m_sendBuffer.size() )
//  {
//    printf( " rank %5i: bufferSize = %6lu, m_sendBuffer.size() = %6lu \n",this->m_rank, bufferSize, m_sendBuffer.size() );
//#if ATK_FOUND
//    SLIC_ERROR("m_sendBuffer.size() isn't what it should be\n");
//#endif
//  }

  return bufferSize;
}

void NeighborCommunication::UnpackBuffer( const std::map<string, sArray1d>& fieldNames)
{
  const char* pbuffer = m_receiveBuffer.data();

  // unpack node contributions to the buffer
  // unpack element contributions to the buffer
  // unpack face contributions to the buffer
  // unpack external face contributions to the buffer
  // unpack virtual face contributions to the buffer

//  const Array1dT<ObjectDataStructureBaseT*>& omap = m_domain->GetObjectDataStructures();
//  std::map< PhysicalDomainT::ObjectDataStructureKeys, sArray1d>::const_iterator it;
//  for(it = fieldNames.begin(); it != fieldNames.end(); ++it)
//  {
//    if(it->first == PhysicalDomainT::FiniteElementElementManager)
//    {
//      for( std::map< ElementManagerT::RegKeyType, ElementRegionT >::iterator ielemReg=m_domain->m_feElementManager.m_ElementRegions.begin() ;
//           ielemReg!=m_domain->m_feElementManager.m_ElementRegions.end() ; ++ielemReg )
//      {
//        const std::string& elementRegionName = ielemReg->first;
//        ElementRegionT& elementRegion = ielemReg->second;
//        const lArray1d* const elementRegionReceiveLocalIndices = stlMapLookupPointer( m_elementRegionsReceiveLocalIndices,
//                                                                                      elementRegionName );
//
//        if( elementRegionReceiveLocalIndices )
//        elementRegion.UnpackFieldsFromBuffer( pbuffer,
//                                              it->second,
//                                              *elementRegionReceiveLocalIndices );
//      }
//    }
//    else
//    {
//      // pack node contributions to the buffer
//      m_domain->GetObjectDataStructure(it->first).UnpackFieldsFromBuffer( pbuffer,
//                                               it->second,
//                                               m_receiveLocalIndices[it->first] );
//    }
//  }
}



void NeighborCommunication::SendReceive(char* inBuffer,
                                        const int inSize,
                                        char* outBuffer,
                                        const int outSize,
                                        MPI_Request& inReq,
                                        MPI_Request& outReq,
                                        const CommRegistry::commID commID )
{
  const int stag = CommRegistry::CommTag(m_rank, m_neighborRank, commID);
  //m_rank * m_size + m_neighborRank + m_size*m_size*commID;
  MPI_Isend( outBuffer, outSize, MPI_CHAR, m_neighborRank, stag, MPI_COMM_WORLD, &outReq);

  const int rtag = CommRegistry::CommTag(m_neighborRank, m_rank, commID);
  //m_neighborRank * m_size + m_rank + m_size*m_size*commID;
  MPI_Irecv( inBuffer, inSize, MPI_CHAR, m_neighborRank, rtag, MPI_COMM_WORLD, &inReq );
}

void NeighborCommunication::SendReceiveSizes(const CommRegistry::commID commID)
{
  this->m_sendSize[commID] = m_sendBuffer.size();
  SendReceive( reinterpret_cast<char*>(&m_receiveSize[commID]), sizeof(int),
               reinterpret_cast<char*>(&m_sendSize[commID]), sizeof(int),
               mpiRecvSizeRequest[commID], mpiSendSizeRequest[commID], commID);
}

void NeighborCommunication::SendReceiveBufferSizes( const CommRegistry::commID commID,
                                                    MPI_Request& mpiSendRequest,
                                                    MPI_Request& mpiRecvRequest )
{
  this->m_sendSize[commID] = m_sendBuffer.size();
  SendReceive( &m_sendSize[commID], 1, mpiSendRequest,
               &m_receiveSize[commID], 1, mpiRecvRequest, commID);
}


void NeighborCommunication::SendReceiveBuffers(const CommRegistry::commID commID)
{
  this->m_receiveBuffer.resize( this->m_receiveSize[commID] );
  SendReceive(m_receiveBuffer.data(), m_receiveBuffer.size(),
              m_sendBuffer.data(), m_sendBuffer.size(),
              mpiRecvBufferRequest[commID], mpiSendBufferRequest[commID], commID);
}

void NeighborCommunication::SendReceiveBuffers(const CommRegistry::commID commID,
                                               MPI_Request& mpiSendRequest,
                                               MPI_Request& mpiRecvRequest)
{
  this->m_receiveBuffer.resize( this->m_receiveSize[commID] );
  SendReceive( m_receiveBuffer.data(), m_receiveBuffer.size(),
               m_sendBuffer.data(), m_sendBuffer.size(),
               mpiRecvRequest, mpiSendRequest, commID);

}

void NeighborCommunication::SendReceiveBuffers(const CommRegistry::commID sizeCommID,
                                               const CommRegistry::commID commID,
                                               MPI_Request& mpiSendRequest,
                                               MPI_Request& mpiRecvRequest)
{
  this->m_receiveBuffer.resize( this->m_receiveSize[sizeCommID] );
  SendReceive( m_receiveBuffer.data(), m_receiveBuffer.size(),
               m_sendBuffer.data(), m_sendBuffer.size(),
               mpiRecvRequest, mpiSendRequest, commID);

}
//
//void NeighborCommunication::WriteSilo( SiloFile& siloFile )
//{
//  siloFile.DBWriteWrapper("m_neighborRank",m_neighborRank);
//  siloFile.DBWriteWrapper("m_rank",m_rank);
//  siloFile.DBWriteWrapper("m_size",m_size);
//  siloFile.DBWriteWrapper("m_rankOfNeighborNeighbors",m_rankOfNeighborNeighbors);
//
//  std::map<std::string, lArray1d> sendLocalIndices;
//  std::map<std::string, lArray1d> receiveLocalIndices;
//
//  for( std::map<PhysicalDomainT::ObjectDataStructureKeys, lArray1d>::const_iterator i=m_sendLocalIndices.begin() ;
//      i!= m_sendLocalIndices.end() ; ++i )
//  {
////    if( !(i->second.empty()) )
//    {
//      sendLocalIndices[ PhysicalDomainT::GetObjectDataStructureName( i->first ) ] = i->second;
//    }
//  }
//
//  for( std::map<PhysicalDomainT::ObjectDataStructureKeys, lArray1d>::const_iterator i=m_receiveLocalIndices.begin() ;
//      i!= m_receiveLocalIndices.end() ; ++i )
//  {
////    if( !(i->second.empty()) )
//    {
//      receiveLocalIndices[ PhysicalDomainT::GetObjectDataStructureName( i->first ) ] = i->second;
//    }
//  }
//
//  siloFile.DBWriteWrapper("sendLocalIndices",sendLocalIndices);
//  siloFile.DBWriteWrapper("m_elementRegionsSendLocalIndices",m_elementRegionsSendLocalIndices);
//  siloFile.DBWriteWrapper("receiveLocalIndices",receiveLocalIndices);
//  siloFile.DBWriteWrapper("m_elementRegionsReceiveLocalIndices",m_elementRegionsReceiveLocalIndices);
//
//  siloFile.DBWriteWrapper("m_sendSize",m_sendSize,CommRegistry::maxComm);
//  siloFile.DBWriteWrapper("m_receiveSize",m_receiveSize,CommRegistry::maxComm);
//
//
//}
//
//void NeighborCommunication::ReadSilo( const SiloFile& siloFile )
//{
//  siloFile.DBReadWrapper("m_neighborRank",m_neighborRank);
//  siloFile.DBReadWrapper("m_rank",m_rank);
//  siloFile.DBReadWrapper("m_size",m_size);
//  siloFile.DBReadWrapper("m_rankOfNeighborNeighbors",m_rankOfNeighborNeighbors);
//
//  std::map<std::string, lArray1d> sendLocalIndices;
//  std::map<std::string, lArray1d> receiveLocalIndices;
//
//
///*
//  for( std::map<PhysicalDomainT::ObjectDataStructureKeys, lArray1d>::const_iterator i=m_sendLocalIndices.begin() ;
//      i!= m_sendLocalIndices.end() ; ++i )
//  {
//    if( !(i->second.empty()) )
//    {
//      sendLocalIndices[ PhysicalDomainT::GetObjectDataStructureName( i->first ) ];
//    }
//  }
//
//  for( std::map<PhysicalDomainT::ObjectDataStructureKeys, lArray1d>::const_iterator i=m_receiveLocalIndices.begin() ;
//      i!= m_receiveLocalIndices.end() ; ++i )
//  {
//    if( !(i->second.empty()) )
//    {
//      receiveLocalIndices[ PhysicalDomainT::GetObjectDataStructureName( i->first ) ];
//    }
//  }
//*/
//
//
//
//
//
//  siloFile.DBReadWrapper("sendLocalIndices",sendLocalIndices);
//  siloFile.DBReadWrapper("m_elementRegionsSendLocalIndices",m_elementRegionsSendLocalIndices);
//  siloFile.DBReadWrapper("receiveLocalIndices",receiveLocalIndices);
//  siloFile.DBReadWrapper("m_elementRegionsReceiveLocalIndices",m_elementRegionsReceiveLocalIndices);
//
//
//  siloFile.DBReadWrapper("m_sendSize",m_sendSize,CommRegistry::maxComm);
//  siloFile.DBReadWrapper("m_receiveSize",m_receiveSize,CommRegistry::maxComm);
//
//
//  for( std::map<std::string, lArray1d>::const_iterator i=sendLocalIndices.begin() ;
//      i!= sendLocalIndices.end() ; ++i )
//  {
//    m_sendLocalIndices[ PhysicalDomainT::GetObjectDataStructureKey(i->first)] = i->second;
//  }
//
//  for( std::map<std::string, lArray1d>::const_iterator i=receiveLocalIndices.begin() ;
//      i!= receiveLocalIndices.end() ; ++i )
//  {
//    m_receiveLocalIndices[ PhysicalDomainT::GetObjectDataStructureKey(i->first)] = i->second;
//  }
//
//
//}
/* Rui Wang, return the size of the neighbor communication to determine whether the neighbor
 * is excessive.
 */
int NeighborCommunication::ReturnNeighborRcvSndSize()
  {
//  int size = m_receiveLocalIndices[PhysicalDomainT::FiniteElementNodeManager].size();
//  size += m_receiveLocalIndices[PhysicalDomainT::FiniteElementEdgeManager].size();
//  size += m_receiveLocalIndices[PhysicalDomainT::FiniteElementFaceManager].size();
//  size += m_sendLocalIndices[PhysicalDomainT::FiniteElementNodeManager].size();
//  size += m_sendLocalIndices[PhysicalDomainT::FiniteElementEdgeManager].size();
//  size += m_sendLocalIndices[PhysicalDomainT::FiniteElementFaceManager].size();
//    return size;
  return 0;
  }
}
