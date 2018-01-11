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
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
// SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
// TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S.
// Department of Energy (DOE). This work was produced at Lawrence Livermore
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National
// Security, LLC nor any of their employees, makes any warranty, express or
//     implied, or assumes any liability or responsibility for the accuracy,
// completeness, or usefulness of any information, apparatus, product, or
//     process disclosed, or represents that its use would not infringe
// privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or
// services by trade name, trademark, manufacturer or otherwise does not
//     necessarily constitute or imply its endorsement, recommendation, or
// favoring by the United States Government or Lawrence Livermore National
// Security,
//     LLC. The views and opinions of authors expressed herein do not
// necessarily state or reflect those of the United States Government or
// Lawrence
//     Livermore National Security, LLC, and shall not be used for advertising
// or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The
// BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @file PartitionBase.cpp
 * @author settgast1
 * @date Mar 16, 2011
 */

#include "PartitionBase.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/NodeManager.hpp"

//#include "IO/BinStream.h"
//#include "ObjectManagers/DomainPartition.h"
#include <limits.h>

#ifdef USE_ATK
#include "slic/slic.hpp"
#endif

namespace geosx
{

using namespace dataRepository;

PartitionBase::PartitionBase( ):
  m_domain(NULL),
  m_hasLocalGhosts(false)
{
  //maxComm
}


PartitionBase::PartitionBase( const unsigned int numPartitions, const unsigned int thisPartiton ):
  m_size(numPartitions),
  m_rank(thisPartiton),
  m_color(0),
  m_numColors(1),
  m_domain(NULL),
  m_t1(0.0),
  m_t2(0.0),
  m_t3(0.0),
  m_t4(0.0),
  m_hasLocalGhosts(false)
{
  //maxComm
}


PartitionBase::~PartitionBase()
{}

/**
 * @brief Call SetDomain on each neighbor
 */
void PartitionBase::SetDomain( DomainPartition * domain )
{
  // set the const pointer "m_domain" by casting away the the const
  // on the address of the pointer, and modifying what the address of
  // the pointer.
  DomainPartition** temp = const_cast<DomainPartition**>(&m_domain);
  *temp = domain;

  for( array<NeighborCommunication>::iterator neighbor=m_neighbors.begin() ; neighbor!=m_neighbors.end() ; ++neighbor )
  {
//    neighbor->SetDomain( domain );
  }
}

///**
// * @brief Call Clear on each neighbor and set up neighbor lists
// */
//void PartitionBase::ResetNeighborLists( DomainPartition * domain,
//                                        const int elementGhostingDepth )
//{
//  for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin() ;
// neighbor!=m_neighbors.end() ; ++neighbor )
//  {
//    neighbor->Clear();
//  }
//  SetUpNeighborLists( domain, 0, elementGhostingDepth );
//}


void PartitionBase::AssignGlobalIndices( DomainPartition * domain )
{
//
//  for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin() ;
// neighbor!=m_neighbors.end() ; ++neighbor )
//  {
//    neighbor->Clear();
//  }
//
//  m_mpiRequest.resize(2 * m_neighbors.size());
//  m_mpiStatus.resize(2 * m_neighbors.size());
//
//  //------------------------------------------
//  //(1) create list of boundary nodes on local process:
//
//
//  // array of global nodes on the boundary by old global index
//  localIndex_array localBoundaryNodesLocalIndex;
//  gArray1d localBoundaryNodesGlobalIndex;
//  gArray1d localBoundaryNodesNewGlobalIndex;
//
//  domain->m_feNodeManager.ConstructListOfBoundaryObjects(
// localBoundaryNodesLocalIndex );
//
//  for( localIndex_array::const_iterator a=localBoundaryNodesLocalIndex.begin()
// ; a!=localBoundaryNodesLocalIndex.end() ; ++a )
//  {
//    localBoundaryNodesGlobalIndex.push_back(
// domain->m_feNodeManager.m_localToGlobalMap[*a] );
//    localBoundaryNodesNewGlobalIndex.push_back( GlobalIndexManager::Index(
// this->m_rank, *a ) );
//  }
//  //(1)
//
//  //------------------------------------------
//  //(2) SENDRECV the sizes of the lists of boundary nodes
//  //    loop over all neighbors and communicate size of boundary nodes and
// faces
//  int i_mpireq = 0;
//  {
//    gArray1d::size_type size = localBoundaryNodesGlobalIndex.size();
//    for (VectorT<NeighborCommunication>::iterator neighbor =
// m_neighbors.begin();
//        neighbor != m_neighbors.end(); ++neighbor)
//    {
//      neighbor->tempNeighborData.neighborBoundaryObjectsSizes.resize(1, 0);
//      neighbor->SendReceive(
//          reinterpret_cast<char*>(neighbor->tempNeighborData.neighborBoundaryObjectsSizes.data()),
//          sizeof(gArray1d::size_type), reinterpret_cast<char*>(&size),
// sizeof(gArray1d::size_type),
//          m_mpiRequest[i_mpireq], m_mpiRequest[i_mpireq + 1]);
//      i_mpireq += 2;
//    }
//    MPI_Waitall(m_mpiRequest.size(), m_mpiRequest.data(), m_mpiStatus.data());
//  }
//  //(2)
//
//
//  m_mpiRequest.resize(2 * 2 * m_neighbors.size());
//  m_mpiStatus.resize(2 * 2 * m_neighbors.size());
//
//  //------------------------------------------
//  //(3) SENDRECV the existing global Indices of the nodes, given by the mesh
//  i_mpireq = 0;
//  for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin() ;
// neighbor!=m_neighbors.end() ; ++neighbor )
//  {
//    //(3a) initialize and resize the "tempNeighborData.neighborNumbers"
// buffers to receives objects
//    neighbor->tempNeighborData.neighborNumbers[DomainPartition::FiniteElementNodeManager].resize(neighbor->tempNeighborData.neighborBoundaryObjectsSizes[0]);
//    neighbor->tempNeighborData.neighborNumbers[DomainPartition::FiniteElementEdgeManager].resize(neighbor->tempNeighborData.neighborBoundaryObjectsSizes[0]);
//
//    //(3b) receive the "tempNeighborData.neighborNumbers" buffers
//    neighbor->SendReceive(
// reinterpret_cast<char*>(neighbor->tempNeighborData.neighborNumbers[DomainPartition::FiniteElementNodeManager].data()),
//                           neighbor->tempNeighborData.neighborNumbers[DomainPartition::FiniteElementNodeManager].size()*sizeof(globalIndex),
//                           reinterpret_cast<char*>(localBoundaryNodesGlobalIndex.data()),
//                           localBoundaryNodesNewGlobalIndex.size()*sizeof(globalIndex),
//                           m_mpiRequest[i_mpireq], m_mpiRequest[i_mpireq+1] );
//
//    neighbor->SendReceive(
// reinterpret_cast<char*>(neighbor->tempNeighborData.neighborNumbers[DomainPartition::FiniteElementEdgeManager].data()),
//                           neighbor->tempNeighborData.neighborNumbers[DomainPartition::FiniteElementEdgeManager].size()*sizeof(globalIndex),
//                           reinterpret_cast<char*>(localBoundaryNodesNewGlobalIndex.data()),
//                           localBoundaryNodesNewGlobalIndex.size()*sizeof(globalIndex),
//                           m_mpiRequest[i_mpireq+2], m_mpiRequest[i_mpireq+3]
// );
//
//    i_mpireq += 4;
//  }
//  MPI_Waitall( m_mpiRequest.size(), m_mpiRequest.data(), m_mpiStatus.data() );
//  //(3)
//
//
//  //------------------------------------------
//  //(4) SET NEW GLOBAL NODE NUMBERS
//  {
//    std::map<globalIndex,globalIndex> globalIndexTransformation;
//
//    for( localIndex a=0 ; a<domain->m_feNodeManager.m_numNodes ; ++a )
//    {
//      globalIndexTransformation.insert( std::make_pair(
// domain->m_feNodeManager.m_localToGlobalMap[a],
// GlobalIndexManager::Index(this->m_rank,a) ) );
//    }
//
//    // now check the global num
//    for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin()
// ; neighbor!=m_neighbors.end() ; ++neighbor )
//    {
//      if( neighbor->NeighborRank() < this->m_rank )
//      {
//        gArray1d& neighborOldGlobal =
// neighbor->tempNeighborData.neighborNumbers[DomainPartition::FiniteElementNodeManager];
//        gArray1d& neighborNewGlobal =
// neighbor->tempNeighborData.neighborNumbers[DomainPartition::FiniteElementEdgeManager];
//        for( gArray1d::size_type a=0 ; a<neighborOldGlobal.size() ; ++a )
//        {
//          std::map<globalIndex,globalIndex>::iterator
// iter_globalIndexTransformation =
// globalIndexTransformation.find(neighborOldGlobal[a]);
//          if( iter_globalIndexTransformation !=
// globalIndexTransformation.end() )
//          {
//            if( iter_globalIndexTransformation->second > neighborNewGlobal[a]
// )
//              iter_globalIndexTransformation->second = neighborNewGlobal[a];
//          }
//        }
//      }
//    }
//
//    for( localIndex a=0 ; a<domain->m_feNodeManager.m_numNodes ; ++a )
//    {
//      domain->m_feNodeManager.m_localToGlobalMap[a] =
// globalIndexTransformation[domain->m_feNodeManager.m_localToGlobalMap[a]];
//    }
//    domain->m_feNodeManager.ResetGlobalToLocalMap();
//  }
//  //(4)
//
//  //(5) Now do the rest of the objects
//  array<DomainPartition::ObjectDataStructureKeys> objectNames;
//  NeighborCommunication::SyncNames(objectNames);
//
//  for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin() ;
// neighbor!=m_neighbors.end() ; ++neighbor )
//  {
//    neighbor->Clear();
//  }
//  AssignGlobalIndices( domain->m_feFaceManager, domain->m_feNodeManager );
//
//  for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin() ;
// neighbor!=m_neighbors.end() ; ++neighbor )
//  {
//    neighbor->Clear();
//  }
//  AssignGlobalIndices( domain->m_feEdgeManager, domain->m_feNodeManager );
//
//
//  // set the global indices of elements
//  for( std::map< ElementManagerT::RegKeyType, ElementRegionT >::iterator
// iter_elemReg=domain->m_feElementManager.m_ElementRegions.begin() ;
// iter_elemReg!=domain->m_feElementManager.m_ElementRegions.end();
// ++iter_elemReg )
//  {
//    ElementRegionT& elemRegion = iter_elemReg->second;
//    for( localIndex a=0 ; a<elemRegion.DataLengths() ; ++a)
//    {
//      elemRegion.m_localToGlobalMap[a] =
// GlobalIndexManager::Index(this->m_rank,a) ;
//    }
//  }
//
//  // set the global indices for the discrete element surface nodes
//  for( localIndex a=0 ; a<domain->m_discreteElementSurfaceNodes.DataLengths()
// ; ++a )
//  {
//    domain->m_discreteElementSurfaceNodes.m_localToGlobalMap[a] =
// GlobalIndexManager::Index(this->m_rank,a);
//  }
//
//  for( localIndex a=0 ; a<domain->m_discreteElementSurfaceFaces.DataLengths()
// ; ++a )
//  {
//    domain->m_discreteElementSurfaceFaces.m_localToGlobalMap[a] =
// GlobalIndexManager::Index(this->m_rank,a);
//  }
//
//  // set the global indices for the discrete elements
//  for( localIndex a=0 ; a<domain->m_discreteElementManager.DataLengths() ; ++a
// )
//  {
//    domain->m_discreteElementManager.m_localToGlobalMap[a] =
// GlobalIndexManager::Index(this->m_rank,a);
//  }
//
//  // set the global indices for the ellipsoidal discrete elements
//  for( localIndex a=0 ;
// a<domain->m_ellipsoidalDiscreteElementManager.DataLengths() ; ++a )
//  {
//    domain->m_ellipsoidalDiscreteElementManager.m_localToGlobalMap[a] =
// GlobalIndexManager::Index(this->m_rank,a);
//  }
//
//#ifdef SRC_EXTERNAL
//  // set the global indices for the fault nodes
//  for( localIndex a=0 ; a<domain->m_faultPatchNodes.DataLengths() ; ++a )
//  {
//    domain->m_faultPatchNodes.m_localToGlobalMap[a] =
// GlobalIndexManager::Index(this->m_rank,a);
//  }
//
//  // set the global indices for the fault elements
//  for( localIndex a=0 ; a<domain->m_faultPatchFaces.DataLengths() ; ++a )
//  {
//    domain->m_faultPatchFaces.m_localToGlobalMap[a] =
// GlobalIndexManager::Index(this->m_rank,a);
//  }
//#endif
}

/**
 *
 * @param object
 * @param compositionObject
 * assign global indices for
 */
void PartitionBase::AssignGlobalIndices( ObjectDataStructureBaseT& object, const ObjectDataStructureBaseT& compositionObject )
{
//
//  // set the global indices as if they were all local to this process
//  for( localIndex a=0 ; a<object.DataLengths() ; ++a )
//  {
//    object.m_localToGlobalMap[a] = GlobalIndexManager::Index(this->m_rank,a) ;
//  }
//
//
//  // get the relation to the composition object used that will be used to
// identify the main object. For example,
//  // a face can be identified by its nodes.
//  array<gArray1d> objectToCompositionObject;
//  object.ExtractMapFromObjectForAssignGlobalObjectNumbers( compositionObject,
// objectToCompositionObject );
//  gArray1d objectToCompositionObjectBuffer;
//
//
//  // now arrange the data from objectToCompositionObject into a map
// "indexByFirstCompositionIndex", such that the key
//  // is the lowest global index of the composition object that make up this
// object. The value of the map is a pair, with the
//  // array being the remaining composition object global indices, and the
// second being the global index of the object
//  // itself.
//  std::map<globalIndex, array<std::pair<gArray1d,localIndex> > >
// indexByFirstCompositionIndex;
//  for( array<gArray1d>::const_iterator a=objectToCompositionObject.begin() ;
// a!=objectToCompositionObject.end() ; ++a )
//  {
//    // global index of the object
//    const globalIndex gIndex = (*a)[0];
//
//    // grab the first global index of the composition objects
//    const globalIndex firstCompositionIndex = (*a)[1];
//
//    // create a temporary to hold the pair
//    std::pair<gArray1d,localIndex> tempComp;
//
//    // fill the array with the remaining composition object global indices
//    tempComp.first.insert( tempComp.first.begin(), a->begin()+2, a->end() );
//
//    // set the second value of the pair to the localIndex of the object.
//    tempComp.second = GlobalIndexManager::OwningIndex(gIndex);
//
//    // push the tempComp onto the map.
//    indexByFirstCompositionIndex[firstCompositionIndex].push_back(tempComp);
//
//  }
//
//
//  // put the map into a buffer
//  for( array<gArray1d>::const_iterator a=objectToCompositionObject.begin() ;
// a!=objectToCompositionObject.end() ; ++a )
//  {
//    objectToCompositionObjectBuffer.push_back( a->size() );
//    for( gArray1d::const_iterator b=a->begin() ; b!=a->end() ; ++b )
//    {
//      objectToCompositionObjectBuffer.push_back( *b );
//    }
//  }
//
//
//  // send the sizes
//  {
//    m_mpiRequest.resize(2 * m_neighbors.size());
//    m_mpiStatus.resize(2 * m_neighbors.size());
//    int i_mpireq = 0;
//    gArray1d::size_type sendsize = objectToCompositionObjectBuffer.size();
//    for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin()
// ; neighbor!=m_neighbors.end() ; ++neighbor )
//    {
//      neighbor->Clear();
//      neighbor->tempNeighborData.neighborBoundaryObjectsSizes.resize(1, 0);
//      neighbor->SendReceive(
// reinterpret_cast<char*>(neighbor->tempNeighborData.neighborBoundaryObjectsSizes.data()),
//                             sizeof(gArray1d::size_type),
//                             reinterpret_cast<char*>(&sendsize),
//                             sizeof(gArray1d::size_type),
//                             m_mpiRequest[i_mpireq], m_mpiRequest[i_mpireq+1]
// );
//      i_mpireq += 2;
//    }
//    MPI_Waitall( m_mpiRequest.size() , m_mpiRequest.data(), m_mpiStatus.data()
// );
//  }
//
//
//  // send/receive the actual data, using
// neighbor->tempNeighborData.neighborNumbers[DomainPartition::FiniteElementNodeManager]
//  // to store the neighbor data
//  {
//    int i_mpireq = 0;
//    for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin()
// ; neighbor!=m_neighbors.end() ; ++neighbor )
//    {
//      neighbor->tempNeighborData.neighborNumbers[DomainPartition::FiniteElementNodeManager].resize(neighbor->tempNeighborData.neighborBoundaryObjectsSizes[0]);
//      neighbor->SendReceive(reinterpret_cast<char*>(neighbor->tempNeighborData.neighborNumbers[DomainPartition::FiniteElementNodeManager].data()),
//                            neighbor->tempNeighborData.neighborNumbers[DomainPartition::FiniteElementNodeManager].size()*sizeof(gArray1d::size_type),
//                            reinterpret_cast<char*>(objectToCompositionObjectBuffer.data()),
//                            objectToCompositionObjectBuffer.size()*sizeof(gArray1d::size_type),
//                            m_mpiRequest[i_mpireq], m_mpiRequest[i_mpireq+1]
// );
//      i_mpireq += 2;
//    }
//    MPI_Waitall( m_mpiRequest.size() , m_mpiRequest.data(), m_mpiStatus.data()
// );
//  }
//
//
//  // unpack the data from
// neighbor->tempNeighborData.neighborNumbers[DomainPartition::FiniteElementNodeManager]
// to
//  // the local arrays
//
//  // object to receive the neighbor data
//  // this baby is and Array (for each neighbor) of maps, with the key of
// lowest composition index, and a value containing
//  // an array containing the std::pairs of the remaining composition indices,
// and the globalIndex of the object.
//  array< std::map<globalIndex, array<std::pair<gArray1d,globalIndex> > > >
// neighborCompositionObjects( this->m_neighbors.size() );
//
//  {
//    // we are going to need a counter for the neighbors
//    int neighborNum = 0;
//    for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin()
// ; neighbor!=m_neighbors.end() ; ++neighbor, ++neighborNum )
//    {
//      // iterate over data that was just received
//      for( gArray1d::const_iterator
// a=neighbor->tempNeighborData.neighborNumbers[DomainPartition::FiniteElementNodeManager].begin()
// ;
//           a!=neighbor->tempNeighborData.neighborNumbers[DomainPartition::FiniteElementNodeManager].end()
// ; ++a )
//      {
//        // the first thing packed was the data size for a given object
//        const gArray1d::size_type dataSize = *(a++);
//
//        // the second thing packed was the globalIndex of that object
//        const globalIndex neighborGlobalIndex = *(a++);
//
//        // the global indices of the composition objects were next. they are
// ordered, so the lowest one is first.
//        const globalIndex firstCompositionIndex = *(a);
//
//
//        // the remaining composition object indices.
//        gArray1d temp;
//        for( gArray1d::size_type b=2 ; b<dataSize ; ++b )
//        {
//          ++a;
//          temp.push_back(*a);
//        }
//
//        // fill neighborCompositionObjects
//        std::pair<gArray1d,globalIndex> tempComp(
// std::make_pair(temp,neighborGlobalIndex) );
//        neighborCompositionObjects[neighborNum][firstCompositionIndex].push_back(tempComp);
//
//      }
//    }
//  }
//
//
//  // now check to see if the global index is valid. We do this by checking the
// contents of neighborCompositionObjects
//  // with indexByFirstCompositionIndex
//  int neighborNum = 0;
//  for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin() ;
// neighbor!=m_neighbors.end() ; ++neighbor, ++neighborNum )
//  {
//    // it only matters if the neighbor rank is lower than this rank
//    if( neighbor->NeighborRank() < this->m_rank )
//    {
//      // now we are going to need to do a pretty tricky loop. Set iterators to
// the beginning of each indexByFirstCompositionIndex,
//      // and neighborCompositionObjects[neighborNum].
//      std::map<globalIndex, array<std::pair<gArray1d,localIndex> >
// >::const_iterator iter_local = indexByFirstCompositionIndex.begin();
//      std::map<globalIndex, array<std::pair<gArray1d,globalIndex> >
// >::const_iterator iter_neighbor =
//  neighborCompositionObjects[neighborNum].begin();
//
//      // now we continue the while loop as long as both of our iterators are
// in range.
//      while( iter_local!=indexByFirstCompositionIndex.end() &&
// iter_neighbor!=neighborCompositionObjects[neighborNum].end() )
//      {
//        // check to see if the map keys (first composition index) are the
// same.
//        if( iter_local->first == iter_neighbor->first )
//        {
//          // first we loop over all local composition arrays (objects with the
// matched key)
//          for( array<std::pair<gArray1d,localIndex> >::const_iterator
// iter_local2 = iter_local->second.begin() ;
//              iter_local2 != iter_local->second.end() ; ++iter_local2 )
//          {
//            // and loop over all of the neighbor composition arrays (objects
// with the matched key)
//            for( array<std::pair<gArray1d,globalIndex> >::const_iterator
// iter_neighbor2 = iter_neighbor->second.begin() ;
//                iter_neighbor2 != iter_neighbor->second.end() ;
// ++iter_neighbor2 )
//            {
//              // now compare the composition arrays
//              if( iter_local2->first.size() == iter_neighbor2->first.size() &&
//                  std::equal( iter_local2->first.begin(),
// iter_local2->first.end(), iter_neighbor2->first.begin() ) )
//              {
//                // they are equal, so we need to overwrite the global index
// for the object
//                if( iter_neighbor2->second <
// object.m_localToGlobalMap[iter_local2->second] )
//                {
//                  object.m_localToGlobalMap[iter_local2->second] =
// iter_neighbor2->second;
//                }
//
//                // we should break out of the iter_local2 loop since we aren't
// going to find another match.
//                break;
//              }
//            }
//          }
//
//          ++iter_local;
//          ++iter_neighbor;
//
//        }
//        else if( iter_local->first < iter_neighbor->first )
//        {
//          ++iter_local;
//        }
//        else if( iter_local->first > iter_neighbor->first )
//        {
//          ++iter_neighbor;
//        }
//
//      }
//    }
//
//  }

}


template< typename T >
void PartitionBase::SendReceive( const array<array<T> >& sendArray, array<array<T> >& recvArray )
{

  if( sendArray.size() != m_neighbors.size() || recvArray.size() != m_neighbors.size() )
  {
#ifdef USE_ATK
    SLIC_ERROR("PartitionBase::SendRecieve: size of arrays do not equal number of neighbors");
#endif
  }


  array<MPI_Request> mpiSendRequest(m_neighbors.size());
  array<MPI_Request> mpiRecvRequest(m_neighbors.size());
  array<MPI_Status>  mpiSendStatus(m_neighbors.size());
  array<MPI_Status>  mpiRecvStatus(m_neighbors.size());


  array<typename array<T>::size_type> recvSize(m_neighbors.size());

  for( localIndex i=0 ; i<m_neighbors.size() ; ++i )
  {
    NeighborCommunication& neighbor = m_neighbors[i];
    typename array<T>::size_type sendSize = sendArray[i].size();

    neighbor.SendReceive( &sendSize, 1, mpiSendRequest[i],
                          &(recvSize[i]), 1, mpiRecvRequest[i] );

  }

  MPI_Waitall( mpiSendRequest.size(), mpiSendRequest.data(), mpiSendStatus.data() );
  MPI_Waitall( mpiRecvRequest.size(), mpiRecvRequest.data(), mpiRecvStatus.data() );

  for( localIndex i=0 ; i<m_neighbors.size() ; ++i )
  {
    NeighborCommunication& neighbor = m_neighbors[i];
    recvArray[i].resize( recvSize[i] );

    neighbor.SendReceive( sendArray[i].data(), sendArray[i].size(), mpiSendRequest[i],
                          recvArray[i].data(), recvArray[i].size(), mpiRecvRequest[i] );

  }
  MPI_Waitall( mpiSendRequest.size(), mpiSendRequest.data(), mpiSendStatus.data() );
  MPI_Waitall( mpiRecvRequest.size(), mpiRecvRequest.data(), mpiRecvStatus.data() );

}



/**
 * @brief Set up the lists on each neighbor
 * @param[in,out] domain Physical domain object
 * @param[in] contactActive Flag whether contact is used
 * @param[in] elementGhostingDepth Depth
 */
void PartitionBase::SetUpNeighborLists( DomainPartition * domain,
                                        const bool contactActive )
{
  for( array<NeighborCommunication>::iterator neighbor=m_neighbors.begin() ; neighbor!=m_neighbors.end() ; ++neighbor )
  {
    neighbor->Clear();
  }

  //------------------------------------------
  //(2) Collect FEM objects that are owned by other processes. This is know from
  // the global indices.

//
//
//  FindMatchedBoundaryIndices( keys::nodeManager,
// domain->GetGroup<ObjectManagerBase>(keys::nodeManager) );
//  FindMatchedBoundaryIndices( DomainPartition::FiniteElementFaceManager,
// domain->m_feFaceManager );
//
//  for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin() ;
// neighbor!=m_neighbors.end() ; ++neighbor )
//  {
//    localIndex_array& matchedBoundaryFaceIndices = stlMapLookup(
// neighbor->tempNeighborData.matchedIndices,
// DomainPartition::FiniteElementFaceManager, "Cannot find
// FiniteElementFaceManager in matchedIndices in
// PartitionBase::SetupNeighborLists" );
//    localIndex_array& matchedBoundaryNodeIndices = stlMapLookup(
// neighbor->tempNeighborData.matchedIndices,
// DomainPartition::FiniteElementNodeManager, "Cannot find
// FiniteElementNodeManager in matchedIndices in
// PartitionBase::SetupNeighborLists" );
//
//    // create an array to contain all faces that have a matched face on a
// neighbor. That means that the face CANNOT
//    // be an external face.
//    domain->m_feFaceManager.m_matchedBoundaryFaces.insert(matchedBoundaryFaceIndices.begin(),
//                                                       matchedBoundaryFaceIndices.end()
// );
//
//    domain->m_feNodeManager.m_matchedBoundaryNodes.insert(matchedBoundaryNodeIndices.begin(),
//                                                       matchedBoundaryNodeIndices.end()
// );
//
//  }
//
//  domain->m_feFaceManager.SetIsExternal();
//  domain->m_feNodeManager.SetIsExternal( &(domain->m_feFaceManager) );
//
//  domain->m_discreteElementSurfaceFaces.m_isExternal = 1;
//  domain->m_discreteElementSurfaceNodes.m_isExternal = 1;
//
//#ifdef SRC_EXTERNAL
//  domain->m_faultPatchFaces.SetIsExternal();
//  domain->m_faultPatchNodes.SetIsExternal( &(domain->m_faultPatchFaces) );
//#endif
//  //------------------------------------------
//  //(5) SENDRECV "GHOST" SIZES
//  //    send the size of the ghost node objects
//  {
//    for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin()
// ; neighbor!=m_neighbors.end() ; ++neighbor )
//    {
//      // pack the element to send as ghosts, and find a list of nodes, edges,
// and faces that must be sent as ghosts.
//      neighbor->FindGhosts( contactActive, m_ghostDepth );
//    }
//
//    CommunicateRequiredObjectIndices();
//
//    for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin()
// ; neighbor!=m_neighbors.end() ; ++neighbor )
//    {
//      neighbor->FindPackGhosts_Step2();
//    }
//
//  }
//
//  for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin() ;
// neighbor!=m_neighbors.end() ; ++neighbor )
//  {
//    // pack the element to send as ghosts, and find a list of nodes, edges,
// and faces that must be sent as ghosts.
//    neighbor->FindPackGhostsDiscreteElement( contactActive, m_ghostDepth );
//#ifdef SRC_EXTERNAL
//    neighbor->FindPackGhostsFaultElement( m_ghostDepth );
//#endif
//  }
//
//
//
//
//
//
//
//
//  {
//    int i_mpireq = 0;
//    for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin()
// ; neighbor!=m_neighbors.end() ; ++neighbor )
//    {
//      //resize the size arrays
//      neighbor->tempNeighborData.sendSizes.resize(neighbor->tempNeighborData.objectsToSend.size(),
// 0);
//      neighbor->tempNeighborData.receiveSizes.resize(neighbor->tempNeighborData.objectsToSend.size(),
// 0);
//
//      //"objectsToSend" filled in NeighborCommunication::FindPackGhosts
//      localIndex index = 0;
//      for(std::map<DomainPartition::ObjectDataStructureKeys,
// bufvector>::const_iterator it =
// neighbor->tempNeighborData.objectsToSend.begin();
//          it != neighbor->tempNeighborData.objectsToSend.end(); ++it, ++index)
//      {
//        neighbor->tempNeighborData.sendSizes[index] = it->second.size();
//      }
//
//      //sendrecv the "receiveSizes" and "sendSizes": note that the sizes are
// calculated the same on all processes (const number of fields)
//      //  so we don't need to send the sizes
//      neighbor->SendReceive(
// reinterpret_cast<char*>(neighbor->tempNeighborData.receiveSizes.data()),
//                             neighbor->tempNeighborData.receiveSizes.size() *
// sizeof(bufvector::size_type),
//                             reinterpret_cast<char*>(neighbor->tempNeighborData.sendSizes.data()),
//                             neighbor->tempNeighborData.sendSizes.size() *
// sizeof(bufvector::size_type),
//                             m_mpiRequest[i_mpireq], m_mpiRequest[i_mpireq+1]
// );
//      i_mpireq += 2;
//    }
//    MPI_Waitall( m_mpiRequest.size() , m_mpiRequest.data(), m_mpiStatus.data()
// );
//
//    //initialize and resize the "tempNeighborData.objectsToReceive" buffers
//    for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin()
// ; neighbor!=m_neighbors.end() ; ++neighbor )
//    {
//      localIndex index = 0;
//      for(std::map<DomainPartition::ObjectDataStructureKeys,
// bufvector>::const_iterator it =
// neighbor->tempNeighborData.objectsToSend.begin();
//          it != neighbor->tempNeighborData.objectsToSend.end(); ++it, ++index)
//      {
//        neighbor->tempNeighborData.objectsToReceive[it->first].resize(
// neighbor->tempNeighborData.receiveSizes[index]);
//      }
//    }
//  }
//
//  array<DomainPartition::ObjectDataStructureKeys> objectNames;
//  NeighborCommunication::SyncNames(objectNames);
//
//  //------------------------------------------
//  // (6) SENDRECV "GHOST" OBJECTS
//  //     send the the ghost objects
//  {
//    DomainPartition::ObjectDataStructureKeys last =
// DomainPartition::Last_ObjectDataStructureNames_Index;
//    localIndex index = 0;
//
//    std::map<DomainPartition::ObjectDataStructureKeys, localIndex_array>
// newIndices;
//    std::map<std::string,localIndex_array> newElementIndices;
//
//    // (6a) do the communication
//    for(array<DomainPartition::ObjectDataStructureKeys>::const_iterator it =
// objectNames.begin() ; it != objectNames.end(); ++it, ++index)
//    {
//      int i_mpireq = 0;
//      for( VectorT<NeighborCommunication>::iterator
// neighbor=m_neighbors.begin() ; neighbor!=m_neighbors.end() ; ++neighbor )
//      {
//        std::map<DomainPartition::ObjectDataStructureKeys,
// bufvector>::iterator itmp =
// neighbor->tempNeighborData.objectsToSend.find(*it);
//        if(itmp == neighbor->tempNeighborData.objectsToSend.end())
//#ifdef USE_ATK
//          SLIC_ERROR("Cannot find name " + toString<int>(*it) + " in
// objectsToSend");
//#endif
//
//        bufvector& send = itmp->second;
//        bufvector& recv = neighbor->tempNeighborData.objectsToReceive[*it];
//
//        neighbor->SendReceive( reinterpret_cast<char*>(recv.data() ),
// recv.size(),
//                               reinterpret_cast<char*>(send.data() ),
// send.size(),
//                               m_mpiRequest[i_mpireq],
// m_mpiRequest[i_mpireq+1] );
//
//        i_mpireq += 2;
//        if(index > 0)
//        {
//          //per usual: handle elements as a special case
//          if(last == DomainPartition::FiniteElementElementManager)
//            neighbor->UnpackGhostElements( newElementIndices );
//          else
//            neighbor->UnpackGhosts(last, newIndices[last]);
//        }
//      }
//
//      last = *it;
//
//      // wait for ghost node objects
//      MPI_Waitall( m_mpiRequest.size() , m_mpiRequest.data(),
// m_mpiStatus.data() );
//    }
//
//    // (6b) sweep up the last to be received
//    for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin()
// ; neighbor!=m_neighbors.end() ; ++neighbor )
//    {
//      //per usual: handle elements as a special case
//      if(last == DomainPartition::FiniteElementElementManager)
//        neighbor->UnpackGhostElements( newElementIndices );
//      else
//        neighbor->UnpackGhosts(last, newIndices[last]);
//    }
//
//
//
//
//
//
//    // (6c) per usual: handle elements as a special case
//    for( std::map<std::string,localIndex_array>::iterator
// i=newElementIndices.begin() ; i!=newElementIndices.end() ; ++i )
//    {
//      localIndex_array& newElems = i->second;
//      // sort the entries
//      std::sort(newElems.begin(),newElems.end());
//      // now remove the duplicates
//      localIndex_array::iterator iend =
// std::unique(newElems.begin(),newElems.end());
//      newElems.resize( iend - newElems.begin() );
//    }
//
//
//
//
//
//
//
//
//
//
//
//
//
//    // (6d) update maps in the anti-dependency direction
//    //      must update the "upward" pointing maps, as the unpacking only
// updated "downward" pointing maps
//    domain->m_feNodeManager.AddToNodeToElementMap( domain->m_feElementManager,
// newElementIndices );
//    domain->m_feFaceManager.AddToFaceToElementMap( domain->m_feElementManager,
// newElementIndices );
//
//    std::map<DomainPartition::ObjectDataStructureKeys,
// localIndex_array>::iterator iArr;
//
//    //Faces
//    {
//      iArr = newIndices.find(DomainPartition::FiniteElementFaceManager);
//      if( iArr != newIndices.end() )
//      {
//        localIndex_array& newFaceIndices = iArr->second;
//        domain->m_feNodeManager.AddToVariableOneToManyFromInverse(
// "nodeToFaceMap", domain->m_feFaceManager.m_toNodesRelation,  newFaceIndices);
//        domain->m_feEdgeManager.AddToVariableOneToManyFromInverse(
// "edgesToFaces", domain->m_feFaceManager.m_toEdgesRelation, newFaceIndices);
//      }
//    }
//
//    //Edges
//    {
//      iArr = newIndices.find(DomainPartition::FiniteElementEdgeManager);
//      if( iArr != newIndices.end() )
//      {
//        localIndex_array& newEdgeIndices = iArr->second;
//        domain->m_feNodeManager.AddToNodeToEdgeMap( domain->m_feEdgeManager,
// newEdgeIndices );
//      }
//    }
//
//
//    //Discrete Elements
//    {
//      iArr = newIndices.find(DomainPartition::DiscreteElementFaceManager);
//      if(iArr != newIndices.end())
//      {
//        localIndex_array& newDEFaceIndices = iArr->second;
//        domain->m_discreteElementSurfaceNodes.AddToNodeToFaceMap(
// domain->m_discreteElementSurfaceFaces, newDEFaceIndices);
//      }
//    }
//  }
//
//  //------------------------------------------
//  // (7) CLEAN UP
//  for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin() ;
// neighbor!=m_neighbors.end() ; ++neighbor )
//  {
//    neighbor->tempNeighborData.clear();
//  }
}


/**
 * @author settgast
 *
 * function to inform partitions of what they need to send to their neighbors
 */
void PartitionBase::CommunicateRequiredObjectIndices()
{

//
//  const DomainPartition::ObjectDataStructureKeys keys[3] = {
// DomainPartition::FiniteElementNodeManager,
//                                                             DomainPartition::FiniteElementEdgeManager,
//                                                             DomainPartition::FiniteElementFaceManager
// };
//
//
//
//
//  //***** Send the list of global indices that are going to be sent to
// neighbor (stored in NeighborData.objectGlobalIndicesToSend),
//  //      and receive in tempNeighborData.objectGlobalIndicesToRecieve
//  int i_mpireq = 0;
//  for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin() ;
// neighbor!=m_neighbors.end() ; ++neighbor )
//  {
//    //resize the size arrays
//    neighbor->tempNeighborData.sendSizes.resize(DomainPartition::numObjectDataStructureNames,
// 0);
//    neighbor->tempNeighborData.receiveSizes.resize(DomainPartition::numObjectDataStructureNames,
// 0);
//
//    // send/recv sizes for the global indices that are required to exist on
// the neighbor
//    for( int i=0 ; i<3 ; ++i )
//    {
//      neighbor->tempNeighborData.sendSizes[keys[i]] =
// neighbor->tempNeighborData.objectGlobalIndicesToSend[keys[i]].size();
//    }
//    neighbor->SendReceive( neighbor->tempNeighborData.sendSizes.data(),
//                           neighbor->tempNeighborData.sendSizes.size(),
//                           m_mpiRequest[i_mpireq],
//                           neighbor->tempNeighborData.receiveSizes.data(),
//                           neighbor->tempNeighborData.receiveSizes.size(),
//                           m_mpiRequest[i_mpireq+1] );
//    i_mpireq += 2;
//  }
//  MPI_Waitall( m_mpiRequest.size() , m_mpiRequest.data(), m_mpiStatus.data()
// );
//
//
//  // now allocate space to receive the global indices that are required.
//  for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin() ;
// neighbor!=m_neighbors.end() ; ++neighbor )
//  {
//    for( int i=0 ; i<3 ; ++i )
//    {
//      neighbor->tempNeighborData.objectGlobalIndicesToRecieve[keys[i]].resize(neighbor->tempNeighborData.receiveSizes[keys[i]]);
//    }
//  }
//
//
//  // send/receive the indices that are required.
//  for( int i=0 ; i<3 ; ++i )
//  {
//
//    i_mpireq = 0;
//    for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin()
// ; neighbor!=m_neighbors.end() ; ++neighbor )
//    {
//      neighbor->SendReceive(
// neighbor->tempNeighborData.objectGlobalIndicesToSend[keys[i]].data(),
//                             neighbor->tempNeighborData.objectGlobalIndicesToSend[keys[i]].size(),
//                             m_mpiRequest[i_mpireq],
//                             neighbor->tempNeighborData.objectGlobalIndicesToRecieve[keys[i]].data(),
//                             neighbor->tempNeighborData.objectGlobalIndicesToRecieve[keys[i]].size(),
//                             m_mpiRequest[i_mpireq+1] );
//      i_mpireq += 2;
//    }
//
//    MPI_Waitall( m_mpiRequest.size() , m_mpiRequest.data(), m_mpiStatus.data()
// );
//
//    // add all required indices into the "requiredObjects" set.
//    gSet requiredObjects;
//    for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin()
// ; neighbor!=m_neighbors.end() ; ++neighbor )
//    {
//      requiredObjects.insert(
// neighbor->tempNeighborData.objectGlobalIndicesToRecieve[keys[i]].begin(),
//                                 neighbor->tempNeighborData.objectGlobalIndicesToRecieve[keys[i]].end()
// );
//    }
//
//    std::map<globalIndex,int> requiredObjectFlag;
//    for( gSet::const_iterator gi=requiredObjects.begin() ;
// gi!=requiredObjects.end() ; ++gi )
//    {
//      requiredObjectFlag[*gi] = 0;
//    }
//    // now we have to determine what neighbor owns the objects we need on this
// process...reuse tempNeighborData.objectGlobalIndicesToSend
//    for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin()
// ; neighbor!=m_neighbors.end() ; ++neighbor )
//    {
//      neighbor->tempNeighborData.objectGlobalIndicesToSend[keys[i]].clear();
//      for( gSet::const_iterator gi=requiredObjects.begin() ;
// gi!=requiredObjects.end() ; ++gi )
//      {
//        if( GlobalIndexManager::OwningRank(*gi) == neighbor->NeighborRank() )
//        {
//          requiredObjectFlag[*gi] = 1;
//          neighbor->tempNeighborData.objectGlobalIndicesToSend[keys[i]].push_back(*gi);
//        }
//      }
//    }
//
//    {
//      int throwFlag = 0;
//      gSet failList;
//      for( std::map<globalIndex,int>::const_iterator
// iter=requiredObjectFlag.begin() ; iter!=requiredObjectFlag.end() ; ++iter )
//      {
//        if( iter->second == 0 )
//        {
//          throwFlag = 1;
//          failList.insert( iter->first );
//        }
//      }
//
//      std::stringstream st;
//      if( throwFlag==1 )
//      {
//        std::string objectType;
//        switch( i )
//        {
//          case(0):
//            objectType = "node";
//            break;
//          case(1):
//            objectType = "edge";
//            break;
//          case(2):
//            objectType = "face";
//            break;
//          default:
//            break;
//        }
//
//        st<< "Error in PartitionBase::CommunicateRequiredObjectIndices()\n
// required "<<objectType<<"/s ";
//        st<<"on rank "<< this->m_rank << " is not owned by a neighbor.\n";
//        for( gSet::const_iterator failed=failList.begin() ;
// failed!=failList.end() ; ++failed )
//        {
//          st<< *failed <<"\n";
//        }
////#ifdef USE_ATK
////        SLIC_ERROR(st.str().c_str());
////#endif
//
//      }
//
//      GPException::MPI_Throw( throwFlag, st.str().c_str() );
//    }
//  }
//
//
//
//
//
//
//  //***** send the global indices that are requred from a neighbor.
// ****************************************************
//
//  // send the sizes
//  i_mpireq = 0;
//  for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin() ;
// neighbor!=m_neighbors.end() ; ++neighbor )
//  {
//    for( int i=0 ; i<3 ; ++i )
//    {
//      neighbor->tempNeighborData.sendSizes[keys[i]] =
// neighbor->tempNeighborData.objectGlobalIndicesToSend[keys[i]].size();
//    }
//    neighbor->SendReceive( neighbor->tempNeighborData.sendSizes.data(),
//                           neighbor->tempNeighborData.sendSizes.size(),
//                           m_mpiRequest[i_mpireq],
//                           neighbor->tempNeighborData.receiveSizes.data(),
//                           neighbor->tempNeighborData.receiveSizes.size(),
//                           m_mpiRequest[i_mpireq+1] );
//    i_mpireq += 2;
//  }
//  MPI_Waitall( m_mpiRequest.size() , m_mpiRequest.data(), m_mpiStatus.data()
// );
//
//  for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin() ;
// neighbor!=m_neighbors.end() ; ++neighbor )
//  {
//    for( int i=0 ; i<3 ; ++i )
//    {
//      neighbor->tempNeighborData.objectGlobalIndicesToRecieve[keys[i]].resize(neighbor->tempNeighborData.receiveSizes[keys[i]]);
//    }
//  }
//
//  // send the lists
//  for( int i=0 ; i<3 ; ++i )
//  {
//    i_mpireq = 0;
//    for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin()
// ; neighbor!=m_neighbors.end() ; ++neighbor )
//    {
//      neighbor->SendReceive(
// neighbor->tempNeighborData.objectGlobalIndicesToSend[keys[i]].data(),
//                             neighbor->tempNeighborData.objectGlobalIndicesToSend[keys[i]].size(),
//                             m_mpiRequest[i_mpireq],
//                             neighbor->tempNeighborData.objectGlobalIndicesToRecieve[keys[i]].data(),
//                             neighbor->tempNeighborData.objectGlobalIndicesToRecieve[keys[i]].size(),
//                             m_mpiRequest[i_mpireq+1] );
//      i_mpireq += 2;
//    }
//
//    MPI_Waitall( m_mpiRequest.size() , m_mpiRequest.data(), m_mpiStatus.data()
// );
//
//  }
//
//  // now we amend the send sets
//  for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin() ;
// neighbor!=m_neighbors.end() ; ++neighbor )
//  {
//    for( int i=0 ; i<3 ; ++i )
//    {
//      lSet& indicesToSend =
// neighbor->tempNeighborData.objectLocalIndicesToSend[keys[i]];
//      const gArray1d& indicesToAdd =
// neighbor->tempNeighborData.objectGlobalIndicesToRecieve[keys[i]];
//
//      const ObjectDataStructureBaseT& object =
// m_domain->GetObjectDataStructure( keys[i] );
//
//      for( gArray1d::const_iterator a=indicesToAdd.begin() ;
// a!=indicesToAdd.end() ; ++a )
//      {
//        if( GlobalIndexManager::OwningRank(*a) == m_rank )
//        {
//          const localIndex li = stlMapLookup( object.m_globalToLocalMap, *a );
//          indicesToSend.insert( li );
//        }
//      }
//    }
//  }


}


/*
 * Finds those boundary indices that are the same on two neighboring processors.
 *
 *
 */
void PartitionBase::FindMatchedBoundaryIndices( string const & key,
                                                const ObjectManagerBase& object )
{
//
//  gArray1d boundaryObjectGlobalIndices;
//  object.ConstructListOfBoundaryObjects( boundaryObjectGlobalIndices );
//
//
//  // send the size of the boundaryObjectGlobalIndices to neighbors
//  {
//    VectorT<MPI_Request> mpiSendRequest( m_neighbors.size() );
//    VectorT<MPI_Request> mpiRecvRequest( m_neighbors.size() );
//    VectorT<MPI_Status>  mpiSendStatus( m_neighbors.size() );
//    VectorT<MPI_Status>  mpiRecvStatus( m_neighbors.size() );
//
//    int i_mpi = 0;
//    for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin()
// ; neighbor!=m_neighbors.end() ; ++neighbor, ++i_mpi )
//    {
//      neighbor->tempNeighborData.sendSizes.resize(DomainPartition::Last_ObjectDataStructureNames_Index+1);
//      neighbor->tempNeighborData.receiveSizes.resize(DomainPartition::Last_ObjectDataStructureNames_Index+1);
//
//      neighbor->tempNeighborData.sendSizes[key] =
// boundaryObjectGlobalIndices.size() ;
//
//      neighbor->SendReceive(
// reinterpret_cast<char*>(&(neighbor->tempNeighborData.receiveSizes[key])),
//                             sizeof(bufvector::size_type),
//                             reinterpret_cast<char*>(&(neighbor->tempNeighborData.sendSizes[key])),
//                             sizeof(bufvector::size_type),
//                             mpiRecvRequest[i_mpi], mpiSendRequest[i_mpi] );
//    }
//    MPI_Waitall( mpiRecvRequest.size() , mpiRecvRequest.data(),
// mpiRecvStatus.data() );
//    MPI_Waitall( mpiSendRequest.size() , mpiSendRequest.data(),
// mpiSendStatus.data() );
//  }
//
//
//  {
//    VectorT<MPI_Request> mpiSendRequest( m_neighbors.size() );
//    VectorT<MPI_Request> mpiRecvRequest( m_neighbors.size() );
//    VectorT<MPI_Status>  mpiSendStatus( m_neighbors.size() );
//    VectorT<MPI_Status>  mpiRecvStatus( m_neighbors.size() );
//
//    int i_mpi = 0;
//    for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin()
// ; neighbor!=m_neighbors.end() ; ++neighbor, ++i_mpi )
//    {
//      neighbor->tempNeighborData.neighborNumbers[key].resize(neighbor->tempNeighborData.receiveSizes[key]);
//
//      neighbor->SendReceive( reinterpret_cast<char*>(
// (neighbor->tempNeighborData.neighborNumbers[key]).data() ),
//                             (neighbor->tempNeighborData.neighborNumbers[key]).size()
// * sizeof(gArray1d::size_type),
//                             reinterpret_cast<char*>(
// boundaryObjectGlobalIndices.data() ),
//                             boundaryObjectGlobalIndices.size() *
// sizeof(gArray1d::size_type),
//                             mpiRecvRequest[i_mpi], mpiSendRequest[i_mpi]  );
//    }
//
//    for( unsigned int count=0 ; count<m_neighbors.size() ; ++count )
//    {
//      int index;
//      MPI_Waitany( mpiRecvRequest.size(), mpiRecvRequest.data(), &index,
// mpiRecvStatus.data() );
//
//      NeighborCommunication& neighbor = this->m_neighbors[index];
//
//      neighbor.tempNeighborData.matchedIndices[key].clear();
//
//      // compare lists
//      gArray1d::const_iterator iter_local =
// boundaryObjectGlobalIndices.begin();
//      gArray1d::const_iterator iter_neighbor =
// neighbor.tempNeighborData.neighborNumbers[key].begin();
//      while( iter_local!=boundaryObjectGlobalIndices.end() &&
// iter_neighbor!=neighbor.tempNeighborData.neighborNumbers[key].end() )
//      {
//        if( *iter_local==*iter_neighbor )
//        {
//          const localIndex& li = stlMapLookup( object.m_globalToLocalMap,
// *iter_local );
//          neighbor.tempNeighborData.matchedIndices[key].push_back(li);
//          ++iter_local;
//          ++iter_neighbor;
//        }
//        else if( *iter_local > *iter_neighbor )
//        {
//          ++iter_neighbor;
//        }
//        else if( *iter_local < *iter_neighbor )
//        {
//          ++iter_local;
//        }
//      }
//
//
//    }
//
//    MPI_Waitall( mpiSendRequest.size() , mpiSendRequest.data(),
// mpiSendStatus.data() );
//
//  }
}



#define METHOD 1


#if METHOD == 3
void PartitionBase::ModifyGhostsAndNeighborLists( const ModifiedObjectLists& modifiedObjects )
{

  // count the number of new objects that are owned by a neighbor, and send that
  // number to the neighbor


  array<MPI_Request> mpiSendSizeRequest( m_neighbors.size() );
  array<MPI_Request> mpiRecvSizeRequest( m_neighbors.size() );
  array<MPI_Status>  mpiSendSizeStatus( m_neighbors.size() );
  array<MPI_Status>  mpiRecvSizeStatus( m_neighbors.size() );

  array<MPI_Request> mpiSendBufferRequest( m_neighbors.size() );
  array<MPI_Request> mpiRecvBufferRequest( m_neighbors.size() );
  array<MPI_Status>  mpiSendBufferStatus( m_neighbors.size() );
  array<MPI_Status>  mpiRecvBufferStatus( m_neighbors.size() );



  //******************************
  // 1) for a given processor, send new objects created on this processor and
  // the connectivities to neighbors
  //    that either are a ghost on this processor(owned by neighbor), and the
  // ones that are a ghost on the neighbor.

  // pack the buffers, and send the size of the buffers
  for( unsigned int neighborIndex=0 ; neighborIndex<m_neighbors.size() ; ++neighborIndex )
  {
    NeighborCommunication& neighbor = m_neighbors[neighborIndex];

    neighbor.ResizeSendBuffer(0);

    // pack the new/modified nodes that are owned by this process
    neighbor.PackTopologyModifications( DomainPartition::FiniteElementNodeManager, modifiedObjects.newNodes,  modifiedObjects.modifiedNodes, true, false );

    // pack the new/modified edges that are owned by the neighbor.
    neighbor.PackTopologyModifications( DomainPartition::FiniteElementEdgeManager, modifiedObjects.newEdges,  modifiedObjects.modifiedEdges, true,  true );
    // pack the new/modified edges that are owned by this process
    neighbor.PackTopologyModifications( DomainPartition::FiniteElementEdgeManager, modifiedObjects.newEdges,  modifiedObjects.modifiedEdges, true,  false );

    // pack the new/modified faces that are owned by the neighbor.
    neighbor.PackTopologyModifications( DomainPartition::FiniteElementFaceManager, modifiedObjects.newFaces,  modifiedObjects.modifiedFaces, true,  true );
    // pack the new/modified faces that are owned by this process
    neighbor.PackTopologyModifications( DomainPartition::FiniteElementFaceManager, modifiedObjects.newFaces,  modifiedObjects.modifiedFaces, true,  false );

    // pack the new/modified elements that are owned by the neighbor.
    neighbor.PackTopologyModifications( DomainPartition::FiniteElementElementManager, std::map< std::string,
                                                                                                lSet >(),  modifiedObjects.modifiedElements, true,  true );
    // pack the new/modified elements that are owned by this process
    neighbor.PackTopologyModifications( DomainPartition::FiniteElementElementManager, std::map< std::string,
                                                                                                lSet >(),  modifiedObjects.modifiedElements, true,  false );

    neighbor.SendReceiveBufferSizes(CommRegistry::genericComm01, mpiSendSizeRequest[neighborIndex], mpiRecvSizeRequest[neighborIndex] );
  }

  // send/recv the buffers
  for( unsigned int count=0 ; count<m_neighbors.size() ; ++count )
  {
    int neighborIndex;
    MPI_Waitany( mpiRecvSizeRequest.size(), mpiRecvSizeRequest.data(), &neighborIndex, mpiRecvSizeStatus.data() );

    NeighborCommunication& neighbor = m_neighbors[neighborIndex];

    neighbor.SendReceiveBuffers( CommRegistry::genericComm01, CommRegistry::genericComm02, mpiSendBufferRequest[neighborIndex],
                                 mpiRecvBufferRequest[neighborIndex] );
  }


  lSet allNewAndModifiedLocalNodes;

  lSet allNewNodes, allModifiedNodes;
  lSet allNewEdges, allModifiedEdges;
  lSet allNewFaces, allModifiedFaces;
  std::map< std::string, lSet> allModifiedElements;
  // unpack the buffers
  for( unsigned int count=0 ; count<m_neighbors.size() ; ++count )
  {
    localIndex_array newLocalNodes, modifiedLocalNodes;
    localIndex_array newGhostNodes, modifiedGhostNodes;

    localIndex_array newLocalEdges, modifiedLocalEdges;
    localIndex_array newGhostEdges, modifiedGhostEdges;

    localIndex_array newLocalFaces, modifiedLocalFaces;
    localIndex_array newGhostFaces, modifiedGhostFaces;

    std::map< std::string, localIndex_array> modifiedElements;



    int neighborIndex;
    MPI_Waitany( mpiRecvBufferRequest.size(), mpiRecvBufferRequest.data(), &neighborIndex, mpiRecvBufferStatus.data() );

    NeighborCommunication& neighbor = this->m_neighbors[neighborIndex];

    const char* pbuffer = neighbor.ReceiveBuffer().data();

    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementNodeManager, pbuffer, newGhostNodes, modifiedGhostNodes, false );
    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementEdgeManager, pbuffer, newLocalEdges, modifiedLocalEdges, true );
    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementEdgeManager, pbuffer, newGhostEdges, modifiedGhostEdges, false );
    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementFaceManager, pbuffer, newLocalFaces, modifiedLocalFaces, true );
    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementFaceManager, pbuffer, newGhostFaces, modifiedGhostFaces, false );
    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementElementManager, pbuffer, modifiedElements );
    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementElementManager, pbuffer, modifiedElements );

    allNewNodes.insert( newLocalNodes.begin(), newLocalNodes.end() );
    allNewNodes.insert( newGhostNodes.begin(), newGhostNodes.end() );

    allModifiedNodes.insert( modifiedLocalNodes.begin(), modifiedLocalNodes.end() );
    allModifiedNodes.insert( modifiedGhostNodes.begin(), modifiedGhostNodes.end() );

    allNewAndModifiedLocalNodes.insert( newLocalNodes.begin(), newLocalNodes.end() );
    allNewAndModifiedLocalNodes.insert( modifiedLocalNodes.begin(), modifiedLocalNodes.end() );

    allNewEdges.insert( newLocalEdges.begin(), newLocalEdges.end() );
    allNewEdges.insert( newGhostEdges.begin(), newGhostEdges.end() );

    allModifiedEdges.insert( modifiedLocalEdges.begin(), modifiedLocalEdges.end() );
    allModifiedEdges.insert( modifiedGhostEdges.begin(), modifiedGhostEdges.end() );

    allNewFaces.insert( newLocalFaces.begin(), newLocalFaces.end() );
    allNewFaces.insert( newGhostFaces.begin(), newGhostFaces.end() );

    allModifiedFaces.insert( modifiedLocalFaces.begin(), modifiedLocalFaces.end() );
    allModifiedFaces.insert( modifiedGhostFaces.begin(), modifiedGhostFaces.end() );

    for( std::map< std::string, localIndex_array>::const_iterator i=modifiedElements.begin() ; i!=modifiedElements.end() ; ++i )
    {
      allModifiedElements[i->first].insert( i->second.begin(), i->second.end() );
    }
  }

  MPI_Waitall( mpiSendSizeRequest.size(), mpiSendSizeRequest.data(), mpiSendSizeStatus.data() );
  MPI_Waitall( mpiSendBufferRequest.size(), mpiSendBufferRequest.data(), mpiSendBufferStatus.data() );


  // must send the new objects that are local to this processor, but created on
  // a neighbor, to the other neighbors.
  for( unsigned int neighborIndex=0 ; neighborIndex<m_neighbors.size() ; ++neighborIndex )
  {
    NeighborCommunication& neighbor = m_neighbors[neighborIndex];

    neighbor.ResizeSendBuffer(0);

    neighbor.PackTopologyModifications( DomainPartition::FiniteElementEdgeManager, allNewEdges,  allModifiedEdges, false, false );
    neighbor.PackTopologyModifications( DomainPartition::FiniteElementFaceManager, allNewFaces,  allModifiedFaces, false, false );
    neighbor.PackTopologyModifications( DomainPartition::FiniteElementElementManager, std::map< std::string, lSet >(),  allModifiedElements, false, false );

    neighbor.SendReceiveBufferSizes(CommRegistry::genericComm01, mpiSendSizeRequest[neighborIndex], mpiRecvSizeRequest[neighborIndex] );
  }

  // send/recv the buffers
  for( unsigned int count=0 ; count<m_neighbors.size() ; ++count )
  {
    int neighborIndex;
    MPI_Waitany( mpiRecvSizeRequest.size(), mpiRecvSizeRequest.data(), &neighborIndex, mpiRecvSizeStatus.data() );

    NeighborCommunication& neighbor = m_neighbors[neighborIndex];

    neighbor.SendReceiveBuffers( CommRegistry::genericComm01, CommRegistry::genericComm02, mpiSendBufferRequest[neighborIndex],
                                 mpiRecvBufferRequest[neighborIndex] );
  }
  // unpack the buffers
  for( unsigned int count=0 ; count<m_neighbors.size() ; ++count )
  {
    localIndex_array newLocalNodes, modifiedLocalNodes;
    localIndex_array newGhostNodes, modifiedGhostNodes;

    localIndex_array newLocalEdges, modifiedLocalEdges;
    localIndex_array newGhostEdges, modifiedGhostEdges;

    localIndex_array newLocalFaces, modifiedLocalFaces;
    localIndex_array newGhostFaces, modifiedGhostFaces;

    std::map< std::string, localIndex_array> modifiedElements;



    int neighborIndex;
    MPI_Waitany( mpiRecvBufferRequest.size(), mpiRecvBufferRequest.data(), &neighborIndex, mpiRecvBufferStatus.data() );

    NeighborCommunication& neighbor = this->m_neighbors[neighborIndex];

    const char* pbuffer = neighbor.ReceiveBuffer().data();

    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementEdgeManager, pbuffer, newGhostEdges, modifiedGhostEdges, false );
    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementFaceManager, pbuffer, newGhostFaces, modifiedGhostFaces, false );
    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementElementManager, pbuffer, modifiedElements );

    allNewNodes.insert( newLocalNodes.begin(), newLocalNodes.end() );
    allNewNodes.insert( newGhostNodes.begin(), newGhostNodes.end() );

    allModifiedNodes.insert( modifiedLocalNodes.begin(), modifiedLocalNodes.end() );
    allModifiedNodes.insert( modifiedGhostNodes.begin(), modifiedGhostNodes.end() );

    allNewAndModifiedLocalNodes.insert( newLocalNodes.begin(), newLocalNodes.end() );
    allNewAndModifiedLocalNodes.insert( modifiedLocalNodes.begin(), modifiedLocalNodes.end() );

    allNewEdges.insert( newLocalEdges.begin(), newLocalEdges.end() );
    allNewEdges.insert( newGhostEdges.begin(), newGhostEdges.end() );

    allModifiedEdges.insert( modifiedLocalEdges.begin(), modifiedLocalEdges.end() );
    allModifiedEdges.insert( modifiedGhostEdges.begin(), modifiedGhostEdges.end() );

    allNewFaces.insert( newLocalFaces.begin(), newLocalFaces.end() );
    allNewFaces.insert( newGhostFaces.begin(), newGhostFaces.end() );

    allModifiedFaces.insert( modifiedLocalFaces.begin(), modifiedLocalFaces.end() );
    allModifiedFaces.insert( modifiedGhostFaces.begin(), modifiedGhostFaces.end() );

    for( std::map< std::string, localIndex_array>::const_iterator i=modifiedElements.begin() ; i!=modifiedElements.end() ; ++i )
    {
      allModifiedElements[i->first].insert( i->second.begin(), i->second.end() );
    }
  }

  MPI_Waitall( mpiSendSizeRequest.size(), mpiSendSizeRequest.data(), mpiSendSizeStatus.data() );


  lSet allReceivedNodes;
  allReceivedNodes.insert( allNewNodes.begin(), allNewNodes.end() );
  allReceivedNodes.insert( allModifiedNodes.begin(), allModifiedNodes.end() );

  allReceivedNodes.erase( allNewAndModifiedLocalNodes.begin(), allNewAndModifiedLocalNodes.end() );

  m_domain->m_feNodeManager.ConnectivityFromGlobalToLocal( allNewAndModifiedLocalNodes,
                                                           allReceivedNodes,
                                                           m_domain->m_feFaceManager.m_globalToLocalMap );


  lSet allReceivedEdges;
  allReceivedEdges.insert( allNewEdges.begin(), allNewEdges.end() );
  allReceivedEdges.insert( allModifiedEdges.begin(), allModifiedEdges.end() );

  m_domain->m_feEdgeManager.ConnectivityFromGlobalToLocal( allReceivedEdges, m_domain->m_feNodeManager.m_globalToLocalMap );

  lSet allReceivedFaces;
  allReceivedFaces.insert( allNewFaces.begin(), allNewFaces.end() );
  allReceivedFaces.insert( allModifiedFaces.begin(), allModifiedFaces.end() );

  m_domain->m_feFaceManager.ConnectivityFromGlobalToLocal( allReceivedFaces,
                                                           m_domain->m_feNodeManager.m_globalToLocalMap,
                                                           m_domain->m_feEdgeManager.m_globalToLocalMap );


  m_domain->m_feElementManager.ConnectivityFromGlobalToLocal( allModifiedElements,
                                                              m_domain->m_feNodeManager.m_globalToLocalMap,
                                                              m_domain->m_feFaceManager.m_globalToLocalMap );


  m_domain->m_feNodeManager.ModifyNodeToEdgeMapFromSplit( m_domain->m_feEdgeManager,
                                                          allNewEdges,
                                                          allModifiedEdges );

  m_domain->m_feFaceManager.ModifyToFaceMapsFromSplit( allNewFaces,
                                                       allModifiedFaces,
                                                       m_domain->m_feNodeManager,
                                                       m_domain->m_feEdgeManager );

  m_domain->m_feElementManager.ModifyToElementMapsFromSplit( allModifiedElements,
                                                             m_domain->m_feNodeManager,
                                                             m_domain->m_feFaceManager );



}


#elif METHOD == 2
void PartitionBase::ModifyGhostsAndNeighborLists( const ModifiedObjectLists& modifiedObjects )
{



  lSet newLocalNodes( modifiedObjects.newNodes );
  lSet modifiedLocalNodes( modifiedObjects.modifiedNodes );

  lSet newLocalEdges( modifiedObjects.newEdges );
  lSet modifiedLocalEdges( modifiedObjects.modifiedEdges );

  lSet newLocalFaces( modifiedObjects.newFaces );
  lSet modifiedLocalFaces( modifiedObjects.modifiedFaces );

  std::map< std::string, lSet> modifiedLocalElements(modifiedObjects.modifiedElements);



  // buffers and MPI objects for sending data about new/modified objects back to
  // the partition that owns them.
  array<bufvector> send_buffer0( m_neighbors.size() );
  array<bufvector::size_type> sendSize0( m_neighbors.size() );

  array<bufvector> recv_buffer0( m_neighbors.size() );
  array<bufvector::size_type> recvSize0( m_neighbors.size() );

  array<MPI_Request> mpiSendSizeRequest0( m_neighbors.size() );
  array<MPI_Request> mpiRecvSizeRequest0( m_neighbors.size() );
  array<MPI_Status>  mpiSendSizeStatus0( m_neighbors.size() );
  array<MPI_Status>  mpiRecvSizeStatus0( m_neighbors.size() );


  array<MPI_Request> mpiSendBufferRequest0( m_neighbors.size() );
  array<MPI_Request> mpiRecvBufferRequest0( m_neighbors.size() );
  array<MPI_Status>  mpiSendBufferStatus0( m_neighbors.size() );
  array<MPI_Status>  mpiRecvBufferStatus0( m_neighbors.size() );

  for( unsigned int neighborIndex=0 ; neighborIndex<m_neighbors.size() ; ++neighborIndex )
  {
    NeighborCommunication& neighbor = m_neighbors[neighborIndex];

//    neighbor.PackNewAndModifiedGhostObjects(
// DomainPartition::FiniteElementNodeManager, newLocalNodes,
//  modifiedLocalNodes, send_buffer0[count] );
    neighbor.PackNewAndModifiedGhostObjects( DomainPartition::FiniteElementEdgeManager, newLocalEdges,  newLocalEdges, send_buffer0[neighborIndex] );
    neighbor.PackNewAndModifiedGhostObjects( DomainPartition::FiniteElementFaceManager, newLocalFaces,  modifiedLocalFaces, send_buffer0[neighborIndex] );
    neighbor.PackNewAndModifiedGhostObjects( DomainPartition::FiniteElementElementManager, modifiedLocalElements, send_buffer0[neighborIndex] );

    neighbor.PackNewAndModifiedLocalObjectsFromThisPartition( DomainPartition::FiniteElementNodeManager, newLocalNodes,  modifiedLocalNodes,
                                                              send_buffer0[neighborIndex] );
    neighbor.PackNewAndModifiedLocalObjectsFromThisPartition( DomainPartition::FiniteElementEdgeManager, newLocalEdges,  newLocalEdges,
                                                              send_buffer0[neighborIndex] );
    neighbor.PackNewAndModifiedLocalObjectsFromThisPartition( DomainPartition::FiniteElementFaceManager, newLocalFaces,  modifiedLocalFaces,
                                                              send_buffer0[neighborIndex] );
    neighbor.PackNewAndModifiedLocalObjectsFromThisPartition( DomainPartition::FiniteElementElementManager, modifiedLocalElements,
                                                              send_buffer0[neighborIndex] );



    sendSize0[neighborIndex] = send_buffer0[neighborIndex].size();
    neighbor.SendReceive( &(sendSize0[neighborIndex]), 1, mpiSendSizeRequest0[neighborIndex],
                          &(recvSize0[neighborIndex]), 1, mpiRecvSizeRequest0[neighborIndex] );

  }

  MPI_Waitall( mpiRecvSizeRequest0.size(), mpiRecvSizeRequest0.data(), mpiRecvSizeStatus0.data() );
  MPI_Waitall( mpiSendSizeRequest0.size(), mpiSendSizeRequest0.data(), mpiSendSizeStatus0.data() );

  // send/recv the buffers
  for( unsigned int neighborIndex=0 ; neighborIndex<m_neighbors.size() ; ++neighborIndex )
  {
    NeighborCommunication& neighbor = m_neighbors[neighborIndex];

    recv_buffer0[neighborIndex].resize(recvSize0[neighborIndex]);
    neighbor.SendReceive( send_buffer0[neighborIndex], mpiSendBufferRequest0[neighborIndex],
                          recv_buffer0[neighborIndex], mpiRecvBufferRequest0[neighborIndex] );
  }



  lSet newLocalNodesFromNeighbor;
  lSet modifiedLocalNodesFromNeighbor;

  lSet newLocalEdgesFromNeighbor;
  lSet modifiedLocalEdgesFromNeighbor;

  lSet newLocalFacesFromNeighbor;
  lSet modifiedLocalFacesFromNeighbor;

  std::map< std::string, lSet> modifiedLocalElementsFromNeighbor;


  lSet newGhostNodesFromNeighbor;
  lSet modifiedGhostNodesFromNeighbor;

  lSet newGhostEdgesFromNeighbor;
  lSet modifiedGhostEdgesFromNeighbor;

  lSet newGhostFacesFromNeighbor;
  lSet modifiedGhostFacesFromNeighbor;

  std::map< std::string, lSet> modifiedGhostElementsFromNeighbor;


  MPI_Waitall( mpiRecvBufferRequest0.size(), mpiRecvBufferRequest0.data(), mpiRecvBufferStatus0.data() );
  MPI_Waitall( mpiSendBufferRequest0.size(), mpiSendBufferRequest0.data(), mpiSendBufferStatus0.data() );


  // unpack the buffers
  for( unsigned int neighborIndex=0 ; neighborIndex<m_neighbors.size() ; ++neighborIndex )
  {
    localIndex_array newNodes, modifiedNodes;
    localIndex_array newEdges, modifiedEdges;
    localIndex_array newFaces, modifiedFaces;
    std::map< std::string, localIndex_array> modifiedElements;

//    int neighborIndex = -1;
//    MPI_Waitany( mpiRecvBufferRequest0.size(), mpiRecvBufferRequest0.data(),
// &neighborIndex, mpiRecvBufferStatus0.data() );

    NeighborCommunication& neighbor = this->m_neighbors[neighborIndex];

    const char* pbuffer = recv_buffer0[neighborIndex].data();

//    neighbor.UnpackNewAndModifiedLocalObjectsFromNeighborGhosts(
// DomainPartition::FiniteElementNodeManager, pbuffer, newNodes, modifiedNodes
// );
    neighbor.UnpackNewAndModifiedLocalObjectsFromNeighborGhosts( DomainPartition::FiniteElementEdgeManager, pbuffer, newEdges, modifiedEdges );
    neighbor.UnpackNewAndModifiedLocalObjectsFromNeighborGhosts( DomainPartition::FiniteElementFaceManager, pbuffer, newFaces, modifiedFaces );
    neighbor.UnpackNewAndModifiedLocalObjectsFromNeighborGhosts( DomainPartition::FiniteElementElementManager, pbuffer, modifiedElements );

    newLocalNodesFromNeighbor.insert( newNodes.begin(), newNodes.end() );
    modifiedLocalNodesFromNeighbor.insert( modifiedNodes.begin(), modifiedNodes.end() );

    newLocalEdgesFromNeighbor.insert( newEdges.begin(), newEdges.end() );
    modifiedLocalEdgesFromNeighbor.insert( modifiedEdges.begin(), modifiedEdges.end() );

    newLocalFacesFromNeighbor.insert( newFaces.begin(), newFaces.end() );
    modifiedLocalFacesFromNeighbor.insert( modifiedFaces.begin(), modifiedFaces.end() );

    for( std::map< std::string, localIndex_array>::const_iterator i=modifiedElements.begin() ; i!=modifiedElements.end() ; ++i )
    {
      modifiedLocalElementsFromNeighbor[i->first].insert( i->second.begin(), i->second.end() );
    }



    neighbor.UnpackNewAndModifiedDirectGhostObjects( DomainPartition::FiniteElementNodeManager, pbuffer, newNodes, modifiedNodes );
    neighbor.UnpackNewAndModifiedDirectGhostObjects( DomainPartition::FiniteElementEdgeManager, pbuffer, newEdges, modifiedEdges );
    neighbor.UnpackNewAndModifiedDirectGhostObjects( DomainPartition::FiniteElementFaceManager, pbuffer, newFaces, modifiedFaces );
    neighbor.UnpackNewAndModifiedGhostObjects( DomainPartition::FiniteElementElementManager, pbuffer, modifiedElements );

    newGhostNodesFromNeighbor.insert( newNodes.begin(), newNodes.end() );
    modifiedGhostNodesFromNeighbor.insert( modifiedNodes.begin(), modifiedNodes.end() );

    newGhostEdgesFromNeighbor.insert( newEdges.begin(), newEdges.end() );
    modifiedGhostEdgesFromNeighbor.insert( modifiedEdges.begin(), modifiedEdges.end() );

    newGhostFacesFromNeighbor.insert( newFaces.begin(), newFaces.end() );
    modifiedGhostFacesFromNeighbor.insert( modifiedFaces.begin(), modifiedFaces.end() );

    for( std::map< std::string, localIndex_array>::const_iterator i=modifiedElements.begin() ; i!=modifiedElements.end() ; ++i )
    {
      modifiedGhostElementsFromNeighbor[i->first].insert( i->second.begin(), i->second.end() );
    }

  }



/*



   // buffers and MPI objects for sending data from new/modified local objects
      to the neighbors.
   array<bufvector> send_buffer1( m_neighbors.size() );
   array<bufvector::size_type> sendSize1( m_neighbors.size() );

   array<bufvector> recv_buffer1( m_neighbors.size() );
   array<bufvector::size_type> recvSize1( m_neighbors.size() );

   array<MPI_Request> mpiSendSizeRequest1( m_neighbors.size() );
   array<MPI_Request> mpiRecvSizeRequest1( m_neighbors.size() );
   array<MPI_Status>  mpiSendSizeStatus1( m_neighbors.size() );
   array<MPI_Status>  mpiRecvSizeStatus1( m_neighbors.size() );


   array<MPI_Request> mpiSendBufferRequest1( m_neighbors.size() );
   array<MPI_Request> mpiRecvBufferRequest1( m_neighbors.size() );
   array<MPI_Status>  mpiSendBufferStatus1( m_neighbors.size() );
   array<MPI_Status>  mpiRecvBufferStatus1( m_neighbors.size() );



   for( unsigned int count=0 ; count<m_neighbors.size() ; ++count )
   {
    NeighborCommunication& neighbor = m_neighbors[count];

    neighbor.PackNewAndModifiedLocalObjectsFromNeighbor(
       DomainPartition::FiniteElementNodeManager, newLocalNodesFromNeighbor,
        modifiedLocalNodesFromNeighbor, send_buffer1[count] );
    neighbor.PackNewAndModifiedLocalObjectsFromNeighbor(
       DomainPartition::FiniteElementEdgeManager, newLocalEdgesFromNeighbor,
        modifiedLocalEdgesFromNeighbor, send_buffer1[count] );
    neighbor.PackNewAndModifiedLocalObjectsFromNeighbor(
       DomainPartition::FiniteElementFaceManager, newLocalFacesFromNeighbor,
        modifiedLocalFacesFromNeighbor, send_buffer1[count] );
    neighbor.PackNewAndModifiedLocalObjectsFromNeighbor(
       DomainPartition::FiniteElementElementManager,
        modifiedLocalElementsFromNeighbor, send_buffer1[count] );

    sendSize1[count] = send_buffer1[count].size();
    neighbor.SendReceive( &(sendSize1[count]), 1, mpiSendSizeRequest1[count],
                          &(recvSize1[count]), 1, mpiRecvSizeRequest1[count] );

   }
   // send/recv the buffers
   for( unsigned int count=0 ; count<m_neighbors.size() ; ++count )
   {
    int neighborIndex;
    MPI_Waitany( mpiRecvSizeRequest1.size(), mpiRecvSizeRequest1.data(),
       &neighborIndex, mpiRecvSizeStatus1.data() );

    NeighborCommunication& neighbor = m_neighbors[neighborIndex];

    recv_buffer0[count].resize(recvSize1[count]);
    neighbor.SendReceive( send_buffer1[count], mpiSendBufferRequest1[count],
       recv_buffer1[count], mpiRecvBufferRequest1[count] );
   }

   // unpack the buffers
   for( unsigned int count=0 ; count<m_neighbors.size() ; ++count )
   {
    localIndex_array newNodes, modifiedNodes;
    localIndex_array newEdges, modifiedEdges;
    localIndex_array newFaces, modifiedFaces;
    std::map< std::string, localIndex_array> modifiedElements;

    int neighborIndex;
    MPI_Waitany( mpiRecvBufferRequest0.size(), mpiRecvBufferRequest0.data(),
       &neighborIndex, mpiRecvBufferStatus0.data() );

    NeighborCommunication& neighbor = this->m_neighbors[neighborIndex];

    const char* pbuffer = recv_buffer0[count].data();

    neighbor.UnpackNewAndModifiedIndirectGhostObjects(
       DomainPartition::FiniteElementNodeManager, pbuffer, newNodes,
       modifiedNodes );
    neighbor.UnpackNewAndModifiedIndirectGhostObjects(
       DomainPartition::FiniteElementEdgeManager, pbuffer, newEdges,
       modifiedEdges );
    neighbor.UnpackNewAndModifiedIndirectGhostObjects(
       DomainPartition::FiniteElementFaceManager, pbuffer, newFaces,
       modifiedFaces );
    neighbor.UnpackNewAndModifiedGhostObjects(
       DomainPartition::FiniteElementElementManager, pbuffer, modifiedElements
       );

    newGhostNodesFromNeighbor.insert( newNodes.begin(), newNodes.end() );
    modifiedGhostNodesFromNeighbor.insert( modifiedNodes.begin(),
       modifiedNodes.end() );

    newGhostEdgesFromNeighbor.insert( newEdges.begin(), newEdges.end() );
    modifiedGhostEdgesFromNeighbor.insert( modifiedEdges.begin(),
       modifiedEdges.end() );

    newGhostFacesFromNeighbor.insert( newFaces.begin(), newFaces.end() );
    modifiedGhostFacesFromNeighbor.insert( modifiedFaces.begin(),
       modifiedFaces.end() );

    for( std::map< std::string, localIndex_array>::const_iterator
       i=modifiedElements.begin() ; i!=modifiedElements.end() ; ++i )
    {
      modifiedGhostElementsFromNeighbor[i->first].insert( i->second.begin(),
         i->second.end() );
    }

   }
 */

  lSet allReceivedEdges;
  allReceivedEdges.insert( newLocalEdgesFromNeighbor.begin(), newLocalEdgesFromNeighbor.end() );
  allReceivedEdges.insert( newGhostEdgesFromNeighbor.begin(), newGhostEdgesFromNeighbor.end() );
  allReceivedEdges.insert( modifiedLocalEdgesFromNeighbor.begin(), modifiedLocalEdgesFromNeighbor.end() );
  allReceivedEdges.insert( modifiedGhostEdgesFromNeighbor.begin(), modifiedGhostEdgesFromNeighbor.end() );

  m_domain->m_feEdgeManager.ConnectivityFromGlobalToLocal( allReceivedEdges, m_domain->m_feNodeManager.m_globalToLocalMap );

  lSet allReceivedFaces;
  allReceivedFaces.insert( newLocalFacesFromNeighbor.begin(), newLocalFacesFromNeighbor.end() );
  allReceivedFaces.insert( newGhostFacesFromNeighbor.begin(), newGhostFacesFromNeighbor.end() );
  allReceivedFaces.insert( modifiedLocalFacesFromNeighbor.begin(), modifiedLocalFacesFromNeighbor.end() );
  allReceivedFaces.insert( modifiedGhostFacesFromNeighbor.begin(), modifiedGhostFacesFromNeighbor.end() );

  m_domain->m_feFaceManager.ConnectivityFromGlobalToLocal( allReceivedFaces,
                                                           m_domain->m_feNodeManager.m_globalToLocalMap,
                                                           m_domain->m_feEdgeManager.m_globalToLocalMap );


  std::map< std::string, lSet> allReceivedElements;
  allReceivedElements.insert( modifiedLocalElementsFromNeighbor.begin(), modifiedLocalElementsFromNeighbor.end() );
  allReceivedElements.insert( modifiedGhostElementsFromNeighbor.begin(), modifiedGhostElementsFromNeighbor.end() );


  m_domain->m_feElementManager.ConnectivityFromGlobalToLocal( allReceivedElements,
                                                              m_domain->m_feNodeManager.m_globalToLocalMap,
                                                              m_domain->m_feFaceManager.m_globalToLocalMap );


  m_domain->m_feNodeManager.ModifyNodeToEdgeMapFromSplit( m_domain->m_feEdgeManager,
                                                          allReceivedEdges,
                                                          lSet() );

  m_domain->m_feFaceManager.ModifyToFaceMapsFromSplit( allReceivedFaces,
                                                       lSet(),
                                                       m_domain->m_feNodeManager,
                                                       m_domain->m_feEdgeManager );

  m_domain->m_feElementManager.ModifyToElementMapsFromSplit( allReceivedElements,
                                                             m_domain->m_feNodeManager,
                                                             m_domain->m_feFaceManager );



}

#elif METHOD == 1
void PartitionBase::ModifyGhostsAndNeighborLists( const ModifiedObjectLists& modifiedObjects )
{
//
//  const bool pack=true;
//  const realT t0=MPI_Wtime();
//
//
//  // count the number of new objects that are owned by a neighbor, and send
// that number to the neighbor
//
//  array<MPI_Request> mpiSendGlobalRequest( m_neighbors.size() );
//  array<MPI_Request> mpiRecvGlobalRequest( m_neighbors.size() );
//  array<MPI_Status>  mpiSendGlobalStatus( m_neighbors.size() );
//  array<MPI_Status>  mpiRecvGlobalStatus( m_neighbors.size() );
//
//  array<MPI_Request> mpiSendSizeRequest( m_neighbors.size() );
//  array<MPI_Request> mpiRecvSizeRequest( m_neighbors.size() );
//  array<MPI_Status>  mpiSendSizeStatus( m_neighbors.size() );
//  array<MPI_Status>  mpiRecvSizeStatus( m_neighbors.size() );
//
//  array<MPI_Request> mpiSendBufferRequest( m_neighbors.size() );
//  array<MPI_Request> mpiRecvBufferRequest( m_neighbors.size() );
//  array<MPI_Status>  mpiSendBufferStatus( m_neighbors.size() );
//  array<MPI_Status>  mpiRecvBufferStatus( m_neighbors.size() );
//
//
//
//  for( unsigned int neighborIndex=0 ; neighborIndex<m_neighbors.size() ;
// ++neighborIndex )
//  {
//    NeighborCommunication& neighbor = m_neighbors[neighborIndex];
//
//    neighbor.tempNeighborData.clear();
//
//    if(pack)
//    neighbor.PackNewGlobalIndexRequests(modifiedObjects);
//    neighbor.SendReceive( neighbor.tempNeighborData.sendGlobalIndexRequests,
// 3, mpiSendGlobalRequest[neighborIndex],
//                          neighbor.tempNeighborData.recvGlobalIndexRequests,
// 3, mpiRecvGlobalRequest[neighborIndex],
//                          CommRegistry::genericComm01 );
//  }
//
//
//  // once we have the number of objects this rank needs to create, we can
// resize() and return the first new
//  // global index to the requesting neighbor.
//  const realT t1=MPI_Wtime();
//
//  array<MPI_Request> mpiSendNewGlobalRequest( m_neighbors.size() );
//  array<MPI_Request> mpiRecvNewGlobalRequest( m_neighbors.size() );
//  array<MPI_Status>  mpiSendNewGlobalStatus( m_neighbors.size() );
//  array<MPI_Status>  mpiRecvNewGlobalStatus( m_neighbors.size() );
//
//  lSet newNodeGlobals;
//  lSet newEdgeGlobals;
//  lSet newFaceGlobals;
//  for( unsigned int count=0 ; count<m_neighbors.size() ; ++count )
//  {
//    int neighborIndex;
//    MPI_Waitany( mpiRecvGlobalRequest.size(), mpiRecvGlobalRequest.data(),
// &neighborIndex, mpiRecvGlobalStatus.data() );
//
//    NeighborCommunication& neighbor = m_neighbors[neighborIndex];
//
//    if(pack)
//    neighbor.ProcessNewGlobalIndexRequests( newNodeGlobals, newEdgeGlobals,
// newFaceGlobals );
//    neighbor.SendReceive( neighbor.tempNeighborData.sendFirstNewGlobalIndices,
// 3, mpiSendNewGlobalRequest[neighborIndex],
//                          neighbor.tempNeighborData.recvFirstNewGlobalIndices,
// 3, mpiRecvNewGlobalRequest[neighborIndex],
//                          CommRegistry::genericComm02 );
//  }
//
//  // now we can assign global indices based on the first recieved by the
// neighbor.
//  for( unsigned int count=0 ; count<m_neighbors.size() ; ++count )
//  {
//    int neighborIndex;
//    MPI_Waitany( mpiRecvNewGlobalRequest.size(),
// mpiRecvNewGlobalRequest.data(), &neighborIndex, mpiRecvNewGlobalStatus.data()
// );
//    NeighborCommunication& neighbor = m_neighbors[neighborIndex];
//    if(pack)
//    neighbor.UnpackNewGlobalIndices(  );
//
//  }
//
//
//
//  MPI_Waitall( mpiSendGlobalRequest.size() , mpiSendGlobalRequest.data(),
// mpiSendGlobalStatus.data() );
//  MPI_Waitall( mpiSendNewGlobalRequest.size() ,
// mpiSendNewGlobalRequest.data(), mpiSendNewGlobalStatus.data() );
//
//  const realT t2=MPI_Wtime();
//  const realT t3=MPI_Wtime();
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//  // first we need to send over the new objects to our neighbors. This will
// consist of all new objects owned by this
//  // partition that are ghosts on the neighbor partition, or are owned on the
// neighbor but a ghost on this partition.
//
//
//
//  //******************************
//  // 1) for a given processor, send new objects created on this processor and
// the connectivities to neighbors
//  //    that either are a ghost on this processor(owned by neighbor), and the
// ones that are a ghost on the neighbor.
//
//  // pack the buffers, and send the size of the buffers
//  for( unsigned int neighborIndex=0 ; neighborIndex<m_neighbors.size() ;
// ++neighborIndex )
//  {
//    NeighborCommunication& neighbor = m_neighbors[neighborIndex];
//
//    neighbor.ResizeSendBuffer(0);
//
//    if(pack)
//    {
//    // pack the new/modified nodes that are owned by this process
//    neighbor.PackTopologyModifications(
// DomainPartition::FiniteElementNodeManager, modifiedObjects.newNodes,
//  modifiedObjects.modifiedNodes, true, false );
//    // pack the new/modified edges that are owned by the neighbor.
//    neighbor.PackTopologyModifications(
// DomainPartition::FiniteElementEdgeManager, modifiedObjects.newEdges,
//  modifiedObjects.modifiedEdges, true,  true );
//    // pack the new/modified edges that are owned by this process
//    neighbor.PackTopologyModifications(
// DomainPartition::FiniteElementEdgeManager, modifiedObjects.newEdges,
//  modifiedObjects.modifiedEdges, true,  false );
//
//    // pack the new/modified faces that are owned by the neighbor.
//    neighbor.PackTopologyModifications(
// DomainPartition::FiniteElementFaceManager, modifiedObjects.newFaces,
//  modifiedObjects.modifiedFaces, true,  true );
//    // pack the new/modified faces that are owned by this process
//    neighbor.PackTopologyModifications(
// DomainPartition::FiniteElementFaceManager, modifiedObjects.newFaces,
//  modifiedObjects.modifiedFaces, true,  false );
//
//    // pack the new/modified elements that are owned by the neighbor.
//    neighbor.PackTopologyModifications(
// DomainPartition::FiniteElementElementManager, std::map< std::string, lSet
// >(),  modifiedObjects.modifiedElements, true,  true );
//    // pack the new/modified elements that are owned by this process
//    neighbor.PackTopologyModifications(
// DomainPartition::FiniteElementElementManager, std::map< std::string, lSet
// >(),  modifiedObjects.modifiedElements, true,  false );
//    }
//    neighbor.SendReceiveBufferSizes(CommRegistry::genericComm01,
// mpiSendSizeRequest[neighborIndex], mpiRecvSizeRequest[neighborIndex] );
//  }
//
//  // send/recv the buffers
//  for( unsigned int count=0 ; count<m_neighbors.size() ; ++count )
//  {
//    int neighborIndex;
//    MPI_Waitany( mpiRecvSizeRequest.size(), mpiRecvSizeRequest.data(),
// &neighborIndex, mpiRecvSizeStatus.data() );
//
//    NeighborCommunication& neighbor = m_neighbors[neighborIndex];
//
//    neighbor.SendReceiveBuffers( CommRegistry::genericComm01,
// CommRegistry::genericComm02, mpiSendBufferRequest[neighborIndex],
// mpiRecvBufferRequest[neighborIndex] );
//  }
//
//
//
//
//
//  lSet allNewAndModifiedLocalNodes;
//
//  lSet allNewNodes, allModifiedNodes;
//  lSet allNewEdges, allModifiedEdges;
//  lSet allNewFaces, allModifiedFaces;
//  std::map< std::string, lSet> allModifiedElements;
//  // unpack the buffers
//  for( unsigned int count=0 ; count<m_neighbors.size() ; ++count )
//  {
//    localIndex_array newLocalNodes, modifiedLocalNodes;
//    localIndex_array newGhostNodes, modifiedGhostNodes;
//
//    localIndex_array newLocalEdges, modifiedLocalEdges;
//    localIndex_array newGhostEdges, modifiedGhostEdges;
//
//    localIndex_array newLocalFaces, modifiedLocalFaces;
//    localIndex_array newGhostFaces, modifiedGhostFaces;
//
//    std::map< std::string, localIndex_array> modifiedElements;
//
//
//
//    int neighborIndex;
//    MPI_Waitany( mpiRecvBufferRequest.size(), mpiRecvBufferRequest.data(),
// &neighborIndex, mpiRecvBufferStatus.data() );
//
//    NeighborCommunication& neighbor = this->m_neighbors[neighborIndex];
//
//    const char* pbuffer = neighbor.ReceiveBuffer().data();
//
//    if(pack)
//    {
//    neighbor.UnpackTopologyModifications(
// DomainPartition::FiniteElementNodeManager, pbuffer, newGhostNodes,
// modifiedGhostNodes, false );
//    neighbor.UnpackTopologyModifications(
// DomainPartition::FiniteElementEdgeManager, pbuffer, newLocalEdges,
// modifiedLocalEdges, true );
//    neighbor.UnpackTopologyModifications(
// DomainPartition::FiniteElementEdgeManager, pbuffer, newGhostEdges,
// modifiedGhostEdges, false );
//    neighbor.UnpackTopologyModifications(
// DomainPartition::FiniteElementFaceManager, pbuffer, newLocalFaces,
// modifiedLocalFaces, true );
//    neighbor.UnpackTopologyModifications(
// DomainPartition::FiniteElementFaceManager, pbuffer, newGhostFaces,
// modifiedGhostFaces, false );
//    neighbor.UnpackTopologyModifications(
// DomainPartition::FiniteElementElementManager, pbuffer, modifiedElements );
//    neighbor.UnpackTopologyModifications(
// DomainPartition::FiniteElementElementManager, pbuffer, modifiedElements );
//    }
//    allNewNodes.insert( newLocalNodes.begin(), newLocalNodes.end() );
//    allNewNodes.insert( newGhostNodes.begin(), newGhostNodes.end() );
//
//    allModifiedNodes.insert( modifiedLocalNodes.begin(),
// modifiedLocalNodes.end() );
//    allModifiedNodes.insert( modifiedGhostNodes.begin(),
// modifiedGhostNodes.end() );
//
//    allNewAndModifiedLocalNodes.insert( newLocalNodes.begin(),
// newLocalNodes.end() );
//    allNewAndModifiedLocalNodes.insert( modifiedLocalNodes.begin(),
// modifiedLocalNodes.end() );
//
//    allNewEdges.insert( newLocalEdges.begin(), newLocalEdges.end() );
//    allNewEdges.insert( newGhostEdges.begin(), newGhostEdges.end() );
//
//    allModifiedEdges.insert( modifiedLocalEdges.begin(),
// modifiedLocalEdges.end() );
//    allModifiedEdges.insert( modifiedGhostEdges.begin(),
// modifiedGhostEdges.end() );
//
//    allNewFaces.insert( newLocalFaces.begin(), newLocalFaces.end() );
//    allNewFaces.insert( newGhostFaces.begin(), newGhostFaces.end() );
//
//    allModifiedFaces.insert( modifiedLocalFaces.begin(),
// modifiedLocalFaces.end() );
//    allModifiedFaces.insert( modifiedGhostFaces.begin(),
// modifiedGhostFaces.end() );
//
//    for( std::map< std::string, localIndex_array>::const_iterator
// i=modifiedElements.begin() ; i!=modifiedElements.end() ; ++i )
//    {
//      allModifiedElements[i->first].insert( i->second.begin(), i->second.end()
// );
//    }
//  }
//
//
//  MPI_Waitall( mpiSendSizeRequest.size() , mpiSendSizeRequest.data(),
// mpiSendSizeStatus.data() );
//  MPI_Waitall( mpiSendBufferRequest.size() , mpiSendBufferRequest.data(),
// mpiSendBufferStatus.data() );
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//  // must send the new objects that are local to this processor, but created
// on a neighbor, to the other neighbors.
//  for( unsigned int neighborIndex=0 ; neighborIndex<m_neighbors.size() ;
// ++neighborIndex )
//  {
//    NeighborCommunication& neighbor = m_neighbors[neighborIndex];
//
//    neighbor.ResizeSendBuffer(0);
//
//    if(pack)
//    {
//    neighbor.PackTopologyModifications(
// DomainPartition::FiniteElementEdgeManager, allNewEdges,  allModifiedEdges,
// false, false );
//    neighbor.PackTopologyModifications(
// DomainPartition::FiniteElementFaceManager, allNewFaces,  allModifiedFaces,
// false, false );
//    neighbor.PackTopologyModifications(
// DomainPartition::FiniteElementElementManager, std::map< std::string, lSet
// >(),  allModifiedElements, false, false );
//    }
//    neighbor.SendReceiveBufferSizes(CommRegistry::genericComm01,
// mpiSendSizeRequest[neighborIndex], mpiRecvSizeRequest[neighborIndex] );
//  }
//
//  // send/recv the buffers
//  for( unsigned int count=0 ; count<m_neighbors.size() ; ++count )
//  {
//    int neighborIndex;
//    MPI_Waitany( mpiRecvSizeRequest.size(), mpiRecvSizeRequest.data(),
// &neighborIndex, mpiRecvSizeStatus.data() );
//
//    NeighborCommunication& neighbor = m_neighbors[neighborIndex];
//
//    neighbor.SendReceiveBuffers( CommRegistry::genericComm01,
// CommRegistry::genericComm02, mpiSendBufferRequest[neighborIndex],
// mpiRecvBufferRequest[neighborIndex] );
//  }
//
//  // unpack the buffers
//  for( unsigned int count=0 ; count<m_neighbors.size() ; ++count )
//  {
//    localIndex_array newLocalNodes, modifiedLocalNodes;
//    localIndex_array newGhostNodes, modifiedGhostNodes;
//
//    localIndex_array newLocalEdges, modifiedLocalEdges;
//    localIndex_array newGhostEdges, modifiedGhostEdges;
//
//    localIndex_array newLocalFaces, modifiedLocalFaces;
//    localIndex_array newGhostFaces, modifiedGhostFaces;
//
//    std::map< std::string, localIndex_array> modifiedElements;
//
//
//
//    int neighborIndex;
//    MPI_Waitany( mpiRecvBufferRequest.size(), mpiRecvBufferRequest.data(),
// &neighborIndex, mpiRecvBufferStatus.data() );
//
//    NeighborCommunication& neighbor = this->m_neighbors[neighborIndex];
//
//    const char* pbuffer = neighbor.ReceiveBuffer().data();
//
//    if(pack)
//    {
//      neighbor.UnpackTopologyModifications(DomainPartition::FiniteElementEdgeManager,
// pbuffer,
//                                           newGhostEdges, modifiedGhostEdges,
// false);
//      neighbor.UnpackTopologyModifications(DomainPartition::FiniteElementFaceManager,
// pbuffer,
//                                           newGhostFaces, modifiedGhostFaces,
// false);
//      neighbor.UnpackTopologyModifications(DomainPartition::FiniteElementElementManager,
// pbuffer,
//                                           modifiedElements);
//    }
//    allNewNodes.insert( newLocalNodes.begin(), newLocalNodes.end() );
//    allNewNodes.insert( newGhostNodes.begin(), newGhostNodes.end() );
//
//    allModifiedNodes.insert( modifiedLocalNodes.begin(),
// modifiedLocalNodes.end() );
//    allModifiedNodes.insert( modifiedGhostNodes.begin(),
// modifiedGhostNodes.end() );
//
//    allNewAndModifiedLocalNodes.insert( newLocalNodes.begin(),
// newLocalNodes.end() );
//    allNewAndModifiedLocalNodes.insert( modifiedLocalNodes.begin(),
// modifiedLocalNodes.end() );
//
//    allNewEdges.insert( newLocalEdges.begin(), newLocalEdges.end() );
//    allNewEdges.insert( newGhostEdges.begin(), newGhostEdges.end() );
//
//    allModifiedEdges.insert( modifiedLocalEdges.begin(),
// modifiedLocalEdges.end() );
//    allModifiedEdges.insert( modifiedGhostEdges.begin(),
// modifiedGhostEdges.end() );
//
//    allNewFaces.insert( newLocalFaces.begin(), newLocalFaces.end() );
//    allNewFaces.insert( newGhostFaces.begin(), newGhostFaces.end() );
//
//    allModifiedFaces.insert( modifiedLocalFaces.begin(),
// modifiedLocalFaces.end() );
//    allModifiedFaces.insert( modifiedGhostFaces.begin(),
// modifiedGhostFaces.end() );
//
//    for( std::map< std::string, localIndex_array>::const_iterator
// i=modifiedElements.begin() ; i!=modifiedElements.end() ; ++i )
//    {
//      allModifiedElements[i->first].insert( i->second.begin(), i->second.end()
// );
//    }
//  }
//
//
//  MPI_Waitall( mpiSendSizeRequest.size() , mpiSendSizeRequest.data(),
// mpiSendSizeStatus.data() );
//  MPI_Waitall( mpiSendBufferRequest.size() , mpiSendBufferRequest.data(),
// mpiSendBufferStatus.data() );
//
//
//
//
//  lSet allReceivedNodes;
//  allReceivedNodes.insert( allNewNodes.begin(), allNewNodes.end() );
//  allReceivedNodes.insert( allModifiedNodes.begin(), allModifiedNodes.end() );
//
//  allReceivedNodes.erase( allNewAndModifiedLocalNodes.begin(),
// allNewAndModifiedLocalNodes.end() );
//
//  m_domain->m_feNodeManager.ConnectivityFromGlobalToLocal(
// allNewAndModifiedLocalNodes,
//                                                           allReceivedNodes,
//                                                           m_domain->m_feFaceManager.m_globalToLocalMap
// );
//
//
//  lSet allReceivedEdges;
//  allReceivedEdges.insert( allNewEdges.begin(), allNewEdges.end() );
//  allReceivedEdges.insert( allModifiedEdges.begin(), allModifiedEdges.end() );
//
//  m_domain->m_feEdgeManager.ConnectivityFromGlobalToLocal( allReceivedEdges,
//                                                           m_domain->m_feNodeManager.m_globalToLocalMap,
//                                                           m_domain->m_feFaceManager.m_globalToLocalMap
// );
//
//  lSet allReceivedFaces;
//  allReceivedFaces.insert( allNewFaces.begin(), allNewFaces.end() );
//  allReceivedFaces.insert( allModifiedFaces.begin(), allModifiedFaces.end() );
//
//  m_domain->m_feFaceManager.ConnectivityFromGlobalToLocal( allReceivedFaces,
//                                                           m_domain->m_feNodeManager.m_globalToLocalMap,
//                                                           m_domain->m_feEdgeManager.m_globalToLocalMap
// );
//
//
//  m_domain->m_feElementManager.ConnectivityFromGlobalToLocal(
// allModifiedElements,
//                                                              m_domain->m_feNodeManager.m_globalToLocalMap,
//                                                              m_domain->m_feFaceManager.m_globalToLocalMap
// );
//
//
//  m_domain->m_feNodeManager.ModifyNodeToEdgeMapFromSplit(
// m_domain->m_feEdgeManager,
//                                                          allNewEdges,
//                                                          allModifiedEdges );
//
//  m_domain->m_feFaceManager.ModifyToFaceMapsFromSplit( allNewFaces,
//                                                       allModifiedFaces,
//                                                       m_domain->m_feNodeManager,
//                                                       m_domain->m_feEdgeManager,
//                                                       m_domain->m_externalFaces
// );
//
//  m_domain->m_feElementManager.ModifyToElementMapsFromSplit(
// allModifiedElements,
//                                                             m_domain->m_feNodeManager,
//                                                             m_domain->m_feFaceManager
// );
//
//  //This is to update the isExternal attribute of tip node/edge that has just
// become tip because of splitting of a ghost node.
//  m_domain->m_feElementManager.UpdateExternalityFromSplit(
// allModifiedElements,
//                                                           m_domain->m_feNodeManager,
//                                                           m_domain->m_feEdgeManager,
//                                                           m_domain->m_feFaceManager);
////
////  m_domain->m_feEdgeManager.UpdateEdgeExternalityFromSplit(
// m_domain->m_feFaceManager,
////                                                            allNewEdges,
////
//                                                            allModifiedEdges);
//
//
//  const realT t4=MPI_Wtime();
//
//
//
//  m_t1 += t1-t0;
//  m_t2 += t2-t1;
//  m_t3 += t3-t2;
//  m_t4 += t4-t0;
//
//
//

}

#elif 0
/**
 * @author settgast
 * @param modifiedObjects
 * This routine will modify the topology maps after a separation event has
 * occurred. This procedure involves:
 *    1) for a given processor, send new objects created on this processor and
 * the connectivities to neighbors
 *       that either are a ghost on this processor(owned by neighbor), and the
 * ones that are a ghost on the neighbor.
 *    2) once the new objects are received for the neighbor, they are unpacked,
 * thus creating new objects on the
 *       current process. Objects owned by this process are assigned a global
 * index.
 *    4) global indices are sent back to the process that created the object.
 *    5) the new objects that are owned by the local process, but created on a
 * neighbor, are then packed and sent to
 *       all the other neighbors that require that object.
 *    6) the upward pointing maps are modified.
 */
void PartitionBase::ModifyGhostsAndNeighborLists( const ModifiedObjectLists& modifiedObjects )
{

  // first we need to send over the new objects to our neighbors. This will
  // consist of all new objects owned by this
  // partition that are ghosts on the neighbor partition, or are owned on the
  // neighbor but a ghost on this partition.

  array<MPI_Request> mpiSendSizeRequest( m_neighbors.size() );
  array<MPI_Request> mpiRecvSizeRequest( m_neighbors.size() );
  array<MPI_Status>  mpiSendSizeStatus( m_neighbors.size() );
  array<MPI_Status>  mpiRecvSizeStatus( m_neighbors.size() );

  array<MPI_Request> mpiSendBufferRequest( m_neighbors.size() );
  array<MPI_Request> mpiRecvBufferRequest( m_neighbors.size() );
  array<MPI_Status>  mpiSendBufferStatus( m_neighbors.size() );
  array<MPI_Status>  mpiRecvBufferStatus( m_neighbors.size() );



  //******************************
  // 1) for a given processor, send new objects created on this processor and
  // the connectivities to neighbors
  //    that either are a ghost on this processor(owned by neighbor), and the
  // ones that are a ghost on the neighbor.

  // pack the buffers, and send the size of the buffers
  for( unsigned int count=0 ; count<m_neighbors.size() ; ++count )
  {
    NeighborCommunication& neighbor = m_neighbors[count];

    neighbor.ResizeSendBuffer(0);

    neighbor.PackTopologyModifications( DomainPartition::FiniteElementNodeManager, modifiedObjects.newNodes,  modifiedObjects.modifiedNodes, true );
    neighbor.PackTopologyModifications( DomainPartition::FiniteElementNodeManager, modifiedObjects.newNodes,  modifiedObjects.modifiedNodes, false );

    neighbor.PackTopologyModifications( DomainPartition::FiniteElementEdgeManager, modifiedObjects.newEdges,  modifiedObjects.modifiedEdges, true );
    neighbor.PackTopologyModifications( DomainPartition::FiniteElementEdgeManager, modifiedObjects.newEdges,  modifiedObjects.modifiedEdges, false );

    neighbor.PackTopologyModifications( DomainPartition::FiniteElementFaceManager, modifiedObjects.newFaces,  modifiedObjects.modifiedFaces, true );
    neighbor.PackTopologyModifications( DomainPartition::FiniteElementFaceManager, modifiedObjects.newFaces,  modifiedObjects.modifiedFaces, false );

    neighbor.PackTopologyModifications( DomainPartition::FiniteElementElementManager, std::map< std::string, lSet >(),  modifiedObjects.modifiedElements,
                                        true );
    neighbor.PackTopologyModifications( DomainPartition::FiniteElementElementManager, std::map< std::string, lSet >(),  modifiedObjects.modifiedElements,
                                        false );

    neighbor.SendReceiveBufferSizes(CommRegistry::genericComm01, mpiSendSizeRequest[count], mpiRecvSizeRequest[count] );
  }

  // send/recv the buffers
  for( unsigned int count=0 ; count<m_neighbors.size() ; ++count )
  {
    int neighborIndex;
    MPI_Waitany( mpiRecvSizeRequest.size(), mpiRecvSizeRequest.data(), &neighborIndex, mpiRecvSizeStatus.data() );

    NeighborCommunication& neighbor = m_neighbors[neighborIndex];

    neighbor.SendReceiveBuffers( CommRegistry::genericComm01, mpiSendBufferRequest[neighborIndex], mpiRecvBufferRequest[neighborIndex] );
  }



  MPI_Waitall( mpiSendSizeRequest.size(), mpiSendSizeRequest.data(), mpiSendSizeStatus.data() );
  MPI_Waitall( mpiSendBufferRequest.size(), mpiSendBufferRequest.data(), mpiSendBufferStatus.data() );



  localIndex_array allNewLocalNodes, allModifiedLocalNodes;
  localIndex_array allNewLocalEdges, allModifiedLocalEdges;
  localIndex_array allNewLocalFaces, allModifiedLocalFaces;
  std::map< std::string, localIndex_array> allModifiedElements;



  // unpack the buffers
  for( unsigned int count=0 ; count<m_neighbors.size() ; ++count )
  {
    localIndex_array newLocalNodes, modifiedLocalNodes;
    localIndex_array newGhostNodes, modifiedGhostNodes;

    localIndex_array newLocalEdges, modifiedLocalEdges;
    localIndex_array newGhostEdges, modifiedGhostEdges;

    localIndex_array newLocalFaces, modifiedLocalFaces;
    localIndex_array newGhostFaces, modifiedGhostFaces;

    std::map< std::string, localIndex_array> modifiedElements;



    int neighborIndex;
    MPI_Waitany( mpiRecvBufferRequest.size(), mpiRecvBufferRequest.data(), &neighborIndex, mpiRecvBufferStatus.data() );

    NeighborCommunication& neighbor = this->m_neighbors[neighborIndex];

    const char* pbuffer = neighbor.ReceiveBuffer().data();

    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementNodeManager, pbuffer, newLocalNodes, modifiedLocalNodes, true );
    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementNodeManager, pbuffer, newGhostNodes, modifiedGhostNodes, false );

    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementEdgeManager, pbuffer, newLocalEdges, modifiedLocalEdges, true );
    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementEdgeManager, pbuffer, newGhostEdges, modifiedGhostEdges, false );

    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementFaceManager, pbuffer, newLocalFaces, modifiedLocalFaces, true );
    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementFaceManager, pbuffer, newGhostFaces, modifiedGhostFaces, false );

    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementElementManager, pbuffer, modifiedElements, true );
    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementElementManager, pbuffer, modifiedElements, false );

    allNewLocalNodes.insert( allNewLocalNodes.end(), newLocalNodes.begin(), newLocalNodes.end() );
    allModifiedLocalNodes.insert( allModifiedLocalNodes.end(), modifiedLocalNodes.begin(), modifiedLocalNodes.end() );
    allNewLocalEdges.insert( allNewLocalEdges.end(), newLocalEdges.begin(), newLocalEdges.end() );
    allModifiedLocalEdges.insert( allModifiedLocalEdges.end(), modifiedLocalEdges.begin(), modifiedLocalEdges.end() );
    allNewLocalFaces.insert( allNewLocalFaces.end(), newLocalFaces.begin(), newLocalFaces.end() );
    allModifiedLocalFaces.insert( allModifiedLocalFaces.end(), modifiedLocalFaces.begin(), modifiedLocalFaces.end() );

    for( std::map< std::string, localIndex_array>::const_iterator i=modifiedElements.begin() ; i!=modifiedElements.end() ; ++i )
    {
      allModifiedElements[i->first].insert( allModifiedElements[i->first].end(), i->second.begin(), i->second.end() );
    }

    // pack the objects that have new global indices into the buffers for the
    // neighbor that created them. This has to
    // be done because the regular unpacking routine will not know where to put
    // them on the neighbor, as there is no
    // valid global number on the neighbor.
    neighbor.PackNewGlobalNumbers(  );
  }



  for( unsigned int count=0 ; count<m_neighbors.size() ; ++count )
  {
    NeighborCommunication& neighbor = m_neighbors[count];


    neighbor.PackTopologyModifications( DomainPartition::FiniteElementNodeManager, allNewLocalNodes,  allModifiedLocalNodes, false );
    neighbor.PackTopologyModifications( DomainPartition::FiniteElementEdgeManager, allNewLocalEdges,  allModifiedLocalEdges, false );
    neighbor.PackTopologyModifications( DomainPartition::FiniteElementFaceManager, allNewLocalFaces,  allModifiedLocalFaces, false );
    neighbor.PackTopologyModifications( DomainPartition::FiniteElementElementManager, std::map< std::string, lSet >(),  modifiedObjects.modifiedElements,
                                        false );


    neighbor.SendReceiveBufferSizes(CommRegistry::genericComm01, mpiSendSizeRequest[count], mpiRecvSizeRequest[count] );


  }
  // send/recv the buffers
  for( unsigned int count=0 ; count<m_neighbors.size() ; ++count )
  {
    int neighborIndex;
    MPI_Waitany( mpiRecvSizeRequest.size(), mpiRecvSizeRequest.data(), &neighborIndex, mpiRecvSizeStatus.data() );

    NeighborCommunication& neighbor = m_neighbors[neighborIndex];

    neighbor.SendReceiveBuffers( CommRegistry::genericComm01, mpiSendBufferRequest[neighborIndex], mpiRecvBufferRequest[neighborIndex] );
  }
  // unpack the buffers
  for( unsigned int count=0 ; count<m_neighbors.size() ; ++count )
  {
    localIndex_array newGhostNodes, modifiedGhostNodes;
    localIndex_array newGhostEdges, modifiedGhostEdges;
    localIndex_array newGhostFaces, modifiedGhostFaces;

    std::map< std::string, localIndex_array> modifiedElements;



    int neighborIndex;
    MPI_Waitany( mpiRecvBufferRequest.size(), mpiRecvBufferRequest.data(), &neighborIndex, mpiRecvBufferStatus.data() );

    NeighborCommunication& neighbor = this->m_neighbors[neighborIndex];

    const char* pbuffer = neighbor.ReceiveBuffer().data();

    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementNodeManager, pbuffer, newGhostNodes, modifiedGhostNodes, false );
    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementEdgeManager, pbuffer, newGhostEdges, modifiedGhostEdges, false );
    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementFaceManager, pbuffer, newGhostFaces, modifiedGhostFaces, false );
    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementElementManager, pbuffer, modifiedElements, false );
  }


  MPI_Waitall( mpiSendSizeRequest.size(), mpiSendSizeRequest.data(), mpiSendSizeStatus.data() );
  MPI_Waitall( mpiSendBufferRequest.size(), mpiSendBufferRequest.data(), mpiSendBufferStatus.data() );



}

#elif 0
void PartitionBase::ModifyGhostsAndNeighborLists( const ModifiedObjectLists& modifiedObjects )
{

  array<MPI_Request> mpiSendSizeRequest( m_neighbors.size() );
  array<MPI_Request> mpiRecvSizeRequest( m_neighbors.size() );
  array<MPI_Status>  mpiSendSizeStatus( m_neighbors.size() );
  array<MPI_Status>  mpiRecvSizeStatus( m_neighbors.size() );

  array<MPI_Request> mpiSendBufferRequest( m_neighbors.size() );
  array<MPI_Request> mpiRecvBufferRequest( m_neighbors.size() );
  array<MPI_Status>  mpiSendBufferStatus( m_neighbors.size() );
  array<MPI_Status>  mpiRecvBufferStatus( m_neighbors.size() );



  //***** Reverse communication to update any changed ghost objects on their
  // parent domains

  // pack the buffers, and send the size of the buffers
  {
    int neighborNum = 0;
    for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin() ; neighbor!=m_neighbors.end() ; ++neighbor, ++neighborNum )
    {
      neighbor->ResizeSendBuffer(0);
      array<bufvector::size_type> sendSizeArray;
      neighbor->PackTopologyModifications( DomainPartition::FiniteElementNodeManager, modifiedObjects.newNodes,  modifiedObjects.modifiedNodes, true );
      neighbor->PackTopologyModifications( DomainPartition::FiniteElementEdgeManager, modifiedObjects.newEdges,  modifiedObjects.modifiedEdges, true );
      neighbor->PackTopologyModifications( DomainPartition::FiniteElementFaceManager, modifiedObjects.newFaces,  modifiedObjects.modifiedFaces, true );
      neighbor->PackTopologyModifications( DomainPartition::FiniteElementElementManager, lSet(),  modifiedObjects.modifiedElements, true );

      neighbor->SendReceiveBufferSizes(CommRegistry::genericComm01, mpiSendSizeRequest[neighborNum], mpiRecvSizeRequest[neighborNum] );
    }
  }


  // send/recv the buffers
  for( unsigned int count=0 ; count<m_neighbors.size() ; ++count )
  {
    int neighborIndex;
    MPI_Waitany( mpiRecvSizeRequest.size(), mpiRecvSizeRequest.data(), &neighborIndex, mpiRecvSizeStatus.data() );

    NeighborCommunication& neighbor = this->m_neighbors[neighborIndex];

    neighbor.SendReceiveBuffers( CommRegistry::genericComm01, mpiSendBufferRequest[neighborIndex], mpiRecvBufferRequest[neighborIndex] );

  }


  lSet localNewNodes, localModifiedNodes;
  lSet localNewEdges, localModifiedEdges;
  lSet localNewFaces, localModifiedFaces;
  std::map< std::string, lSet> localModifiedElements;

  // unpack the buffers
  for( unsigned int count=0 ; count<m_neighbors.size() ; ++count )
  {
    localIndex_array newNodes, modifiedNodes;
    localIndex_array newEdges, modifiedEdges;
    localIndex_array newFaces, modifiedFaces;
    std::map< std::string, localIndex_array> modifiedElements;

    int neighborIndex;
    MPI_Waitany( mpiRecvBufferRequest.size(), mpiRecvBufferRequest.data(), &neighborIndex, mpiRecvBufferStatus.data() );

    NeighborCommunication& neighbor = this->m_neighbors[neighborIndex];

    const char* pbuffer = neighbor.ReceiveBuffer().data();

    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementNodeManager, pbuffer, newNodes, modifiedNodes, true );
    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementEdgeManager, pbuffer, newEdges, modifiedEdges, true );
    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementFaceManager, pbuffer, newFaces, modifiedFaces, true );
    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementElementManager, pbuffer, modifiedElements, true );

    localNewNodes.insert( newNodes.begin(), newNodes.end() );
    localModifiedNodes.insert( modifiedNodes.begin(), modifiedNodes.end() );

    localNewEdges.insert( newEdges.begin(), newEdges.end() );
    localModifiedEdges.insert( modifiedEdges.begin(), modifiedEdges.end() );

    localNewFaces.insert( newFaces.begin(), newFaces.end() );
    localModifiedFaces.insert( modifiedFaces.begin(), modifiedFaces.end() );


    for( std::map< std::string, localIndex_array>::const_iterator i=modifiedElements.begin() ; i!=modifiedElements.end() ; ++i )
    {
      localModifiedElements[i->first].insert( i->second.begin(), i->second.end() );
    }
  }

  m_domain->m_feNodeManager.ModifyNodeToEdgeMapFromSplit( m_domain->m_feEdgeManager,
                                                          localNewEdges,
                                                          localModifiedEdges );

  m_domain->m_feFaceManager.ModifyToFaceMapsFromSplit( localNewFaces,
                                                       localModifiedFaces,
                                                       m_domain->m_feNodeManager,
                                                       m_domain->m_feEdgeManager );

  m_domain->m_feElementManager.ModifyToElementMapsFromSplit( localModifiedElements,
                                                             m_domain->m_feNodeManager,
                                                             m_domain->m_feFaceManager );

  MPI_Waitall( mpiSendSizeRequest.size(), mpiSendSizeRequest.data(), mpiSendSizeStatus.data() );
  MPI_Waitall( mpiSendBufferRequest.size(), mpiSendBufferRequest.data(), mpiSendBufferStatus.data() );



  //***** Now push from the owning domain to all ghost on the

  localNewNodes.insert( modifiedObjects.newNodes.begin(), modifiedObjects.newNodes.end() );
  localModifiedNodes.insert( modifiedObjects.modifiedNodes.begin(), modifiedObjects.modifiedNodes.end() );

  localNewEdges.insert( modifiedObjects.newEdges.begin(), modifiedObjects.newEdges.end() );
  localModifiedEdges.insert( modifiedObjects.modifiedEdges.begin(), modifiedObjects.modifiedEdges.end() );

  localNewFaces.insert( modifiedObjects.newFaces.begin(), modifiedObjects.newFaces.end() );
  localModifiedFaces.insert( modifiedObjects.modifiedFaces.begin(), modifiedObjects.modifiedFaces.end() );


  for( std::map< std::string, lSet>::const_iterator i=modifiedObjects.modifiedElements.begin() ; i!=modifiedObjects.modifiedElements.end() ; ++i )
  {
    localModifiedElements[i->first].insert( i->second.begin(), i->second.end() );
  }



  // pack the buffers, and send the size of the buffers
  {
    int neighborNum = 0;
    for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin() ; neighbor!=m_neighbors.end() ; ++neighbor, ++neighborNum )
    {
      neighbor->ResizeSendBuffer(0);
      array<bufvector::size_type> sendSizeArray;
      neighbor->PackTopologyModifications( DomainPartition::FiniteElementNodeManager, localNewNodes,  localModifiedNodes, false );
      neighbor->PackTopologyModifications( DomainPartition::FiniteElementEdgeManager, localNewEdges,  localModifiedEdges, false );
      neighbor->PackTopologyModifications( DomainPartition::FiniteElementFaceManager, localNewFaces,  localModifiedFaces, false );
      neighbor->PackTopologyModifications( DomainPartition::FiniteElementElementManager, lSet(),  localModifiedElements, false );

      neighbor->SendReceiveBufferSizes(CommRegistry::genericComm01, mpiSendSizeRequest[neighborNum], mpiRecvSizeRequest[neighborNum] );
    }
  }


  // send/recv the buffers
  for( unsigned int count=0 ; count<m_neighbors.size() ; ++count )
  {
    int neighborIndex;
    MPI_Waitany( mpiRecvSizeRequest.size(), mpiRecvSizeRequest.data(), &neighborIndex, mpiRecvSizeStatus.data() );

    NeighborCommunication& neighbor = this->m_neighbors[neighborIndex];

    neighbor.SendReceiveBuffers( CommRegistry::genericComm01, mpiSendBufferRequest[neighborIndex], mpiRecvBufferRequest[neighborIndex] );

  }


  lSet ghostNewNodes, ghostModifiedNodes;
  lSet ghostNewEdges, ghostModifiedEdges;
  lSet ghostNewFaces, ghostModifiedFaces;
  std::map< std::string, lSet> ghostModifiedElements;


  // unpack the buffers
  for( unsigned int count=0 ; count<m_neighbors.size() ; ++count )
  {
    localIndex_array newNodes, modifiedNodes;
    localIndex_array newEdges, modifiedEdges;
    localIndex_array newFaces, modifiedFaces;
    std::map< std::string, localIndex_array> modifiedElements;

    int neighborIndex;
    MPI_Waitany( mpiRecvBufferRequest.size(), mpiRecvBufferRequest.data(), &neighborIndex, mpiRecvBufferStatus.data() );

    NeighborCommunication& neighbor = this->m_neighbors[neighborIndex];

    const char* pbuffer = neighbor.ReceiveBuffer().data();

    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementNodeManager, pbuffer, newNodes, modifiedNodes, false );
    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementEdgeManager, pbuffer, newEdges, modifiedEdges, false );
    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementFaceManager, pbuffer, newFaces, modifiedFaces, false );
    neighbor.UnpackTopologyModifications( DomainPartition::FiniteElementElementManager, pbuffer, modifiedElements, false );

    ghostNewNodes.insert( newNodes.begin(), newNodes.end() );
    ghostModifiedNodes.insert( modifiedNodes.begin(), modifiedNodes.end() );

    ghostNewEdges.insert( newEdges.begin(), newEdges.end() );
    ghostModifiedEdges.insert( modifiedEdges.begin(), modifiedEdges.end() );

    ghostNewFaces.insert( newFaces.begin(), newFaces.end() );
    ghostModifiedFaces.insert( modifiedFaces.begin(), modifiedFaces.end() );


    for( std::map< std::string, localIndex_array>::const_iterator i=modifiedElements.begin() ; i!=modifiedElements.end() ; ++i )
    {
      ghostModifiedElements[i->first].insert( i->second.begin(), i->second.end() );
    }
  }


  m_domain->m_feNodeManager.ModifyNodeToEdgeMapFromSplit( m_domain->m_feEdgeManager,
                                                          ghostNewEdges,
                                                          ghostModifiedEdges );

  m_domain->m_feFaceManager.ModifyToFaceMapsFromSplit( ghostNewFaces,
                                                       ghostModifiedFaces,
                                                       m_domain->m_feNodeManager,
                                                       m_domain->m_feEdgeManager );

  m_domain->m_feElementManager.ModifyToElementMapsFromSplit( ghostModifiedElements,
                                                             m_domain->m_feNodeManager,
                                                             m_domain->m_feFaceManager );

  MPI_Waitall( mpiSendSizeRequest.size(), mpiSendSizeRequest.data(), mpiSendSizeStatus.data() );
  MPI_Waitall( mpiSendBufferRequest.size(), mpiSendBufferRequest.data(), mpiSendBufferStatus.data() );



}

#endif


/**
 *
 * @param fieldNames
 * @param commID
 */
void PartitionBase::SetBufferSizes( const std::map<string, array<string> >& fieldNames,
                                    const CommRegistry::commID commID  )
{
  // get buffer sizes, and send/receive sizes
  for( array<NeighborCommunication>::iterator neighbor=m_neighbors.begin() ; neighbor!=m_neighbors.end() ; ++neighbor )
  {

    neighbor->GetPackedBufferSize( fieldNames,
                                   commID );

    neighbor->SendReceiveSizes(commID);
  }

  for( array<NeighborCommunication>::iterator neighbor=m_neighbors.begin() ; neighbor!=m_neighbors.end() ; ++neighbor )
  {
    neighbor->MPI_Wait_RecvSizeRequest(commID);
    neighbor->MPI_Wait_SendSizeRequest(commID);
  }
}

void PartitionBase::SynchronizeFields( const std::map<string, array<string> >& fieldNames,
                                       const CommRegistry::commID commID  )
{
//
//  if(m_hasLocalGhosts){
//    // synchronize local ghosts prior to external communication
//    std::map<DomainPartition::ObjectDataStructureKeys,
// array<string>>::const_iterator it  = fieldNames.begin();
//    std::map<DomainPartition::ObjectDataStructureKeys,
// array<string>>::const_iterator iend  = fieldNames.end();
//    for(;it!=iend;++it){
//
//      if(it->first == DomainPartition::FiniteElementElementManager)
//      {
//        for( std::map< ElementManagerT::RegKeyType, ElementRegionT >::iterator
// ielemReg=m_domain->m_feElementManager.m_ElementRegions.begin() ;
//             ielemReg!=m_domain->m_feElementManager.m_ElementRegions.end() ;
// ++ielemReg )
//        {
//          ElementRegionT& elementRegion = ielemReg->second;
//
//          const std::string& regionName = ielemReg->first;
//
//          const localIndex_array& source =
// m_elementRegionsLocalGhostSources[regionName];
//          const localIndex_array& target =
// m_elementRegionsLocalGhosts[regionName];
//
//          for( array<string>::size_type i =0; i < it->second.size(); ++i){
//            const std::string& fieldName = it->second[i];
//            FieldType fieldType = elementRegion.GetFieldType(fieldName);
//
//            if(fieldType != FieldInfo::numFieldTypes)
//              elementRegion.CopyFieldSubset( fieldType, fieldName,  source,
// target);
//
//          }
//
//        }
//      } else {
//        ObjectDataStructureBaseT& object =
// m_domain->GetObjectDataStructure(it->first);
//
//        const localIndex_array& source = m_localGhostSources[it->first];
//        const localIndex_array& target = m_localGhosts[it->first];
//
//        for( array<string>::size_type i =0; i < it->second.size(); ++i){
//          const std::string& fieldName = it->second[i];
//          FieldType fieldType = object.GetFieldType(fieldName);
//
//          if(fieldType != FieldInfo::numFieldTypes)
//            object.CopyFieldSubset( fieldType, fieldName,  source, target);
//
//        }
//
//      }
//    }
//  }
//
//  SetBufferSizes(fieldNames, commID);
//


#if 0
  // send and receive buffers
  for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin() ; neighbor!=m_neighbors.end() ; ++neighbor )
  {
    neighbor->PackBuffer( fieldNames, commID  );
    neighbor->SendReceiveBuffers( commID);
    neighbor->m_unpacked = false;
  }

  // unpack buffers
  for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin() ; neighbor!=m_neighbors.end() ; ++neighbor )
  {
    neighbor->MPI_Wait_RecvBufferRequest(commID);
    //neighbor->UnpackBuffer( fieldNames );
  }

  for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.begin() ; neighbor!=m_neighbors.end() ; ++neighbor )
  {
    neighbor->MPI_Wait_SendBufferRequest(commID);
  }
#else

  array<MPI_Request> mpiSendBufferRequest( m_neighbors.size() );
  array<MPI_Request> mpiRecvBufferRequest( m_neighbors.size() );
  array<MPI_Status>  mpiSendBufferStatus( m_neighbors.size() );
  array<MPI_Status>  mpiRecvBufferStatus( m_neighbors.size() );

  // send and receive buffers
  for( int neighborIndex=0 ; neighborIndex<m_neighbors.size() ; ++neighborIndex )
  {
    NeighborCommunication& neighbor = m_neighbors[neighborIndex];
    neighbor.PackBuffer( fieldNames, commID  );
    neighbor.SendReceiveBuffers( commID, mpiSendBufferRequest[neighborIndex], mpiRecvBufferRequest[neighborIndex] );

  }


  // unpack the buffers
  for( int count=0 ; count<m_neighbors.size() ; ++count )
  {
    int neighborIndex;
    MPI_Waitany( mpiRecvBufferRequest.size(), mpiRecvBufferRequest.data(), &neighborIndex, mpiRecvBufferStatus.data() );

    NeighborCommunication& neighbor = this->m_neighbors[neighborIndex];
    neighbor.UnpackBuffer( fieldNames );
  }

  MPI_Waitall( mpiSendBufferRequest.size(), mpiSendBufferRequest.data(), mpiSendBufferStatus.data() );


#endif

}


/**
 *
 * Runs through exchange indicies (send recv) and if it is a recv index it sets
 * the object's ghost rank.
 *
 */
void PartitionBase::SetGhostArrays( DomainPartition * domain )
{
  //(1) initialize ghost arrays
//  std::map<DomainPartition::ObjectDataStructureKeys,
// array<Field<FieldInfo::ghostRank>::Type>*> ghostRank;
//  {
//    const localIndex n = NeighborCommunication::NumberOfSyncNames();
//    array<DomainPartition::ObjectDataStructureKeys> objectNames;
//    NeighborCommunication::SyncNames(objectNames);
//
//    for(localIndex i = 0; i < n; ++i)
//    {
//      if(objectNames[i] != DomainPartition::FiniteElementElementManager)
//      {
//        ghostRank[objectNames[i]] =
// &(domain->GetObjectDataStructure(objectNames[i]).GetFieldData<FieldInfo::ghostRank>());
//        *ghostRank[objectNames[i]] = INT_MIN;
//      }
//      else
//      {
//        ghostRank[DomainPartition::FiniteElementElementManager] = NULL;
//        //--per usual: elements are treated as a special case
//        for( std::map< ElementManagerT::RegKeyType, ElementRegionT >::iterator
// iregion=domain->m_feElementManager.m_ElementRegions.begin() ;
//             iregion!=domain->m_feElementManager.m_ElementRegions.end() ;
// ++iregion )
//        {
//          ElementRegionT& elemRegion = iregion->second;
//          array<Field<FieldInfo::ghostRank>::Type>& ghostRankCurr =
// elemRegion.GetFieldData<FieldInfo::ghostRank>();
//          ghostRankCurr = INT_MIN;
//        }
//      }
//    }
//  }
//
//  //(2) set values in individual arrays
//  std::map<DomainPartition::ObjectDataStructureKeys,
// array<Field<FieldInfo::ghostRank>::Type>*>::iterator it;
//  for(it = ghostRank.begin(); it != ghostRank.end(); ++it)
//  {
//    if(it->first != DomainPartition::FiniteElementElementManager)
//    {
//      for( VectorT<NeighborCommunication>::iterator
// neighbor=m_neighbors.begin() ; neighbor!=m_neighbors.end() ; ++neighbor )
//      {
//        const int neighborRank = neighbor->NeighborRank();
//        array<Field<FieldInfo::ghostRank>::Type>& ghostRankCurr =
// *ghostRank[it->first];
//
//        const localIndex_array& receiveLocalIndices =
// neighbor->ReceiveLocalIndices(it->first);
//        for( localIndex_array::const_iterator i=receiveLocalIndices.begin() ;
// i!=receiveLocalIndices.end() ; ++i )
//        {
//          ghostRankCurr[*i] = neighborRank;
//        }
//
//        const localIndex_array& sendLocalIndices =
// neighbor->SendLocalIndices(it->first);
//        for( localIndex_array::const_iterator i=sendLocalIndices.begin() ;
// i!=sendLocalIndices.end() ; ++i )
//        {
//          ghostRankCurr[*i] = -1;
//        }
//      }
//    }
//    else
//    {
//      for( VectorT<NeighborCommunication>::iterator
// neighbor=m_neighbors.begin() ; neighbor!=m_neighbors.end() ; ++neighbor )
//      {
//        //--per usual: elements are treated as a special case
//        const int neighborRank = neighbor->NeighborRank();
//        const std::map<std::string,localIndex_array>&
// elementRegionsReceiveLocalIndices =
// neighbor->ElementRegionsReceiveLocalIndices();
//        for( std::map< ElementManagerT::RegKeyType, ElementRegionT >::iterator
// iregion=domain->m_feElementManager.m_ElementRegions.begin() ;
//             iregion!=domain->m_feElementManager.m_ElementRegions.end() ;
// ++iregion )
//        {
//          const std::string& elemRegionName = iregion->first;
//          ElementRegionT& elemRegion = iregion->second;
//          array<Field<FieldInfo::ghostRank>::Type>& ghostRankCurr =
// elemRegion.GetFieldData<FieldInfo::ghostRank>();
//
//
//          const localIndex_array* const receiveLocalIndices =
// stlMapLookupPointer( elementRegionsReceiveLocalIndices,
//                                                                           elemRegionName
// );
//
//          //const localIndex_array& receiveLocalIndices =
// neighbor->ElementRegionReceiveLocalIndices(elemRegionName);
//          if( receiveLocalIndices )
//          for( localIndex_array::const_iterator i=receiveLocalIndices->begin()
// ; i!=receiveLocalIndices->end() ; ++i )
//          {
//            ghostRankCurr[*i] = neighborRank;
//          }
//
//          const localIndex_array* const sendLocalIndices =
// stlMapLookupPointer( neighbor->ElementRegionsSendLocalIndices() ,
//                                                                             elemRegionName
// );
//          if( sendLocalIndices )
//          for( localIndex_array::const_iterator i=sendLocalIndices->begin() ;
// i!=sendLocalIndices->end() ; ++i )
//          {
//            ghostRankCurr[*i] = -1;
//          }
//        }
//      }
//    }
//  }
}

void PartitionBase::SetRankOfNeighborNeighbors()
{
  array<integer> ranks;
  array<array<integer> > neighborRanks(m_neighbors.size());

  for( array<NeighborCommunication>::const_iterator neighbor=m_neighbors.begin() ; neighbor!=m_neighbors.end() ; ++neighbor )
  {
    ranks.push_back( neighbor->NeighborRank() );
  }

  array<integer>::size_type sendSize = ranks.size();
  array<array<integer>::size_type> recvSize(m_neighbors.size());
  array<MPI_Request> mpiSendRequest(m_neighbors.size());
  array<MPI_Request> mpiRecvRequest(m_neighbors.size());
  array<MPI_Status>  mpiSendStatus(m_neighbors.size());
  array<MPI_Status>  mpiRecvStatus(m_neighbors.size());

  for( int i=0 ; i<m_neighbors.size() ; ++i )
  {
    NeighborCommunication& neighbor = m_neighbors[i];
    neighbor.SendReceive( &sendSize, 1, mpiSendRequest[i],
                          &(recvSize[i]), 1, mpiRecvRequest[i] );
  }
  MPI_Waitall( mpiSendRequest.size(), mpiSendRequest.data(), mpiSendStatus.data() );
  MPI_Waitall( mpiRecvRequest.size(), mpiRecvRequest.data(), mpiRecvStatus.data() );

  for( int i=0 ; i<m_neighbors.size() ; ++i )
  {
    NeighborCommunication& neighbor = m_neighbors[i];
    neighborRanks[i].resize( recvSize[i] );

    neighbor.SendReceive( ranks.data(), ranks.size(), mpiSendRequest[i],
                          neighborRanks[i].data(), neighborRanks[i].size(), mpiRecvRequest[i] );
  }
  MPI_Waitall( mpiSendRequest.size(), mpiSendRequest.data(), mpiSendStatus.data() );
  MPI_Waitall( mpiRecvRequest.size(), mpiRecvRequest.data(), mpiRecvStatus.data() );

  for( int i=0 ; i<m_neighbors.size() ; ++i )
  {
    NeighborCommunication& neighbor = m_neighbors[i];
    neighbor.SetRankOfNeighborNeighbors( neighborRanks[i] );
  }

}

//void PartitionBase::WriteSilo( SiloFile& siloFile )
//{
//  std::string subDirectory =   "PartitionManager";
//  DBMkDir( siloFile.m_dbFilePtr, subDirectory.c_str() );
//  DBSetDir(siloFile.m_dbFilePtr, subDirectory.c_str());
//
//  siloFile.DBWriteWrapper("m_size",m_size);
//  siloFile.DBWriteWrapper("m_rank",m_rank);
//  siloFile.DBWriteWrapper("m_contactGhostMin",m_contactGhostMin);
//  siloFile.DBWriteWrapper("m_contactGhostMax",m_contactGhostMax);
//  siloFile.DBWriteWrapper("m_color",m_color);
//  siloFile.DBWriteWrapper("m_numColors",m_numColors);
//
//  siloFile.DBWriteWrapper("numNeighbors",m_neighbors.size() );
//
//  for( VectorT<NeighborCommunication>::size_type i=0 ; i<m_neighbors.size() ;
// ++i )
//  {
//
//    char dirName[200] = { 0 };
//    sprintf(dirName, "neighbor_%02lu", i);
//    DBMkDir( siloFile.m_dbFilePtr, dirName );
//    DBSetDir(siloFile.m_dbFilePtr, dirName );
//
//    m_neighbors[i].WriteSilo( siloFile );
//
//    DBSetDir(siloFile.m_dbFilePtr, "..");
//
//  }
//
//
//  WriteSiloDerived( siloFile );
//
//  DBSetDir(siloFile.m_dbFilePtr, "..");
//
//}
//
//void PartitionBase::ReadSilo( const SiloFile& siloFile )
//{
//
//  if( DBSetDir(siloFile.m_dbFilePtr, "PartitionManager" ) != -1 )
//  {
//
//    siloFile.DBReadWrapper("m_size",m_size);
//    siloFile.DBReadWrapper("m_rank",m_rank);
//    siloFile.DBReadWrapper("m_contactGhostMin",m_contactGhostMin);
//    siloFile.DBReadWrapper("m_contactGhostMax",m_contactGhostMax);
//    siloFile.DBReadWrapper("m_color",m_color);
//    siloFile.DBReadWrapper("m_numColors",m_numColors);
//
//
//    VectorT<NeighborCommunication>::size_type numNeighbors;
//    siloFile.DBReadWrapper("numNeighbors",numNeighbors );
//
//    m_neighbors.resize(numNeighbors);
//    for( VectorT<NeighborCommunication>::size_type i=0 ; i<numNeighbors ; ++i
// )
//    {
//
//      char dirName[200] = { 0 };
//      sprintf(dirName, "neighbor_%02lu", i);
//      DBSetDir(siloFile.m_dbFilePtr, dirName );
//
//      m_neighbors[i].ReadSilo( siloFile );
//
//      DBSetDir(siloFile.m_dbFilePtr, "..");
//
//    }
//
//
//
//    ReadSiloDerived( siloFile );
//
//    DBSetDir(siloFile.m_dbFilePtr, "..");
//  }
//
//}

//Delete the neighbors that do not communicate
void PartitionBase::DeleteExcessNeighbors()
{
  for( array<NeighborCommunication>::iterator neighbor = m_neighbors.end()-1 ; neighbor!=m_neighbors.begin()-1 ; --neighbor )
  {
//    std::cout << m_rank << ":" << neighbor-> ReturnNeighborRank() << ": " <<
// neighbor->ReturnNeighborRcvSndSize() << std::endl;
    if (neighbor->ReturnNeighborRcvSndSize() == 0)
    {
      neighbor->Clear();
      //m_partition.m_neighbors.erase(neighbor);
      m_neighbors.erase(neighbor);
    }
  }
//  for( VectorT<NeighborCommunication>::iterator neighbor=m_neighbors.end()-1 ;
// neighbor!=m_neighbors.begin()-1 ; --neighbor )
//    {
//      std::cout << m_rank << ":" << neighbor-> ReturnNeighborRank() << ": " <<
// neighbor->ReturnNeighborRcvSndSize() << std::endl;
//    }
}

void PartitionBase::GraphBasedColoring()
//Use the so-called greedy coloring method for metis-based partitions
{
  //First collect the partition graph to rank 0

  if (m_rank == 0)
    std::cout<<"Coloring partitions ... ";
  array<integer> localNeighborList(1);
  localNeighborList = 0;
  for( array<NeighborCommunication>::iterator neighbor = m_neighbors.begin() ; neighbor!=m_neighbors.end() ; ++neighbor )
  {
    localNeighborList.push_back(neighbor->NeighborRank());
  }
  // localNeighborList.insert(localNeighborList.begin(),
  // localNeighborList.size());
  localNeighborList[0] = localNeighborList.size() - 1;

  int localNumNeighbors = localNeighborList.size();
  int maxLocalNumNeighbors = 0;
  MPI_Allreduce(&localNumNeighbors, &maxLocalNumNeighbors, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  array<integer> allNeighborList;
  allNeighborList.resize(maxLocalNumNeighbors * m_size);

  for (localIndex i = 0 ; i < maxLocalNumNeighbors - localNeighborList.size() ; ++i)
  {
    localNeighborList.push_back(-1);
  }

  MPI_Gather( localNeighborList.data(), maxLocalNumNeighbors, MPI_INT, allNeighborList.data(), maxLocalNumNeighbors, MPI_INT, 0, MPI_COMM_WORLD);

  array<integer> colorByRank(m_size);

  if (m_rank == 0)
  {
    array<array<integer> > listNeighbors;
    listNeighbors.resize(m_size);

    int count = 0;
    int COUNT = 0;
    localIndex myrank = -1;
    for (array<integer>::const_iterator it = allNeighborList.begin() ; it != allNeighborList.end() ; ++it)
    {
      if (COUNT == 0)
      {
        count = *it;
        COUNT = maxLocalNumNeighbors;
        myrank++;
      }
      else if (COUNT > 0 && count > 0 && *it >= 0)
      {
        listNeighbors[myrank].push_back(*it);
        count--;
      }
      COUNT--;
    }
    // Now we loop through and assign colors

    colorByRank = -1;
    array<integer> countColor(m_size);
    m_numColors = -1;

    for (localIndex rank = 0 ; rank < listNeighbors.size() ; ++rank)
    {
      countColor = 0;
      for (localIndex i = 0 ; i < listNeighbors[rank].size() ; ++i)
      {
        if (colorByRank[listNeighbors[rank][i]] >= 0)
          countColor[colorByRank[listNeighbors[rank][i]]]++;
      }

      array<integer>::const_iterator it = countColor.begin();
      colorByRank[rank] = 0;
      while (*it > 0)
      {
        ++it;
        ++colorByRank[rank];
      }

      m_numColors = std::max(m_numColors, colorByRank[rank]);
    }

    m_numColors++;
    std::cout<< "Done. " << m_numColors << " colors was used." << std::endl;
    std::cout<< "Rank: color - list of neighbors" << std::endl;
    for (localIndex rank = 0 ; rank < listNeighbors.size() ; ++rank)
    {
      std::cout<< rank << ": " << colorByRank[rank] << " - ";
      for (array<integer>::const_iterator it = listNeighbors[rank].begin() ; it != listNeighbors[rank].end() ; ++it)
      {
        std::cout << *it << " ,";
      }
      std::cout << std::endl;
    }


    //Sanity check
    for (localIndex rank = 0 ; rank < listNeighbors.size() ; ++rank)
    {
      countColor = 0;
      for (localIndex i = 0 ; i < listNeighbors[rank].size() ; ++i)
      {
        if (colorByRank[listNeighbors[rank][i]] ==  colorByRank[rank])
        {
#ifdef USE_ATK
          SLIC_ERROR("ERROR: Two neighbors were assigned the same color.");
#endif
        }
      }
    }

  }

  MPI_Bcast(colorByRank.data(), m_size, MPI_INT, 0, MPI_COMM_WORLD);
  m_color = colorByRank[m_rank];

  MPI_Bcast(&m_numColors, 1, MPI_INT, 0, MPI_COMM_WORLD);
}
}
