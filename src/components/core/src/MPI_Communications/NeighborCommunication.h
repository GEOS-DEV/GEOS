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
 * @file NeighborCommunication.h
 * @author settgast1
 * @date Mar 9, 2011
 */

#ifndef NEIGHBORCOMMUNICATION_H_
#define NEIGHBORCOMMUNICATION_H_
#include "common/DataTypes.hpp"
//#include "../ObjectManagers/PhysicalDomainT.h"
//#include "../Utilities/Utilities.h"
#include "Communication.h"

class PhysicalDomainT;
class ObjectDataStructureBaseT;
class FaceManagerT;
class NodeManager;
class DiscreteElementManagerT;
class EllipsoidalDiscreteElementManagerT;
class oBinStream;
class iBinStream;

struct TempNeighborData
{
//
  //sizes to receive
  Array1dT<gArray1d::size_type> neighborBoundaryObjectsSizes;
  Array1dT<bufvector::size_type> receiveSizes;
  Array1dT<bufvector::size_type> sendSizes;

  //node, edge, face (3)
  std::map<PhysicalDomainT::ObjectDataStructureKeys, gArray1d> neighborNumbers;
  std::map<PhysicalDomainT::ObjectDataStructureKeys, gArray1d> matchedNumbers;
  std::map<PhysicalDomainT::ObjectDataStructureKeys, lArray1d> matchedIndices;

  //node, face (2)
  std::map<PhysicalDomainT::ObjectDataStructureKeys, lSet> indicesInRange;

  lSet objectLocalIndicesToSend[PhysicalDomainT::numObjectDataStructureNames] ;
  gArray1d objectGlobalIndicesToSend[PhysicalDomainT::numObjectDataStructureNames] ;
  gArray1d objectGlobalIndicesToRecieve[PhysicalDomainT::numObjectDataStructureNames] ;

  //node, edge, face, element (4)
  std::map<PhysicalDomainT::ObjectDataStructureKeys, bufvector> objectsToSend;
  std::map<PhysicalDomainT::ObjectDataStructureKeys, bufvector> objectsToReceive;

  std::map<PhysicalDomainT::ObjectDataStructureKeys, lArray1d > objectWithNewGlobalIndex;
  std::map<PhysicalDomainT::ObjectDataStructureKeys, lArray1d > objectWithUnassignedGlobal;

  gArray1d::size_type sendGlobalIndexRequests[3];
  gArray1d::size_type recvGlobalIndexRequests[3];
  globalIndex sendFirstNewGlobalIndices[3];
  globalIndex recvFirstNewGlobalIndices[3];





  void clear()
  {
    //clear sizes
    neighborBoundaryObjectsSizes.clear();
    receiveSizes.clear();
    sendSizes.clear();

    //node, edge, face (3)
    ClearStlMapValues(neighborNumbers);
    ClearStlMapValues(matchedNumbers);
    ClearStlMapValues(matchedIndices);

    //node, face (2)
    ClearStlMapValues(indicesInRange);

    //node, edge, face, element (4)
    ClearStlMapValues(objectsToSend);
    ClearStlMapValues(objectsToReceive);

    ClearStlMapValues(objectWithNewGlobalIndex);
    ClearStlMapValues(objectWithUnassignedGlobal);


    for( int i=0 ; i<3 ; ++i )
    {
      sendGlobalIndexRequests[i] = 0;
      recvGlobalIndexRequests[i] = 0;
      sendFirstNewGlobalIndices[i] = GLOBALINDEX_MAX;
      recvFirstNewGlobalIndices[i] = GLOBALINDEX_MAX;
    }

  }
};


struct ModifiedObjectLists
{
  lSet newNodes;
  lSet newEdges;
  lSet newFaces;
  lSet modifiedNodes;
  lSet modifiedEdges;
  lSet modifiedFaces;
  std::map< std::string, lSet > modifiedElements;

};


/**
 *
 */
class NeighborCommunication
{
public:
  NeighborCommunication();

  NeighborCommunication( const NeighborCommunication& init );

  NeighborCommunication& operator=( const NeighborCommunication& init );

  virtual ~NeighborCommunication();

//  void SetUpExternalFaceNeigborList();

  void Initialize(const int neighborRank, const int rank, const int size );

  void Clear();

  void ResizeSendBuffer( const bufvector::size_type newsize )
  {
    this->m_sendBuffer.resize( newsize );
  }

  void SetDomain( PhysicalDomainT& domain );

  bufvector::size_type PackBuffer( const std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d>& fieldNames,
                                   const CommRegistry::commID commID,
                                   const bool doBufferPacking = 1 );

  bufvector::size_type GetPackedBufferSize( const std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d>& fieldNames,
                                            const CommRegistry::commID commID )
  {

    const bufvector::size_type bufferSize = PackBuffer( fieldNames,
                                                        commID, 0 );
    m_sendBuffer.resize( bufferSize );
    return bufferSize;
  }

  void UnpackBuffer( const std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d>& fieldNames);

  void UnpackGhostElements( std::map<std::string,lArray1d>& newElementIndices );
  int UnpackGhosts(const PhysicalDomainT::ObjectDataStructureKeys name, lArray1d& newIndices );

  void SendReceiveSizes(const CommRegistry::commID commID = CommRegistry::genericComm01  );

  void SendReceiveBufferSizes( const CommRegistry::commID commID,
                               MPI_Request& mpiSendRequest,
                               MPI_Request& mpiRecvRequest );

  void SendReceiveBuffers(const CommRegistry::commID commID = CommRegistry::genericComm01);

  void SendReceiveBuffers(const CommRegistry::commID commID,
                          MPI_Request& mpiSendRequest,
                          MPI_Request& mpiRecvRequest);

  void SendReceiveBuffers(const CommRegistry::commID sizeCommID,
                          const CommRegistry::commID commID,
                          MPI_Request& mpiSendRequest,
                          MPI_Request& mpiRecvRequest);


  void CommunicatePackedObjectBufferSizes(  );



  void DetermineMatchedBoundaryObject( const ObjectDataStructureBaseT& object,
                                       const PhysicalDomainT::ObjectDataStructureKeys name,
                                       const gArray1d& localObjectNumbers);

  void FindPackGhostsDiscreteElement( const bool contactActive, const int depth = 1 );

#ifdef SRC_EXTERNAL
  void FindPackGhostsFaultElement( const int depth = 1 );
#endif

  void FindGhosts( const bool contactActive,
                   const int depth = 1 );

  void FindPackGhosts_Step2();




  const lArray1d& GetSendLocalIndices(const PhysicalDomainT::ObjectDataStructureKeys key){
    return m_sendLocalIndices[key];
  };

  const lArray1d& GetElementRegionSendLocalIndices(const std::string& regionName){
    return  m_elementRegionsSendLocalIndices[regionName];
  };

  const std::map< std::string, lArray1d>& GetElementRegionSendLocalIndices(){
    return  m_elementRegionsSendLocalIndices;
  };




  template< typename TYPE >
  void PackTopologyModifications( const PhysicalDomainT::ObjectDataStructureKeys key,
                                  const TYPE& newObjects,
                                  const TYPE& modifiedObjects,
                                  const bool packConnectivityToGlobal,
                                  const bool reverseOp );




  void UnpackTopologyModifications( const PhysicalDomainT::ObjectDataStructureKeys key,
                                    const char*& pbuffer,
                                    lArray1d& newIndices,
                                    lArray1d& modifiedIndices,
                                    const bool reverseOp  );

  void UnpackTopologyModifications( const PhysicalDomainT::ObjectDataStructureKeys key,
                                    const char*& pbuffer,
                                    std::map< std::string, lArray1d>& modifiedIndices );


  void PackNewGlobalIndexRequests( const ModifiedObjectLists& modifiedObjects );

  void UnpackNewGlobalIndices(  );

  void ProcessNewGlobalIndexRequests( lSet& newNodeGlobals, lSet& newEdgeGlobals, lSet& newFaceGlobals );

  TempNeighborData tempNeighborData;

  template< typename T >
  void SendReceive(T* sendBuffer,
                   const int sendSize,
                   MPI_Request& sendReq,
                   T* recvBuffer,
                   const int recvSize,
                   MPI_Request& recvReq,
                   const CommRegistry::commID commID = CommRegistry::genericComm01 )
  {
    SendReceive( reinterpret_cast<char*>(recvBuffer), recvSize*sizeof(T),
                 reinterpret_cast<char*>(sendBuffer), sendSize*sizeof(T),
                 recvReq, sendReq, commID );
  }

  void SendReceive( bufvector& sendBuffer,
                    MPI_Request& sendReq,
                    bufvector& recvBuffer,
                    MPI_Request& recvReq,
                    const CommRegistry::commID commID = CommRegistry::genericComm01 )
  {
    SendReceive( reinterpret_cast<char*>(recvBuffer.data()), recvBuffer.size()*sizeof(bufvector::size_type),
                 reinterpret_cast<char*>(sendBuffer.data()), sendBuffer.size()*sizeof(bufvector::size_type),
                 recvReq, sendReq, commID );
  }


  void SendReceive(char* inBuffer,
                   const int inSize,
                   char* outBuffer,
                   const int outSize,
                   MPI_Request& inReq,
                   MPI_Request& outReq,
                   const CommRegistry::commID commID = CommRegistry::genericComm01 );

  void WriteSilo( SiloFile& siloFile );

  void ReadSilo( const SiloFile& siloFile );

  int ReturnNeighborRank()
  {
    return m_neighborRank;
  }
  int ReturnNeighborRcvSndSize();





private:

  int m_neighborRank;
  std::map<PhysicalDomainT::ObjectDataStructureKeys, lArray1d> m_receiveLocalIndices;

  void SetAllOwned(const gArray1d& localToGlobal, lArray1d& indices) const;

  int UnpackGhosts(const PhysicalDomainT::ObjectDataStructureKeys name);

  PhysicalDomainT* const m_domain;
  int m_rank;
  int m_size;
  iSet m_rankOfNeighborNeighbors;

  //NOTE: element regions are currently being dealt with using a special structure ... everything else with general one

  std::map<PhysicalDomainT::ObjectDataStructureKeys, lArray1d> m_sendLocalIndices;
  std::map< std::string, lArray1d> m_elementRegionsSendLocalIndices;

  std::map< std::string, lArray1d> m_elementRegionsReceiveLocalIndices;

  bufvector m_sendBuffer;
  int m_sendSize[CommRegistry::maxComm];
  bufvector m_receiveBuffer;
  int m_receiveSize[CommRegistry::maxComm];

  MPI_Request mpiSendSizeRequest[CommRegistry::maxComm];
  MPI_Request mpiRecvSizeRequest[CommRegistry::maxComm];
  MPI_Status  mpiSendSizeStatus[CommRegistry::maxComm];
  MPI_Status  mpiRecvSizeStatus[CommRegistry::maxComm];
  MPI_Request mpiSendBufferRequest[CommRegistry::maxComm];
  MPI_Request mpiRecvBufferRequest[CommRegistry::maxComm];
  MPI_Status  mpiSendBufferStatus[CommRegistry::maxComm];
  MPI_Status  mpiRecvBufferStatus[CommRegistry::maxComm];


public:

  void ClearSendBuffer() { m_sendBuffer.clear(); }

  const bufvector& ReceiveBuffer() const { return this->m_receiveBuffer; }

  void MPI_Wait_SendSizeRequest( const CommRegistry::commID commID );
  void MPI_Wait_RecvSizeRequest( const CommRegistry::commID commID );
  void MPI_Wait_SendBufferRequest( const CommRegistry::commID commID );
  void MPI_Wait_RecvBufferRequest( const CommRegistry::commID commID );

  bool MPI_Test_SendSizeRequest( const CommRegistry::commID commID );
  bool MPI_Test_RecvSizeRequest( const CommRegistry::commID commID );
  bool MPI_Test_SendBufferRequest( const CommRegistry::commID commID );
  bool MPI_Test_RecvBufferRequest( const CommRegistry::commID commID );

  void SetRankOfNeighborNeighbors( const iArray1d& ranks )
  {

    m_rankOfNeighborNeighbors.clear();
    m_rankOfNeighborNeighbors.insert( ranks.begin(), ranks.end() );
  }

  int Rank() const {return m_rank;}
  int NeighborRank() const {return m_neighborRank;}

  inline static localIndex NumberOfSyncNames() { return 8;}
  inline static PhysicalDomainT::ObjectDataStructureKeys SyncName(const localIndex i)
  {
    switch(i)
    {
      case 0:
        return PhysicalDomainT::FiniteElementNodeManager;
      case 1:
        return PhysicalDomainT::FiniteElementEdgeManager;
      case 2:
        return PhysicalDomainT::FiniteElementFaceManager;
      case 3:
        return PhysicalDomainT::FiniteElementElementManager;
      case 4:
        return PhysicalDomainT::DiscreteElementNodeManager;
      case 5:
        return PhysicalDomainT::DiscreteElementFaceManager;
      case 6:
        return PhysicalDomainT::DiscreteElementManager;
      case 7:
        return PhysicalDomainT::EllipsoidalDiscreteElementManager;
      default:
        throw GPException("Unrecognized index for SyncName");
    }
  }
  inline static void SyncNames(Array1dT<PhysicalDomainT::ObjectDataStructureKeys>& syncNames)
  {
    syncNames.resize(8);
    syncNames[0] = PhysicalDomainT::FiniteElementNodeManager;
    syncNames[1] = PhysicalDomainT::FiniteElementEdgeManager;
    syncNames[2] = PhysicalDomainT::FiniteElementFaceManager;
    syncNames[3] = PhysicalDomainT::FiniteElementElementManager;
    syncNames[4] = PhysicalDomainT::DiscreteElementNodeManager;
    syncNames[5] = PhysicalDomainT::DiscreteElementFaceManager;
    syncNames[6] = PhysicalDomainT::DiscreteElementManager;
    syncNames[7] = PhysicalDomainT::EllipsoidalDiscreteElementManager;
  }

  const lArray1d& ElementRegionsReceiveLocalIndices( const std::string& name ) const  {return stlMapLookup( m_elementRegionsReceiveLocalIndices, name );}
  const std::map< std::string, lArray1d>& ElementRegionsReceiveLocalIndices() const  {return m_elementRegionsReceiveLocalIndices;}

  const lArray1d& ElementRegionsSendLocalIndices( const std::string& name ) const  {return stlMapLookup( m_elementRegionsSendLocalIndices, name );}
  const std::map< std::string, lArray1d>& ElementRegionsSendLocalIndices() const  {return m_elementRegionsSendLocalIndices;}

  const lArray1d& ReceiveLocalIndices( const PhysicalDomainT::ObjectDataStructureKeys name) const
  {
    std::map<PhysicalDomainT::ObjectDataStructureKeys, lArray1d>::const_iterator it = this->m_receiveLocalIndices.find(name);
    if(it == this->m_receiveLocalIndices.end())
      throw GPException("Failed to find name " + toString<int>(name) + " in m_receiveLocalIndices (ReceiveLocalIndices function)");
    return it->second;
  }

  const lArray1d& SendLocalIndices( const PhysicalDomainT::ObjectDataStructureKeys name) const
  {
    std::map<PhysicalDomainT::ObjectDataStructureKeys, lArray1d>::const_iterator it = this->m_sendLocalIndices.find(name);
    if(it == this->m_sendLocalIndices.end())
      throw GPException("Failed to find name " + toString<int>(name) + " in m_sendLocalIndices (SendLocalIndices function)");
    return it->second;
  }
};

inline void NeighborCommunication::MPI_Wait_SendSizeRequest( const CommRegistry::commID commID )
{  MPI_Wait( &(mpiSendSizeRequest[commID]), &(mpiSendSizeStatus[commID]) ); }

inline void NeighborCommunication::MPI_Wait_RecvSizeRequest( const CommRegistry::commID commID )
{  MPI_Wait( &(mpiRecvSizeRequest[commID]), &(mpiRecvSizeStatus[commID]) ); }

inline void NeighborCommunication::MPI_Wait_SendBufferRequest( const CommRegistry::commID commID )
{
  if( mpiSendBufferRequest[commID] != MPI_REQUEST_NULL )
    MPI_Wait( &(mpiSendBufferRequest[commID]), &(mpiSendBufferStatus[commID]) );
}

inline void NeighborCommunication::MPI_Wait_RecvBufferRequest( const CommRegistry::commID commID )
{  MPI_Wait( &(mpiRecvBufferRequest[commID]), &(mpiRecvBufferStatus[commID]) ); }


inline bool NeighborCommunication::MPI_Test_SendSizeRequest( const CommRegistry::commID commID )
{
  int flag;
  MPI_Test( &(mpiSendSizeRequest[commID]), &flag, &(mpiSendSizeStatus[commID]) );
  return flag;
}

inline bool NeighborCommunication::MPI_Test_RecvSizeRequest( const CommRegistry::commID commID )
{
  int flag;
  MPI_Test( &(mpiRecvSizeRequest[commID]), &flag, &(mpiRecvSizeStatus[commID]) );
  return flag;
}

inline bool NeighborCommunication::MPI_Test_SendBufferRequest( const CommRegistry::commID commID )
{
  int flag;
  MPI_Test( &(mpiSendBufferRequest[commID]), &flag, &(mpiSendBufferStatus[commID]) );
  return flag;
}

inline bool NeighborCommunication::MPI_Test_RecvBufferRequest( const CommRegistry::commID commID )
{
  int flag;
  MPI_Test( &(mpiRecvBufferRequest[commID]), &flag, &(mpiRecvBufferStatus[commID]) );
  return flag;
}


#endif /* NEIGHBORCOMMUNICATION_H_ */
