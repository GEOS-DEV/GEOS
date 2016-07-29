/*
 * BifroestBase.h
 *
 *  Created on: Aug 15, 2012
 *      Author: johnson346
 */

#ifndef BIFROESTBASE_H_
#define BIFROESTBASE_H_

#include "Common/Common.h"
#include "Common/typedefs.h"
#include "MPI_Communications/Communication.h"
#include <mpi.h>

enum TempBifroestSizeFields
{
  NumberOfNodes,
  NumberOfFaces,
  NumberOfFacePairs,
  SenderRank,
  NumberOfFields
};

enum TempBifroestRequests
{
  RequestNodes,
  RequestFaces,
  RequestFacePairs,
  RequestReturnedData,
  NumberOfRequests
};

struct TempBifroestNodeSendData
{
  //int rankOfOriginatingProcess;
  globalIndex nodeIndex;
  R1Tensor x, dx;
};

struct TempBifroestFaceSendData
{
  //int rankOfOriginatingProcess;
  globalIndex faceIndex;
  R1Tensor normal;
  int numberOfNodes;
  globalIndex nodeIndices[4];
};

struct TempBifroestFacePairSendData
{
  //int rankOfOriginatingProcess;
  globalIndex faceIndex1;
  globalIndex faceIndex2;
};

struct TempBifroestNodeReturnedData
{
  //int rankOfOriginatingProcess;
  globalIndex nodeIndex;
  R1Tensor dforce;
};

struct TempBifroest
{
  //neighbor rank
  int neighborRank;
  bool receiveFlag;

  //data buffers
  iArray1d sizes;
  Array1dT<TempBifroestNodeSendData> nodeData;
  Array1dT<TempBifroestFaceSendData> faceData;
  Array1dT<TempBifroestFacePairSendData> facePairData;

  //returned data buffer
  Array1dT<TempBifroestNodeReturnedData> returnedNodeData;

  Array1dT<MPI_Request> requests;

  void Clear()
  {
    receiveFlag = false;
    nodeData.clear();
    faceData.clear();
    facePairData.clear();

    if(sizes.size() != NumberOfFields)
      sizes.resize(NumberOfFields, -1);
    sizes[NumberOfNodes] = 0;
    sizes[NumberOfFaces] = 0;
    sizes[NumberOfFacePairs] = 0;
    if(requests.size() != NumberOfRequests)
      requests.resize(NumberOfRequests);

    returnedNodeData.clear();
  };

  void ResetReceiveFlag() { receiveFlag = false; };
};

struct DelegateProcessingToAnother
{
  TempBifroest data;

  void SetSizesAndSend(const int m_rank)
  {
    data.sizes.resize(NumberOfFields, 0);
    data.sizes[NumberOfNodes] = data.nodeData.size();
    data.sizes[NumberOfFaces] = data.faceData.size();
    data.sizes[NumberOfFacePairs] = data.facePairData.size();
    data.sizes[SenderRank] = m_rank;
    data.receiveFlag = false;
    data.returnedNodeData.resize(data.sizes[NumberOfNodes]);
    data.requests.resize(NumberOfRequests);

    MPI_Isend(data.sizes.data(), 4, MPI_INTEGER, data.neighborRank,
              CommRegistry::CommTag(m_rank, data.neighborRank, CommRegistry::processShareSendSize),
              MPI_COMM_WORLD, &data.requests[RequestNodes]);
  };

  void SendPackage()
  {
    //NODES
    if(data.sizes[NumberOfNodes] > 0)
    {
      const int tag = CommRegistry::CommTag(data.sizes[SenderRank], data.neighborRank, CommRegistry::processShareSendNodeData);
      const size_t bufferSize = data.sizes[NumberOfNodes] * sizeof(TempBifroestNodeSendData);
      //std::cout << tag << ": (package) being sent on " << data.sizes[SenderRank] << std::endl;
      MPI_Isend( data.nodeData.data(), bufferSize, MPI_CHAR, data.neighborRank, tag, MPI_COMM_WORLD, &data.requests[RequestNodes]);
    }
    //FACES
    if(data.sizes[NumberOfFaces] > 0)
    {
      const int tag = CommRegistry::CommTag(data.sizes[SenderRank], data.neighborRank, CommRegistry::processShareSendFaceData);
      const size_t bufferSize = data.sizes[NumberOfFaces] * sizeof(TempBifroestFaceSendData);
      MPI_Isend( data.faceData.data(), bufferSize, MPI_CHAR, data.neighborRank, tag, MPI_COMM_WORLD, &data.requests[RequestFaces]);
    }
    //FACEPAIRS
    if(data.sizes[NumberOfFacePairs] > 0)
    {
      const int tag = CommRegistry::CommTag(data.sizes[SenderRank], data.neighborRank, CommRegistry::processShareSendFacePairData);
      const size_t bufferSize = data.sizes[NumberOfFacePairs] * sizeof(TempBifroestFacePairSendData);
      MPI_Isend( data.facePairData.data(), bufferSize, MPI_CHAR, data.neighborRank, tag, MPI_COMM_WORLD, &data.requests[RequestFacePairs]);
    }
  };

  void ReceiveReturnedData()
  {
    if(data.sizes[NumberOfNodes] == 0)
      return;
    const size_t bufferSize = data.sizes[NumberOfNodes] * sizeof(TempBifroestNodeReturnedData);
    const int tag = CommRegistry::CommTag(data.neighborRank, data.sizes[SenderRank], CommRegistry::processShareReturn);
    //std::cout << tag << ": (returned data) being received on " << data.sizes[SenderRank] << std::endl;
    MPI_Irecv( data.returnedNodeData.data(), bufferSize, MPI_CHAR, data.neighborRank, tag, MPI_COMM_WORLD, &data.requests[RequestReturnedData] );
  };

  bool ReturnedDataReadyToProcess()
  {
    if(data.sizes[NumberOfNodes] == 0)
      return true;
    int flag;
    MPI_Status status;
    MPI_Test(&data.requests[RequestReturnedData], &flag, &status);
    return flag == 1;
  };
};

struct ProcessingDelegatedToMe
{
  TempBifroest data;

  void Allocate(const int m_rank, const int neighborRank, const iArray1d& tmpSizes)
  {
    data.sizes = tmpSizes;
    data.sizes[SenderRank] = m_rank;

    data.neighborRank = neighborRank;
    data.nodeData.resize(data.sizes[NumberOfNodes]);
    data.faceData.resize(data.sizes[NumberOfFaces]);
    data.facePairData.resize(data.sizes[NumberOfFacePairs]);
    data.returnedNodeData.resize(data.sizes[NumberOfNodes]);

    data.requests.resize(NumberOfRequests);
    data.receiveFlag = false;
  };

  void ReceivePackage()
  {
    //NODES
    if(data.sizes[NumberOfNodes] > 0)
    {
      const int tag = CommRegistry::CommTag(data.neighborRank, data.sizes[SenderRank], CommRegistry::processShareSendNodeData);
      //std::cout << tag << ": (package) being received on " << data.sizes[SenderRank] << std::endl;

      const size_t bufferSize = data.sizes[NumberOfNodes] * sizeof(TempBifroestNodeSendData);
      MPI_Irecv( data.nodeData.data(), bufferSize, MPI_CHAR, data.neighborRank, tag, MPI_COMM_WORLD, &data.requests[RequestNodes]);
    }
    //FACES
    if(data.sizes[NumberOfFaces] > 0)
    {
      const int tag = CommRegistry::CommTag(data.neighborRank, data.sizes[SenderRank], CommRegistry::processShareSendFaceData);
      const size_t bufferSize = data.sizes[NumberOfFaces] * sizeof(TempBifroestFaceSendData);
      MPI_Irecv( data.faceData.data(), bufferSize, MPI_CHAR, data.neighborRank, tag, MPI_COMM_WORLD, &data.requests[RequestFaces]);
    }
    //FACEPAIRS
    if(data.sizes[NumberOfFacePairs] > 0)
    {
      const int tag = CommRegistry::CommTag(data.neighborRank, data.sizes[SenderRank], CommRegistry::processShareSendFacePairData);
      const size_t bufferSize = data.sizes[NumberOfFacePairs] * sizeof(TempBifroestFacePairSendData);
      MPI_Irecv( data.facePairData.data(), bufferSize, MPI_CHAR, data.neighborRank, tag, MPI_COMM_WORLD, &data.requests[RequestFacePairs]);
    }
  };

  bool PackageReadyToProcess()
  {
    int flag;
    MPI_Status status;
    if(data.receiveFlag)
      return true;
    bool ok = true;
    if(data.sizes[NumberOfNodes] > 0)
    {
      MPI_Test(&data.requests[RequestNodes], &flag, &status);
      if(flag != 1)
        ok = false;
    }
//    const int tag = CommRegistry::CommTag(data.neighborRank, data.sizes[SenderRank], CommRegistry::processShareSendNodeData);
    //std::cout << ":: source " << status.MPI_SOURCE << ": " << status.MPI_TAG << " (should be " << tag << ") returning " << flag << std::endl;
    if(data.sizes[NumberOfFaces] > 0)
    {
      MPI_Test(&data.requests[RequestFaces], &flag, &status);
      if(flag != 1)
        ok = false;
    }
    if(data.sizes[NumberOfFacePairs] > 0)
    {
      MPI_Test(&data.requests[RequestFacePairs], &flag, &status);
      if(flag != 1)
        ok = false;
    }
    return ok;
  };

  void SendReturnedData()
  {
    if(data.sizes[NumberOfNodes] == 0)
      return;
    const size_t bufferSize = data.sizes[NumberOfNodes] * sizeof(TempBifroestNodeReturnedData);
    const int tag = CommRegistry::CommTag(data.sizes[SenderRank], data.neighborRank, CommRegistry::processShareReturn);
    //std::cout << tag << ": (returned data) being sent on " << data.sizes[SenderRank] << std::endl;
    MPI_Isend( data.returnedNodeData.data(), bufferSize, MPI_CHAR, data.neighborRank, tag, MPI_COMM_WORLD, &data.requests[RequestReturnedData] );
  };
};


class BifroestBase
{
public:
  BifroestBase();
  virtual ~BifroestBase();
  void Initialize();

  void AddShare(const int destination, const TempBifroestNodeSendData share);
  void AddShare(const int destination, const TempBifroestFaceSendData share);
  void AddShare(const int destination, const TempBifroestFacePairSendData share);

  virtual void CoordinateShares();

  void Synchronize();

  void Clear();

protected:
  virtual void ProcessPackage(ProcessingDelegatedToMe& current) = 0;
  virtual void ProcessReturnedData(DelegateProcessingToAnother& current) = 0;

  int m_rank, m_size;
  std::map<int, DelegateProcessingToAnother> toProcessData;
  std::map<int, ProcessingDelegatedToMe> fromProcessData;
public:
  inline int Rank() const { return m_rank; };
  inline int Size() const { return m_size; };
};

#endif /* BIFROESTBASE_H_ */
