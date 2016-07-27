/*
 * BifroestBase.cpp
 *
 *  Created on: Aug 15, 2012
 *      Author: johnson346
 */

#include "BifroestBase.h"

BifroestBase::BifroestBase()
{
}

BifroestBase::~BifroestBase()
{
}

void BifroestBase::Initialize()
{
  MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &m_size);
}

void BifroestBase::Clear()
{
  for(std::map<int, DelegateProcessingToAnother>::iterator iter = this->toProcessData.begin();
      iter != this->toProcessData.end(); ++iter)
  {
    iter->second.data.Clear();
  }
  for(std::map<int, ProcessingDelegatedToMe>::iterator iter = this->fromProcessData.begin();
      iter != this->fromProcessData.end(); ++iter)
  {
    iter->second.data.Clear();
  }
}

void BifroestBase::AddShare(const int destination, const TempBifroestNodeSendData share)
{
  if(destination==m_rank)
    return;
  std::map<int, DelegateProcessingToAnother>::iterator iter = toProcessData.find( destination );
  if( iter==toProcessData.end()  )
  {
    toProcessData[destination].data.neighborRank = destination;
    toProcessData[destination].data.nodeData.push_back(share);
  }
  else
  {
    iter->second.data.nodeData.push_back(share);
  }
}
void BifroestBase::AddShare(const int destination, const TempBifroestFaceSendData share)
{
  if(destination==m_rank)
    return;
  std::map<int, DelegateProcessingToAnother>::iterator iter = toProcessData.find( destination );
  if( iter==toProcessData.end()  )
  {
    toProcessData[destination].data.faceData.push_back(share);
  }
  else
  {
    iter->second.data.faceData.push_back(share);
  }
}
void BifroestBase::AddShare(const int destination, const TempBifroestFacePairSendData share)
{
  if(destination==m_rank)
    return;
  std::map<int, DelegateProcessingToAnother>::iterator iter = toProcessData.find( destination );
  if( iter==toProcessData.end()  )
  {
    toProcessData[destination].data.facePairData.push_back(share);
  }
  else
  {
    iter->second.data.facePairData.push_back(share);
  }
}

void BifroestBase::CoordinateShares()
{
  //construct the map of those processes I am sending to ... also send out my size and info to those processes
  iArray1d sendMap(m_size, 0);
  {
    localIndex i = 0;
    for(std::map<int, DelegateProcessingToAnother>::iterator iter = this->toProcessData.begin();
        iter != this->toProcessData.end(); ++iter, ++i)
    {
      iter->second.SetSizesAndSend(m_rank);
      sendMap[iter->first] = 1;
    }
  }

  //get the expected number of requests to send me data
  int m_numberRequests = 0;
  {
    iArray1d receiveArray(m_size,0);
    MPI_Allreduce(sendMap.data(), receiveArray.data(), m_size,
                  MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
    m_numberRequests = receiveArray[m_rank];
  }

  //receive the size and info for each request to send me data
  Array1dT<iArray1d> requestedSizes(m_numberRequests);
  Array1dT<MPI_Request> requests(m_numberRequests);
  {
    for(int i = 0; i < m_numberRequests; i++)
    {
      requestedSizes[i].resize(NumberOfFields, 0);
      MPI_Irecv(requestedSizes[i].data(), NumberOfFields, MPI_INTEGER, MPI_ANY_SOURCE,
                MPI_ANY_TAG, MPI_COMM_WORLD, &requests[i]);
    }
  }

  //wait until all size/info packets have been received
  MPI_Barrier(MPI_COMM_WORLD);
  //std::cout << this->m_rank << ": done reducing requests" << std::endl;

  //allocate data to the requestors
  for(int i = 0; i < m_numberRequests; i++)
  {
    const int neighbor_rank = requestedSizes[i][SenderRank];
    fromProcessData[neighbor_rank].Allocate(m_rank, neighbor_rank, requestedSizes[i]);
  }
}

void BifroestBase::Synchronize()
{
  //send data to those I requested to send data to
  for(std::map<int, DelegateProcessingToAnother>::iterator iter = this->toProcessData.begin();
      iter != this->toProcessData.end(); ++iter)
  {
    iter->second.data.ResetReceiveFlag();
    iter->second.SendPackage();
    iter->second.ReceiveReturnedData();
  }

  //receive data from those that requested to send me data
  for(std::map<int, ProcessingDelegatedToMe>::iterator iter = this->fromProcessData.begin();
      iter != this->fromProcessData.end(); ++iter)
  {
    iter->second.data.ResetReceiveFlag();
    iter->second.ReceivePackage();
  }

  //process returned data and calculate data to return
  bool stillWaiting = true;
  while(stillWaiting)
  {
    stillWaiting = false;

    //send data to those I requested to send data to
    for(std::map<int, DelegateProcessingToAnother>::iterator iter = this->toProcessData.begin();
        iter != this->toProcessData.end(); ++iter)
    {
      const bool ok = iter->second.ReturnedDataReadyToProcess();
      if(ok && !iter->second.data.receiveFlag)
      {
        ProcessReturnedData(iter->second);
        iter->second.data.receiveFlag = true;
      }
      stillWaiting = stillWaiting || !iter->second.data.receiveFlag;
    }

    //receive data from those that requested to send me data
    for(std::map<int, ProcessingDelegatedToMe>::iterator iter = this->fromProcessData.begin();
        iter != this->fromProcessData.end(); ++iter)
    {
      const bool ok = iter->second.PackageReadyToProcess();
      if(ok && !iter->second.data.receiveFlag)
      {
        ProcessPackage(iter->second);
        iter->second.data.receiveFlag = true;
      }
      stillWaiting = stillWaiting || !iter->second.data.receiveFlag;
    }
  }

  //wait until all size/info packets have been received
  MPI_Barrier(MPI_COMM_WORLD);
}
