/*
 * DistributedBifroest.cpp
 *
 *  Created on: Aug 15, 2012
 *      Author: johnson346
 */

#include "DistributedBifroest.h"

DistributedBifroest::DistributedBifroest()
{
  // TODO Auto-generated constructor stub

}

DistributedBifroest::~DistributedBifroest()
{
  // TODO Auto-generated destructor stub
}

int DistributedBifroest::Status(const int numberOfContacts) const
{
  return (numberOfContacts < this->floorMean) ?
      (numberOfContacts - this->floorMean)  : ((numberOfContacts > this->ceilingMean) ?
          (numberOfContacts - this->ceilingMax) : 0);
}

//void DistributedBifroest::CoordinateShares()
//{
//  int rank, size;
//  {
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//    allCurrentContacts.resize(size);
//    //Send this->m_currentContacts to all processes and get theirs in allCurrentContacts
//  }
//
//  int contactSum = 0;
//  for(int i = 0; i < size; ++i)
//    contactSum += allCurrentContacts[i];
//  this->floorMean = contactSum / size;
//  this->ceilingMean = this->floorMean + (size*this->floorMean == contactSum) ? 0 : 1;
//
//  this->supplies.reserve(size);
//  this->demands.reserve(size);
//  for(int i = 0; i < size; ++i)
//  {
//    int status = Status(allCurrentContacts[i]);
//    if(i == rank || status == 0)
//      continue;
//    std::vector<RankCurrentContactsTuple>& v = status < 0 ? this->demands : this->supplies;
//    v.resize(v.size()+1);
//    v.back().status = status;
//    v.back().rank = i;
//  }
//  this->status = Status(allCurrentContacts[rank]);
//}

void DistributedBifroest::CoordinateShares()
{
  //want to match demands with supplies such that as few connections as possible take place
  //std::sort(this->demands.begin(), this->demands.end(), DistributedBifroest::CompareDecreasing());
  //std::sort(this->supplies.begin(), this->supplies.end(), DistributedBifroest::CompareDecreasing());
  //do marching matching and return the max connections per node as well as the total connections
  //fixme add marching match

  //std::reverse(this->supplies.begin(), this->supplies.end());
  //do marching matching and return the max connections per node as well as the total connections
  //fixme add marching match

  //do clustering then marching match and return the max connections per node as well as the total connections
  //fixme add clustering then marching match

  //KEEP THE BEST PERFORMER!!
}
