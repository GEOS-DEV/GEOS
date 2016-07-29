/*
 * TestBifroest.cpp
 *
 *  Created on: Aug 29, 2012
 *      Author: johnson346
 */

#include "TestBifroest.h"

TestBifroest::TestBifroest()
{
  // TODO Auto-generated constructor stub

}

TestBifroest::~TestBifroest()
{
  // TODO Auto-generated destructor stub
}

void TestBifroest::ProcessPackage(ProcessingDelegatedToMe& current)
{
  std::cout << this->m_rank << ": received package data from " << current.data.neighborRank << " of size " << current.data.sizes[NumberOfNodes] << std::endl;
  current.SendReturnedData();
}

void TestBifroest::ProcessReturnedData(DelegateProcessingToAnother& current)
{
  std::cout << this->m_rank << ": received return data from " << current.data.neighborRank << " of size " << current.data.sizes[NumberOfNodes] << std::endl;
}


