#include "TestBifroest.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#if GPAC_MPI
#include <mpi.h>
#endif


void fill(const int rank, const int size, Array1dT<R1Tensor>& tensors)
{
  const realT xdim = 1.0 / size;
  const realT x0 = xdim * rank;
  for(Array1dT<R1Tensor>::size_type i = 0; i < tensors.size(); i++)
  {
    tensors[i](0) = x0 + xdim * rand()/(RAND_MAX + 1.0);
    tensors[i](1) = rand()/(RAND_MAX + 1.0);
    tensors[i](2) = rand()/(RAND_MAX + 1.0);
  }
}

void set_shares(const Array1dT<R1Tensor>& tensors, TestBifroest& share)
{
  const realT ydim = 1.0 / share.Size();
  TempBifroestNodeSendData s;
  for(Array1dT<R1Tensor>::size_type i = 0; i < tensors.size(); i++)
  {
    int bifroestRank = static_cast<int>(tensors[i](1) / ydim);
    if(bifroestRank >= share.Size())
      bifroestRank = share.Size() - 1;
    if(bifroestRank == share.Rank())
      continue;
    s.nodeIndex = i + tensors.size() * share.Rank();
    s.x = tensors[i];
    share.AddShare(bifroestRank, s);
  }
}

int main(int argc, char** argv)
{
#if GPAC_MPI
  MPI_Init(&argc, &argv);
#endif
  //Read parameters from the first file
  if(argc != 2)
  {
    std::cout << "usage: <number of points per process>\n" << std::endl;
    throw GPException("usage: <number of points per process>\n");
  }

  //Instantiate process share class
  TestBifroest share;
  share.Initialize();

  //Set tensors
  Array1dT<R1Tensor> tensors(atoi(argv[1]));
  fill(share.Rank(), share.Size(), tensors);
  set_shares(tensors, share);

  //Determine/distribute the bifroest communication topology
  share.CoordinateShares();

  //Communicate data, do something with it, return
  share.Synchronize();

#if GPAC_MPI
  MPI_Finalize();
#endif
  return 0;
}
