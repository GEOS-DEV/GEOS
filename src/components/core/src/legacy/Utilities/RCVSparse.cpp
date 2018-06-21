// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.

/// Very simple sparse matrix data structure

#include <iostream>
#include <fstream>
#include "RCVSparse.h"

// sort functions
inline bool rowColumnOrder(const rcv& A, const rcv& B){  return (A.r < B.r)|| ( (A.r==B.r) && (A.c < B.c) );}
inline bool columnRowOrder(const rcv& A, const rcv& B){  return (A.c < B.c)|| ( (A.c==B.c) && (A.r < B.r) );}

/// Sum multiple entries with same row and column in K
/// and reduce the vector size
///
/// K = the sparse array stored in row, column value format
void ConsolidateSparseArray(array<rcv>& K){
  sort (K.begin(), K.end(), rowColumnOrder);
  localIndex i = 0;
  localIndex ii = 1;

  while(ii < K.size())
  {
    if(  K[i].c == K[ii].c && K[i].r == K[ii].r )
    {
      K[i].v +=K[ii].v;
    }
    else
    {
      ++i;
      if(i != ii)
      {
        K[i] = K[ii];
      }
    }
    ++ii;
  }
  K.resize(i+1);
}

void WriteSparseArrayToFile(std::string filename, const array<rcv>& K){
  std::ofstream outputfile(filename.c_str());
  for(unsigned int i =0 ; i < K.size() ; ++i)
  {
    outputfile <<  K[i].r << " " <<  K[i].c << " "  <<  K[i].v << "\n";
  }
  outputfile.close();
}
