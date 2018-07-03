/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifndef RCVSPARSE_H_
#define RCVSPARSE_H_

/**
 * @file RCVSparse.h
 * @author walsh24
 */

#include <string>
#include "Common/Common.h"

/// Struct for simple row,column,value sparse matrix representation
struct rcv
{
  rcv(localIndex rr, localIndex cc, realT vv): r(rr),c(cc),v(vv){};
  rcv(): r(0),c(0),v(0){};
  localIndex r;
  localIndex c;
  realT v;
};

bool rowColumnOrder(const rcv& A, const rcv& B);
bool columnRowOrder(const rcv& A, const rcv& B);

void ConsolidateSparseArray(array<rcv>& K);

void WriteSparseArrayToFile(std::string filename, const array<rcv>& K);

#endif
