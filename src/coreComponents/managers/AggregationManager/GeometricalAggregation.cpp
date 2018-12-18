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
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file GeometricalAggregation.cpp
 */


#include "GeometricalAggregation.hpp"
#include <metis.h>

namespace geosx
{
void Aggregate(SparsityPattern const & fineSparsityPattern,SparsityPattern & coarseSparsityPattern, int totalNumberOfAggregates)
{

  /// STEP 1 : Conversion of fineSparsityPattern into a datastructure compatible with METIS
  
  // Four next variables stands for the METIS data structure to represent the connectivity graph of the fine scale
  idx_t nbOfFineCells;
  idx_t ncon;
  idx_t * xadj = nullptr; // (for the draft, it will be written using fineSparsityPattern)
  idx_t * adjncy = nullptr; // (for the draft, it will be written using fineSparsityPattern)

  idx_t * vwgt = nullptr;
  idx_t * vsize = nullptr;
  idx_t * adjwgt = nullptr;

  idx_t nparts = totalNumberOfAggregates;

  real_t * tpwgts = nullptr;
  real_t * ubvec = nullptr;
  idx_t * options = nullptr; // Maybe there will be some option to set.

  idx_t objval; 
  idx_t * part = nullptr;  //This vector will contain the map between the new coarse cells and the fine cells
  
  /// STEP 2 : Call to METIS
  
  if( totalNumberOfAggregates < 8 )
  {
    METIS_PartGraphRecursive(&nbOfFineCells, &ncon, xadj, adjncy, vwgt, vsize, adjwgt, &nparts, tpwgts, ubvec, options, &objval, part);
  }
  else
  {
    METIS_PartGraphKway(&nbOfFineCells, &ncon, xadj, adjncy, vwgt, vsize, adjwgt, &nparts, tpwgts, ubvec, options, &objval, part);
  }

  /// STEP 3 : Write the new Sparsity Pattern
  
}
}
