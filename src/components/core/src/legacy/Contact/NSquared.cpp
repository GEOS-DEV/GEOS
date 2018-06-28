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

/************************************************************
* @file NSquared.cpp
* @date Nov 11, 2011
* @author Scott Johnson
*
* @brief Order N^2 spatial sorting class
************************************************************/

#include "NSquared.h"

#include "SpatialSorterFactory.h"

namespace SpatialSorting
{

NSquared::NSquared()
{
  // TODO Auto-generated constructor stub

}

NSquared::~NSquared()
{
  // TODO Auto-generated destructor stub
}

bool NSquared::Update(const array<real64>& radii,
                      const array<R1Tensor>& x,
                      const lSet& toResort,
                      array<lArray1d>& neighborList,
                      array<lSet>& neighborListInverse,
                      const int* const excludeFromSorting)
{
  return Sort(radii, x, neighborList, neighborListInverse, excludeFromSorting);
}

/**
 * @brief Brief
 * @author Scott Johnson
 * Description
 * @param[in] radii Radii vector
 * @param[in] centers Centers vector
 * @param[out] neighborList Neighbor
 */
bool NSquared::Sort(const array<real64>& radii,
                    const array<R1Tensor>& centers,
                    array<lArray1d>& neighborList,
                    array<lSet>& neighborListInverse,
                    const int* const excludeFromContact)
{
  localIndex num = radii.size();

  //here, we at least limit by bounding sphere interaction (2)
  R1Tensor dd(0.);
  for (localIndex kf0 = 0 ; kf0 < num ; ++kf0)
  {
    neighborList[kf0].clear();
    neighborListInverse[kf0].clear();
  }

  //iterate through first
  if(excludeFromContact != 0)
  {
    for (localIndex kf0 = 0 ; kf0 < num ; ++kf0)
    {
      if(excludeFromContact[kf0]>0)
        continue;

      //iterate through second
      for (localIndex kf1 = kf0 + 1 ; kf1 < num ; ++kf1)
      {
        if(excludeFromContact[kf1]>0)
          continue;
        SpatialSorterBase::AddIfClose(kf0, kf1, radii, centers, neighborList, neighborListInverse);
      }
    }
  }
  else
  {
    for (localIndex kf0 = 0 ; kf0 < num ; ++kf0)
      for (localIndex kf1 = kf0 + 1 ; kf1 < num ; ++kf1)
        SpatialSorterBase::AddIfClose(kf0, kf1, radii, centers, neighborList, neighborListInverse);
  }
  return true;
}  //Sort
}

/// Register spatial sorter in the spatial sorter factory
REGISTER_SPATIALSORTER( NSquared )
