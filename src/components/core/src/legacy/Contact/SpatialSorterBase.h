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
* @file SpatialSorter.h
* @date Nov 11, 2011
* @author Scott Johnson
*
* @brief Abstract base spatial sorting class
************************************************************/

#ifndef SPATIALSORTER_H_
#define SPATIALSORTER_H_

#include "math/TensorT/R1TensorT.h"
#include "Common/Common.h"
#include "Common/intrinsic_typedefs.h"
#include "Common/typedefs.h"
#include "ArrayT/ArrayT.h"

namespace SpatialSorting
{

class SpatialSorterBase
{
public:
  SpatialSorterBase();
  virtual ~SpatialSorterBase();

  inline static bool Close(const realT r0, const R1Tensor& x0, const realT r1, const R1Tensor& x1)
  {
    R1Tensor dd(x0);
    dd -= x1;
    realT sumr2 = r0;
    sumr2 += r1;
    sumr2 *= sumr2;
    return Dot(dd,dd) < sumr2;
  }

  /**
   * @author Scott Johnson
   * Determines whether two spherical potentials interfere
   * @param[in] kf0 First index in the lists
   * @param[in] kf1 Second index in the lists
   * @param[in] radii List of radii of the potentials
   * @param[in] centers List of centers of the potentials
   * @return Whether the potentials interfere
   */
  inline static bool Close(const localIndex kf0, const localIndex kf1, const array<real64>& radii, const array<R1Tensor>& centers)
  {
    return Close(radii[kf0], centers[kf0], radii[kf1], centers[kf1]);
  }

  static void Remove(const lSet& toRemove,
                     array<lArray1d>& neighborList,
                     array<lSet>& neighborListInverse);

  static void RemoveDuplicates(array<lArray1d>& neighborList);

  static void insert(const localIndex i, const localIndex count,
                     array<lArray1d>& neighborList,
                     array<lSet>& neighborListInverse);

  static void Add(const localIndex kf0, const localIndex kf1,
                  array<lArray1d>& neighborList,
                  array<lSet>& neighborListInverse);

  static void AddIfClose(const localIndex kf0, const localIndex kf1,
                         const array<real64>& radii, const array<R1Tensor>& centers,
                         array<lArray1d>& neighborList,
                         array<lSet>& neighborListInverse);

  static void OffsetIndexing(const localIndex start, const localIndex offset,
                             array<lArray1d>& neighborList,
                             array<lSet>& neighborListInverse);


  virtual bool Sort(const array<real64>& radii,
                    const array<R1Tensor>& x,
                    array<lArray1d>& neighborList,
                    array<lSet>& neighborListInverse,
                    const int* const excludeFromSorting = 0) = 0;

  virtual bool Update(const array<real64>& radii,
                      const array<R1Tensor>& x,
                      const lSet& toResort,
                      array<lArray1d>& neighborList,
                      array<lSet>& neighborListInverse,
                      const int* const excludeFromSorting = 0) = 0;


  virtual void Clear();
};
}

#endif /* SPATIALSORTER_H_ */
