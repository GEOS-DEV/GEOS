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
* @file CellVerlet.h
* @date Nov 11, 2011
* @author Scott Johnson
*
* @brief Cell-Verlet spatial sorting class
************************************************************/

#ifndef CELLVERLET_H_
#define CELLVERLET_H_

#include "SpatialSorterBase.h"

namespace SpatialSorting
{

class CellVerlet : public SpatialSorterBase
{

public:
  CellVerlet();
  virtual ~CellVerlet();

  static std::string SpatialSorterName() { return "CellVerlet"; }

  virtual bool Sort(const array<real64>& radii, const array<R1Tensor>& x,
                    array<lArray1d>& neighborList, array<lSet>& neighborListInverse,
                    const int* const excludeFromSorting = 0);

  virtual bool Update(const array<real64>& radii, const array<R1Tensor>& x, const lSet& toResort,
                      array<lArray1d>& neighborList, array<lSet>& neighborListInverse,
                      const int* const excludeFromSorting = 0);

  virtual void Clear();

private:

  void RemoveToCheck(const std::map<localIndex, array<integer> >& toCheckFurther);
  void Add(const array<real64>& radii, const array<R1Tensor>& x, const localIndex ixfc0, const array<integer>& index,
           array<lArray1d>& neighborList, array<lSet>& neighborListInverse,
           const int* const excludeFromSorting);

  bool UpdateMinMaxDimension(const array<real64>& radii,
                             const array<R1Tensor>& x,
                             const int* const excludeFromContact);

  bool UpdateMinMaxDimension(const array<real64>& radii,
                             const array<R1Tensor>& x,
                             const lSet& toResort);

  bool SortSub(const bool reset, const array<real64>& radii, const array<R1Tensor>& x,
               array<lArray1d>& neighborList, array<lSet>& neighborListInverse,
               const int* const excludeFromSorting);

  bool Remove(const localIndex, const array<integer>& guess);

  realT binFct;
  Array3dT< lArray1d > bins;
  Array3dT< array<lArray1d*> > neighborBins;
  R1Tensor xmin, xmax;
  realT dx;

  inline void BinIndices(const R1Tensor& x, array<integer>& bin)
  {
    R1Tensor xx(x);
    xx -= this->xmin;   //adjust by the minimum
    xx *= 1.0 / this->dx;
    for (unsigned int i = 0 ; i < nsdof ; i++)
      bin[i] = static_cast<int>(xx(i));
  };
};

}

#endif /* CELLVERLET_H_ */
