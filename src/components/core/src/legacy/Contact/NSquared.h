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
* @file NSquared.h
* @date Nov 11, 2011
* @author Scott Johnson
*
* @brief Order N^2 spatial sorting class
************************************************************/

#ifndef NSQUARED_H_
#define NSQUARED_H_

#include "SpatialSorterBase.h"

namespace SpatialSorting
{

class NSquared : public SpatialSorterBase
{
public:
  NSquared();
  virtual ~NSquared();

  static std::string SpatialSorterName() { return "N2"; }

  virtual bool Sort(const array<real64>& radii,
                    const array<R1Tensor>& x,
                    array<lArray1d>& neighborList,
                    array<lSet>& neighborListInverse,
                    const int* const excludeFromSorting = 0);

  virtual bool Update(const array<real64>& radii,
                      const array<R1Tensor>& x,
                      const lSet& toResort,
                      array<lArray1d>& neighborList,
                      array<lSet>& neighborListInverse,
                      const int* const excludeFromSorting = 0);
};
}

#endif /* NSQUARED_H_ */
