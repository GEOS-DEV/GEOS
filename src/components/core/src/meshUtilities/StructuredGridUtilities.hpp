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

#ifndef STRUCTURED_GRID_UTILITIES_H
#define STRUCTURED_GRID_UTILITIES_H

/**
 * @file StructuredGridUtilities.h
 * @author white230
 */

//#include "legacy/Common/Common.h"
#include <cassert>

namespace StructuredGrid
{
/*!
 *  Given n, compute n^d, where d is the spatial dimension.
 */

template <int dim>
int dimpower(int n);

template <> inline int dimpower<1>(int n) { return n; }
template <> inline int dimpower<2>(int n) { return n*n; }
template <> inline int dimpower<3>(int n) { return n*n*n; }

/*!
 * Given a lexographical index N, map back to the original
 * i,j,k indices of the point. The first variation here assumes
 * a uniform number of points nnx in all coordinate directions.
 */

template <int dim>
void map_index(const int index,
               const int nnx,
               std::vector<int> &indices);

template <>
inline
void map_index<1>(const int index,
                  const int nnx,
                  std::vector<int> &indices)
{
  assert(index < nnx);
  indices[0] = index;
}

template <>
inline
void map_index<2>(const int index,
                  const int nnx,
                  std::vector<int> &indices)
{
  assert(index < nnx*nnx);
  indices[0] = index % nnx;
  indices[1] = index / nnx;
}

template <>
inline
void map_index<3>(const int index,
                  const int nnx,
                  std::vector<int> &indices)
{
  assert(index < nnx*nnx*nnx);
  indices[0] = index % nnx;
  indices[1] = (index / nnx) % nnx;
  indices[2] = index / (nnx*nnx);
}

} // end namespace

#endif
