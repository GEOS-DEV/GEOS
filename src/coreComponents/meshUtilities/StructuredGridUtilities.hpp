/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef STRUCTURED_GRID_UTILITIES_H
#define STRUCTURED_GRID_UTILITIES_H

/**
 * @file StructuredGridUtilities.hpp
 */

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
  GEOSX_ASSERT_GT(nnx, index);
  GEOSX_DEBUG_VAR(nnx);
  indices[0] = index;
}

template <>
inline
void map_index<2>(const int index,
                  const int nnx,
                  std::vector<int> &indices)
{
  GEOSX_ASSERT_GT(nnx*nnx, index);
  indices[0] = index % nnx;
  indices[1] = index / nnx;
}

template <>
inline
void map_index<3>(const int index,
                  const int nnx,
                  std::vector<int> &indices)
{
  GEOSX_ASSERT_GT(nnx*nnx*nnx, index);
  indices[0] = index % nnx;
  indices[1] = (index / nnx) % nnx;
  indices[2] = index / (nnx*nnx);
}

} // end namespace

#endif
