/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MultiVectorOperations.hpp
 *
 * Operations on blocks of vectors that run in the memory space of the linear algebra backend. A block of dot
 * products is computed by one kernel and returned to the host as a single small matrix, instead of one
 * reduction and one device synchronization for each pair of vectors. A linear combination of many vectors is
 * done in one fused kernel instead of one axpy for each vector. Only the small coefficient matrices cross the
 * host/device boundary.
 */

#ifndef GEOS_LINEARALGEBRA_UTILITIES_MULTIVECTOROPERATIONS_HPP_
#define GEOS_LINEARALGEBRA_UTILITIES_MULTIVECTOROPERATIONS_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "common/MpiWrapper.hpp"

#include <vector>

namespace geos
{

namespace multiVectorOperations
{

namespace internal
{

/// Number of rows handled sequentially by one thread in the block dot product kernel
constexpr localIndex dotChunkSize = 1024;

/// Number of vectors handled by one kernel. The views of these vectors are passed to the kernel by value, so that
/// no pointer table needs to be uploaded to the device.
constexpr localIndex tileSize = 16;

/// Views of up to tileSize vectors, passed by value to a kernel
struct ViewTile
{
  /// The views
  arrayView1d< real64 const > views[tileSize];
};

/**
 * @brief Collect the views of a tile of vectors.
 * @tparam VECTOR type of the vectors
 * @param[in] vectors the vectors
 * @param[in] begin index of the first vector of the tile
 * @param[in] count number of vectors of the tile, at most tileSize
 * @return the views
 */
template< typename VECTOR >
ViewTile collectTile( std::vector< VECTOR const * > const & vectors, localIndex const begin, localIndex const count )
{
  ViewTile tile;
  for( localIndex k = 0; k < count; ++k )
  {
    tile.views[k] = vectors[begin + k]->values();
  }
  return tile;
}

} // namespace internal

/**
 * @brief Convert pointers to non-const vectors to pointers to const vectors.
 * @tparam VECTOR type of the vectors
 * @param[in] vectors the pointers
 * @return the pointers to const vectors
 */
template< typename VECTOR >
std::vector< VECTOR const * > constPointers( std::vector< VECTOR * > const & vectors )
{
  return std::vector< VECTOR const * >( vectors.begin(), vectors.end() );
}

/**
 * @brief Compute the matrix of dot products between two blocks of vectors.
 * @tparam VECTOR type of the vectors
 * @param[in] X first block, nx vectors
 * @param[in] Y second block, ny vectors
 * @param[out] result the nx x ny matrix with entries X[i] . Y[j], on the host
 *
 * The products are computed on the device, one kernel for each tile of at most 16 x 16 pairs of vectors. They
 * are reduced over the rows in a second kernel, and copied to the host in a single transfer of nx * ny values.
 * One reduction over the processes completes them. Nothing is uploaded to the device.
 */
template< typename VECTOR >
void dots( std::vector< VECTOR const * > const & X,
           std::vector< VECTOR const * > const & Y,
           array2d< real64 > & result )
{
  localIndex const nx = LvArray::integerConversion< localIndex >( X.size() );
  localIndex const ny = LvArray::integerConversion< localIndex >( Y.size() );
  result.resize( nx, ny );
  if( nx == 0 || ny == 0 )
  {
    return;
  }

  localIndex const n = X[0]->localSize();
  localIndex const chunkSize = internal::dotChunkSize;
  localIndex const numChunks = LvArray::math::max( ( n + chunkSize - 1 ) / chunkSize, localIndex( 1 ) );
  localIndex constexpr tile = internal::tileSize;

  // The scratch arrays are allocated directly in the memory space of the kernels: allocating them on the host
  // would upload their content to the device at every call
  array1d< real64 > local;
  local.resizeWithoutInitializationOrDestruction( parallelDeviceMemorySpace, nx * ny );
  arrayView1d< real64 > const localView = local.toView();
  array1d< real64 > partial;
  partial.resizeWithoutInitializationOrDestruction( parallelDeviceMemorySpace, tile * tile * numChunks );
  arrayView1d< real64 > const partialView = partial.toView();
  arrayView1d< real64 const > const partialConst = partial.toViewConst();

  for( localIndex x0 = 0; x0 < nx; x0 += tile )
  {
    localIndex const tx = LvArray::math::min( tile, nx - x0 );
    internal::ViewTile const xs = internal::collectTile( X, x0, tx );
    for( localIndex y0 = 0; y0 < ny; y0 += tile )
    {
      localIndex const ty = LvArray::math::min( tile, ny - y0 );
      internal::ViewTile const ys = internal::collectTile( Y, y0, ty );

      forAll< parallelDevicePolicy<> >( tx * ty * numChunks, [=] GEOS_HOST_DEVICE ( localIndex const w )
      {
        localIndex const pair = w / numChunks;
        localIndex const chunk = w - pair * numChunks;
        localIndex const i = pair / ty;
        localIndex const j = pair - i * ty;
        localIndex const begin = chunk * chunkSize;
        localIndex const end = LvArray::math::min( n, begin + chunkSize );
        real64 sum = 0.0;
        for( localIndex r = begin; r < end; ++r )
        {
          sum += xs.views[i][r] * ys.views[j][r];
        }
        partialView[w] = sum;
      } );

      forAll< parallelDevicePolicy<> >( tx * ty, [=] GEOS_HOST_DEVICE ( localIndex const pair )
      {
        real64 sum = 0.0;
        for( localIndex c = 0; c < numChunks; ++c )
        {
          sum += partialConst[pair * numChunks + c];
        }
        localIndex const i = pair / ty;
        localIndex const j = pair - i * ty;
        localView[( x0 + i ) * ny + y0 + j] = sum;
      } );
    }
  }

  // One transfer of nx * ny values to the host, and one reduction over the processes
  array1d< real64 > global( nx * ny );
  local.move( hostMemorySpace, true );
  MpiWrapper::allReduce( local, global, MpiWrapper::Reduction::Sum, X[0]->comm() );
  for( localIndex i = 0; i < nx; ++i )
  {
    for( localIndex j = 0; j < ny; ++j )
    {
      result( i, j ) = global[i * ny + j];
    }
  }
}

/**
 * @brief Linear combination of a block of vectors: out = ( accumulate ? out : 0 ) + sum_j coefficients[j] V[j].
 * @tparam VECTOR type of the vectors
 * @param[in] V the vectors, which must not include @p out
 * @param[in] coefficients one coefficient for each vector
 * @param[in,out] out the result
 * @param[in] accumulate whether to add to the current content of @p out
 *
 * One fused kernel reads a tile of at most 16 vectors. The views and the coefficients are passed to the kernel
 * by value, so nothing is uploaded to the device.
 */
template< typename VECTOR >
void combine( std::vector< VECTOR const * > const & V,
              std::vector< real64 > const & coefficients,
              VECTOR & out,
              bool const accumulate )
{
  localIndex const nv = LvArray::integerConversion< localIndex >( V.size() );
  GEOS_ERROR_IF( coefficients.size() != V.size(), "multiVectorOperations::combine: one coefficient per vector is required" );
  if( nv == 0 )
  {
    if( !accumulate )
    {
      out.zero();
    }
    return;
  }

  localIndex const n = out.localSize();
  localIndex constexpr tile = internal::tileSize;
  arrayView1d< real64 > const o = out.open();
  for( localIndex v0 = 0; v0 < nv; v0 += tile )
  {
    localIndex const count = LvArray::math::min( tile, nv - v0 );
    internal::ViewTile const vs = internal::collectTile( V, v0, count );
    real64 weights[internal::tileSize] = {};
    for( localIndex k = 0; k < count; ++k )
    {
      weights[k] = coefficients[v0 + k];
    }
    bool const add = accumulate || v0 > 0;
    forAll< parallelDevicePolicy<> >( n, [=] GEOS_HOST_DEVICE ( localIndex const r )
    {
      real64 sum = add ? o[r] : 0.0;
      for( localIndex k = 0; k < count; ++k )
      {
        sum += weights[k] * vs.views[k][r];
      }
      o[r] = sum;
    } );
  }
  out.close();
}

} // namespace multiVectorOperations

} // namespace geos

#endif /* GEOS_LINEARALGEBRA_UTILITIES_MULTIVECTOROPERATIONS_HPP_ */
