/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/*
 * LAIHelperFunctions.hpp
 */


#ifndef GEOSX_LINEARALGEBRA_UTILITIES_LAIHELPERFUNCTIONS_HPP_
#define GEOSX_LINEARALGEBRA_UTILITIES_LAIHELPERFUNCTIONS_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/ElementRegionManager.hpp"

namespace geosx
{
namespace LAIHelperFunctions
{

/**
 * @brief Create a permutation matrix for a given nodal variable.
 * @param[in]  nodeManager       the node manager
 * @param[in]  nRows             number of local rows in the matrix
 * @param[in]  nCols             number of local columns in the matrix
 * @param[in]  nDofPerNode       number of degrees-of-freedom per node
 * @param[in]  dofKey            DofManager key used to access dof index array
 * @param[out] permutationMatrix the target matrix
 */
void CreatePermutationMatrix( NodeManager const * const nodeManager,
                              localIndex const nRows,
                              localIndex const nCols,
                              int const nDofPerNode,
                              string const dofKey,
                              ParallelMatrix & permutationMatrix );

/**
 * @brief Create a permutation matrix for a given nodal variable.
 * @param[in]  elemManager       the element region manager
 * @param[in]  nRows             number of local rows in the matrix
 * @param[in]  nCols             number of local columns in the matrix
 * @param[in]  nDofPerNode       number of degrees-of-freedom per node
 * @param[in]  dofKey            DofManager key used to access dof index array
 * @param[out] permutationMatrix the target matrix
 */
void CreatePermutationMatrix( ElementRegionManager const * const elemManager,
                              localIndex const nRows,
                              localIndex const nCols,
                              int const nDofPerNode,
                              string const DofKey,
                              ParallelMatrix & permutationMatrix );

/**
 * @brief Permute a vector.
 * @param[in] vector            the source vector
 * @param[in] permutationMatrix the permutation matrix
 * @return the permuted vector
 */
ParallelVector PermuteVector( ParallelVector const & vector,
                              ParallelMatrix const & permutationMatrix );

/**
 * @brief Permute rows and columns of a square matrix.
 * @param[in] matrix            the source matrix
 * @param[in] permutationMatrix permutation matrix
 * @return the permuted matrix
 */
ParallelMatrix PermuteMatrix( ParallelMatrix const & matrix,
                              ParallelMatrix const & permutationMatrix );

/**
 * Permute rows and columns of a rectangular matrix
 * @param[in] matrix                 the source matrix
 * @param[in] permutationMatrixLeft  left permutation matrix
 * @param[in] permutationMatrixRight right permutation matrix
 * @return the permuted matrix
 */
ParallelMatrix PermuteMatrix( ParallelMatrix const & matrix,
                              ParallelMatrix const & permutationMatrixLeft,
                              ParallelMatrix const & permutationMatrixRight );

void PrintPermutedVector( ParallelVector const & vector,
                          ParallelMatrix const & permuationMatrix,
                          std::ostream & os );


void PrintPermutedMatrix( ParallelMatrix const & matrix,
                          ParallelMatrix const & permutationMatrix,
                          std::ostream & os );

void PrintPermutedMatrix( ParallelMatrix const & matrix,
                          ParallelMatrix const & permutationMatrixLeft,
                          ParallelMatrix const & permutationMatrixRight,
                          std::ostream & os );

/**
 * @brief Apply a separate component approximation (filter) to a matrix.
 * @tparam MATRIX the type of matrices
 * @param src         the source matrix
 * @param dst         the target (filtered) matrix
 * @param dofsPerNode number of degrees-of-freedom per node
 */
template< typename MATRIX >
void SeparateComponentFilter( MATRIX const & src,
                              MATRIX & dst,
                              const localIndex dofsPerNode )
{
  GEOSX_ERROR_IF( dofsPerNode < 2, "Function requires dofsPerNode > 1" );

  const localIndex localRows  = src.numLocalRows();
  const localIndex maxEntries = src.maxRowLength();
  const localIndex maxDstEntries = maxEntries / dofsPerNode;


#if 0
  dst.createWithLocalSize( localRows, maxEntries, MPI_COMM_WORLD );
  dst.open();
  array1d< real64 > srcValues;
  array1d< real64 > dstValues( maxDstEntries );

  array1d< globalIndex > srcIndices;
  array1d< globalIndex > dstIndices( maxDstEntries );

#if defined(OVERRIDE_CREATE)
#if defined(GEOSX_USE_CUDA)
  srcIndices.move( LvArray::MemorySpace::GPU, false );
  srcValues.move( LvArray::MemorySpace::GPU, false );
#endif
#endif

  for( globalIndex row=src.ilower(); row<src.iupper(); ++row )
  {
    const globalIndex rowComponent = row % dofsPerNode;
    const localIndex rowLength = src.globalRowLength( row );
    srcIndices.resize( rowLength );
    srcValues.resize( rowLength );

    src.getRowCopy( row, srcIndices, srcValues );
    globalIndex const * indicesPtr = srcIndices.data();
    real64 const * valuesPtr = srcValues.data();

//    forAll< parallelDevicePolicy<> >( 1, [=] GEOSX_DEVICE ( localIndex const )
//    {
//      printf( "srcIndices = { " ); for( localIndex i=0 ; i<rowLength ; ++i ) { printf( "%4d, ",indicesPtr[i] ); }  printf( " }\n" );
//      printf( "srcValues  = { " ); for( localIndex i=0 ; i<rowLength ; ++i ) { printf( "%4.1g, ",valuesPtr[i] ); }  printf( " }\n" );
//    });

#if defined(OVERRIDE_CREATE)
#if defined(GEOSX_USE_CUDA)
  srcIndices.move( LvArray::MemorySpace::CPU, false );
  srcValues.move( LvArray::MemorySpace::CPU, false );
  dstIndices.move( LvArray::MemorySpace::CPU, false );
  dstValues.move( LvArray::MemorySpace::CPU, false );
#endif
#endif

//  printf( "srcIndices = { " ); for( localIndex i=0 ; i<rowLength ; ++i ) { printf( "%4d, ",srcIndices[i] ); }  printf( " }\n" );
//  printf( "srcValues  = { " ); for( localIndex i=0 ; i<rowLength ; ++i ) { printf( "%4.1g, ",srcValues[i] ); }  printf( " }\n" );

    localIndex k=0;
    for( localIndex col=0; col<rowLength; ++col )
    {
      const globalIndex colComponent = srcIndices[col] % dofsPerNode;
      if( rowComponent == colComponent )
      {
        dstValues[k] = srcValues[col];
        dstIndices[k] = srcIndices[col];
        ++k;
      }
    }
#if defined(OVERRIDE_CREATE)
#if defined(GEOSX_USE_CUDA)
  dstIndices.move( LvArray::MemorySpace::GPU, false );
  dstValues.move( LvArray::MemorySpace::GPU, false );
#endif
#endif
  indicesPtr = dstIndices.data();
  valuesPtr = dstValues.data();

//  printf( "k = %d \n", k );
//  forAll< parallelDevicePolicy<> >( 1, [=] GEOSX_DEVICE ( localIndex const )
//  {
//    printf( "dstIndices = { " ); for( localIndex i=0 ; i<k ; ++i ) { printf( "%6d, ",indicesPtr[i] ); }  printf( " }\n" );
//    printf( "dstValues  = { " ); for( localIndex i=0 ; i<k ; ++i ) { printf( "%6.1g, ",valuesPtr[i] ); }  printf( " }\n" );
//  });

    dst.insert( row, dstIndices.data(), dstValues.data(), k );
  }
  dst.close();

#else

  array2d< globalIndex > srcIndices( localRows, maxEntries);;
  array2d< real64 > srcValues( localRows, maxEntries);

  CRSMatrix< real64 > tempMat;
  tempMat.resize( localRows, src.numLocalCols(), maxDstEntries );

  for( globalIndex row=src.ilower(); row<src.iupper(); ++row )
  {
    const globalIndex rowComponent = row % dofsPerNode;
    const localIndex rowLength = src.globalRowLength( row );

    src.getRowCopy( row, srcIndices[row], srcValues[row] );

    for( localIndex col=0; col<rowLength; ++col )
    {
      const globalIndex colComponent = srcIndices(row,col) % dofsPerNode;
      if( rowComponent == colComponent )
      {
        tempMat.insertNonZero( row, col, srcValues(row,col) );
      }
    }
  }

  dst.create( tempMat.toViewConst(), MPI_COMM_GEOSX );

#endif
}

} // LAIHelperFunctions namespace

} // geosx namespace

#endif /*GEOSX_LINEARALGEBRA_UTILITIES_LAIHELPERFUNCTIONS_HPP_*/
