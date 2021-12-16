/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HypreKernels.cpp
 */

#include "HypreKernels.hpp"

namespace geosx
{
namespace hypre
{

HYPRE_BigInt const * getOffdColumnMap( hypre_ParCSRMatrix const * const mat )
{
  return ( hypre_CSRMatrixMemoryLocation( hypre_ParCSRMatrixDiag( mat ) ) == HYPRE_MEMORY_HOST )
         ? hypre_ParCSRMatrixColMapOffd( mat )
         : hypre_ParCSRMatrixDeviceColMapOffd( mat );
}

void scaleMatrixValues( hypre_CSRMatrix * const mat,
                        real64 const factor )
{
  GEOSX_ASSERT( mat != nullptr );
  if( hypre_CSRMatrixNumCols( mat ) == 0 )
  {
    return;
  }
  if( isEqual( factor, 1.0 ) )
  {
    return;
  }
  HYPRE_Real * const va = hypre_CSRMatrixData( mat );
  forAll< execPolicy >( hypre_CSRMatrixNumNonzeros( mat ), [=] GEOSX_HYPRE_DEVICE ( HYPRE_Int const i )
  {
    va[i] *= factor;
  } );
}

void scaleMatrixRows( hypre_CSRMatrix * const mat,
                      hypre_Vector const * const vec )
{
  GEOSX_ASSERT( mat != nullptr );
  if( hypre_CSRMatrixNumCols( mat ) == 0 )
  {
    return;
  }

  CSRData< false > const csr{ mat };
  HYPRE_Real const * const scalingFactors = hypre_VectorData( vec );

  forAll< execPolicy >( csr.nrow, [=] GEOSX_HYPRE_DEVICE ( HYPRE_Int const localRow )
  {
    real64 const factor = scalingFactors[localRow];
    if( !isEqual( factor, 1.0 ) )
    {
      for( HYPRE_Int j = csr.rowptr[localRow]; j < csr.rowptr[localRow + 1]; ++j )
      {
        csr.values[j] *= factor;
      }
    }
  } );
}

void clampMatrixEntries( hypre_CSRMatrix * const mat,
                         real64 const lo,
                         real64 const hi,
                         bool const skip_diag )
{
  GEOSX_ASSERT( mat != nullptr );
  if( hypre_CSRMatrixNumCols( mat ) == 0 )
  {
    return;
  }

  CSRData< false > const csr{ mat };
  forAll< execPolicy >( csr.nrow, [=] GEOSX_HYPRE_DEVICE ( HYPRE_Int const localRow )
  {
    // Hypre stores diagonal element at the beginning of each row, we assume it's always present
    for( HYPRE_Int k = csr.rowptr[localRow] + skip_diag; k < csr.rowptr[localRow+1]; ++k )
    {
      csr.values[k] = LvArray::math::min( hi, LvArray::math::max( lo, csr.values[k] ) );
    }
  } );
}

real64 computeMaxNorm( hypre_CSRMatrix const * const mat )
{
  GEOSX_ASSERT( mat != nullptr );
  if( hypre_CSRMatrixNumCols( mat ) == 0 )
  {
    return 0;
  }

  HYPRE_Real const * const va = hypre_CSRMatrixData( mat );
  RAJA::ReduceMax< ReducePolicy< execPolicy >, real64 > maxAbsElement( 0.0 );
  forAll< execPolicy >( hypre_CSRMatrixNumNonzeros( mat ), [=] GEOSX_HYPRE_DEVICE ( HYPRE_Int const k )
  {
    maxAbsElement.max( LvArray::math::abs( va[k] ) );
  } );
  return maxAbsElement.get();
}

real64 computeMaxNorm( hypre_CSRMatrix const * const mat,
                       arrayView1d< globalIndex const > const & rowIndices,
                       globalIndex const firstLocalRow )
{
  GEOSX_ASSERT( mat != nullptr );
  if( hypre_CSRMatrixNumCols( mat ) == 0 )
  {
    return 0;
  }

  CSRData< true > const csr{ mat };
  HYPRE_Int const numRows = hypre_CSRMatrixNumRows( mat );
  GEOSX_DEBUG_VAR( numRows );

  RAJA::ReduceMax< ReducePolicy< execPolicy >, real64 > maxAbsElement( 0.0 );
  forAll< execPolicy >( rowIndices.size(), [=] GEOSX_HYPRE_DEVICE ( localIndex const i )
  {
    localIndex const localRow = rowIndices[i] - firstLocalRow;
    GEOSX_ASSERT( 0 <= localRow && localRow < numRows );
    for( HYPRE_Int j = csr.rowptr[localRow]; j < csr.rowptr[localRow + 1]; ++j )
    {
      maxAbsElement.max( LvArray::math::abs( csr.values[j] ) );
    }
  } );
  return maxAbsElement.get();
}

} // namespace hypre
} // namespace geosx
