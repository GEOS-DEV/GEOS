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

/**
 * @file SuperluUtils.cpp
 */

#include "SuperluUtils.hpp"
#include "common/Stopwatch.hpp"

namespace geosx
{

void ConvertToSuperMatrix( HypreMatrix const & matrix,
                           array1d< globalIndex > & rowPtr,
                           array1d< globalIndex > & cols,
                           array1d< real64 > & vals,
                           SuperMatrix & SLUDMat )
{
  localIndex const numLocalRows = matrix.numLocalRows();
  rowPtr.resize( numLocalRows+1 );
  cols.reserve( matrix.numLocalNonzeros() );
  vals.reserve( matrix.numLocalNonzeros() );
  rowPtr[0] = 0;
  for( localIndex i = 0; i < numLocalRows; ++i )
  {
    localIndex const nonZeros = matrix.localRowLength( i );
    array1d< globalIndex > colIndices( nonZeros );
    array1d< real64 > values( nonZeros );
    matrix.getRowCopy( matrix.getGlobalRowID( i ), colIndices, values );
    cols.insert( rowPtr[i], colIndices.begin(), colIndices.end() );
    vals.insert( rowPtr[i], values.begin(), values.end() );
    rowPtr[i+1] = rowPtr[i] + nonZeros;
  }

  dCreate_CompRowLoc_Matrix_dist( &SLUDMat,
                                  toSuperlu_intT( matrix.numGlobalRows() ),
                                  toSuperlu_intT( matrix.numGlobalRows() ),
                                  toSuperlu_intT( matrix.numLocalNonzeros() ),
                                  toSuperlu_intT( numLocalRows ),
                                  toSuperlu_intT( matrix.ilower() ),
                                  vals.data(),
                                  toSuperlu_intT( cols.data() ),
                                  toSuperlu_intT( rowPtr.data() ),
                                  SLU_NR_loc,
                                  SLU_D,
                                  SLU_GE );
}

int SolveSuperMatrix( SuperMatrix & SLUDMat,
                      HypreVector const & b,
                      HypreVector & x,
                      MPI_Comm const & comm,
                      superlu_dist_options_t & options,
                      integer const & logLevel,
                      real64 & timeFact,
                      real64 & timeSolve )
{
  Stopwatch watch;

  int_t const m = SLUDMat.nrow;
  int_t const n = SLUDMat.ncol;

  // Initialize ScalePermstruct.
  dScalePermstruct_t ScalePermstruct;
  dScalePermstructInit( m, n, &ScalePermstruct );

  // Initialize LUstruct.
  dLUstruct_t LUstruct;
  dLUstructInit( n, &LUstruct );

  // Initialize the statistics variables.
  SuperLUStat_t stat;
  PStatInit( &stat );

  // Create process grid.
  int const num_procs = MpiWrapper::Comm_size( comm );
  int pcols = 1;
  int prows = 1;
  while( prows*pcols <= num_procs )
  {
    ++prows;
  }
  --prows;
  pcols = num_procs/prows;
  while( prows*pcols != num_procs )
  {
    prows -= 1;
    pcols = num_procs/prows;
  }
  gridinfo_t grid;
  superlu_gridinit( comm, prows, pcols, &grid );

  // Call the linear equation solver.
  int nrhs = 0;
  int const ldb = b.localSize();
  dSOLVEstruct_t SOLVEstruct;
  int info = 0;

  options.Fact = DOFACT;
  pdgssvx( &options,
           &SLUDMat,
           &ScalePermstruct,
           NULL,
           ldb,
           nrhs,
           &grid,
           &LUstruct,
           &SOLVEstruct,
           NULL,
           &stat,
           &info );

  timeFact = watch.elapsedTime();
  watch.zero();

  if( info == 0 )
  {
    options.Fact = FACTORED;
    nrhs = 1;
    array1d< real64 > berr( nrhs );
    x.copy( b );
    pdgssvx( &options,
             &SLUDMat,
             &ScalePermstruct,
             x.extractLocalVector(),
             ldb,
             nrhs,
             &grid,
             &LUstruct,
             &SOLVEstruct,
             berr.data(),
             &stat,
             &info );
  }

  timeSolve = watch.elapsedTime();

  if( logLevel > 0 )
  {
    // Print the statistics.
    PStatPrint( &options, &stat, &grid );
  }

  // Deallocate other SuperLU data structures
  dScalePermstructFree( &ScalePermstruct );
  dDestroy_LU( n, &grid, &LUstruct );
  dLUstructFree( &LUstruct );
  PStatFree( &stat );
  superlu_gridexit( &grid );
  if( options.SolveInitialized )
  {
    dSolveFinalize( &options, &SOLVEstruct );
  }

  return info;
}

void DestroySuperMatrix( SuperMatrix & SLUDMat )
{
  SUPERLU_FREE( SLUDMat.Store );
}

colperm_t const & getColPermType( string const & value )
{
  static std::map< string, colperm_t > const optionMap =
  {
    { "none", NATURAL },
    { "MMD_At+A", MMD_AT_PLUS_A },
    { "MMD_AtA", MMD_ATA },
    { "metis", METIS_AT_PLUS_A },
    { "parmetis", PARMETIS },
  };

  GEOSX_LAI_ASSERT_MSG( optionMap.count( value ) > 0, "Unsupported SuperLU_Dist columns permutation option: " << value );
  return optionMap.at( value );
}

rowperm_t const & getRowPermType( string const & value )
{
  static std::map< string, rowperm_t > const optionMap =
  {
    { "none", NOROWPERM },
    { "mc64", LargeDiag_MC64 },
    { "awpm", LargeDiag_AWPM },
  };

  GEOSX_LAI_ASSERT_MSG( optionMap.count( value ) > 0, "Unsupported SuperLU_Dist rows permutation option: " << value );
  return optionMap.at( value );
}

}
