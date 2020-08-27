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
 * @file HypreSuperlu.cpp
 */

#include "HypreSuperlu.hpp"
#include "common/Stopwatch.hpp"

#include "HYPRE.h"
#include "_hypre_parcsr_mv.h"

namespace geosx
{

// Check matching requirements on index/value types between GEOSX and SuperLU_Dist

static_assert( sizeof( int_t ) == sizeof( globalIndex ),
               "SuperLU_Dist int_t and geosx::globalIndex must have the same size" );

static_assert( std::is_signed< int_t >::value == std::is_signed< globalIndex >::value,
               "SuperLU_Dist int_t and geoex::globalIndex must both be signed or unsigned" );

static_assert( std::is_same< double, real64 >::value,
               "SuperLU_Dist real and geosx::real64 must be the same type" );

void ConvertToSuperMatrix( HypreMatrix const & matrix,
                           SuperLU_DistData & SLUDData )
{
  // Merge diag and offd into one matrix (global ids)
  SLUDData.localStrip = hypre_MergeDiagAndOffd( matrix.unwrapped() );

  HYPRE_Int const * const hypreI = hypre_CSRMatrixI( SLUDData.localStrip );
  SLUDData.rowPtr = new int_t[matrix.numLocalRows()+1];
  for( localIndex i = 0; i <= matrix.numLocalRows(); ++i )
  {
    SLUDData.rowPtr[i] = LvArray::integerConversion< int_t >( hypreI[i] );
  }

  dCreate_CompRowLoc_Matrix_dist( &SLUDData.mat,
                                  toSuperlu_intT( matrix.numGlobalRows() ),
                                  toSuperlu_intT( matrix.numGlobalRows() ),
                                  toSuperlu_intT( matrix.numLocalNonzeros() ),
                                  toSuperlu_intT( matrix.numLocalRows() ),
                                  toSuperlu_intT( matrix.ilower() ),
                                  hypre_CSRMatrixData( SLUDData.localStrip ),
                                  toSuperlu_intT( hypre_CSRMatrixBigJ( SLUDData.localStrip ) ),
                                  SLUDData.rowPtr,
                                  SLU_NR_loc,
                                  SLU_D,
                                  SLU_GE );
}

void SuperLU_DistCreate( HypreMatrix const & matrix,
                         LinearSolverParameters const & params,
                         SuperLU_DistData & SLUDData )
{
  // Initialize options.
  set_default_options_dist( &SLUDData.options );
  if( params.logLevel > 1 )
  {
    SLUDData.options.PrintStat = YES;
  }
  else
  {
    SLUDData.options.PrintStat = NO;
  }

  if( params.direct.equilibrate )
  {
    SLUDData.options.Equil = YES;
  }
  else
  {
    SLUDData.options.Equil = NO;
  }
  SLUDData.options.ColPerm = HypreGetColPermType( params.direct.colPerm );
  SLUDData.options.RowPerm = HypreGetRowPermType( params.direct.rowPerm );
  if( params.direct.replaceTinyPivot )
  {
    SLUDData.options.ReplaceTinyPivot = YES;
  }
  else
  {
    SLUDData.options.ReplaceTinyPivot = NO;
  }
  if( params.direct.iterativeRefine )
  {
    SLUDData.options.IterRefine = SLU_DOUBLE;
  }
  else
  {
    SLUDData.options.IterRefine = NOREFINE;
  }

  if( params.logLevel > 0 )
  {
    print_sp_ienv_dist( &SLUDData.options );
    print_options_dist( &SLUDData.options );
  }

  // Convert matrix from Hypre to SuperLU_Dist format
  ConvertToSuperMatrix( matrix, SLUDData );

  // Save communicator
  SLUDData.comm = matrix.getComm();
}

int SuperLU_DistSetup( SuperLU_DistData & SLUDData,
                       real64 & time )
{
  Stopwatch watch;

  int_t const m = SLUDData.mat.nrow;
  int_t const n = SLUDData.mat.ncol;

  // Initialize ScalePermstruct.
  dScalePermstructInit( m, n, &SLUDData.ScalePermstruct );

  // Initialize LUstruct.
  dLUstructInit( n, &SLUDData.LUstruct );

  // Initialize the statistics variables.
  PStatInit( &SLUDData.stat );

  // Create process grid.
  int const num_procs = MpiWrapper::Comm_size( SLUDData.comm );
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
  superlu_gridinit( SLUDData.comm, prows, pcols, &SLUDData.grid );

  // Call the linear equation solver to factorize the matrix.
  int const nrhs = 0;
  int info = 0;

  SLUDData.options.Fact = DOFACT;
  pdgssvx( &SLUDData.options,
           &SLUDData.mat,
           &SLUDData.ScalePermstruct,
           NULL,
           n,
           nrhs,
           &SLUDData.grid,
           &SLUDData.LUstruct,
           &SLUDData.SOLVEstruct,
           NULL,
           &SLUDData.stat,
           &info );

  time = watch.elapsedTime();

  if( SLUDData.options.PrintStat == YES )
  {
    // Print the statistics.
    PStatPrint( &SLUDData.options, &SLUDData.stat, &SLUDData.grid );
  }

  return info;
}

int SuperLU_DistSolve( SuperLU_DistData & SLUDData,
                       HypreVector const & b,
                       HypreVector & x,
                       real64 & time )
{
  Stopwatch watch;

  x.copy( b );

  // Call the linear equation solver to solve the matrix.
  int const nrhs = 1;
  int const ldb = b.localSize();
  array1d< real64 > berr( nrhs );
  int info = 0;

  SLUDData.options.Fact = FACTORED;
  pdgssvx( &SLUDData.options,
           &SLUDData.mat,
           &SLUDData.ScalePermstruct,
           x.extractLocalVector(),
           ldb,
           nrhs,
           &SLUDData.grid,
           &SLUDData.LUstruct,
           &SLUDData.SOLVEstruct,
           berr.data(),
           &SLUDData.stat,
           &info );

  time = watch.elapsedTime();

  if( SLUDData.options.PrintStat == YES )
  {
    // Print the statistics.
    PStatPrint( &SLUDData.options, &SLUDData.stat, &SLUDData.grid );
  }

  return info;
}

void SuperLU_DistDestroy( SuperLU_DistData & SLUDData )
{
  // Deallocate other SuperLU data structures
  dScalePermstructFree( &SLUDData.ScalePermstruct );
  dDestroy_LU( SLUDData.mat.nrow, &SLUDData.grid, &SLUDData.LUstruct );
  dLUstructFree( &SLUDData.LUstruct );
  PStatFree( &SLUDData.stat );
  superlu_gridexit( &SLUDData.grid );
  if( SLUDData.options.SolveInitialized )
  {
    dSolveFinalize( &SLUDData.options, &SLUDData.SOLVEstruct );
  }

  // From HYPRE SuperLU_Dist interfaces (superlu.c)
  // SuperLU frees assigned data, so set them to null before
  // calling hypre_CSRMatrixdestroy on localStrip to avoid memory errors.
  hypre_CSRMatrixI( SLUDData.localStrip ) = NULL;
  hypre_CSRMatrixData( SLUDData.localStrip ) = NULL;
  hypre_CSRMatrixBigJ( SLUDData.localStrip ) = NULL;
  hypre_CSRMatrixDestroy( SLUDData.localStrip );

  Destroy_CompRowLoc_Matrix_dist( &SLUDData.mat );
}

colperm_t const & HypreGetColPermType( string const & value )
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

rowperm_t const & HypreGetRowPermType( string const & value )
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
