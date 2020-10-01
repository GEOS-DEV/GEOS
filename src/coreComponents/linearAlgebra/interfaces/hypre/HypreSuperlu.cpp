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
               "SuperLU_Dist int_t and geosx::globalIndex must both be signed or unsigned" );

static_assert( std::is_same< double, real64 >::value,
               "SuperLU_Dist real and geosx::real64 must be the same type" );

namespace
{

/**
 * @brief Convert GEOSX globalIndex value to SLUD int_t
 * @param index the input value
 * @return the converted value
 */
inline int_t toSuperlu_intT( globalIndex const index )
{
  return LvArray::integerConversion< int_t >( index );
}

/**
 * @brief Converts a non-const array from GEOSX globalIndex type to SLUD int_t
 * @param[in] index the input array
 * @return the converted array
 */
inline int_t * toSuperlu_intT( globalIndex * const index )
{
  return reinterpret_cast< int_t * >( index );
}

/**
 * @brief Converts from GEOSX to SuperLU_Dist columns permutation option
 * @param[in] value the GEOSX option
 * @return the SuperLU_Dist option
 */
colperm_t const & getColPermType( LinearSolverParameters::Direct::ColPerm const & value )
{
  static std::map< LinearSolverParameters::Direct::ColPerm, colperm_t > const OPTION_MAP =
  {
    { LinearSolverParameters::Direct::ColPerm::none, NATURAL },
    { LinearSolverParameters::Direct::ColPerm::MMD_AtplusA, MMD_AT_PLUS_A },
    { LinearSolverParameters::Direct::ColPerm::MMD_AtA, MMD_ATA },
    { LinearSolverParameters::Direct::ColPerm::colAMD, COLAMD },
    { LinearSolverParameters::Direct::ColPerm::metis, METIS_AT_PLUS_A },
    { LinearSolverParameters::Direct::ColPerm::parmetis, PARMETIS },
  };

  GEOSX_LAI_ASSERT_MSG( OPTION_MAP.count( value ) > 0, "Unsupported SuperLU_Dist columns permutation option: " << value );
  return OPTION_MAP.at( value );
}

/**
 * @brief Converts from GEOSX to SuperLU_Dist rows permutation option
 * @param[in] value the GEOSX option
 * @return the SuperLU_Dist option
 */
rowperm_t const & getRowPermType( LinearSolverParameters::Direct::RowPerm const & value )
{
  static std::map< LinearSolverParameters::Direct::RowPerm, rowperm_t > const OPTION_MAP =
  {
    { LinearSolverParameters::Direct::RowPerm::none, NOROWPERM },
    { LinearSolverParameters::Direct::RowPerm::mc64, LargeDiag_MC64 },
  };

  GEOSX_LAI_ASSERT_MSG( OPTION_MAP.count( value ) > 0, "Unsupported SuperLU_Dist rows permutation option: " << value );
  return OPTION_MAP.at( value );
}

/**
 * @brief Converts a matrix from Hypre to SuperLU_Dist format
 * @param[in] matrix the HypreMatrix object
 * @param[out] SLUDData the structure containing the matrix in SuperLU_Dist format
 */
void ConvertToSuperMatrix( HypreMatrix const & matrix,
                           SuperLU_DistData & sludData )
{
  // Merge diag and offd into one matrix (global ids)
  sludData.localStrip = hypre_MergeDiagAndOffd( matrix.unwrapped() );

  HYPRE_Int const * const hypreI = hypre_CSRMatrixI( sludData.localStrip );
  sludData.rowPtr = new int_t[matrix.numLocalRows()+1];
  for( localIndex i = 0; i <= matrix.numLocalRows(); ++i )
  {
    sludData.rowPtr[i] = LvArray::integerConversion< int_t >( hypreI[i] );
  }

  dCreate_CompRowLoc_Matrix_dist( &sludData.mat,
                                  toSuperlu_intT( matrix.numGlobalRows() ),
                                  toSuperlu_intT( matrix.numGlobalRows() ),
                                  toSuperlu_intT( matrix.numLocalNonzeros() ),
                                  toSuperlu_intT( matrix.numLocalRows() ),
                                  toSuperlu_intT( matrix.ilower() ),
                                  hypre_CSRMatrixData( sludData.localStrip ),
                                  toSuperlu_intT( hypre_CSRMatrixBigJ( sludData.localStrip ) ),
                                  sludData.rowPtr,
                                  SLU_NR_loc,
                                  SLU_D,
                                  SLU_GE );
}
}

void SuperLU_DistCreate( HypreMatrix const & matrix,
                         LinearSolverParameters const & params,
                         SuperLU_DistData & sludData )
{
  // Initialize options.
  set_default_options_dist( &sludData.options );
  sludData.options.PrintStat = params.logLevel > 1 ? YES : NO;
  sludData.options.Equil = params.direct.equilibrate ? YES : NO;
  sludData.options.ColPerm = getColPermType( params.direct.colPerm );
  sludData.options.RowPerm = getRowPermType( params.direct.rowPerm );
  sludData.options.ParSymbFact = params.direct.colPerm == LinearSolverParameters::Direct::ColPerm::parmetis ? YES : NO;
  sludData.options.ReplaceTinyPivot = params.direct.replaceTinyPivot ? YES : NO;
  sludData.options.IterRefine = params.direct.iterativeRefine ? SLU_DOUBLE : NOREFINE;

  if( params.logLevel > 0 )
  {
    print_sp_ienv_dist( &sludData.options );
    print_options_dist( &sludData.options );
  }

  // Convert matrix from Hypre to SuperLU_Dist format
  ConvertToSuperMatrix( matrix, sludData );

  // Save communicator
  sludData.comm = matrix.getComm();
}

int SuperLU_DistSetup( SuperLU_DistData & sludData,
                       real64 & time )
{
  Stopwatch watch;

  int_t const m = sludData.mat.nrow;
  int_t const n = sludData.mat.ncol;

  // Initialize ScalePermstruct.
  dScalePermstructInit( m, n, &sludData.ScalePermstruct );

  // Initialize LUstruct.
  dLUstructInit( n, &sludData.LUstruct );

  // Initialize the statistics variables.
  PStatInit( &sludData.stat );

  // Create process grid: the target is to have the process grid as square as possible
  int const num_procs = MpiWrapper::Comm_size( sludData.comm );
  int prows = static_cast< int >( std::sqrt( num_procs ) );
  while( num_procs % prows )
  {
    --prows;
  }
  int pcols = num_procs/prows;
  std::tie( prows, pcols ) = std::minmax( prows, pcols );

  superlu_gridinit( sludData.comm, prows, pcols, &sludData.grid );

  // Call the linear equation solver to factorize the matrix.
  int const nrhs = 0;
  int info = 0;

  sludData.options.Fact = DOFACT;
  pdgssvx( &sludData.options,
           &sludData.mat,
           &sludData.ScalePermstruct,
           NULL,
           n,
           nrhs,
           &sludData.grid,
           &sludData.LUstruct,
           &sludData.SOLVEstruct,
           NULL,
           &sludData.stat,
           &info );

  time = watch.elapsedTime();

  if( sludData.options.PrintStat == YES )
  {
    // Print the statistics.
    PStatPrint( &sludData.options, &sludData.stat, &sludData.grid );
  }

  return info;
}

int SuperLU_DistSolve( SuperLU_DistData & sludData,
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

  sludData.options.Fact = FACTORED;
  pdgssvx( &sludData.options,
           &sludData.mat,
           &sludData.ScalePermstruct,
           x.extractLocalVector(),
           ldb,
           nrhs,
           &sludData.grid,
           &sludData.LUstruct,
           &sludData.SOLVEstruct,
           berr.data(),
           &sludData.stat,
           &info );

  time = watch.elapsedTime();

  // Check for nan or inf
  if( std::isnan( berr[0] ) || std::isinf( berr[0] ) )
  {
    info = 1;
  }

  if( sludData.options.PrintStat == YES )
  {
    // Print the statistics.
    PStatPrint( &sludData.options, &sludData.stat, &sludData.grid );
  }

  return info;
}

void SuperLU_DistDestroy( SuperLU_DistData & sludData )
{
  // Deallocate other SuperLU data structures
  dScalePermstructFree( &sludData.ScalePermstruct );
  dDestroy_LU( sludData.mat.nrow, &sludData.grid, &sludData.LUstruct );
  dLUstructFree( &sludData.LUstruct );
  PStatFree( &sludData.stat );
  superlu_gridexit( &sludData.grid );
  if( sludData.options.SolveInitialized )
  {
    dSolveFinalize( &sludData.options, &sludData.SOLVEstruct );
  }

  // From HYPRE SuperLU_Dist interfaces (superlu.c)
  // SuperLU frees assigned data, so set them to null before
  // calling hypre_CSRMatrixdestroy on localStrip to avoid memory errors.
  hypre_CSRMatrixI( sludData.localStrip ) = NULL;
  hypre_CSRMatrixData( sludData.localStrip ) = NULL;
  hypre_CSRMatrixBigJ( sludData.localStrip ) = NULL;
  hypre_CSRMatrixDestroy( sludData.localStrip );

  Destroy_CompRowLoc_Matrix_dist( &sludData.mat );
}

}
