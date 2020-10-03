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

//static_assert( sizeof( int_t ) == sizeof( globalIndex ),
//               "SuperLU_Dist int_t and geosx::globalIndex must have the same size" );

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
  static std::map< LinearSolverParameters::Direct::ColPerm, colperm_t > const optionMap =
  {
    { LinearSolverParameters::Direct::ColPerm::none, NATURAL },
    { LinearSolverParameters::Direct::ColPerm::MMD_AtplusA, MMD_AT_PLUS_A },
    { LinearSolverParameters::Direct::ColPerm::MMD_AtA, MMD_ATA },
    { LinearSolverParameters::Direct::ColPerm::colAMD, COLAMD },
    { LinearSolverParameters::Direct::ColPerm::metis, METIS_AT_PLUS_A },
    { LinearSolverParameters::Direct::ColPerm::parmetis, PARMETIS },
  };

  GEOSX_LAI_ASSERT_MSG( optionMap.count( value ) > 0, "Unsupported SuperLU_Dist columns permutation option: " << value );
  return optionMap.at( value );
}

/**
 * @brief Converts from GEOSX to SuperLU_Dist rows permutation option
 * @param[in] value the GEOSX option
 * @return the SuperLU_Dist option
 */
rowperm_t const & getRowPermType( LinearSolverParameters::Direct::RowPerm const & value )
{
  static std::map< LinearSolverParameters::Direct::RowPerm, rowperm_t > const optionMap =
  {
    { LinearSolverParameters::Direct::RowPerm::none, NOROWPERM },
    { LinearSolverParameters::Direct::RowPerm::mc64, LargeDiag_MC64 },
  };

  GEOSX_LAI_ASSERT_MSG( optionMap.count( value ) > 0, "Unsupported SuperLU_Dist rows permutation option: " << value );
  return optionMap.at( value );
}

/**
 * @brief Converts a matrix from Hypre to SuperLU_Dist format
 * @param[in] matrix the HypreMatrix object
 * @param[out] SLUDData the structure containing the matrix in SuperLU_Dist format
 */
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
}

void SuperLU_DistCreate( HypreMatrix const & matrix,
                         LinearSolverParameters const & params,
                         SuperLU_DistData & SLUDData )
{
  // Initialize options.
  set_default_options_dist( &SLUDData.options );
  SLUDData.options.PrintStat = params.logLevel > 1 ? YES : NO;
  SLUDData.options.Equil = params.direct.equilibrate ? YES : NO;
  SLUDData.options.ColPerm = getColPermType( params.direct.colPerm );
  SLUDData.options.RowPerm = getRowPermType( params.direct.rowPerm );
  SLUDData.options.ParSymbFact = params.direct.colPerm == LinearSolverParameters::Direct::ColPerm::parmetis ? YES : NO;
  SLUDData.options.ReplaceTinyPivot = params.direct.replaceTinyPivot ? YES : NO;
  SLUDData.options.IterRefine = params.direct.iterativeRefine ? SLU_DOUBLE : NOREFINE;

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

  // Create process grid: the target is to have the process grid as square as possible
  int const num_procs = MpiWrapper::Comm_size( SLUDData.comm );
  int prows = static_cast< int >( std::sqrt( num_procs ) );
  while( num_procs % prows )
  {
    --prows;
  }
  int pcols = num_procs/prows;
  std::tie( prows, pcols ) = std::minmax( prows, pcols );

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

  // Check for nan or inf
  if( std::isnan( berr[0] ) || std::isinf( berr[0] ) )
  {
    info = 1;
  }

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

}
