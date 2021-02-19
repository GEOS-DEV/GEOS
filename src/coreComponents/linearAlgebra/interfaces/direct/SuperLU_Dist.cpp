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
 * @file SuperLU_Dist.cpp
 */

#include "SuperLU_Dist.hpp"
#include "linearAlgebra/common.hpp"
#include "mpiCommunications/MpiWrapper.hpp"
#include "common/Stopwatch.hpp"

namespace geosx
{

// Check matching requirements on index/value types between GEOSX and SuperLU_Dist

//static_assert( sizeof( int_t ) == sizeof( globalIndex ),
//               "SuperLU_Dist int_t and geosx::globalIndex must have the same size" );
//
//static_assert( std::is_signed< int_t >::value == std::is_signed< globalIndex >::value,
//               "SuperLU_Dist int_t and geosx::globalIndex must both be signed or unsigned" );

static_assert( std::is_same< double, real64 >::value,
               "SuperLU_Dist real and geosx::real64 must be the same type" );

namespace
{

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
}

SuperLU_Dist::SuperLU_Dist():
  m_condEst( -1 ),
  m_setupTime( 0 ),
  m_solveTime( 0 )
{}

SuperLU_Dist::SuperLU_Dist( LinearSolverParameters const & params ):
  m_condEst( -1 ),
  m_setupTime( 0 ),
  m_solveTime( 0 )
{
  create( params );
}

SuperLU_Dist::~SuperLU_Dist()
{
  destroy();
}

void SuperLU_Dist::create( LinearSolverParameters const & params )
{
  // Initialize options.
  set_default_options_dist( &m_options );
  m_options.PrintStat = params.logLevel > 1 ? YES : NO;
  m_options.Equil = params.direct.equilibrate ? YES : NO;
  m_options.ColPerm = getColPermType( params.direct.colPerm );
  m_options.RowPerm = getRowPermType( params.direct.rowPerm );
  m_options.ParSymbFact = params.direct.colPerm == LinearSolverParameters::Direct::ColPerm::parmetis ? YES : NO;
  m_options.ReplaceTinyPivot = params.direct.replaceTinyPivot ? YES : NO;
  m_options.IterRefine = params.direct.iterativeRefine ? SLU_DOUBLE : NOREFINE;

  if( params.logLevel > 0 )
  {
    print_sp_ienv_dist( &m_options );
    print_options_dist( &m_options );
  }

  // Save parameters
  m_params = params;
}

int SuperLU_Dist::setup()
{
  Stopwatch watch;

  int_t const m = m_mat.nrow;
  int_t const n = m_mat.ncol;

  // Initialize ScalePermstruct.
  dScalePermstructInit( m, n, &m_ScalePermstruct );

  // Initialize LUstruct.
  dLUstructInit( n, &m_LUstruct );

  // Initialize the statistics variables.
  PStatInit( &m_stat );

  // Create process grid: the target is to have the process grid as square as possible
  int const num_procs = MpiWrapper::commSize( m_comm );
  int prows = static_cast< int >( std::sqrt( num_procs ) );
  while( num_procs % prows )
  {
    --prows;
  }
  int pcols = num_procs / prows;
  std::tie( prows, pcols ) = std::minmax( prows, pcols );

  superlu_gridinit( m_comm, prows, pcols, &m_grid );

  // Call the linear equation solver to factorize the matrix.
  int const nrhs = 0;
  int info = 0;

  m_options.Fact = DOFACT;
  pdgssvx( &m_options,
           &m_mat,
           &m_ScalePermstruct,
           NULL,
           n,
           nrhs,
           &m_grid,
           &m_LUstruct,
           &m_SOLVEstruct,
           NULL,
           &m_stat,
           &info );

  m_setupTime = watch.elapsedTime();

  if( m_options.PrintStat == YES )
  {
    // Print the statistics.
    PStatPrint( &m_options, &m_stat, &m_grid );
  }

  return info;
}

int SuperLU_Dist::solve( real64 const * b, real64 * x )
{
  Stopwatch watch;

  // Call the linear equation solver to solve the matrix.
  int const nrhs = 1;
  array1d< real64 > berr( nrhs );
  int info = 0;

  // Copy rhs in solution vector (SuperLU_Dist works in place)
  std::copy( b, b+m_numLocalRows, x );

  m_options.Fact = FACTORED;
  pdgssvx( &m_options,
           &m_mat,
           &m_ScalePermstruct,
           x,
           m_numLocalRows,
           nrhs,
           &m_grid,
           &m_LUstruct,
           &m_SOLVEstruct,
           berr.data(),
           &m_stat,
           &info );

  m_solveTime = watch.elapsedTime();

  // Check for nan or inf
  if( std::isnan( berr[0] ) || std::isinf( berr[0] ) )
  {
    info = 1;
  }

  if( m_options.PrintStat == YES )
  {
    // Print the statistics.
    PStatPrint( &m_options, &m_stat, &m_grid );
  }

  return info;
}

real64 SuperLU_Dist::condEst()
{
  array1d< real64 > diagU( m_mat.nrow );
  pdGetDiagU( m_mat.nrow,
              &m_LUstruct,
              &m_grid,
              diagU.data() );

  if( m_options.Equil == YES )
  {
    real64 const * R = m_ScalePermstruct.R;
    real64 const * C = m_ScalePermstruct.C;
    int_t const * perm_c = m_ScalePermstruct.perm_c;

    real64 minU = std::abs( diagU[perm_c[0]] / C[0] );
    real64 maxU = std::abs( diagU[perm_c[0]] / C[0] );
    real64 minL = std::abs( 1.0 / R[0] );
    real64 maxL = std::abs( 1.0 / R[0] );

    for( globalIndex i = 1; i < LvArray::integerConversion< globalIndex >( m_mat.nrow ); ++i )
    {
      real64 const u = std::abs( diagU[perm_c[i]] / C[i] );
      minU = ( u < minU ) ? u : minU;
      maxU = ( u > maxU ) ? u : maxU;
      real64 const l = std::abs( 1.0 / R[i] );
      minL = ( l < minL ) ? l : minL;
      maxL = ( l > maxL ) ? l : maxL;
    }
    m_condEst = ( maxU / minU ) * ( maxL / minL );
  }
  else
  {
    real64 minU = std::abs( diagU[0] );
    real64 maxU = std::abs( diagU[0] );
    for( globalIndex i = 1; i < LvArray::integerConversion< globalIndex >( m_mat.nrow ); ++i )
    {
      minU = ( std::abs( diagU[i] ) < minU ) ? std::abs( diagU[i] ) : minU;
      maxU = ( std::abs( diagU[i] ) > maxU ) ? std::abs( diagU[i] ) : maxU;
    }
    m_condEst = ( maxU / minU );
  }
  return m_condEst;
}

real64 SuperLU_Dist::relativeTolerance()
{
  if( m_condEst < 0 )
  {
    condEst();
  }
  return m_condEst * m_precisionTolerance;
}

void SuperLU_Dist::destroy()
{
  // Deallocate other SuperLU data structures
  dScalePermstructFree( &m_ScalePermstruct );
  dDestroy_LU( m_mat.nrow, &m_grid, &m_LUstruct );
  dLUstructFree( &m_LUstruct );
  PStatFree( &m_stat );
  superlu_gridexit( &m_grid );
  if( m_options.SolveInitialized )
  {
    dSolveFinalize( &m_options, &m_SOLVEstruct );
  }
  Destroy_CompRowLoc_Matrix_dist( &m_mat );
}

void SuperLU_Dist::setNumGlobalRows( int_t const numGlobalRows )
{
  m_numGlobalRows = numGlobalRows;
}

int_t SuperLU_Dist::numGlobalRows() const
{
  return m_numGlobalRows;
}

void SuperLU_Dist::setNumGlobalCols( int_t const numGlobalCols )
{
  m_numGlobalCols = numGlobalCols;
}

int_t SuperLU_Dist::numGlobalCols() const
{
  return m_numGlobalCols;
}

int_t SuperLU_Dist::numLocalRows() const
{
  return m_numLocalRows;
}

void SuperLU_Dist::setComm( MPI_Comm const comm )
{
  m_comm = comm;
}

MPI_Comm SuperLU_Dist::getComm() const
{
  return m_comm;
}

void SuperLU_Dist::resize( localIndex const numLocalRows, localIndex const numLocalNonzeros )
{
  m_numLocalRows = numLocalRows;
  m_numLocalNonzeros = numLocalNonzeros;
  // This will be deleted by Destroy_CompRowLoc_Matrix_dist.
  // No need for explicit call to delete!!!
  m_rowPtr = intMalloc_dist( numLocalRows+1 );
  m_colIndices = intMalloc_dist( numLocalNonzeros );
  m_values = doubleMalloc_dist( numLocalNonzeros );
}

void SuperLU_Dist::createSuperMatrix( globalIndex const ilower )
{
  dCreate_CompRowLoc_Matrix_dist( &m_mat,
                                  toSuperLU_intT( m_numGlobalRows ),
                                  toSuperLU_intT( m_numGlobalRows ),
                                  toSuperLU_intT( m_numLocalNonzeros ),
                                  toSuperLU_intT( m_numLocalRows ),
                                  toSuperLU_intT( ilower ),
                                  m_values,
                                  m_colIndices,
                                  m_rowPtr,
                                  SLU_NR_loc,
                                  SLU_D,
                                  SLU_GE );
}

int_t * SuperLU_Dist::rowPtr()
{
  return m_rowPtr;
}

int_t * SuperLU_Dist::colIndices()
{
  return m_colIndices;
}

real64 * SuperLU_Dist::values()
{
  return m_values;
}

real64 SuperLU_Dist::setupTime() const
{
  return m_setupTime;
}

real64 SuperLU_Dist::solveTime() const
{
  return m_solveTime;
}

LinearSolverParameters SuperLU_Dist::getParameters() const
{
  return m_params;
}

real64 SuperLU_Dist::precisionTolerance() const
{
  return m_precisionTolerance;
}

}
