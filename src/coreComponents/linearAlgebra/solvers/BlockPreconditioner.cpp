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
 * @file BlockPreconditioner.cpp
 */

#include "BlockPreconditioner.hpp"
#include "codingUtilities/Utilities.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/utilities/BlockVectorWrapper.hpp"

#include <numeric>

namespace geos
{

template< typename LAI >
BlockPreconditioner< LAI >::BlockPreconditioner( LinearSolverParameters::Block params )
  : Base(),
  m_params( params ),
  m_matBlocks( 2, 2 ),
  m_solvers{},
  m_scaling{ 1.0, 1.0 },
  m_rhs( 2 ),
  m_sol( 2 )
{}

template< typename LAI >
void BlockPreconditioner< LAI >::reinitialize( Matrix const & mat )
{
  MPI_Comm const & comm = mat.comm();

  if( m_blockDofs[1].empty() )
  {
    GEOS_LAI_ASSERT( mat.dofManager() != nullptr );
    m_blockDofs[1] = mat.dofManager()->filterDofs( m_blockDofs[0] );
  }

  for( localIndex i = 0; i < 2; ++i )
  {
    if( m_prolongators[i] == nullptr )
    {
      GEOS_LAI_ASSERT( mat.dofManager() != nullptr );
      mat.dofManager()->makeRestrictor( m_blockDofs[i], comm, true, m_prolongatorsOwned[i] );
      m_prolongators[i] = &m_prolongatorsOwned[i];
    }
    m_rhs( i ).create( m_prolongators[i]->numLocalCols(), comm );
    m_sol( i ).create( m_prolongators[i]->numLocalCols(), comm );
  }
}

template< typename LAI >
void BlockPreconditioner< LAI >::setupBlock( localIndex const blockIndex,
                                             std::vector< DofManager::SubComponent > blockDofs,
                                             std::unique_ptr< PreconditionerBase< LAI > > solver,
                                             real64 const scaling )
{
  GEOS_LAI_ASSERT_GT( 2, blockIndex );
  GEOS_LAI_ASSERT( solver );
  GEOS_LAI_ASSERT( !blockDofs.empty() );
  GEOS_LAI_ASSERT_GT( scaling, 0.0 );

  m_blockDofs[blockIndex] = std::move( blockDofs );
  m_solversOwned[blockIndex] = std::move( solver );
  m_solvers[blockIndex] = m_solversOwned[blockIndex].get();
  m_scaling[blockIndex] = scaling;
}

template< typename LAI >
void BlockPreconditioner< LAI >::setupBlock( localIndex const blockIndex,
                                             std::vector< DofManager::SubComponent > blockDofs,
                                             PreconditionerBase< LAI > * const solver,
                                             real64 const scaling )
{
  GEOS_LAI_ASSERT_GT( 2, blockIndex );
  GEOS_LAI_ASSERT( solver );
  GEOS_LAI_ASSERT( !blockDofs.empty() );
  GEOS_LAI_ASSERT_GT( scaling, 0.0 );

  m_blockDofs[blockIndex] = std::move( blockDofs );
  m_solversOwned[blockIndex].reset();
  m_solvers[blockIndex] = solver;
  m_scaling[blockIndex] = scaling;
}

template< typename LAI >
void BlockPreconditioner< LAI >::setProlongation( localIndex const blockIndex,
                                                  Matrix const & P )
{
  GEOS_LAI_ASSERT_GT( 2, blockIndex );

  m_prolongatorsOwned[blockIndex].reset();
  m_prolongators[blockIndex] = &P;
}

template< typename LAI >
void BlockPreconditioner< LAI >::applyBlockScaling()
{
  if( m_params.scaling != LinearSolverParameters::Block::Scaling::None )
  {
    if( m_params.scaling == LinearSolverParameters::Block::Scaling::FrobeniusNorm )
    {
      real64 const norms[2] = { m_matBlocks( 0, 0 ).normFrobenius(), m_matBlocks( 1, 1 ).normFrobenius() };
      m_scaling[0] = std::min( norms[1] / norms[0], 1.0 );
      m_scaling[1] = std::min( norms[0] / norms[1], 1.0 );
    }

    for( localIndex i = 0; i < 2; ++i )
    {
      for( localIndex j = 0; j < 2; ++j )
      {
        m_matBlocks( i, j ).scale( m_scaling[i] );
      }
    }
  }
}

template< typename LAI >
void BlockPreconditioner< LAI >::computeSchurComplement()
{
  switch( m_params.schurType )
  {
    case LinearSolverParameters::Block::SchurType::None:
    {
      // nothing to do
      break;
    }
    case LinearSolverParameters::Block::SchurType::FirstBlockDiagonal:
    {
      m_matBlocks( 0, 0 ).extractDiagonal( m_rhs( 0 ) );
      m_rhs( 0 ).reciprocal();
      m_matBlocks( 0, 1 ).leftScale( m_rhs( 0 ) );
      Matrix mat11;
      m_matBlocks( 1, 0 ).multiply( m_matBlocks( 0, 1 ), mat11 );
      m_matBlocks( 1, 1 ).addEntries( mat11, MatrixPatternOp::Restrict, -1.0 );
      // Restore original scaling
      m_rhs( 0 ).reciprocal();
      m_matBlocks( 0, 1 ).leftScale( m_rhs( 0 ) );
      break;
    }
    case LinearSolverParameters::Block::SchurType::RowsumDiagonalProbing:
    {
      m_sol( 1 ).set( -1.0 );
      m_matBlocks( 0, 1 ).apply( m_sol( 1 ), m_rhs( 0 ) );
      m_solvers[0]->apply( m_rhs( 0 ), m_sol( 0 ) );
      m_matBlocks( 1, 0 ).apply( m_sol( 0 ), m_rhs( 1 ) );
      m_matBlocks( 1, 1 ).addDiagonal( m_rhs( 1 ), 1.0 );
      break;
    }
    case LinearSolverParameters::Block::SchurType::FirstBlockUserDefined:
    {
      Matrix const & prec00 = m_solvers[0]->preconditionerMatrix();
      Matrix mat11;
      prec00.multiplyRAP( m_matBlocks( 1, 0 ), m_matBlocks( 0, 1 ), mat11 );
      m_matBlocks( 1, 1 ).addEntries( mat11, MatrixPatternOp::Extend, -1.0 );
      break;
    }
    default:
    {
      GEOS_ERROR( "BlockPreconditioner: unsupported Schur complement option" );
    }
  }
}

template< typename LAI >
void BlockPreconditioner< LAI >::setup( Matrix const & mat )
{
  // Check that user has set block solvers
  GEOS_LAI_ASSERT( m_solvers[0] != nullptr );
  GEOS_LAI_ASSERT( m_solvers[1] != nullptr );

  // Compare old sizes vs new matris sizes.
  // A change in size indicates a new matrix structure.
  // This is done before Base::compute() since it overwrites old sizes.
  bool const newSize = !this->ready() ||
                       mat.numGlobalRows() != this->numGlobalRows() ||
                       mat.numGlobalCols() != this->numGlobalRows();

  Base::setup( mat );

  // If the matrix size/structure has changed, need to resize internal LA objects and recompute restrictors.
  if( newSize )
  {
    reinitialize( mat );
  }

  // Extract diagonal blocks
  mat.multiplyPtAP( *m_prolongators[0], m_matBlocks( 0, 0 ) );
  mat.multiplyPtAP( *m_prolongators[1], m_matBlocks( 1, 1 ) );

  // HACK: a coupled DofManager is technically not compatible with diagonal blocks.
  // We can create "reduced" managers using DofManager::setupFrom(), but it can be a waste of time,
  // Instead, we expect block solvers to not rely on global information and instead query by field.
  m_matBlocks( 0, 0 ).setDofManager( mat.dofManager() );
  m_matBlocks( 1, 1 ).setDofManager( mat.dofManager() );

  // Extract off-diagonal blocks only if used
  if( m_params.schurType != LinearSolverParameters::Block::SchurType::None
      || m_params.shape != LinearSolverParameters::Block::Shape::Diagonal )
  {
    mat.multiplyPtAP( *m_prolongators[0], *m_prolongators[1], m_matBlocks( 0, 1 ) );
    mat.multiplyPtAP( *m_prolongators[1], *m_prolongators[0], m_matBlocks( 1, 0 ) );
  }

  applyBlockScaling();
  m_solvers[0]->setup( m_matBlocks( 0, 0 ) );
  computeSchurComplement();
  m_solvers[1]->setup( m_matBlocks( 1, 1 ) );
}

template< typename LAI >
void BlockPreconditioner< LAI >::apply( Vector const & src,
                                        Vector & dst ) const
{
  using Shape = LinearSolverParameters::Block::Shape;

  m_prolongators[0]->applyTranspose( src, m_rhs( 0 ) );
  m_prolongators[1]->applyTranspose( src, m_rhs( 1 ) );

  for( localIndex i = 0; i < 2; ++i )
  {
    m_rhs( i ).scale( m_scaling[i] );
  }

  if( m_params.shape == Shape::LowerUpperTriangular || m_params.shape == Shape::LowerTriangular )
  {
    // Solve the (0,0)-block and update (1,1)-block rhs
    m_solvers[0]->apply( m_rhs( 0 ), m_sol( 0 ) );
    m_matBlocks( 1, 0 ).residual( m_sol( 0 ), m_rhs( 1 ), m_rhs( 1 ) );
  }

  // Solve the (1,1)-block modified via Schur complement
  m_solvers[1]->apply( m_rhs( 1 ), m_sol( 1 ) );

  if( m_params.shape != Shape::LowerTriangular )
  {
    if( m_params.shape != Shape::Diagonal )
    {
      // Update the 0-block rhs
      m_matBlocks( 0, 1 ).residual( m_sol( 1 ), m_rhs( 0 ), m_rhs( 0 ) );
    }

    // Solve the (0,0) block with the current rhs
    m_solvers[0]->apply( m_rhs( 0 ), m_sol( 0 ) );
  }

  // Combine block solutions into global solution vector
  m_prolongators[0]->apply( m_sol( 0 ), dst );
  m_prolongators[1]->gemv( 1.0, m_sol( 1 ), 1.0, dst );
}

template< typename LAI >
void BlockPreconditioner< LAI >::clear()
{
  Base::clear();
  for( localIndex i = 0; i < 2; ++i )
  {
    m_prolongatorsOwned[i].reset();
    m_prolongators[i] = nullptr;
    m_solvers[i]->clear();
    m_rhs( i ).reset();
    m_sol( i ).reset();
    for( localIndex j = 0; j < 2; ++j )
    {
      m_matBlocks( i, j ).reset();
    }
  }
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class BlockPreconditioner< TrilinosInterface >;
#endif

#ifdef GEOSX_USE_HYPRE
template class BlockPreconditioner< HypreInterface >;
#endif

#ifdef GEOSX_USE_PETSC
template class BlockPreconditioner< PetscInterface >;
#endif

}
