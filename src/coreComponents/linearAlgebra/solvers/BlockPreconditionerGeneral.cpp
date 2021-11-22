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

/**
 * @file BlockPreconditionerGeneral.cpp
 */

#include "BlockPreconditionerGeneral.hpp"
#include "codingUtilities/Utilities.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/utilities/BlockVectorWrapper.hpp"

namespace geosx
{

template< typename LAI >
BlockPreconditionerGeneral< LAI >::BlockPreconditionerGeneral( localIndex const numBlocks,
                                                               BlockShapeOption const shapeOption,
                                                               std::vector< SchurComplementOption > const & schurOption )
  : Base(),
  m_numBlocks( numBlocks ),
  m_shapeOption( shapeOption ),
  m_matBlocks( numBlocks, numBlocks ),
  m_solvers{},
  m_rhs( 2*numBlocks ),
  m_sol( 2*numBlocks )
{
  m_schurOption.resize( numBlocks );
  m_blockDofs.resize( numBlocks );
  m_restrictors.resize( numBlocks );
  m_prolongators.resize( numBlocks );
  m_solvers.resize( numBlocks );

  for( localIndex iBlock = 0; iBlock < m_numBlocks-1; ++iBlock )
  {
    m_schurOption[iBlock] = schurOption[iBlock];
  }
}

template< typename LAI >
BlockPreconditionerGeneral< LAI >::~BlockPreconditionerGeneral() = default;

template< typename LAI >
void BlockPreconditionerGeneral< LAI >::reinitialize( Matrix const & mat, DofManager const & dofManager, localIndex const iBlock )
{
  MPI_Comm const & comm = mat.comm();

  std::vector< DofManager::SubComponent > dofsList;
  for( localIndex i = 0; i <= iBlock; ++i )
  {
    for( std::size_t j = 0; j < m_blockDofs[i].size(); ++j )
    {
      dofsList.push_back( m_blockDofs[i][j] );
    }
  }

  std::vector< DofManager::SubComponent > const blockDofs = dofManager.filterDofs( dofsList );

  dofManager.makeRestrictor( m_blockDofs[iBlock], comm, false, m_restrictors[iBlock].current );
  dofManager.makeRestrictor( m_blockDofs[iBlock], comm, true, m_prolongators[iBlock].current );
  dofManager.makeRestrictor( blockDofs, comm, false, m_restrictors[iBlock].next );
  dofManager.makeRestrictor( blockDofs, comm, true, m_prolongators[iBlock].next );

  for( localIndex i = 0; i < iBlock; ++i )
  {
    Matrix prolongatorCurrent( m_prolongators[iBlock].current );
    m_restrictors[i].next.multiply( prolongatorCurrent, m_prolongators[iBlock].current );
    Matrix prolongatorNext( m_prolongators[iBlock].next );
    m_restrictors[i].next.multiply( prolongatorNext, m_prolongators[iBlock].next );

    Matrix restrictorCurrent( m_restrictors[iBlock].current );
    restrictorCurrent.multiply( m_prolongators[i].next, m_restrictors[iBlock].current );
    Matrix restrictorNext( m_restrictors[iBlock].next );
    restrictorNext.multiply( m_prolongators[i].next, m_restrictors[iBlock].next );
  }

  m_rhs( 2*iBlock ).create( m_restrictors[iBlock].current.numLocalRows(), comm );
  m_rhs( 2*iBlock+1 ).create( m_restrictors[iBlock].next.numLocalRows(), comm );
  m_sol( 2*iBlock ).create( m_restrictors[iBlock].current.numLocalRows(), comm );
  m_sol( 2*iBlock+1 ).create( m_restrictors[iBlock].next.numLocalRows(), comm );
}

template< typename LAI >
void BlockPreconditionerGeneral< LAI >::setupBlock( localIndex const blockIndex,
                                                    std::vector< DofManager::SubComponent > blockDofs,
                                                    std::unique_ptr< PreconditionerBase< LAI > > solver )
{
  GEOSX_LAI_ASSERT_GT( m_numBlocks, blockIndex );
  GEOSX_LAI_ASSERT( solver );
  GEOSX_LAI_ASSERT( !blockDofs.empty() );

  m_blockDofs[blockIndex] = std::move( blockDofs );
  m_solvers[blockIndex] = std::move( solver );
}

template< typename LAI >
void BlockPreconditionerGeneral< LAI >::computeSchurComplement( localIndex const iBlock )
{
  switch( m_schurOption[iBlock] )
  {
    case SchurComplementOption::None:
    {
      // nothing to do
      break;
    }
    case SchurComplementOption::Diagonal:
    {
      m_matBlocks( iBlock, iBlock ).extractDiagonal( m_rhs( 2*iBlock ) );
      m_rhs( 2*iBlock ).reciprocal();
      m_matBlocks( iBlock, iBlock+1 ).leftScale( m_rhs( 2*iBlock ) );
      Matrix mat11;
      m_matBlocks( iBlock+1, iBlock ).multiply( m_matBlocks( iBlock, iBlock+1 ), mat11 );
      m_matBlocks( iBlock+1, iBlock+1 ).addEntries( mat11, MatrixPatternOp::Restrict, -1.0 );
      // Restore original scaling
      m_rhs( 2*iBlock ).reciprocal();
      m_matBlocks( iBlock, iBlock+1 ).leftScale( m_rhs( 2*iBlock ) );
      break;
    }
    case SchurComplementOption::RowsumDiagonalProbing:
    {
      m_sol( 2*(iBlock+1) ).set( -1.0 );
      m_matBlocks( iBlock, iBlock+1 ).apply( m_sol( 2*(iBlock+1) ), m_rhs( 2*iBlock ) );
      m_solvers[iBlock]->apply( m_rhs( 2*iBlock ), m_sol( 2*iBlock ) );
      m_matBlocks( iBlock+1, iBlock ).apply( m_sol( 2*iBlock ), m_rhs( 2*(iBlock+1) ) );
      m_matBlocks( iBlock+1, iBlock+1 ).addDiagonal( m_rhs( 2*(iBlock+1) ), 1.0 );
      break;
    }
    case SchurComplementOption::UserDefined:
    {
      Matrix const & prec00 = m_solvers[iBlock]->preconditionerMatrix();
      Matrix mat11;
      prec00.multiplyRAP( m_matBlocks( iBlock+1, iBlock ), m_matBlocks( iBlock, iBlock+1 ), mat11 );
      m_matBlocks( iBlock+1, iBlock+1 ).addEntries( mat11, MatrixPatternOp::Extend, -1.0 );
      break;
    }
    default:
    {
      GEOSX_ERROR( "BlockPreconditioner: unsupported Schur complement option" );
    }
  }
}

template< typename LAI >
void BlockPreconditionerGeneral< LAI >::setup( Matrix const & mat )
{
  for( localIndex iBlock = 0; iBlock < m_numBlocks; ++iBlock )
  {
    GEOSX_LAI_ASSERT( m_solvers[iBlock] != nullptr );
  }

  // Compare old sizes vs new matris sizes.
  // A change in size indicates a new matrix structure.
  // This is done before Base::compute() since it overwrites old sizes.
  bool const newSize = !this->ready() ||
                       mat.numGlobalRows() != this->numGlobalRows() ||
                       mat.numGlobalCols() != this->numGlobalRows();

  Base::setup( mat );

  if( newSize )
  {
    reinitialize( mat, *mat.dofManager(), 0 );
  }
  mat.multiplyPtAP( m_prolongators[0].current, m_matBlocks( 0, 0 ) );
  mat.multiplyPtAP( m_prolongators[0].next, m_matBlocks( 1, 1 ) );
  mat.multiplyRAP( m_restrictors[0].current, m_prolongators[0].next, m_matBlocks( 0, 1 ) );
  mat.multiplyRAP( m_restrictors[0].next, m_prolongators[0].current, m_matBlocks( 1, 0 ) );

  m_solvers[0]->setup( m_matBlocks( 0, 0 ) );
  computeSchurComplement( 0 );

  // Compute all Schur complements
  for( localIndex iBlock = 1; iBlock < m_numBlocks-1; ++iBlock )
  {
    // If the matrix size/structure has changed, need to resize internal LA objects and recompute restrictors.
    if( newSize )
    {
      reinitialize( mat, *mat.dofManager(), iBlock );
    }
    Matrix localMat( m_matBlocks( iBlock, iBlock ) );
    localMat.multiplyPtAP( m_prolongators[iBlock].current, m_matBlocks( iBlock, iBlock ) );
    localMat.multiplyPtAP( m_prolongators[iBlock].next, m_matBlocks( iBlock+1, iBlock+1 ) );

    localMat.multiplyRAP( m_restrictors[iBlock].current, m_prolongators[iBlock].next, m_matBlocks( iBlock, iBlock+1 ) );
    localMat.multiplyRAP( m_restrictors[iBlock].next, m_prolongators[iBlock].current, m_matBlocks( iBlock+1, iBlock ) );

    m_solvers[iBlock]->setup( m_matBlocks( iBlock, iBlock ) );
    computeSchurComplement( iBlock );
  }

  // Last Schur complement
  if( newSize )
  {
    reinitialize( mat, *mat.dofManager(), m_numBlocks-1 );
  }
  m_solvers[m_numBlocks-1]->setup( m_matBlocks( m_numBlocks-1, m_numBlocks-1 ) );
}

template< typename LAI >
void BlockPreconditionerGeneral< LAI >::apply( Vector const & src,
                                               Vector & dst ) const
{
  m_restrictors[0].current.apply( src, m_rhs( 0 ) );
  m_restrictors[0].next.apply( src, m_rhs( 1 ) );

  // Perform a predictor step by solving (0,0) block and subtracting from 1-block rhs
  if( m_shapeOption == BlockShapeOption::LowerUpperTriangular )
  {
    m_solvers[0]->apply( m_rhs( 0 ), m_sol( 0 ) );
    m_matBlocks( 1, 0 ).residual( m_sol( 0 ), m_rhs( 1 ), m_rhs( 1 ) );
  }

  // Solve U^{-1} for all blocks
  for( localIndex iBlock = 1; iBlock < m_numBlocks-1; ++iBlock )
  {
    m_restrictors[iBlock].current.apply( m_rhs( 2*iBlock-1 ), m_rhs( 2*iBlock ) );
    m_restrictors[iBlock].next.apply( m_rhs( 2*iBlock-1 ), m_rhs( 2*iBlock+1 ) );

    // Perform a predictor step by solving (i,i) block and subtracting from (i+1)-th block rhs
    if( m_shapeOption == BlockShapeOption::LowerUpperTriangular )
    {
      m_solvers[iBlock]->apply( m_rhs( 2*iBlock ), m_sol( 2*iBlock ) );
      m_matBlocks( iBlock+1, iBlock ).residual( m_sol( 2*iBlock ), m_rhs( 2*iBlock+1 ), m_rhs( 2*iBlock+1 ) );
    }
  }

  // Last Schur complement
  m_solvers[m_numBlocks-1]->apply( m_rhs( 2*(m_numBlocks-1)-1 ), m_sol( 2*(m_numBlocks-1)-1 ) );

  // Solve L^{-1} for all blocks
  for( localIndex iBlock = m_numBlocks-2; iBlock > 0; --iBlock )
  {
    // Update the i-th block rhs
    if( m_shapeOption == BlockShapeOption::LowerUpperTriangular )
    {
      m_matBlocks( iBlock, iBlock+1 ).residual( m_sol( 2*iBlock+1 ), m_rhs( 2*iBlock ), m_rhs( 2*iBlock ) );
      m_solvers[iBlock]->apply( m_rhs( 2*iBlock ), m_sol( 2*iBlock ) );
    }

    // Combine block solutions into global solution vector
    m_prolongators[iBlock].current.apply( m_sol( 2*iBlock ), m_sol( 2*(iBlock-1)+1 ) );
    m_prolongators[iBlock].next.gemv( 1.0, m_sol( 2*iBlock+1 ), 1.0, m_sol( 2*(iBlock-1)+1 ) );
  }

  // Update the 0-block rhs
  if( m_shapeOption != BlockShapeOption::Diagonal )
  {
    m_matBlocks( 0, 1 ).residual( m_sol( 1 ), m_rhs( 0 ), m_rhs( 0 ) );
  }

  // Solve the (0,0) block with the current rhs
  m_solvers[0]->apply( m_rhs( 0 ), m_sol( 0 ) );

  // Combine block solutions into global solution vector
  m_prolongators[0].current.apply( m_sol( 0 ), dst );
  m_prolongators[0].next.gemv( 1.0, m_sol( 1 ), 1.0, dst );
}

template< typename LAI >
void BlockPreconditionerGeneral< LAI >::clear()
{
  Base::clear();
  for( localIndex iBlock = 0; iBlock < m_numBlocks; ++iBlock )
  {
    m_restrictors[iBlock].current.reset();
    m_restrictors[iBlock].next.reset();
    m_prolongators[iBlock].current.reset();
    m_prolongators[iBlock].next.reset();
    m_solvers[iBlock]->clear();
    m_rhs( 2*iBlock ).reset();
    m_rhs( 2*iBlock+1 ).reset();
    m_sol( 2*iBlock ).reset();
    m_sol( 2*iBlock+1 ).reset();
    for( localIndex jBlock = 0; jBlock < m_numBlocks; ++jBlock )
    {
      m_matBlocks( iBlock, jBlock ).reset();
    }
  }
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class BlockPreconditionerGeneral< TrilinosInterface >;
#endif

#ifdef GEOSX_USE_HYPRE
template class BlockPreconditionerGeneral< HypreInterface >;
#endif

#ifdef GEOSX_USE_PETSC
template class BlockPreconditionerGeneral< PetscInterface >;
#endif

}
