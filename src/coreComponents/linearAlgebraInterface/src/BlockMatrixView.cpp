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
 * @file BlockMatrixView.cpp
 */

#include "BlockMatrixView.hpp"

namespace geosx
{

// BEGIN_RST_NARRATIVE BlockMatrixView.rst
// ==============================
// Block Matrix View
// ==============================
// This class contains the block matrix view object. Each block matrix owns
// pointers to a variable number of ParallelMatrices and is templated on the
// LA interface.

// ----------------------------
// Constructors
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Empty constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
template< typename LAI >
BlockMatrixView<LAI>::BlockMatrixView()
{}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Sized constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Constructor with a size (number of rows and columns)
template< typename LAI >
BlockMatrixView<LAI>::BlockMatrixView( localIndex const nRows,
                                       localIndex const nCols )
{
  m_matrices.resize( nRows, nCols );
}

// ----------------------------
// Linear Algebra
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Multiply
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Apply the block matrix to a block vector and compute the result.
template< typename LAI >
void BlockMatrixView<LAI>::multiply( BlockVectorView<LAI> const &x,
                                     BlockVectorView<LAI> &b ) const
{
  for( localIndex row = 0 ; row < m_matrices.size( 0 ) ; row++ )
  {
    b.block( row ).zero();
    ParallelVector temp( b.block( row ));
    for( localIndex col = 0 ; col < m_matrices.size( 1 ) ; col++ )
    {
      if( m_matrices[row][col] != nullptr )
      {
        m_matrices[row][col]->multiply( x.block( col ), temp );
        b.block( row ).axpy( 1.0, temp );
      }
    }
  }
}


// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Residual
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute the residual r = b - Ax.
template< typename LAI >
void BlockMatrixView<LAI>::residual( BlockVectorView<LAI> const &x,
                                     BlockVectorView<LAI> const &b,
                                     BlockVectorView<LAI> &r ) const
{
  for( localIndex row = 0 ; row < m_matrices.size( 0 ) ; row++ )
  {
    r.block( row ).copy( b.block( row ));
    ParallelVector temp( b.block( row ));
    for( localIndex col = 0 ; col < m_matrices.size( 1 ) ; col++ )
    {
      if( m_matrices[row][col] != nullptr )
      {
        m_matrices[row][col]->multiply( x.block( col ), temp );
        r.block( row ).axpy( -1.0, temp );
      }
    }
  }
}


// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Scale
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Scale the matrix with the factor <tt>factor</tt>.
template< typename LAI >
void BlockMatrixView<LAI>::scale( real64 const factor )
{
  for( localIndex row = 0 ; row < m_matrices.size( 0 ) ; row++ )
  {
    for( localIndex col = 0 ; col < m_matrices.size( 1 ) ; col++ )
    {
      if( m_matrices[row][col] != nullptr )
      {
        m_matrices[row][col]->scale( factor );
      }
    }
  }
}


// ----------------------------
// Accessors
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get block
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
template< typename LAI >
typename LAI::ParallelMatrix & BlockMatrixView<LAI>::block( localIndex const blockRowIndex,
                                                            localIndex const blockColIndex ) const
{
  return *m_matrices[blockRowIndex][blockColIndex];
}

// ----------------------------
// Setters
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Set block
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
template< typename LAI >
void BlockMatrixView<LAI>::set( localIndex const blockRowIndex,
                                localIndex const blockColIndex,
                                typename LAI::ParallelMatrix &matrix )
{
  m_matrices[blockRowIndex][blockColIndex] = &matrix;
}

// -----------------------
// Explicit Instantiations
// -----------------------
template class BlockMatrixView<TrilinosInterface>;


} // end geosx namespace
