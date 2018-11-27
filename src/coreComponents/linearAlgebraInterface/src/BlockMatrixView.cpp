/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
BlockMatrixView<LAI>::BlockMatrixView( typename LAI::lid const nRows,
                                       typename LAI::lid const nCols )
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
  for( typename LAI::lid row = 0 ; row < m_matrices.size( 0 ) ; row++ )
  {
    b.block(row).zero();
    ParallelVector temp(b.block(row));
    for( typename LAI::lid col = 0 ; col < m_matrices.size( 1 ) ; col++ )
    {
      if( m_matrices[row][col] != nullptr )
      {
        m_matrices[row][col]->multiply( x.block(col), temp );
        b.block(row).axpy(1.0,temp);
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
  for( typename LAI::lid row = 0 ; row < m_matrices.size( 0 ) ; row++ )
  {
    r.block(row).copy(b.block(row));
    ParallelVector temp(b.block(row));
    for( typename LAI::lid col = 0 ; col < m_matrices.size( 1 ) ; col++ )
    {
      if( m_matrices[row][col] != nullptr )
      {
        m_matrices[row][col]->multiply( x.block( col ), temp );
        r.block(row).axpy(-1.0,temp);
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
  for( typename LAI::lid row = 0 ; row < m_matrices.size( 0 ) ; row++ )
  {
    for( typename LAI::lid col = 0 ; col < m_matrices.size( 1 ) ; col++ )
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
typename LAI::ParallelMatrix & BlockMatrixView<LAI>::block( typename LAI::lid const blockRowIndex,
                                                            typename LAI::lid const blockColIndex ) const
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
void BlockMatrixView<LAI>::set( typename LAI::lid const blockRowIndex,
                                typename LAI::lid const blockColIndex,
                                typename LAI::ParallelMatrix &matrix )
{
  m_matrices[blockRowIndex][blockColIndex] = &matrix;
}

// -----------------------
// Explicit Instantiations
// -----------------------
template class BlockMatrixView<TrilinosInterface>;


} // end geosx namespace

