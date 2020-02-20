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

#include "BlockOperatorView.hpp"

#include "linearAlgebra/interfaces/InterfaceTypes.hpp"

namespace geosx
{

// BEGIN_RST_NARRATIVE BlockOperatorView.rst
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
template< typename VECTOR, typename OPERATOR >
BlockOperatorView< VECTOR, OPERATOR >::BlockOperatorView()
{ }

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Sized constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Constructor with a size (number of rows and columns)
template< typename VECTOR, typename OPERATOR >
BlockOperatorView< VECTOR, OPERATOR >::BlockOperatorView( localIndex const nRows,
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
template< typename VECTOR, typename OPERATOR >
void BlockOperatorView< VECTOR, OPERATOR >::multiply( BlockVectorView< VECTOR > const & x,
                                                      BlockVectorView< VECTOR > & b ) const
{
  for( localIndex row = 0; row < m_matrices.size( 0 ); row++ )
  {
    b.block( row ).zero();
    VECTOR temp( b.block( row ) );
    for( localIndex col = 0; col < m_matrices.size( 1 ); col++ )
    {
      if( m_matrices[row][col] != nullptr )
      {
        m_matrices[row][col]->multiply( x.block( col ), temp );
        b.block( row ).axpy( 1.0, temp );
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
template< typename VECTOR, typename OPERATOR >
OPERATOR & BlockOperatorView< VECTOR, OPERATOR >::block( localIndex const blockRowIndex,
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
template< typename VECTOR, typename OPERATOR >
void BlockOperatorView< VECTOR, OPERATOR >::set( localIndex const blockRowIndex,
                                                 localIndex const blockColIndex,
                                                 OPERATOR & matrix )
{
  m_matrices[blockRowIndex][blockColIndex] = &matrix;
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class BlockOperatorView<TrilinosInterface::ParallelVector>;
#endif

#ifdef GEOSX_USE_HYPRE
template class BlockOperatorView<HypreInterface::ParallelVector>;
#endif

#ifdef GEOSX_USE_PETSC
template class BlockOperatorView<PetscInterface::ParallelVector>;
#endif

} // end geosx namespace
