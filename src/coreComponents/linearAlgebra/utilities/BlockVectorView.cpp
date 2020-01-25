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
 * @file BlockVectorView.cpp
 */

#include "BlockVectorView.hpp"

#include "linearAlgebra/interfaces/InterfaceTypes.hpp"

namespace geosx
{

// BEGIN_RST_NARRATIVE BlockVectorView.rst
// ==============================
// Block Vector View
// ==============================
// This class contains the block vector view object. Each block vector owns
// pointers to a variable number of ParallelVectors and is templated on the
// LA interface.

// ----------------------------
// Constructors
// ----------------------------


// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Size constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Constructor with a size (number of vectors)
template< typename VECTOR >
BlockVectorView<VECTOR>::BlockVectorView( localIndex const nBlocks )
{
  m_vectors.resize( nBlocks );
}

template< typename VECTOR >
BlockVectorView<VECTOR>::BlockVectorView( BlockVectorView< VECTOR > const & x ) = default;

template< typename VECTOR >
BlockVectorView<VECTOR>::BlockVectorView( BlockVectorView< VECTOR > && x ) = default;

template< typename VECTOR >
BlockVectorView< VECTOR > & BlockVectorView< VECTOR >::operator=( BlockVectorView< VECTOR > const & x ) = default;

template< typename VECTOR >
BlockVectorView< VECTOR > & BlockVectorView< VECTOR >::operator=( BlockVectorView< VECTOR > && x ) = default;

// ----------------------------
// Setters
// ----------------------------

template< typename VECTOR >
void BlockVectorView< VECTOR >::copy( BlockVectorView< VECTOR > const & src )
{
  GEOSX_ERROR_IF( blockSize() != src.blockSize(), "Incompatible block size" );
  for( localIndex i = 0 ; i < m_vectors.size() ; i++ )
  {
    m_vectors[i]->copy( *src.m_vectors[i] );
  }
}

// ----------------------------
// Linear Algebra
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Scale
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Scale the whole block vector by factor.
template< typename VECTOR >
void BlockVectorView<VECTOR>::scale( real64 const factor )
{
  for( localIndex i = 0 ; i < m_vectors.size() ; i++ )
  {
    m_vectors[i]->scale( factor );
  }
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Zero
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Zero out the whole block vector.
template< typename VECTOR >
void BlockVectorView<VECTOR>::zero( )
{
  for( localIndex i = 0 ; i < m_vectors.size() ; i++ )
  {
    m_vectors[i]->zero();
  }
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Rand
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Randomize the whole block vector.
template< typename VECTOR >
void BlockVectorView<VECTOR>::rand( unsigned const seed )
{
  for( localIndex i = 0 ; i < m_vectors.size() ; i++ )
  {
    m_vectors[i]->rand( seed );
  }
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Dot
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Perform a dot product.
template< typename VECTOR >
real64 BlockVectorView<VECTOR>::dot( BlockVectorView<VECTOR> const & src ) const
{
  real64 accum = 0;
  for( localIndex i = 0 ; i < m_vectors.size() ; i++ )
  {
    accum += m_vectors[i]->dot( src.block( i ));
  }
  return accum;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 2-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute the 2 norm of the whole block vector.
template< typename VECTOR >
real64 BlockVectorView<VECTOR>::norm2() const
{
  real64 accum = 0;
  for( localIndex i = 0 ; i < m_vectors.size() ; i++ )
  {
    real64 const temp = m_vectors[i]->norm2();
    accum = accum + temp*temp;
  }
  return std::sqrt( accum );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute the inf norm of the whole block vector.
template< typename VECTOR >
real64 BlockVectorView<VECTOR>::normInf() const
{
  real64 result = 0.0;
  for( localIndex i = 0 ; i < m_vectors.size() ; i++ )
  {
    result = std::fmax( result, m_vectors[i]->normInf() );
  }
  return result;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Axpy
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update the whole block vector with y = alpha*x + y
template< typename VECTOR >
void BlockVectorView<VECTOR>::axpy( real64 const alpha,
                                    BlockVectorView<VECTOR> const & x )
{
  for( localIndex i = 0 ; i < m_vectors.size() ; i++ )
  {
    m_vectors[i]->axpy( alpha, x.block( i ) );
  }
}


// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Axpby
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update the whole block vector with y = alpha*x + beta*y
template< typename VECTOR >
void BlockVectorView<VECTOR>::axpby( real64 const alpha,
                                     BlockVectorView<VECTOR> const &x,
                                     real64 const beta )
{
  for( localIndex i = 0 ; i < m_vectors.size() ; i++ )
  {
    m_vectors[i]->axpby( alpha, x.block( i ), beta );
  }
}


// ----------------------------
// Accessor
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the number of blocks
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for block size (number of blocks)
template< typename VECTOR >
localIndex BlockVectorView<VECTOR>::blockSize() const
{
  return m_vectors.size();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the global number of elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for global size across blocks
template< typename VECTOR >
globalIndex BlockVectorView<VECTOR>::globalSize() const
{
  globalIndex size = 0;
  for( localIndex i = 0 ; i < m_vectors.size() ; i++ )
  {
    size += m_vectors[i]->globalSize();
  }
  return size;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get block.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
template< typename VECTOR >
VECTOR & BlockVectorView<VECTOR>::block( localIndex blockIndex ) const
{
  return *m_vectors[blockIndex];
}

template< typename VECTOR >
void BlockVectorView< VECTOR >::print( std::ostream & os ) const
{
  for( localIndex i = 0 ; i < m_vectors.size() ; i++ )
  {
    os << "Block " << i << " of " << m_vectors.size() << ":" << std::endl;
    os << "=============" << std::endl;
    m_vectors[i]->print( os );
  }
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class BlockVectorView<TrilinosInterface::ParallelVector>;
#endif

#ifdef GEOSX_USE_HYPRE
//template class BlockVectorView<HypreInterface::ParallelVector>;
#endif

#ifdef GEOSX_USE_PETSC
template class BlockVectorView<PetscInterface::ParallelVector>;
#endif

} //end geosx namespace

// END_RST_NARRATIVE
