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
 * @file BlockVectorView.cpp
 */

#include "BlockVectorView.hpp"

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
// Empty constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
template< typename LAI >
BlockVectorView<LAI>::BlockVectorView()
{}


// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Size constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Constructor with a size (number of vectors)
template< typename LAI >
BlockVectorView<LAI>::BlockVectorView( typename LAI::lid const nBlocks )
{
  m_vectors.resize( nBlocks );
}


// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// (Shallow) copy constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// DISABLED: The view usage can be non-intuitive, and most users will
// probably expect they are creating a deep copy.  People can use the
// explicit shallowCopy() function below if really needed.

/*
template< typename LAI >
BlockVectorView<LAI>::BlockVectorView( BlockVectorView<LAI> const &src )
{
  shallowCopy(src);
}
*/

// ----------------------------
// Setters
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Shallow copy
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Copy the input block vector pointers
template< typename LAI >
void BlockVectorView<LAI>::shallowCopy( BlockVectorView<LAI> const &src )
{
  m_vectors.resize( src.blockSize()  );
  for( typename LAI::lid i = 0 ; i < blockSize() ; i++ )
    m_vectors[i] = &(src.block(i));
}


// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Setter
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Setter for block.
template< typename LAI >
void BlockVectorView<LAI>::set( typename LAI::lid const blockIndex,
                                ParallelVector &vector )
{
  m_vectors[blockIndex] = &vector;
}


// ----------------------------
// Linear Algebra
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Scale
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Scale the whole block vector by factor.
template< typename LAI >
void BlockVectorView<LAI>::scale( real64 const factor )
{
  for( typename LAI::lid i = 0 ; i < blockSize() ; i++ ) 
    m_vectors[i]->scale( factor );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Dot
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Perform a dot product.
template< typename LAI >
real64 BlockVectorView<LAI>::dot( BlockVectorView<LAI> const &src) const
{
  real64 accum = 0;
  for( typename LAI::lid i = 0 ; i < blockSize() ; i++ )
  {
    real64 temp = m_vectors[i]->dot(src.block(i));
    accum = accum + temp;
  }
  return accum;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 2-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute the 2 norm of the whole block vector.
template< typename LAI >
real64 BlockVectorView<LAI>::norm2() const
{
  real64 accum = 0;
  for( typename LAI::lid i = 0 ; i < blockSize() ; i++ )
  {
    real64 temp = m_vectors[i]->norm2();
    accum = accum + temp*temp;
  }
  return std::sqrt( accum );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute the inf norm of the whole block vector.
template< typename LAI >
real64 BlockVectorView<LAI>::normInf() const
{
  real64 temp = 0;
  real64 result = 0;
  for( typename LAI::lid i = 0 ; i < blockSize() ; i++ )
  {
    temp = m_vectors[i]->normInf();
    if( temp > result )
      result = temp;
  }
  return result;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Axpy
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update the whole block vector with y = alpha*x + y
template< typename LAI >
void BlockVectorView<LAI>::axpy( real64 const alpha,
                                 BlockVectorView<LAI> const &x )
{
  for( typename LAI::lid i = 0 ; i < blockSize() ; i++ )
    m_vectors[i]->axpby( alpha, x.block(i), 1. );
}


// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Axpby
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update the whole block vector with y = alpha*x + beta*y
template< typename LAI >
void BlockVectorView<LAI>::axpby( real64 const alpha,
                                  BlockVectorView<LAI> const &x,
                                  real64 const beta )
{
  for( typename LAI::lid i = 0 ; i < blockSize() ; i++ )
    m_vectors[i]->axpby( alpha, x.block(i), beta );
}


// ----------------------------
// Accessor
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the number of blocks
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for block size (number of blocks)
template< typename LAI >
typename LAI::lid BlockVectorView<LAI>::blockSize() const
{
  return static_cast<typename LAI::lid>(m_vectors.size()); //TODO integer_conversion
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the global number of elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for global size across blocks
template< typename LAI >
typename LAI::gid BlockVectorView<LAI>::globalSize(  ) const
{
  typename LAI::gid size = 0;
  for( typename LAI::lid i = 0 ; i < blockSize() ; i++ )
    size = m_vectors[i]->globalSize() + size;
  return size;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get block.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
template< typename LAI >
typename LAI::ParallelVector & BlockVectorView<LAI>::block(typename LAI::lid blockIndex ) const
{
  return *m_vectors[blockIndex];
}


// -----------------------
// Explicit Instantiations
// -----------------------
template class BlockVectorView<TrilinosInterface>;

} //end geosx namespace

// END_RST_NARRATIVE


