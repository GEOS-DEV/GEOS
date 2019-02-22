/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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
 * @file BlasVector.cpp
 */

// Include the corresponding header file.
#include "BlasVector.hpp"

// Put everything under the geosx namespace.
namespace geosx
{

//----------------------------------------------Constructor/destructor methods---
BlasVector::BlasVector()
{

}

BlasVector::BlasVector( localIndex length )
{
  this->resize( length);
}

BlasVector::BlasVector( BlasVector const & src )
:
    m_size( src.m_size ),
    m_values( src.m_values )
{

}

BlasVector::~BlasVector()
{

}

//-----------------------------------------------------Shaping/sizing methods---
void BlasVector::resize( localIndex length )
{
  GEOS_ASSERT_MSG( length > 0, "Vector size must be > 0" );
  m_size = length;
  m_values.resizeDefault( length );
  this->zero();
}

void BlasVector::zero()
{
  m_values = 0;
}

BlasVector &BlasVector::operator=(double value)
{
   m_values = value;
   return *this;
}

void BlasVector::permute( array1d<int> permVector,
                          const bool forwardPermutation)
{
  // --- check that permutation vector has correct size
  GEOS_ASSERT_MSG( permVector.size() == m_size,
                   "Permutation vector not consistent with vector size" );
  LAPACKE_dlapmr( LAPACK_COL_MAJOR,
                  forwardPermutation,
                  integer_conversion<lapack_int>( m_size ),
                  1,
                  m_values.data(),
                  integer_conversion<lapack_int>( m_size ),
                  permVector.data());

  return;
}

//-------------------------------------------------------Mathematical methods---

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 1-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the one-norm of the vector.
real64 BlasVector::norm1() const
{
  return cblas_dasum( integer_conversion<int>(m_size),
                      m_values.data(),
                      1);
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 2-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the two-norm of the vector.
real64 BlasVector::norm2() const
{
  return cblas_dnrm2( integer_conversion<int>(m_size),
                      m_values.data(),
                      1);
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Infinity-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the infinity-norm of the vector.
real64 BlasVector::normInf() const
{
  real64 max = 0.0;
  for (localIndex i = 0; i < m_size; i++)
  {
     max = std::max(std::abs(m_values[i]), max);
  }
  return max;
}

// matrix-matrix sum (optional scaling)
void BlasVector::vectorAdd( BlasVector const & Vec,
                            real64 const scalarVec )
{

  GEOS_ASSERT_MSG( Vec.getSize() == m_size,
                   "Vector dimensions not compatible for sum" );

  cblas_daxpy( static_cast<int>( m_values.size() ),
               scalarVec,
               Vec.m_values.data(),
               1,
               m_values.data(),
               1 );
  return;
}

// in-place scalar-vector product
void BlasVector::scale(real64 scalarThis)
{
  // Call to BLAS using CBLAS interface
  cblas_dscal(integer_conversion<int>(m_values.size()),
              scalarThis,
              m_values.data(),
              1);
  return;
}

// Dot product with the vector vec.
real64 BlasVector::dot( BlasVector const &vec )
{
  GEOS_ASSERT_MSG( vec.getSize() == m_size,
                   "Vector dimensions not compatible for dot product" );

  // Call to BLAS using CBLAS interface
  return cblas_ddot(integer_conversion<int>(m_values.size()),
             m_values.data(),
             1,
             vec.m_values.data(),
             1);
}

// in-place scalar-vector product
void BlasVector::copy( BlasVector const &vec )
{
  // Call to BLAS using CBLAS interface
  cblas_dcopy(integer_conversion<int>(m_values.size()),
              vec.m_values.data(),
              1,
              m_values.data(),
              1);
  return;
}

//
//------------------------------------------------------Data Accessor methods---

// Returns number of matrix rows.
localIndex BlasVector::getSize() const
{
  return m_size;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get unwrapped pointer
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get pointer to raw Epetra object, with const and non-const versions.
array1d<real64> const * BlasVector::getValues() const
{
  return &m_values;
}

array1d<real64>* BlasVector::getValues()
{
  return &m_values;
}

//----------------------------------------------------------------I/O methods---
// matrix nice output
void BlasVector::print()
{
  for( localIndex i = 0 ; i < m_size ; ++i )
  {
    printf( "%10.2e ", m_values[i] );
    printf( "\n" );
  }
}

} // end geosx namespace

