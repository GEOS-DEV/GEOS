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

BlasVector::BlasVector( BlasVector const & vec )
:
    m_size( vec.m_size ),
    m_values( vec.m_values )
{

}

BlasVector::~BlasVector()
{

}

//-------------------------------------------Shaping/sizing/permuting methods---
void BlasVector::resize( localIndex length )
{
  GEOS_ASSERT_MSG( length > 0, "Vector size must be > 0" );
  m_size = length;
  m_values.resizeDefault( length );
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

void BlasVector::rand()
{
  std::default_random_engine generator;
  std::uniform_real_distribution<real64>  distribution(0.0, 0.1);

  for (localIndex i = 0; i < m_size; ++i)
    m_values[i] = distribution(generator);

  return;
}

void BlasVector::rand(real64 const rangeFrom,
                      real64 const rangeTo)
{
  std::default_random_engine generator;
  std::uniform_real_distribution<real64>  distribution(rangeFrom, rangeTo);

  for (localIndex i = 0; i < m_size; ++i)
    m_values[i] = distribution(generator);

  return;
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

real64 BlasVector::norm1() const
{
  return cblas_dasum( integer_conversion<int>(m_size),
                      m_values.data(),
                      1);
}

real64 BlasVector::norm2() const
{
  return cblas_dnrm2( integer_conversion<int>(m_size),
                      m_values.data(),
                      1);
}

real64 BlasVector::normInf() const
{
  int ind = cblas_idamax( integer_conversion<int>(m_size),
                          m_values.data(),
                          1);
  return std::abs(m_values[ind]);
}

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

void BlasVector::scale(real64 scalarThis)
{
  // Call to BLAS using CBLAS interface
  cblas_dscal(integer_conversion<int>(m_values.size()),
              scalarThis,
              m_values.data(),
              1);
  return;
}

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

void BlasVector::copy( BlasVector const &vec )
{
  GEOS_ASSERT_MSG( vec.getSize() == m_size,
                   "Vector dimensions not compatible for copying" );

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

localIndex BlasVector::getSize() const
{
  return m_size;
}

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

