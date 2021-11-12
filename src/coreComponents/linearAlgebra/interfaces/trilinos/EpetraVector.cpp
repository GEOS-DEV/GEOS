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
 * @file EpetraVector.cpp
 */

#include "EpetraVector.hpp"

#include "codingUtilities/Utilities.hpp"
#include "linearAlgebra/interfaces/trilinos/EpetraUtils.hpp"

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <EpetraExt_MultiVectorOut.h>

namespace geosx
{

EpetraVector::EpetraVector()
  : VectorBase(),
    m_vec{}
{}

EpetraVector::EpetraVector( EpetraVector const & src )
  : EpetraVector()
{
  *this = src;
}

EpetraVector::EpetraVector( EpetraVector && src ) noexcept
  : EpetraVector()
{
  *this = std::move( src );
}

EpetraVector & EpetraVector::operator=( EpetraVector const & src )
{
  if( &src != this )
  {
    reset();
    if( src.ready() )
    {
      m_vec = std::make_unique< Epetra_Vector >( *src.m_vec );
    }
  }
  return *this;
}

EpetraVector & EpetraVector::operator=( EpetraVector && src ) noexcept
{
  if( &src != this )
  {
    std::swap( m_vec, src.m_vec );
    std::swap( m_closed, src.m_closed );
  }
  return *this;
}

EpetraVector::~EpetraVector() = default;

bool EpetraVector::created() const
{
  return bool( m_vec);
}

void EpetraVector::create( localIndex const localSize,
                           MPI_Comm const & MPI_PARAM( comm ) )
{
  VectorBase::create( localSize, comm );
  Epetra_Map const map( LvArray::integerConversion< long long >( -1 ),
                        LvArray::integerConversion< int >( localSize ),
                        0,
                        trilinos::EpetraComm( MPI_PARAM( comm ) ) );
  m_data.move( LvArray::MemorySpace::host, false );
  m_vec = std::make_unique< Epetra_Vector >( View, map, m_data.data() );
}

void EpetraVector::set( real64 value )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( m_vec->PutScalar( value ) );
  touch();
}

void EpetraVector::rand( unsigned const seed )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( m_vec->SetSeed( seed ) );
  GEOSX_LAI_CHECK_ERROR( m_vec->Random() );
  touch();
}

void EpetraVector::close()
{
  GEOSX_LAI_ASSERT( !closed() );
  m_data.move( LvArray::MemorySpace::host, false );
  m_closed = true;
}

void EpetraVector::touch()
{
  GEOSX_LAI_ASSERT( ready() );
  m_data.registerTouch( LvArray::MemorySpace::host );
}

void EpetraVector::reset()
{
  VectorBase::reset();
  m_vec.reset();
}

void EpetraVector::scale( real64 const scalingFactor )
{
  GEOSX_LAI_ASSERT( ready() );

  if( !isEqual( scalingFactor, 1.0 ) )
  {
    GEOSX_LAI_CHECK_ERROR( m_vec->Scale( scalingFactor ) );
    touch();
  }
}

void EpetraVector::reciprocal()
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( m_vec->Reciprocal( *m_vec ) );
  touch();
}

real64 EpetraVector::dot( EpetraVector const & vec ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( vec.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), vec.globalSize() );

  real64 tmp;
  GEOSX_LAI_CHECK_ERROR( m_vec->Dot( vec.unwrapped(), &tmp ) );
  return tmp;
}

void EpetraVector::copy( EpetraVector const & x )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), x.globalSize() );

  GEOSX_LAI_CHECK_ERROR( m_vec->Update( 1., x.unwrapped(), 0. ) );
  touch();
}

void EpetraVector::axpy( real64 const alpha,
                         EpetraVector const & x )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), x.globalSize() );

  GEOSX_LAI_CHECK_ERROR( m_vec->Update( alpha, x.unwrapped(), 1. ) );
  touch();
}

void EpetraVector::axpby( real64 const alpha,
                          EpetraVector const & x,
                          real64 const beta )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), x.globalSize() );

  GEOSX_LAI_CHECK_ERROR( m_vec->Update( alpha, x.unwrapped(), beta ) );
  touch();
}

void EpetraVector::pointwiseProduct( EpetraVector const & x,
                                     EpetraVector & y ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  GEOSX_LAI_ASSERT( y.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), x.globalSize() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), y.globalSize() );

  GEOSX_LAI_CHECK_ERROR( ( y.unwrapped() ).Multiply( 1.0, unwrapped(), x.unwrapped(), 0.0 ) );
  y.touch();
}

real64 EpetraVector::norm1() const
{
  GEOSX_LAI_ASSERT( ready() );

  real64 tmp;
  m_vec->Norm1( &tmp );
  return tmp;
}

real64 EpetraVector::norm2() const
{
  GEOSX_LAI_ASSERT( ready() );
  real64 tmp;
  m_vec->Norm2( &tmp );
  return tmp;
}

real64 EpetraVector::normInf() const
{
  GEOSX_LAI_ASSERT( ready() );
  real64 tmp;
  m_vec->NormInf( &tmp );
  return tmp;
}

void EpetraVector::print( std::ostream & os ) const
{
  GEOSX_LAI_ASSERT( ready() );
  m_vec->Print( os );
}

void EpetraVector::write( string const & filename,
                          LAIOutputFormat const format ) const
{
  GEOSX_LAI_ASSERT( ready() );
  switch( format )
  {
    case LAIOutputFormat::MATLAB_ASCII:
    {
      GEOSX_LAI_CHECK_ERROR( EpetraExt::MultiVectorToMatlabFile( filename.c_str(), *m_vec ) );
      break;
    }
    case LAIOutputFormat::MATRIX_MARKET:
    {
      GEOSX_LAI_CHECK_ERROR( EpetraExt::MultiVectorToMatrixMarketFile( filename.c_str(), *m_vec ) );
      break;
    }
    default:
      GEOSX_ERROR( "Unsupported vector output format" );
  }
}

Epetra_Vector const & EpetraVector::unwrapped() const
{
  GEOSX_LAI_ASSERT( created() );
  return *m_vec;
}

Epetra_Vector & EpetraVector::unwrapped()
{
  GEOSX_LAI_ASSERT( created() );
  return *m_vec;
}

globalIndex EpetraVector::globalSize() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_vec->GlobalLength64();
}

localIndex EpetraVector::localSize() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_vec->MyLength();
}

globalIndex EpetraVector::ilower() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_vec->Map().MinMyGID64();
}

globalIndex EpetraVector::iupper() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_vec->Map().MaxMyGID64() + 1;
}

MPI_Comm EpetraVector::getComm() const
{
  GEOSX_LAI_ASSERT( created() );
#ifdef GEOSX_USE_MPI
  return dynamicCast< Epetra_MpiComm const & >( m_vec->Map().Comm() ).Comm();
#else
  return MPI_COMM_GEOSX;
#endif
}

} // end geosx
