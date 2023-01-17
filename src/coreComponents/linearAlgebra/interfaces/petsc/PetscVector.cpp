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
 * @file PetscVector.cpp
 */

#include "PetscVector.hpp"

#include "codingUtilities/Utilities.hpp"
#include "linearAlgebra/interfaces/petsc/PetscUtils.hpp"

#include <petscvec.h>

namespace geosx
{

PetscVector::PetscVector()
  : VectorBase(),
  m_vec{}
{}

PetscVector::PetscVector( PetscVector const & src )
  : PetscVector()
{
  *this = src;
}

PetscVector::PetscVector( PetscVector && src ) noexcept
  : PetscVector()
{
  *this = std::move( src );
}

PetscVector & PetscVector::operator=( PetscVector const & src )
{
  if( &src != this )
  {
    reset();
    if( src.created() )
    {
      create( src.localSize(), src.comm() );
      copy( src );
    }
  }
  return *this;
}

PetscVector & PetscVector::operator=( PetscVector && src ) noexcept
{
  if( &src != this )
  {
    m_vec = src.m_vec;
    src.m_vec = nullptr;
    VectorBase::operator=( std::move( src ) );
  }
  return *this;
}

PetscVector::~PetscVector()
{
  reset();
}

void PetscVector::reset()
{
  VectorBase::reset();
  if( m_vec != nullptr )
  {
    GEOSX_LAI_CHECK_ERROR( VecDestroy( &m_vec ) );
    m_vec = nullptr;
  }
}

void PetscVector::create( localIndex const localSize, MPI_Comm const & comm )
{
  VectorBase::create( localSize, comm );
  m_values.move( hostMemorySpace, false );
  GEOSX_LAI_CHECK_ERROR( VecCreateMPIWithArray( comm, 1, localSize, PETSC_DETERMINE, m_values.data(), &m_vec ) );
}

bool PetscVector::created() const
{
  return m_vec != nullptr;
}

void PetscVector::set( real64 const value )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( VecSet( m_vec, value ) );
  touch();
}

void PetscVector::rand( unsigned const seed )
{
  GEOSX_LAI_ASSERT( ready() );

  // create random context
  PetscRandom ran;
  GEOSX_LAI_CHECK_ERROR( PetscRandomCreate( PETSC_COMM_WORLD, &ran ) );
  GEOSX_LAI_CHECK_ERROR( PetscRandomSetInterval( ran, -1.0, 1.0 ) );

  // set random seed
  GEOSX_LAI_CHECK_ERROR( PetscRandomSetSeed( ran, seed ) );
  GEOSX_LAI_CHECK_ERROR( PetscRandomSeed( ran ) );

  // create random vector
  GEOSX_LAI_CHECK_ERROR( VecSetRandom( m_vec, ran ) );
  GEOSX_LAI_CHECK_ERROR( PetscRandomDestroy( &ran ) );

  touch();
}

void PetscVector::close()
{
  GEOSX_LAI_ASSERT( !closed() );
  m_values.move( hostMemorySpace, false );
  m_closed = true;
  GEOSX_LAI_CHECK_ERROR( VecAssemblyBegin( m_vec ) );
  GEOSX_LAI_CHECK_ERROR( VecAssemblyEnd( m_vec ) );
}

void PetscVector::touch()
{
  GEOSX_LAI_ASSERT( ready() );
  m_values.registerTouch( hostMemorySpace );
}

void PetscVector::scale( real64 const scalingFactor )
{
  GEOSX_LAI_ASSERT( ready() );

  if( !isEqual( scalingFactor, 1.0 ) )
  {
    GEOSX_LAI_CHECK_ERROR( VecScale( m_vec, scalingFactor ) );
    touch();
  }
}

void PetscVector::reciprocal()
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( VecReciprocal( m_vec ) );
  touch();
}

real64 PetscVector::dot( PetscVector const & vec ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( vec.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), vec.globalSize() );

  real64 dot;
  GEOSX_LAI_CHECK_ERROR( VecDot( m_vec, vec.m_vec, &dot ) );
  return dot;
}

void PetscVector::copy( PetscVector const & x )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), x.globalSize() );

  GEOSX_LAI_CHECK_ERROR( VecSet( m_vec, 0 ) );
  GEOSX_LAI_CHECK_ERROR( VecAXPY( m_vec, 1.0, x.m_vec ) );
  touch();
}

void PetscVector::axpy( real64 const alpha, PetscVector const & x )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), x.globalSize() );

  if( &x != this )
  {
    GEOSX_LAI_CHECK_ERROR( VecAXPY( m_vec, alpha, x.m_vec ) );
    touch();
  }
  else
  {
    scale( 1.0 + alpha );
  }
}

void PetscVector::axpby( real64 const alpha,
                         PetscVector const & x,
                         real64 const beta )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), x.globalSize() );

  if( &x != this )
  {
    GEOSX_LAI_CHECK_ERROR( VecAXPBY( m_vec, alpha, beta, x.m_vec ) );
    touch();
  }
  else
  {
    scale( alpha + beta );
  }
}

void PetscVector::pointwiseProduct( PetscVector const & x,
                                    PetscVector & y ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  GEOSX_LAI_ASSERT( y.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), x.globalSize() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), y.globalSize() );

  GEOSX_LAI_CHECK_ERROR( VecPointwiseMult( y.m_vec, m_vec, x.m_vec ) );
  y.touch();
}

real64 PetscVector::norm1() const
{
  GEOSX_LAI_ASSERT( ready() );
  real64 result;
  GEOSX_LAI_CHECK_ERROR( VecNorm( m_vec, NORM_1, &result ) );
  return result;
}

real64 PetscVector::norm2() const
{
  GEOSX_LAI_ASSERT( ready() );
  real64 result;
  GEOSX_LAI_CHECK_ERROR( VecNorm( m_vec, NORM_2, &result ) );
  return result;
}

real64 PetscVector::normInf() const
{
  GEOSX_LAI_ASSERT( ready() );
  real64 result;
  GEOSX_LAI_CHECK_ERROR( VecNorm( m_vec, NORM_INFINITY, &result ) );
  return result;
}

globalIndex PetscVector::globalSize() const
{
  GEOSX_LAI_ASSERT( created() );
  PetscInt size;
  GEOSX_LAI_CHECK_ERROR( VecGetSize( m_vec, &size ) );
  return LvArray::integerConversion< globalIndex >( size );
}

localIndex PetscVector::localSize() const
{
  GEOSX_LAI_ASSERT( created() );
  PetscInt size;
  GEOSX_LAI_CHECK_ERROR( VecGetLocalSize( m_vec, &size ) );
  return LvArray::integerConversion< localIndex >( size );
}

globalIndex PetscVector::ilower() const
{
  GEOSX_LAI_ASSERT( created() );
  PetscInt low, high;
  GEOSX_LAI_CHECK_ERROR( VecGetOwnershipRange( m_vec, &low, &high ) );
  return low;
}

globalIndex PetscVector::iupper() const
{
  GEOSX_LAI_ASSERT( created() );
  PetscInt low, high;
  GEOSX_LAI_CHECK_ERROR( VecGetOwnershipRange( m_vec, &low, &high ) );
  return high;
}

void PetscVector::print( std::ostream & os ) const
{
  GEOSX_LAI_ASSERT( ready() );
  if( &os == &std::cout )
  {
    GEOSX_LAI_CHECK_ERROR( VecView( m_vec, PETSC_VIEWER_STDOUT_WORLD ) );
  }
  else if( &os == &std::cerr )
  {
    GEOSX_LAI_CHECK_ERROR( VecView( m_vec, PETSC_VIEWER_STDERR_WORLD ) );
  }
  else
  {
    GEOSX_ERROR( "Output to a generic stream not implemented" );
  }
}

void PetscVector::write( string const & filename,
                         LAIOutputFormat const format ) const
{
  GEOSX_LAI_ASSERT( ready() );
  PetscViewerFormat petscFormat = PETSC_VIEWER_DEFAULT;

  bool useMatrixMarket = false;
  switch( format )
  {
    case LAIOutputFormat::NATIVE_ASCII:
    {
      petscFormat = PETSC_VIEWER_ASCII_COMMON;
      break;
    }
    case LAIOutputFormat::MATLAB_ASCII:
    {
      petscFormat = PETSC_VIEWER_ASCII_MATLAB;
      break;
    }
    case LAIOutputFormat::MATRIX_MARKET:
    {
      useMatrixMarket = true;
      break;
    }
    default:
    {
      GEOSX_ERROR( "Unsupported vector output format" );
    }
  }

  if( !useMatrixMarket )
  {
    PetscViewer viewer;
    GEOSX_LAI_CHECK_ERROR( PetscViewerASCIIOpen( comm(), filename.c_str(), &viewer ) );
    GEOSX_LAI_CHECK_ERROR( PetscViewerPushFormat( viewer, petscFormat ) );
    GEOSX_LAI_CHECK_ERROR( VecView( m_vec, viewer ) );
    GEOSX_LAI_CHECK_ERROR( PetscViewerDestroy( &viewer ) );
  }
  else
  {
    VecScatter scatter;
    Vec globalVec;
    GEOSX_LAI_CHECK_ERROR( VecScatterCreateToAll( m_vec, &scatter, &globalVec ) );
    GEOSX_LAI_CHECK_ERROR( VecScatterBegin( scatter, m_vec, globalVec, INSERT_VALUES, SCATTER_FORWARD ) );
    GEOSX_LAI_CHECK_ERROR( VecScatterEnd( scatter, m_vec, globalVec, INSERT_VALUES, SCATTER_FORWARD ) );
    if( MpiWrapper::commRank( comm() ) == 0 )
    {
      PetscScalar *v;
      GEOSX_LAI_CHECK_ERROR( VecGetArray( globalVec, &v ) );

      FILE * fp = std::fopen( filename.c_str(), "w" );
      fprintf( fp, "%s", "%%MatrixMarket matrix array real general\n" );
      fprintf( fp, "%lld %d\n", globalSize(), 1 );
      for( globalIndex i = 0; i < globalSize(); i++ )
      {
        fprintf( fp, "%.16e\n", v[i] );
      }
      std::fclose( fp );

      GEOSX_LAI_CHECK_ERROR( VecRestoreArray( globalVec, &v ) );
    }
    GEOSX_LAI_CHECK_ERROR( VecScatterDestroy( &scatter ) );
    GEOSX_LAI_CHECK_ERROR( VecDestroy( &globalVec ) );
  }
}

Vec const & PetscVector::unwrapped() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_vec;
}

Vec & PetscVector::unwrapped()
{
  GEOSX_LAI_ASSERT( created() );
  return m_vec;
}

MPI_Comm PetscVector::comm() const
{
  GEOSX_LAI_ASSERT( created() );
  MPI_Comm comm;
  GEOSX_LAI_CHECK_ERROR( PetscObjectGetComm( reinterpret_cast< PetscObject >( m_vec ), &comm ) );
  return comm;
}

} // end geosx
