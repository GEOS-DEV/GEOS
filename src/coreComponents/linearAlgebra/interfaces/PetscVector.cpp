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
 * @file PetscVector.cpp
 */

#include "PetscVector.hpp"
#include "linearAlgebra/interfaces/PetscUtils.hpp"
#include <petscvec.h>

namespace geosx
{

// Check matching requirements on index/value types between GEOSX and PETSc

static_assert( sizeof( PetscInt ) == sizeof( globalIndex ),
               "PetscInt and geosx::globalIndex must have the same size" );

static_assert( std::is_signed< PetscInt >::value == std::is_signed< globalIndex >::value,
               "PetscInt and geoex::globalIndex must both be signed or unsigned");

static_assert( std::is_same< PetscScalar, real64 >::value,
               "PetscScalar and geosx::real64 must be the same type" );

PetscVector::PetscVector()
: VectorBase(),
  m_vec{}
{}

PetscVector::PetscVector( PetscVector const & vec )
: PetscVector()
{
  VecDuplicate( vec.m_vec, &m_vec );
  VecCopy( vec.m_vec, m_vec );
}

PetscVector::PetscVector( PetscVector && src ) noexcept
: PetscVector()
{
  std::swap( m_vec, src.m_vec );
}

PetscVector & PetscVector::operator=( PetscVector const & src )
{
  GEOSX_LAI_ASSERT( &src != this );
  GEOSX_LAI_ASSERT( src.created() );
  VecDuplicate( src.m_vec, &m_vec );
  VecCopy( src.m_vec, m_vec );
  return *this;
}

PetscVector & PetscVector::operator=( PetscVector && src ) noexcept
{
  GEOSX_LAI_ASSERT( &src != this );
  std::swap( m_vec, src.m_vec );
  return *this;
}

PetscVector::~PetscVector()
{
  reset();
}

void PetscVector::reset()
{
  GEOSX_LAI_VECTOR_STATUS( closed() );
  if( m_vec != nullptr )
  {
    VecDestroy( &m_vec );
    m_vec = nullptr;
  }
}

void PetscVector::createWithLocalSize( localIndex const localSize, MPI_Comm const & comm )
{
  GEOSX_LAI_VECTOR_STATUS( closed() );
  GEOSX_LAI_ASSERT_GE( localSize, 0 );
  reset();
  VecCreate( comm, &m_vec );
  VecSetType( m_vec, VECMPI );
  VecSetSizes( m_vec, localSize, PETSC_DETERMINE);
}

void PetscVector::createWithGlobalSize( globalIndex const globalSize, MPI_Comm const & comm )
{
  GEOSX_LAI_VECTOR_STATUS( closed() );
  GEOSX_LAI_ASSERT_GE( globalSize, 0 );
  reset();
  VecCreate( comm, &m_vec );
  VecSetType( m_vec, VECMPI );
  VecSetSizes( m_vec, PETSC_DECIDE, globalSize);
}

void PetscVector::create( arraySlice1d<real64 const> const & localValues, MPI_Comm const & comm )
{
  GEOSX_LAI_VECTOR_STATUS( closed() );
  reset();
  PetscInt const size = localValues.size();
  PetscScalar * values;

  VecCreate( comm, &m_vec );
  VecSetType( m_vec, VECMPI );
  VecSetSizes( m_vec, size, PETSC_DETERMINE );
  VecGetArray( m_vec, &values );

  // set vector values
  for (int i = 0; i < size; i++)
  {
    values[i] = localValues[i];
  }

  VecRestoreArray( m_vec, &values );
}

bool PetscVector::created() const
{
  return m_vec != nullptr;
}

void PetscVector::set( globalIndex const globalRow,
                       real64 const value )
{
  GEOSX_LAI_VECTOR_STATUS( !closed() );
  VecSetValue( m_vec, globalRow, value, INSERT_VALUES );
}

void PetscVector::add( globalIndex const globalRow,
                       real64 const value )
{
  GEOSX_LAI_VECTOR_STATUS( !closed() );
  VecSetValue( m_vec, globalRow, value, ADD_VALUES );
}

void PetscVector::set( globalIndex const * globalIndices,
                       real64 const * values,
                       localIndex size )
{
  GEOSX_LAI_VECTOR_STATUS( !closed() );
  VecSetValues( m_vec, size, toPetscInt( globalIndices), values, INSERT_VALUES );
}

void PetscVector::add( globalIndex const * globalIndices,
                       real64 const * values,
                       localIndex size )
{
  GEOSX_LAI_VECTOR_STATUS( !closed() );
  VecSetValues( m_vec, size, toPetscInt( globalIndices), values, ADD_VALUES );
}

void PetscVector::set( arraySlice1d<globalIndex const> const & globalIndices,
                       arraySlice1d<real64 const> const & values )
{
  GEOSX_LAI_VECTOR_STATUS( !closed() );
  VecSetValues( m_vec, values.size(), toPetscInt( globalIndices.data()), values.data(), INSERT_VALUES );
}

void PetscVector::add( arraySlice1d<globalIndex const> const & globalIndices,
                       arraySlice1d<real64 const> const & values )
{
  GEOSX_LAI_VECTOR_STATUS( !closed() );
  VecSetValues( m_vec, values.size(), toPetscInt( globalIndices.data()), values.data(), ADD_VALUES );
}

void PetscVector::set( real64 const value )
{
  GEOSX_LAI_VECTOR_STATUS( ready() );
  VecSet( m_vec, value );
}

void PetscVector::zero()
{
  set( 0.0 );
}

void PetscVector::rand( unsigned const seed )
{
  GEOSX_LAI_VECTOR_STATUS( ready() );

  // create random context
  PetscRandom ran;
  PetscRandomCreate( PETSC_COMM_WORLD, &ran );

  // set random seed
  PetscRandomSetSeed( ran, seed );
  PetscRandomSeed( ran );

  // create random vector
  VecSetRandom( m_vec, ran );
  PetscRandomDestroy( &ran );
}

void PetscVector::open()
{
  GEOSX_LAI_VECTOR_STATUS( ready() );
  m_closed = false;
}

void PetscVector::close()
{
  GEOSX_LAI_VECTOR_STATUS( !closed() );
  // assemble the vector after setting values
  VecAssemblyBegin( m_vec );
  VecAssemblyEnd( m_vec );
  m_closed = true;
}

void PetscVector::scale( real64 const scalingFactor )
{
  GEOSX_LAI_VECTOR_STATUS( ready() );
  VecScale( m_vec, scalingFactor );
}

real64 PetscVector::dot( PetscVector const & vec ) const
{
  GEOSX_LAI_VECTOR_STATUS( ready() );
  GEOSX_LAI_ASSERT( vec.ready() );
  real64 dot;
  VecDot( m_vec, vec.m_vec, &dot );
  return dot;
}

void PetscVector::copy( PetscVector const & x )
{
  GEOSX_LAI_VECTOR_STATUS( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  //VecCopy( x.m_vec, m_vec );
  VecSet( m_vec, 0 );
  VecAXPY( m_vec, 1.0, x.m_vec );
}

void PetscVector::axpy( real64 const alpha, PetscVector const & x )
{
  GEOSX_LAI_VECTOR_STATUS( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  VecAXPY( m_vec, alpha, x.m_vec );
}

void PetscVector::axpby( real64 const alpha,
                         PetscVector const & x,
                         real64 const beta )
{
  GEOSX_LAI_VECTOR_STATUS( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  VecScale( m_vec, beta );
  VecAXPY( m_vec, alpha, x.m_vec );
}

real64 PetscVector::norm1() const
{
  GEOSX_LAI_VECTOR_STATUS( ready() );
  real64 result;
  VecNorm( m_vec, NORM_1, &result );
  return result;
}

real64 PetscVector::norm2() const
{
  GEOSX_LAI_VECTOR_STATUS( ready() );
  real64 result;
  VecNorm( m_vec, NORM_2, &result );
  return result;
}

real64 PetscVector::normInf() const
{
  GEOSX_LAI_VECTOR_STATUS( ready() );
  real64 result;
  VecNorm( m_vec, NORM_INFINITY, &result );
  return result;
}

void PetscVector::print( std::ostream & os ) const
{
  GEOSX_LAI_VECTOR_STATUS( ready() );
  if( &os == &std::cout )
  {
    VecView( m_vec, PETSC_VIEWER_STDOUT_WORLD );
  }
  else if( &os == &std::cerr )
  {
    VecView( m_vec, PETSC_VIEWER_STDERR_WORLD );
  }
  else
  {
    GEOSX_ERROR( "Output to a generic stream not implemented" );
  }
}

void PetscVector::write( string const & filename,
                         LAIOutputFormat const format ) const
{
  GEOSX_LAI_VECTOR_STATUS( ready() );
  PetscViewer viewer;
  PetscViewerASCIIOpen( getComm(), filename.c_str(), &viewer );
  PetscViewerFormat petscFormat = PETSC_VIEWER_DEFAULT;

  switch( format )
  {
    case LAIOutputFormat::NATIVE_ASCII:
      petscFormat = PETSC_VIEWER_ASCII_COMMON;
      break;
    case LAIOutputFormat::MATLAB_ASCII:
      petscFormat = PETSC_VIEWER_ASCII_MATLAB;
      break;
    case LAIOutputFormat::MATRIX_MARKET:
      petscFormat = PETSC_VIEWER_ASCII_MATRIXMARKET;
      break;
    default:
      GEOSX_ERROR( "Unsupported vector output format" );
  }

  PetscViewerPushFormat( viewer, petscFormat );
  VecView( m_vec, viewer );
  PetscViewerDestroy( &viewer );
}

real64 PetscVector::get( globalIndex globalRow ) const
{
  GEOSX_LAI_VECTOR_STATUS( ready() );
  real64 value;
  VecGetValues( m_vec, 1, toPetscInt( &globalRow ), &value );
  return value;
}

void PetscVector::get( arraySlice1d<globalIndex const> const & globalIndices,
                       arraySlice1d<real64> const & values ) const
{
  GEOSX_LAI_VECTOR_STATUS( ready() );
  VecGetValues( m_vec, globalIndices.size(), toPetscInt( globalIndices.data() ), values.data() );
}

Vec const & PetscVector::unwrapped() const
{
  GEOSX_LAI_VECTOR_STATUS( created() );
  return m_vec;
}

Vec & PetscVector::unwrapped()
{
  GEOSX_LAI_VECTOR_STATUS( created() );
  return m_vec;
}

MPI_Comm PetscVector::getComm() const
{
  GEOSX_LAI_VECTOR_STATUS( created() );
  MPI_Comm comm;
  PetscObjectGetComm( reinterpret_cast<PetscObject>( m_vec ), &comm );
  return comm;
}

globalIndex PetscVector::globalSize() const
{
  GEOSX_LAI_VECTOR_STATUS( created() );
  PetscInt size;
  VecGetSize( m_vec, &size );
  return size;
}

localIndex PetscVector::localSize() const
{
  GEOSX_LAI_VECTOR_STATUS( created() );
  PetscInt size;
  VecGetLocalSize( m_vec, &size );
  return size;
}

localIndex PetscVector::getLocalRowID( globalIndex const index ) const
{
  GEOSX_LAI_VECTOR_STATUS( created() );
  PetscInt low, high;
  VecGetOwnershipRange( m_vec, &low, &high );
  return (index >= low && index < high) ? integer_conversion< localIndex >( index - low ) : -1;
}

globalIndex PetscVector::getGlobalRowID( localIndex const localRow ) const
{
  GEOSX_LAI_VECTOR_STATUS( created() );
  GEOSX_LAI_ASSERT_GE( localRow, 0 );
  GEOSX_LAI_ASSERT_GT( localSize(), localRow );
  return ilower() + localRow;
}

real64 const * PetscVector::extractLocalVector() const
{
  GEOSX_LAI_VECTOR_STATUS( ready() );
  PetscScalar * avec;
  VecGetArray( m_vec, &avec );
  real64 const * const localVector = avec;
  VecRestoreArray( m_vec, &avec );
  return localVector;
}

real64 * PetscVector::extractLocalVector()
{
  GEOSX_LAI_VECTOR_STATUS( ready() );
  PetscScalar * avec;
  VecGetArray( m_vec, &avec );
  real64 * const localVector = avec;
  VecRestoreArray( m_vec, &avec );
  return localVector;
}

globalIndex PetscVector::ilower() const
{
  GEOSX_LAI_VECTOR_STATUS( created() );
  PetscInt low, high;
  VecGetOwnershipRange( m_vec, &low, &high );
  return low;
}

globalIndex PetscVector::iupper() const
{
  GEOSX_LAI_VECTOR_STATUS( created() );
  PetscInt low, high;
  VecGetOwnershipRange( m_vec, &low, &high );
  return high;
}

} // end geosx

