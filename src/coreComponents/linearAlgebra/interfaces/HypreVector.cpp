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
 * @file HypreVector.cpp
 */

// BEGIN_RST_NARRATIVE HypreVector.rst
// ==============================
// Hypre-based Vector Object
// ==============================
// ... .
// Include the corresponding header file.
#include "HypreVector.hpp"

// Include required Hypre headers
#include "HYPRE.h"
#include "_hypre_IJ_mv.h"
#include "_hypre_parcsr_mv.h"

// Put everything under the geosx namespace.
namespace geosx
{

inline HYPRE_Int * toHYPRE_Int( globalIndex * const index )
{
  return reinterpret_cast<HYPRE_Int*>(index);
}

inline HYPRE_Int const * toHYPRE_Int( globalIndex const * const index )
{
  return reinterpret_cast<HYPRE_Int const*>(index);
}

// ----------------------------
// Constructors
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Empty constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Construct as an empty vector

HypreVector::HypreVector()
{
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Copy constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a pointer to the raw vector.  The data from the input vector is
// copied to a new memory location. Checks if the vector to be copied is empty.

HypreVector::HypreVector( HypreVector const &src )
{
  GEOS_ERROR_IF( src.unwrappedPointer() == nullptr,
                 "source vector appears to be empty" );
  // Note: every vector is created initialized to 0 and then 'closed',
  //       i.e. made ready to use

  this->create( src );

}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Destructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
HypreVector::~HypreVector()
{
  if( m_ij_vector )
  {
    HYPRE_IJVectorDestroy( m_ij_vector );
  }
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create with known size
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// There are two variants of this function.  In the first, the user knows
// the local size and wants the global size to be the sum of each processor's
// contributions.  In the second, the user knows the global size and wants a
// near-even distribution of elements across processors. All processors
// get the same number of elements, except proc 0 which gets any remainder
// elements necessary when the number of processors does not divide evenly
// into the vector length.
//
// Vector are created initialized to 0.

void HypreVector::create( HypreVector const & src )
{
  GEOS_ERROR_IF( src.unwrappedPointer() == nullptr,
                 "source vector appears to be empty" );
  // Note: every vector is created initialized to 0 and then 'closed',
  //       i.e. made ready to use

  MPI_Comm comm = hypre_IJVectorComm( *src.unwrappedPointer() );
  HYPRE_Int jlower, jupper;
  HYPRE_Int objectType;

  HYPRE_IJVectorGetLocalRange( *src.unwrappedPointer(),
                               &jlower,
                               &jupper );
  HYPRE_IJVectorGetObjectType( *src.unwrappedPointer(),
                               &objectType );

  HYPRE_IJVectorCreate( comm,
                        jlower,
                        jupper,
                        &m_ij_vector );
  HYPRE_IJVectorSetObjectType( m_ij_vector,
                               objectType );
  HYPRE_IJVectorInitialize( m_ij_vector );

  HYPRE_IJVectorAssemble( m_ij_vector );
  HYPRE_IJVectorGetObject( m_ij_vector,
                           (void **) &m_par_vector );

  this->copy( src );
}

void HypreVector::createWithLocalSize( localIndex const localSize,
                                       MPI_Comm const & comm )
{
  GEOS_ERROR_IF( localSize < 1,
                 "local size is lower than 1" );

  int this_mpi_process;
  int n_mpi_process;
  MPI_Comm_rank( comm, &this_mpi_process );
  MPI_Comm_size( comm, &n_mpi_process );

  array1d<int> localSizeArray( n_mpi_process );
  int tmp = integer_conversion<int>( localSize );

  MPI_Allgather( &tmp,
                 1,
                 MPI_INT,
                 localSizeArray.data(),
                 1,
                 MPI_INT,
                 comm );

  HYPRE_Int jLower, jUpper;

  jLower = 0;
  for( int i = 0 ; i < this_mpi_process ; ++i )
  {
    jLower += integer_conversion<HYPRE_Int>( localSizeArray[i] );
  }
  jUpper = jLower + integer_conversion<HYPRE_Int>( localSize ) - 1;

  HYPRE_IJVectorCreate( comm,
                        jLower,
                        jUpper,
                        &m_ij_vector );
  HYPRE_IJVectorSetObjectType( m_ij_vector,
                               HYPRE_PARCSR );
  HYPRE_IJVectorInitialize( m_ij_vector );

  hypre_IJVectorZeroValues( m_ij_vector );

  HYPRE_IJVectorAssemble( m_ij_vector );
  HYPRE_IJVectorGetObject( m_ij_vector,
                           (void **) &m_par_vector );

}

void HypreVector::createWithGlobalSize( globalIndex const globalSize,
                                        MPI_Comm const & comm )
{
  int this_mpi_process;
  int n_mpi_process;
  MPI_Comm_rank( comm, &this_mpi_process );
  MPI_Comm_size( comm, &n_mpi_process );

  HYPRE_Int localSize = integer_conversion<HYPRE_Int>( globalSize )
                      / integer_conversion<HYPRE_Int>( n_mpi_process );
  HYPRE_Int residual = integer_conversion<HYPRE_Int>( globalSize )
                     % integer_conversion<HYPRE_Int>( n_mpi_process );

  GEOS_ERROR_IF( localSize < 1,
                 "local size is lower than 1: less that one processor per component" );

  HYPRE_Int jLower;
  HYPRE_Int jUpper;
  HYPRE_Int rank = integer_conversion<HYPRE_Int>( this_mpi_process );

  if( this_mpi_process == 0 )
  {
    jLower = 0;
    jUpper = localSize + residual - 1;
  }
  else
  {
    jLower = rank * localSize + residual;
    jUpper = jLower + localSize - 1;
  }

  HYPRE_IJVectorCreate( comm,
                        jLower,
                        jUpper,
                        &m_ij_vector );
  HYPRE_IJVectorSetObjectType( m_ij_vector, HYPRE_PARCSR );
  HYPRE_IJVectorInitialize( m_ij_vector );

  hypre_IJVectorZeroValues( m_ij_vector );

  HYPRE_IJVectorAssemble( m_ij_vector );
  HYPRE_IJVectorGetObject( m_ij_vector, (void **) &m_par_vector );

}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create from array
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a vector from local array data.  The global vector contains
// local arrays stitched together.

//TODO: add integer_conversion

void HypreVector::create( array1d<real64> const & localValues,
                          MPI_Comm const & comm )
{
  HYPRE_Int localSize = integer_conversion<HYPRE_Int>( localValues.size() );

  GEOS_ERROR_IF( localSize < 1, "local size is lower than 1" );

  int this_mpi_process;
  int n_mpi_process;
  MPI_Comm_rank( comm, &this_mpi_process );
  MPI_Comm_size( comm, &n_mpi_process );

  array1d<int> localSizeArray( n_mpi_process );
  int tmp = integer_conversion<int>( localSize );

  MPI_Allgather( &tmp,
                 1,
                 MPI_INT,
                 localSizeArray.data(),
                 1,
                 MPI_INT,
                 comm );

  HYPRE_Int jLower, jUpper;

  jLower = 0;
  for( int i = 0 ; i < this_mpi_process ; ++i )
  {
    jLower += integer_conversion<HYPRE_Int>( localSizeArray[i] );
  }
  jUpper = jLower + integer_conversion<HYPRE_Int>( localSize ) - 1;

  HYPRE_IJVectorCreate( comm,
                        jLower,
                        jUpper,
                        &m_ij_vector );
  HYPRE_IJVectorSetObjectType( m_ij_vector, HYPRE_PARCSR );
  HYPRE_IJVectorInitialize( m_ij_vector );

  HYPRE_IJVectorAssemble( m_ij_vector );
  HYPRE_IJVectorGetObject( m_ij_vector, (void **) &m_par_vector );

  double * local_data = this->extractLocalVector();

  for( localIndex i = 0 ; i < localValues.size() ; ++i )
    local_data[i] = localValues[i];

}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Add/Set value(s)
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Add/set entries in the vector.

// single element options

void HypreVector::set( globalIndex const globalRow,
                       real64 const value )
{
  GEOS_ERROR_IF( m_ij_vector == nullptr,
                 "vector appears to be empty (not created)" );
  GEOS_ASSERT_MSG( this->ilower() <= globalRow &&
                   globalRow < this->iupper(),
                   "HypreVector, it is not possible to set values on other processors");
  HYPRE_IJVectorSetValues( m_ij_vector,
                           1,
                           &globalRow,
                           &value );
}

void HypreVector::add( globalIndex const globalRow,
                       real64 const value )
{
  GEOS_ERROR_IF( m_ij_vector == nullptr,
                 "vector appears to be empty (not created)" );
  HYPRE_IJVectorAddToValues( m_ij_vector,
                             1,
                             &globalRow,
                             &value );
}

// n-element, c-style options

void HypreVector::set( globalIndex const * globalIndices,
                       real64 const * values,
                       localIndex size )
{
  GEOS_ERROR_IF( m_ij_vector == nullptr,
                 "vector appears to be empty (not created)" );
  GEOS_ASSERT_MSG( this->ilower() <= *std::min_element(globalIndices, globalIndices + size) &&
                   *std::max_element(globalIndices, globalIndices + size) < this->iupper(),
                   "HypreVector, it is not possible to set values on other processors");
  HYPRE_IJVectorSetValues( m_ij_vector,
                           integer_conversion<HYPRE_Int>( size ),
                           toHYPRE_Int( globalIndices ),
                           values );
}

void HypreVector::add( globalIndex const * globalIndices,
                       real64 const * values,
                       localIndex size )
{
  GEOS_ERROR_IF( m_ij_vector == nullptr,
                 "vector appears to be empty (not created)" );
  HYPRE_IJVectorAddToValues( m_ij_vector,
                             integer_conversion<HYPRE_Int>( size ),
                             toHYPRE_Int( globalIndices ),
                             values );
}

// n-element, array1d options

void HypreVector::set( array1d<globalIndex> const & globalIndices,
                       array1d<real64> const & values )
{
  GEOS_ERROR_IF( m_ij_vector == nullptr,
                 "vector appears to be empty (not created)" );
  GEOS_ASSERT_MSG( this->ilower() <= *std::min_element(globalIndices.data(),
                                                       globalIndices.data() + globalIndices.size() ) &&
                   *std::max_element(globalIndices.data(),
                                     globalIndices.data() + globalIndices.size()) < this->iupper(),
                   "HypreVector, it is not possible to set values on other processors");
  HYPRE_IJVectorSetValues( m_ij_vector,
                           integer_conversion<HYPRE_Int>( values.size() ),
                           toHYPRE_Int( globalIndices.data() ),
                           values.data() );
}

void HypreVector::add( array1d<globalIndex> const & globalIndices,
                       array1d<real64> const & values )
{
  GEOS_ERROR_IF( m_ij_vector == nullptr,
                 "vector appears to be empty (not created)" );
  HYPRE_IJVectorAddToValues( m_ij_vector,
                             integer_conversion<HYPRE_Int>( values.size() ),
                             toHYPRE_Int( globalIndices.data() ),
                             values.data() );
}

//additional options:

void HypreVector::set( real64 value )
{
  GEOS_ERROR_IF( m_ij_vector == nullptr,
                 "vector appears to be empty (not created)" );
  HYPRE_ParVectorSetConstantValues( m_par_vector, value );
}

void HypreVector::zero()
{
  GEOS_ERROR_IF( m_ij_vector == nullptr,
                 "vector appears to be empty (not created)" );
  hypre_IJVectorZeroValues( m_ij_vector );
}

void HypreVector::rand()
{
  GEOS_ERROR_IF( m_ij_vector == nullptr,
                 "vector appears to be empty (not created)" );
  HYPRE_ParVectorSetRandomValues( m_par_vector, 1984 );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Open / close
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void HypreVector::open()
{
  HYPRE_IJVectorInitialize( m_ij_vector );
}

void HypreVector::close()
{
  HYPRE_IJVectorAssemble( m_ij_vector );
}

// ---------------------------------------------------------
// Linear Algebra
// ---------------------------------------------------------
// The following functions support basic linear algebra ops

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Scale
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Multiply all elements by scalingFactor.

void HypreVector::scale( real64 const scalingFactor )
{
  GEOS_ERROR_IF( m_ij_vector == nullptr,
                 "vector appears to be empty (not created)" );
  HYPRE_ParVectorScale( scalingFactor, m_par_vector );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Dot
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Dot product with the vector vec.

real64 HypreVector::dot( HypreVector const &vec )
{
  GEOS_ERROR_IF( m_ij_vector == nullptr,
                 "vector this appears to be empty (not created)" );
  GEOS_ERROR_IF( m_ij_vector == nullptr,
                 "vector vec appears to be empty (not created)" );

  return hypre_ParVectorInnerProd( m_par_vector,
                                   *vec.getHypreParVectorPointer() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Copy
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update vector as this = x.

void HypreVector::copy( HypreVector const &x )
{
  GEOS_ERROR_IF( m_ij_vector == nullptr,
                 "destination vector appears to be empty (not created)" );
  GEOS_ERROR_IF( *x.unwrappedPointer() == nullptr,
                 "source vector appears to be empty (not created)" );
  // TODO: add dimension checks?
  HYPRE_ParVectorCopy( *x.getHypreParVectorPointer(), m_par_vector );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Axpy
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update vector as this = alpha*x + this.

void HypreVector::axpy( real64 const alpha,
                        HypreVector const &x )
{
  GEOS_ERROR_IF( m_ij_vector == nullptr,
                 "destination vector appears to be empty (not created)" );
  GEOS_ERROR_IF( *x.unwrappedPointer() == nullptr,
                 "source vector appears to be empty (not created)" );
  HYPRE_ParVectorAxpy( alpha, *x.getHypreParVectorPointer(), m_par_vector );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Axpby
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update vector as this = alpha*x + beta*this.

void HypreVector::axpby( real64 const alpha,
                         HypreVector const &x,
                         real64 const beta )
{
  GEOS_ERROR_IF( m_ij_vector == nullptr,
                 "destination vector appears to be empty (not created)" );
  GEOS_ERROR_IF( *x.unwrappedPointer() == nullptr,
                 "source vector appears to be empty (not created)" );
  this->scale( beta );
  this->axpy( alpha, x );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 1-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 1-norm of the vector.

real64 HypreVector::norm1() const
{
  GEOS_ERROR_IF( m_ij_vector == nullptr,
                 "vector appears to be empty (not created)" );

  MPI_Comm comm = hypre_IJVectorComm( m_ij_vector );

  real64 const * local_data = this->extractLocalVector();
//  double * local_data = hypre_VectorData( hypre_ParVectorLocalVector (m_par_vector) );

  real64 norm1 = 0;
  real64 loc_norm1 = 0;
  for( int i = 0 ; i < integer_conversion<int>( this->localSize() ) ; ++i )
    loc_norm1 += std::fabs( local_data[i] );

  MPI_Allreduce( &loc_norm1,
                 &norm1,
                 1,
                 MPI_DOUBLE,
                 MPI_SUM,
                 comm );

  return norm1;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 2-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 2-norm of the vector.

real64 HypreVector::norm2() const
{
  GEOS_ERROR_IF( m_ij_vector == nullptr,
                 "vector appears to be empty (not created)" );

  return std::sqrt( hypre_ParVectorInnerProd( m_par_vector,
                                              m_par_vector ) );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm of the vector.

real64 HypreVector::normInf() const
{
  MPI_Comm comm = hypre_IJVectorComm( m_ij_vector );

  real64 const * local_data = this->extractLocalVector();

  real64 normInf = 0;
  real64 loc_normInf = std::fabs( local_data[0] );
  for( int i = 1 ; i < integer_conversion<int>( this->localSize() ) ; ++i )
    loc_normInf = std::max( loc_normInf, std::fabs( local_data[i] ) );


  MPI_Allreduce( &loc_normInf,
                 &normInf,
                 1,
                 MPI_DOUBLE,
                 MPI_MAX,
                 comm );

  return normInf;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print vector to the std::cout in Trilinos format.

void HypreVector::print() const
{
  GEOS_ERROR_IF( m_ij_vector == nullptr,
                 "Vector appears to be empty (not created)" );
  GEOS_ERROR( "not yet implemented" );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Write to file in HYPRE format
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void HypreVector::write( string const & filename ) const
                         {
  GEOS_ERROR_IF( m_ij_vector == nullptr,
                 "vector appears to be empty (not created)" );
  HYPRE_IJVectorPrint( m_ij_vector, filename.c_str() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get value
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get element globalRow

real64 HypreVector::get( globalIndex globalRow ) const
                         {
  GEOS_ERROR_IF( m_ij_vector == nullptr,
                 "vector appears to be empty (not created)" );
  real64 value;
  HYPRE_IJVectorGetValues( m_ij_vector, 1, &globalRow, &value );
  return value;
}

void HypreVector::get( array1d<globalIndex> const & globalIndices,
                       array1d<real64> & values ) const
                       {
  GEOS_ERROR_IF( m_ij_vector == nullptr,
                 "vector appears to be empty (not created)" );
  GEOS_ERROR_IF( globalIndices.size() > values.size(),
                 "IJ_vector, more globalIndices required than available values" );

  HYPRE_Int ierr;
  ierr = HYPRE_IJVectorGetValues( m_ij_vector,
                                  globalIndices.size(),
                                  globalIndices.data(),
                                  values.data() );
  GEOS_ERROR_IF( ierr != 0,
                 "Error getting IJVector values - error code: " +
                     std::to_string( ierr ) );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get unwrapped pointer
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get pointer to raw HYPRE_IJVector object, with const and non-const versions.

HYPRE_IJVector const * HypreVector::unwrappedPointer() const
{
  return &m_ij_vector;
}

HYPRE_IJVector* HypreVector::unwrappedPointer()
{
  return &m_ij_vector;
}

HYPRE_ParVector const * HypreVector::getHypreParVectorPointer() const
{
  return &m_par_vector;
}

HYPRE_ParVector* HypreVector::getHypreParVectorPointer()
{
  return &m_par_vector;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the number of global elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Return the global size of the vector (total number of elements).

globalIndex HypreVector::globalSize() const
{
  GEOS_ERROR_IF( m_ij_vector == nullptr,
                 "vector appears to be empty (not created)" );
  return hypre_IJVectorGlobalNumRows( m_ij_vector );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the number of local elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Return the local size of the vector (total number of local elements).
localIndex HypreVector::localSize() const
{
  GEOS_ERROR_IF( m_ij_vector == nullptr,
                 "vector appears to be empty (not created)" );
  return hypre_ParVectorActualLocalSize( m_par_vector );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// getLocalRowID
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Map a global row index to local row index
localIndex HypreVector::getLocalRowID( globalIndex const index ) const
{
  return integer_conversion< localIndex >( index - this->ilower() ) ;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// getGlobalRowID
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Map a local row index to global row index
globalIndex HypreVector::getGlobalRowID( localIndex const index ) const
{
  return integer_conversion< globalIndex >( index ) +
         integer_conversion< globalIndex >( this->ilower() ) ;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// extractLocalVector
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Extract a view of the local portion of the array
real64 const * HypreVector::extractLocalVector() const
{
  return hypre_VectorData( hypre_ParVectorLocalVector (m_par_vector) );
}

real64 * HypreVector::extractLocalVector()
{
  return hypre_VectorData( hypre_ParVectorLocalVector (m_par_vector) );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// ilower
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the index of the first global row owned by that processor.
globalIndex HypreVector::ilower() const
{
  return integer_conversion< globalIndex >( hypre_ParVectorFirstIndex(m_par_vector) );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// iupper
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// eturns the next index after last global row owned by that processor.
globalIndex HypreVector::iupper() const
{
  return hypre_ParVectorLastIndex(m_par_vector) + 1;
}


} // end geosx

// END_RST_NARRATIVE
