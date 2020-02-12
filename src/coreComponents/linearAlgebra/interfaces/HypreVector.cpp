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

#include "HypreVector.hpp"

// Include required Hypre headers
#include "HYPRE.h"
#include "_hypre_IJ_mv.h"
#include "_hypre_parcsr_mv.h"
#include "HypreUtils.hpp"

// Put everything under the geosx namespace.
namespace geosx
{

// Helper function that performs the following sequence of IJVEctor
// call: Create, SetObjectType, Initialize.
static void initialize( MPI_Comm const & comm,
                        HYPRE_Int const & jlower,
                        HYPRE_Int const & jupper,
                        HYPRE_IJVector & ij_vector )
{
  HYPRE_IJVectorCreate( comm,
                        jlower,
                        jupper,
                        &ij_vector );
  HYPRE_IJVectorSetObjectType( ij_vector,
                               HYPRE_PARCSR );
  HYPRE_IJVectorInitialize( ij_vector );
}

// Helper function that performs the following sequence of IJVEctor
// call: Assemble, GetObject.
static void finalize( HYPRE_IJVector & ij_vector,
                      HYPRE_ParVector & par_vector )
{
  HYPRE_IJVectorAssemble( ij_vector );
  HYPRE_IJVectorGetObject( ij_vector,
                           (void * *) &par_vector );
}

// ----------------------------
// Constructors
// ----------------------------

// Empty constructor
HypreVector::HypreVector()
{}

// Copy constructor
HypreVector::HypreVector( HypreVector const & src )
{
  GEOSX_ASSERT_MSG( src.m_ij_vector != nullptr,
                    "source vector appears to be empty (not created)" );
  this->create( src );
}

// Destructor
HypreVector::~HypreVector()
{
  if( m_ij_vector )
  {
    HYPRE_IJVectorDestroy( m_ij_vector );
    m_ij_vector = nullptr;
  }
}

// Create from HypreVector
void HypreVector::create( HypreVector const & src )
{
  GEOSX_ASSERT_MSG( src.m_ij_vector != nullptr,
                    "source vector appears to be empty (not created)" );

  HYPRE_BigInt jlower, jupper;
  HYPRE_IJVectorGetLocalRange( src.m_ij_vector,
                               &jlower,
                               &jupper );

  initialize( hypre_IJVectorComm( src.m_ij_vector ),
              jlower,
              jupper,
              m_ij_vector );

  finalize( m_ij_vector,
            m_par_vector );

  this->copy( src );
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
// Note: every vector is created initialized to 0 and then 'closed',
//       i.e. made ready to use

void HypreVector::createWithLocalSize( localIndex const localSize,
                                       MPI_Comm const & comm )
{
  GEOSX_ASSERT_MSG( localSize >= 0,
                    "local size is lower than 0" );

  HYPRE_Int const jlower = MpiWrapper::PrefixSum< HYPRE_Int >( integer_conversion< HYPRE_Int >( localSize ) ).first;
  HYPRE_Int const jupper = jlower + integer_conversion< HYPRE_Int >( localSize ) - 1;

  initialize( comm,
              jlower,
              jupper,
              m_ij_vector );

  hypre_IJVectorZeroValues( m_ij_vector );

  finalize( m_ij_vector,
            m_par_vector );
}

void HypreVector::createWithGlobalSize( globalIndex const globalSize,
                                        MPI_Comm const & comm )
{
  int this_mpi_process = MpiWrapper::Comm_rank( comm );
  int n_mpi_process = MpiWrapper::Comm_size( comm );

  HYPRE_Int localSize = integer_conversion< HYPRE_Int >( globalSize )
                        / integer_conversion< HYPRE_Int >( n_mpi_process );
  HYPRE_Int residual = integer_conversion< HYPRE_Int >( globalSize )
                       % integer_conversion< HYPRE_Int >( n_mpi_process );

  GEOSX_ASSERT_MSG( localSize >= 0,
                    "local size is lower than 0" );

  HYPRE_Int jlower;
  HYPRE_Int jupper;
  HYPRE_Int rank = integer_conversion< HYPRE_Int >( this_mpi_process );

  if( this_mpi_process == 0 )
  {
    jlower = 0;
    jupper = localSize + residual - 1;
  }
  else
  {
    jlower = rank * localSize + residual;
    jupper = jlower + localSize - 1;
  }

  initialize( comm,
              jlower,
              jupper,
              m_ij_vector );

  hypre_IJVectorZeroValues( m_ij_vector );

  finalize( m_ij_vector,
            m_par_vector );

}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create from array
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a vector from local array data.  The global vector contains
// local arrays stitched together.

void HypreVector::create( array1d< real64 > const & localValues,
                          MPI_Comm const & comm )
{
  HYPRE_Int localSize = integer_conversion< HYPRE_Int >( localValues.size() );

  GEOSX_ASSERT_MSG( localSize >= 0,
                    "local size is lower than 0" );

  HYPRE_Int const jlower = MpiWrapper::PrefixSum< HYPRE_Int >( localSize ).first;
  HYPRE_Int const jupper = jlower + integer_conversion< HYPRE_Int >( localSize ) - 1;

  initialize( comm,
              jlower,
              jupper,
              m_ij_vector );

  finalize( m_ij_vector,
            m_par_vector );

  HYPRE_Real * local_data = this->extractLocalVector();

  for( localIndex i = 0 ; i < localValues.size() ; ++i )
  {
    local_data[i] = static_cast< HYPRE_Real >( localValues[i] );
  }

}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Add/Set value(s)
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Add/set entries in the vector.

// single element options

void HypreVector::set( globalIndex const globalRow,
                       real64 const value )
{
  GEOSX_ASSERT_MSG( this->getLocalRowID( globalRow ) >= 0,
                    "HypreVector, it is not possible to set values on other processors" );
  HYPRE_Int ierr;
  ierr = HYPRE_IJVectorSetValues( m_ij_vector,
                                  1,
                                  &globalRow,
                                  toHYPRE_Real( &value ) );
  GEOSX_UNUSED_VAR( ierr );
  GEOSX_ASSERT_MSG( ierr == 0,
                    "Error setting values HypreVector - error code: " +
                    std::to_string( ierr ) );
}

void HypreVector::add( globalIndex const globalRow,
                       real64 const value )
{
  HYPRE_Int ierr;
  ierr = HYPRE_IJVectorAddToValues( m_ij_vector,
                                    1,
                                    &globalRow,
                                    toHYPRE_Real( &value ) );
  GEOSX_UNUSED_VAR( ierr );
  GEOSX_ASSERT_MSG( ierr == 0,
                    "Error adding values HypreVector - error code: " +
                    std::to_string( ierr ) );
}

// n-element, c-style options

void HypreVector::set( globalIndex const * globalIndices,
                       real64 const * values,
                       localIndex size )
{
  GEOSX_ASSERT_MSG( this->getLocalRowID( *std::min_element( globalIndices, globalIndices + size ) ) >= 0 &&
                    this->getLocalRowID( *std::max_element( globalIndices, globalIndices + size ) ) >= 0,
                    "HypreVector, it is not possible to set values on other processors" );
  HYPRE_Int ierr;
  ierr = HYPRE_IJVectorSetValues( m_ij_vector,
                                  integer_conversion< HYPRE_Int >( size ),
                                  toHYPRE_BigInt( globalIndices ),
                                  toHYPRE_Real( values ) );
  GEOSX_UNUSED_VAR( ierr );
  GEOSX_ASSERT_MSG( ierr == 0,
                    "Error setting values HypreVector - error code: " +
                    std::to_string( ierr ) );
}

void HypreVector::add( globalIndex const * globalIndices,
                       real64 const * values,
                       localIndex size )
{
  HYPRE_Int ierr;
  ierr = HYPRE_IJVectorAddToValues( m_ij_vector,
                                    integer_conversion< HYPRE_Int >( size ),
                                    toHYPRE_BigInt( globalIndices ),
                                    toHYPRE_Real( values ) );
  GEOSX_UNUSED_VAR( ierr );
  GEOSX_ASSERT_MSG( ierr == 0,
                    "Error adding values HypreVector - error code: " +
                    std::to_string( ierr ) );
}

// n-element, array1d options

void HypreVector::set( array1d< globalIndex > const & globalIndices,
                       array1d< real64 > const & values )
{

  GEOSX_ASSERT_MSG( this->ilower() <= *std::min_element( globalIndices.begin(), globalIndices.end() ) &&
                    *std::max_element( globalIndices.begin(), globalIndices.end()) < this->iupper(),
                    "HypreVector, it is not possible to set values on other processors" );

  HYPRE_Int ierr;
  ierr = HYPRE_IJVectorSetValues( m_ij_vector,
                                  integer_conversion< HYPRE_Int >( values.size() ),
                                  toHYPRE_BigInt( globalIndices.data() ),
                                  toHYPRE_Real( values.data() ) );
  GEOSX_UNUSED_VAR( ierr );
  GEOSX_ASSERT_MSG( ierr == 0,
                    "Error setting values HypreVector - error code: " +
                    std::to_string( ierr ) );
}

void HypreVector::add( array1d< globalIndex > const & globalIndices,
                       array1d< real64 > const & values )
{
  HYPRE_Int ierr;
  ierr = HYPRE_IJVectorAddToValues( m_ij_vector,
                                    integer_conversion< HYPRE_Int >( values.size() ),
                                    toHYPRE_BigInt( globalIndices.data() ),
                                    toHYPRE_Real( values.data() ) );
  GEOSX_UNUSED_VAR( ierr );
  GEOSX_ASSERT_MSG( ierr == 0,
                    "Error adding values HypreVector - error code: " +
                    std::to_string( ierr ) );
}

//additional options:

void HypreVector::set( real64 value )
{
  GEOSX_ERROR_IF( m_ij_vector == nullptr,
                  "vector appears to be empty (not created)" );
  HYPRE_ParVectorSetConstantValues( m_par_vector,
                                    static_cast< HYPRE_Real >( value ) );
}

void HypreVector::zero()
{
  GEOSX_ERROR_IF( m_ij_vector == nullptr,
                  "vector appears to be empty (not created)" );
  hypre_IJVectorZeroValues( m_ij_vector );
}

void HypreVector::rand()
{
  GEOSX_ERROR_IF( m_ij_vector == nullptr,
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
  GEOSX_ASSERT_MSG( m_ij_vector != nullptr,
                    "vector appears to be empty (not created)" );
  HYPRE_Int ierr;
  ierr = HYPRE_ParVectorScale( static_cast< HYPRE_Real >( scalingFactor ),
                               m_par_vector );
  GEOSX_UNUSED_VAR( ierr );
  GEOSX_ASSERT_MSG( ierr == 0,
                    "Error scaling HypreVector - error code: " +
                    std::to_string( ierr ) );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Dot
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Dot product with the vector vec.

real64 HypreVector::dot( HypreVector const & vec ) const
{
  GEOSX_ASSERT_MSG( m_ij_vector != nullptr,
                    "vector this appears to be empty (not created)" );
  GEOSX_ASSERT_MSG( vec.m_ij_vector != nullptr,
                    "vector vec appears to be empty (not created)" );
  GEOSX_ASSERT_MSG( this->localSize() == vec.localSize(),
                    "HypreVector lengths not compatible for dot product" );
  HYPRE_Int ierr;
  HYPRE_Real result;
  ierr = HYPRE_ParVectorInnerProd( m_par_vector,
                                   vec.m_par_vector,
                                   &result );
  GEOSX_UNUSED_VAR( ierr );
  GEOSX_ASSERT_MSG( ierr == 0,
                    "Error dot product HypreVector - error code: " +
                    std::to_string( ierr ) );
  return static_cast< real64 >( result );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Copy
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update vector as this = x.

void HypreVector::copy( HypreVector const & x )
{
  GEOSX_ASSERT_MSG( m_ij_vector != nullptr,
                    "destination vector appears to be empty (not created)" );
  GEOSX_ASSERT_MSG( x.m_ij_vector != nullptr,
                    "source vector appears to be empty (not created)" );
  GEOSX_ASSERT_MSG( this->localSize() == x.localSize(),
                    "HypreVector lengths not compatible for copying" );
  HYPRE_Int ierr;
  ierr = HYPRE_ParVectorCopy( x.m_par_vector,
                              m_par_vector );
  GEOSX_UNUSED_VAR( ierr );
  GEOSX_ASSERT_MSG( ierr == 0,
                    "Error copying HypreVector - error code: " +
                    std::to_string( ierr ) );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Axpy
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update vector as this = alpha*x + this.

void HypreVector::axpy( real64 const alpha,
                        HypreVector const & x )
{
  GEOSX_ASSERT_MSG( m_ij_vector != nullptr,
                    "this vector appears to be empty (not created)" );
  GEOSX_ASSERT_MSG( x.m_ij_vector != nullptr,
                    "source vector appears to be empty (not created)" );
  GEOSX_ASSERT_MSG( this->localSize() == x.localSize(),
                    "HypreVector lengths not compatible for axpy operation" );
  HYPRE_Int ierr;
  ierr = HYPRE_ParVectorAxpy( static_cast< HYPRE_Real >( alpha ),
                              x.m_par_vector,
                              m_par_vector );
  GEOSX_UNUSED_VAR( ierr );
  GEOSX_ASSERT_MSG( ierr == 0,
                    "Error axpy operation on HypreVectors - error code: " +
                    std::to_string( ierr ) );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Axpby
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update vector as this = alpha*x + beta*this.

void HypreVector::axpby( real64 const alpha,
                         HypreVector const & x,
                         real64 const beta )
{
  this->scale( static_cast< HYPRE_Real >( beta ) );
  this->axpy( static_cast< HYPRE_Real >( alpha ), x );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 1-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 1-norm of the vector.

real64 HypreVector::norm1() const
{
  GEOSX_ASSERT_MSG( m_ij_vector != nullptr,
                    "this vector appears to be empty (not created)" );

  real64 const * local_data = this->extractLocalVector();

  real64 norm1 = 0.0;
  real64 loc_norm1 = 0.0;
  for( HYPRE_Int i = 0 ; i < this->localSize() ; ++i )
    loc_norm1 += std::fabs( local_data[i] );

  MpiWrapper::allReduce( &loc_norm1,
                         &norm1,
                         1,
                         MPI_SUM,
                         hypre_IJVectorComm( m_ij_vector ) );

  return norm1;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 2-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 2-norm of the vector.

real64 HypreVector::norm2() const
{
  GEOSX_ASSERT_MSG( m_ij_vector != nullptr,
                    "this vector appears to be empty (not created)" );
  return std::sqrt( this->dot( *this ) );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm of the vector.

real64 HypreVector::normInf() const
{
  GEOSX_ASSERT_MSG( m_ij_vector != nullptr,
                    "this vector appears to be empty (not created)" );

  real64 const * local_data = this->extractLocalVector();

  real64 normInf = 0.0;
  real64 loc_normInf = 0.0;
  for( HYPRE_Int i = 0 ; i < this->localSize() ; ++i )
    loc_normInf = std::max( loc_normInf, std::fabs( local_data[i] ) );

  MpiWrapper::allReduce( &loc_normInf,
                         &normInf,
                         1,
                         MPI_MAX,
                         hypre_IJVectorComm( m_ij_vector ) );

  return normInf;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print vector to the std::cout in Trilinos format.

void HypreVector::print( std::ostream & os ) const
{
  GEOSX_ERROR_IF( m_ij_vector == nullptr,
                  "matrix appears to be empty (not created) or not finalized" );
  if( MpiWrapper::Comm_rank( hypre_IJMatrixComm( m_ij_vector ) ) == 0 )
  {
    os << "Hypre interface: no output on screen available/n";
    os << "                 use write method";
  }
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Write to file in HYPRE format
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void HypreVector::write( string const & filename,
                         bool const mtxFormat ) const
{
  GEOSX_ERROR_IF( m_ij_vector == nullptr,
                  "vector appears to be empty (not created)" );
  if( mtxFormat )
  {
    std::cout << "MatrixMarket not available for HypreMtrix, default used\n";
  }
  HYPRE_IJVectorPrint( m_ij_vector, filename.c_str() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get value
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get element globalRow

real64 HypreVector::get( globalIndex globalRow ) const
{
  GEOSX_ERROR_IF( m_ij_vector == nullptr,
                  "vector appears to be empty (not created)" );
  GEOSX_ASSERT_MSG( this->ilower() <= globalRow &&
                    globalRow < this->iupper(),
                    "HypreVector, globalRow not in local range" );

  HYPRE_Real value;
  HYPRE_Int ierr;
  GEOSX_UNUSED_VAR( ierr );
  ierr = HYPRE_IJVectorGetValues( m_ij_vector,
                                  1,
                                  toHYPRE_BigInt( &globalRow ),
                                  &value );
  GEOSX_ASSERT_MSG( ierr == 0,
                    "Error getting IJVector values - error code: " +
                    std::to_string( ierr ) );

  return static_cast< real64 >( value );
}

void HypreVector::get( array1d< globalIndex > const & globalIndices,
                       array1d< real64 > & values ) const
{
  GEOSX_ERROR_IF( m_ij_vector == nullptr,
                  "vector appears to be empty (not created)" );
  GEOSX_ERROR_IF( globalIndices.size() > values.size(),
                  "IJ_vector, more globalIndices required than available values" );

  HYPRE_Int ierr;
  GEOSX_UNUSED_VAR( ierr );
  ierr = HYPRE_IJVectorGetValues( m_ij_vector,
                                  integer_conversion< HYPRE_Int >( globalIndices.size() ),
                                  toHYPRE_BigInt( globalIndices.data() ),
                                  toHYPRE_Real( values.data() ) );
  GEOSX_ASSERT_MSG( ierr == 0,
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

HYPRE_IJVector * HypreVector::unwrappedPointer()
{
  return &m_ij_vector;
}

HYPRE_ParVector const * HypreVector::getHypreParVectorPointer() const
{
  return &m_par_vector;
}

HYPRE_ParVector * HypreVector::getHypreParVectorPointer()
{
  return &m_par_vector;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the number of global elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Return the global size of the vector (total number of elements).

globalIndex HypreVector::globalSize() const
{
  GEOSX_ERROR_IF( m_ij_vector == nullptr,
                  "vector appears to be empty (not created)" );
  return hypre_IJVectorGlobalNumRows( m_ij_vector );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the number of local elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Return the local size of the vector (total number of local elements).
localIndex HypreVector::localSize() const
{
  GEOSX_ERROR_IF( m_ij_vector == nullptr,
                  "vector appears to be empty (not created)" );
  return hypre_ParVectorActualLocalSize( m_par_vector );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// getLocalRowID
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Map a global row index to local row index
localIndex HypreVector::getLocalRowID( globalIndex const index ) const
{
  if( index < this->iupper() )
  {
    return std::max( integer_conversion< localIndex >( -1 ),
                     integer_conversion< localIndex >( index - this->ilower() ) );
  }
  else
  {
    return integer_conversion< localIndex >( -1 );
  }

}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// getGlobalRowID
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Map a local row index to global row index
globalIndex HypreVector::getGlobalRowID( localIndex const index ) const
{
  return integer_conversion< globalIndex >( index ) +
         integer_conversion< globalIndex >( this->ilower() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// extractLocalVector
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Extract a view of the local portion of the array
real64 const * HypreVector::extractLocalVector() const
{
  return toGEOSX_real64( hypre_VectorData( hypre_ParVectorLocalVector ( m_par_vector ) ) );
}

real64 * HypreVector::extractLocalVector()
{
  return toGEOSX_real64( hypre_VectorData( hypre_ParVectorLocalVector ( m_par_vector ) ) );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// ilower
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the index of the first global row owned by that processor.
globalIndex HypreVector::ilower() const
{
  return integer_conversion< globalIndex >( hypre_ParVectorFirstIndex( m_par_vector ) );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// iupper
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// eturns the next index after last global row owned by that processor.
globalIndex HypreVector::iupper() const
{
  return hypre_ParVectorLastIndex( m_par_vector ) + 1;
}

std::ostream & operator<<( std::ostream & os,
                           HypreVector const & vec )
{
  vec.print( os );
  return os;
}

} // end geosx

// END_RST_NARRATIVE
