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
 * @file HypreVector.cpp
 */

#include "HypreVector.hpp"

#include "codingUtilities/Utilities.hpp"
#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"

#include "HYPRE.h"
#include "_hypre_IJ_mv.h"
#include "_hypre_parcsr_mv.h"

#include <iomanip>

namespace geosx
{

// Helper function that performs the following sequence of IJVEctor
// call: Create, SetObjectType, Initialize.
static void initialize( MPI_Comm const & comm,
                        HYPRE_BigInt const & jlower,
                        HYPRE_BigInt const & jupper,
                        HYPRE_IJVector & ij_vector )
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorCreate( comm, jlower, jupper, &ij_vector ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorSetObjectType( ij_vector, HYPRE_PARCSR ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorInitialize( ij_vector ) );
}

// Helper function that performs the following sequence of IJVEctor
// call: Assemble, GetObject.
static void finalize( HYPRE_IJVector & ij_vector,
                      HYPRE_ParVector & par_vector )
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorAssemble( ij_vector ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorGetObject( ij_vector, (void * *) &par_vector ) );
}

HypreVector::HypreVector()
  : VectorBase(),
  m_ij_vector{},
  m_par_vector{}
{}

// Copy constructor
HypreVector::HypreVector( HypreVector const & src )
  : HypreVector()
{
  *this = src;
}

// Move constructor
HypreVector::HypreVector( HypreVector && src ) noexcept
  : HypreVector()
{
  *this = std::move( src );
}

HypreVector & HypreVector::operator=( HypreVector const & src )
{
  if( &src != this )
  {
    reset();
    if( src.created() )
    {
      HYPRE_BigInt jlower, jupper;
      GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorGetLocalRange( src.m_ij_vector, &jlower, &jupper ) );
      initialize( src.getComm(), jlower, jupper, m_ij_vector );
      finalize( m_ij_vector, m_par_vector );
      if( src.ready() )
      {
        copy( src );
      }
    }
  }
  return *this;
}

HypreVector & HypreVector::operator=( HypreVector && src ) noexcept
{
  if( &src != this )
  {
    std::swap( m_ij_vector, src.m_ij_vector );
    std::swap( m_par_vector, src.m_par_vector );
    std::swap( m_closed, src.m_closed );
  }
  return *this;
}

void HypreVector::reset()
{
  VectorBase::reset();
  if( m_ij_vector )
  {
    HYPRE_IJVectorDestroy( m_ij_vector );
    m_ij_vector = nullptr;
  }
}

HypreVector::~HypreVector()
{
  reset();
}

void HypreVector::createWithLocalSize( localIndex const localSize,
                                       MPI_Comm const & comm )
{
  GEOSX_LAI_ASSERT( closed() );
  GEOSX_LAI_ASSERT_GE( localSize, 0 );

  reset();

  HYPRE_BigInt const jlower = MpiWrapper::prefixSum< HYPRE_BigInt >( LvArray::integerConversion< HYPRE_BigInt >( localSize ), comm );
  HYPRE_BigInt const jupper = jlower + LvArray::integerConversion< HYPRE_BigInt >( localSize ) - 1;

  initialize( comm, jlower, jupper, m_ij_vector );
  GEOSX_LAI_CHECK_ERROR( hypre_IJVectorZeroValues( m_ij_vector ) );
  finalize( m_ij_vector, m_par_vector );
}

void HypreVector::createWithGlobalSize( globalIndex const globalSize,
                                        MPI_Comm const & comm )
{
  GEOSX_LAI_ASSERT( closed() );
  GEOSX_LAI_ASSERT_GE( globalSize, 0 );

  reset();

  HYPRE_Int const rank  = LvArray::integerConversion< HYPRE_Int >( MpiWrapper::commRank( comm ) );
  HYPRE_Int const nproc = LvArray::integerConversion< HYPRE_Int >( MpiWrapper::commSize( comm ) );

  HYPRE_Int const localSize = LvArray::integerConversion< HYPRE_Int >( globalSize / nproc );
  HYPRE_Int const residual  = LvArray::integerConversion< HYPRE_Int >( globalSize % nproc );

  HYPRE_BigInt const ilower = LvArray::integerConversion< HYPRE_BigInt >( rank * localSize + ( rank == 0 ? 0 : residual ) );
  HYPRE_BigInt const iupper = LvArray::integerConversion< HYPRE_BigInt >( ilower + localSize + ( rank == 0 ? residual : 0 ) - 1 );

  initialize( comm, ilower, iupper, m_ij_vector );
  GEOSX_LAI_CHECK_ERROR( hypre_IJVectorZeroValues( m_ij_vector ) );
  finalize( m_ij_vector, m_par_vector );
}

void HypreVector::create( arrayView1d< real64 const > const & localValues,
                          MPI_Comm const & comm )
{
  GEOSX_LAI_ASSERT( closed() );

  HYPRE_BigInt const localSize = LvArray::integerConversion< HYPRE_BigInt >( localValues.size() );
  HYPRE_BigInt const jlower = MpiWrapper::prefixSum< HYPRE_BigInt >( localSize, comm );
  HYPRE_BigInt const jupper = jlower + localSize - 1;

  // In case the vector was already created, we reset it to prevent any memory leak...
  reset();
  // ... then we can continue with the standard creation process.
  initialize( comm, jlower, jupper, m_ij_vector );
  finalize( m_ij_vector, m_par_vector );

  HYPRE_Real * const local_data = extractLocalVector();

  forAll< hypre::execPolicy >( localValues.size(), [=] GEOSX_HYPRE_HOST_DEVICE ( localIndex const i )
  {
    local_data[i] = localValues[i];
  } );

}

bool HypreVector::created() const
{
  return m_ij_vector != nullptr && m_par_vector != nullptr;
}

void HypreVector::set( globalIndex const globalRow,
                       real64 const value )
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_ASSERT_GE( globalRow, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), globalRow );


#if defined(GEOSX_USE_HYPRE_CUDA)
  array1d< globalIndex > globalRowDevice( 1 );
  array1d< real64 > valueDevice( 1 );
  globalRowDevice[0] = globalRow;
  valueDevice[0] = value;
  globalRowDevice.move( LvArray::MemorySpace::cuda, false );
  valueDevice.move( LvArray::MemorySpace::cuda, false );
  globalIndex const * const pGlobalRow = globalRowDevice.data();
  real64 const * const pValue = valueDevice.data();
#else
  globalIndex const * const pGlobalRow = &globalRow;
  real64 const * const pValue = &value;
#endif
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorSetValues( m_ij_vector, 1, pGlobalRow, pValue ) );
}

void HypreVector::add( globalIndex const globalRow,
                       real64 const value )
{
  GEOSX_LAI_ASSERT( !closed() );
#if defined(GEOSX_USE_HYPRE_CUDA)
  array1d< globalIndex > globalRowDevice( 1 );
  array1d< real64 > valueDevice( 1 );
  globalRowDevice[0] = globalRow;
  valueDevice[0] = value;
  globalRowDevice.move( LvArray::MemorySpace::cuda, false );
  valueDevice.move( LvArray::MemorySpace::cuda, false );
  globalIndex const * const pGlobalRow = globalRowDevice.data();
  real64 const * const pValue = valueDevice.data();
#else
  globalIndex const * const pGlobalRow = &globalRow;
  real64 const * const pValue = &value;
#endif
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorAddToValues( m_ij_vector, 1, pGlobalRow, pValue ) );
}

void HypreVector::set( globalIndex const * globalIndices,
                       real64 const * values,
                       localIndex size )
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_ASSERT_GE( *std::min_element( globalIndices, globalIndices + size ), ilower() );
  GEOSX_LAI_ASSERT_GE( iupper(), getLocalRowID( *std::max_element( globalIndices, globalIndices + size ) ) );
//  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorSetValues( m_ij_vector,
//                                                  LvArray::integerConversion< HYPRE_Int >( size ),
//                                                  hypre::toHypreBigInt( globalIndices ),
//                                                  values ) );

  HYPRE_Real * const local_data = hypre_VectorData( hypre_ParVectorLocalVector ( m_par_vector ) );

  forAll< hypre::execPolicy >( size,
                               [=] GEOSX_HYPRE_HOST_DEVICE ( localIndex const i )
  {
    local_data[i] = values[i];
  } );


}

void HypreVector::add( globalIndex const * globalIndices,
                       real64 const * values,
                       localIndex size )
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorAddToValues( m_ij_vector,
                                                    LvArray::integerConversion< HYPRE_Int >( size ),
                                                    hypre::toHypreBigInt( globalIndices ),
                                                    values ) );
}

void HypreVector::set( arraySlice1d< globalIndex const > const & globalIndices,
                       arraySlice1d< real64 const > const & values )
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_ASSERT_GE( *std::min_element( globalIndices.dataIfContiguous(),
                                          globalIndices.dataIfContiguous() + globalIndices.size() ), ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), *std::max_element( globalIndices.dataIfContiguous(),
                                                    globalIndices.dataIfContiguous() + globalIndices.size() ) );

  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorSetValues( m_ij_vector,
                                                  LvArray::integerConversion< HYPRE_Int >( values.size() ),
                                                  hypre::toHypreBigInt( globalIndices.dataIfContiguous() ),
                                                  values.dataIfContiguous() ) );
}

void HypreVector::add( arraySlice1d< globalIndex const > const & globalIndices,
                       arraySlice1d< real64 const > const & values )
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorAddToValues( m_ij_vector,
                                                    LvArray::integerConversion< HYPRE_Int >( values.size() ),
                                                    hypre::toHypreBigInt( globalIndices.dataIfContiguous() ),
                                                    values.dataIfContiguous() ) );
}

void HypreVector::set( real64 value )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParVectorSetConstantValues( m_par_vector, value ) );
}

void HypreVector::zero()
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( hypre_IJVectorZeroValues( m_ij_vector ) );
}

void HypreVector::rand( unsigned const seed )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParVectorSetRandomValues( m_par_vector, seed ) );
}

void HypreVector::open()
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorInitialize( m_ij_vector ) );
  m_closed = false;
}

void HypreVector::close()
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorAssemble( m_ij_vector ) );
  m_closed = true;
}

void HypreVector::scale( real64 const scalingFactor )
{
  GEOSX_LAI_ASSERT( ready() );
  if( !isEqual( scalingFactor, 1.0 ) )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_ParVectorScale( scalingFactor, m_par_vector ) );
  }
}

void HypreVector::reciprocal()
{
  GEOSX_LAI_ASSERT( ready() );
  real64 * const values = extractLocalVector();
  forAll< hypre::execPolicy >( localSize(), [=] GEOSX_HYPRE_HOST_DEVICE ( localIndex const i )
  {
    values[i] = 1.0 / values[i];
  } );
}

real64 HypreVector::dot( HypreVector const & vec ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( vec.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), vec.globalSize() );

  HYPRE_Real result;
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParVectorInnerProd( m_par_vector, vec.m_par_vector, &result ) );
  return result;
}

void HypreVector::copy( HypreVector const & x )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), x.globalSize() );

  GEOSX_LAI_CHECK_ERROR( HYPRE_ParVectorCopy( x.m_par_vector, m_par_vector ) );
}

void HypreVector::axpy( real64 const alpha,
                        HypreVector const & x )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), x.globalSize() );

  if( !isEqual( alpha, 0.0 ) )
  {
    if( &x != this )
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_ParVectorAxpy( alpha, x.m_par_vector, m_par_vector ) );
    }
    else
    {
      scale( 1.0 + alpha );
    }
  }
}

void HypreVector::axpby( real64 const alpha,
                         HypreVector const & x,
                         real64 const beta )
{
  if( &x != this )
  {
    scale( beta );
    axpy( alpha, x );
  }
  else
  {
    scale( alpha + beta );
  }
}

void HypreVector::pointwiseProduct( HypreVector const & x,
                                    HypreVector & y ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  GEOSX_LAI_ASSERT( y.ready() );
  GEOSX_LAI_ASSERT_EQ( localSize(), x.localSize() );
  GEOSX_LAI_ASSERT_EQ( localSize(), y.localSize() );

  real64 const * const data = extractLocalVector();
  real64 const * const x_data = x.extractLocalVector();
  real64 * const y_data = y.extractLocalVector();
  forAll< hypre::execPolicy >( localSize(), [=] GEOSX_HYPRE_HOST_DEVICE ( localIndex const i )
  {
    y_data[i] = data[i] * x_data[i];
  } );
}

real64 HypreVector::norm1() const
{
  GEOSX_LAI_ASSERT( ready() );

  real64 const * const values = extractLocalVector();
  RAJA::ReduceSum< ReducePolicy< hypre::execPolicy >, real64 > localNorm( 0.0 );
  forAll< hypre::execPolicy >( localSize(), [=] GEOSX_HYPRE_HOST_DEVICE ( localIndex const i )
  {
    localNorm += fabs( values[i] );
  } );
  return MpiWrapper::sum( localNorm.get(), getComm() );
}

real64 HypreVector::norm2() const
{
  GEOSX_LAI_ASSERT( ready() );
  return std::sqrt( dot( *this ) );
}

real64 HypreVector::normInf() const
{
  GEOSX_LAI_ASSERT( ready() );

  real64 const * const values = extractLocalVector();
  RAJA::ReduceMax< ReducePolicy< hypre::execPolicy >, real64 > localNorm( 0.0 );
  forAll< hypre::execPolicy >( localSize(), [=] GEOSX_HYPRE_HOST_DEVICE ( localIndex const i )
  {
    localNorm.max( fabs( values[i] ) );
  } );
  return MpiWrapper::max( localNorm.get(), getComm() );
}

void HypreVector::print( std::ostream & os ) const
{
  GEOSX_LAI_ASSERT( ready() );

  int const myRank = MpiWrapper::commRank( getComm() );
  int const numProcs = MpiWrapper::commSize( getComm() );
  char str[77];

  constexpr char const lineFormat[] = "{:>11}{:>18}{:>28.16e}\n";
  constexpr char const headFormat[] = "{:>11}{:>18}{:>28}\n";

  if( myRank == 0 )
  {
    GEOSX_FMT_TO( str, sizeof( str ), headFormat, "MPI_Process", "GlobalRowID", "Value" );
    os << str;
  }

  for( int rank = 0; rank < numProcs; ++rank )
  {
    MpiWrapper::barrier( getComm() );
    if( rank == myRank )
    {
      real64 const * const local_data = extractLocalVector();
      globalIndex const firstRowID = ilower();
      for( localIndex i = 0; i < localSize(); ++i )
      {
        GEOSX_FMT_TO( str, sizeof( str ), lineFormat,
                      rank,
                      firstRowID + i,
                      local_data[i] );
        os << str;
      }
    }
  }
}

void HypreVector::write( string const & filename,
                         LAIOutputFormat const format ) const
{
  GEOSX_LAI_ASSERT( ready() );
  switch( format )
  {
    case LAIOutputFormat::NATIVE_ASCII:
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorPrint( m_ij_vector, filename.c_str() ) );
      break;
    }
    case LAIOutputFormat::MATRIX_MARKET:
    {
      if( globalSize() == 0 )
      {
        if( MpiWrapper::commRank( getComm() ) == 0 )
        {
          FILE * fp = std::fopen( filename.c_str(), "w" );
          hypre_fprintf( fp, "%s", "%%MatrixMarket matrix array real general\n" );
          hypre_fprintf( fp, "%d %d\n", 0, 1 );
          std::fclose( fp );
        }
      }
      else
      {
        // Copy distributed parVector in a local vector on every process
        // with at least one component
        // Warning: works for a parVector that is smaller than 2^31-1
        hypre_Vector *vector;
        vector = hypre_ParVectorToVectorAll( m_par_vector );

        // Identify the smallest process where vector exists
        int myID = MpiWrapper::commRank( getComm() );
        if( vector == 0 )
        {
          myID = MpiWrapper::commSize( getComm() );
        }
        int printID = MpiWrapper::min( myID, getComm() );

        // Write to file vector
        if( MpiWrapper::commRank( getComm() ) == printID )
        {
          FILE * fp = std::fopen( filename.c_str(), "w" );
          HYPRE_Real * data = hypre_VectorData( vector );
          HYPRE_Int size    = hypre_VectorSize( vector );

          hypre_fprintf( fp, "%s", "%%MatrixMarket matrix array real general\n" );
          hypre_fprintf( fp, "%d %d\n", size, 1 );

          for( HYPRE_Int i = 0; i < size; i++ )
          {
            hypre_fprintf( fp, "%.16e\n", data[i] );
          }

          std::fclose( fp );
        }

        // Destroy vector
        if( vector )
        {
          GEOSX_LAI_CHECK_ERROR( hypre_SeqVectorDestroy( vector ) );
        }
      }
      break;
    }
    default:
      GEOSX_ERROR( "Unsupported vector output format" );
  }
}

real64 HypreVector::get( globalIndex globalRow ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT_GE( globalRow, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), globalRow );

  HYPRE_Real value;
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorGetValues( m_ij_vector,
                                                  1,
                                                  hypre::toHypreBigInt( &globalRow ),
                                                  &value ) );
  return value;
}

void HypreVector::get( arraySlice1d< globalIndex const > const & globalIndices,
                       arraySlice1d< real64 > const & values ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT_GE( values.size(), globalIndices.size() );
  GEOSX_LAI_ASSERT_GE( *std::min_element( globalIndices.dataIfContiguous(), globalIndices.dataIfContiguous() + globalIndices.size() ), ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), *std::max_element( globalIndices.dataIfContiguous(),
                                                    globalIndices.dataIfContiguous() + globalIndices.size() ) );

  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorGetValues( m_ij_vector,
                                                  LvArray::integerConversion< HYPRE_Int >( globalIndices.size() ),
                                                  hypre::toHypreBigInt( globalIndices.dataIfContiguous() ),
                                                  values.dataIfContiguous() ) );
}

HYPRE_ParVector const & HypreVector::unwrapped() const
{
  return m_par_vector;
}

HYPRE_IJVector const & HypreVector::unwrappedIJ() const
{
  return m_ij_vector;
}

globalIndex HypreVector::globalSize() const
{
  GEOSX_LAI_ASSERT( created() );
  return hypre_IJVectorGlobalNumRows( m_ij_vector );
}

localIndex HypreVector::localSize() const
{
  GEOSX_LAI_ASSERT( created() );
  return hypre_ParVectorActualLocalSize( m_par_vector );
}

localIndex HypreVector::getLocalRowID( globalIndex const index ) const
{
  GEOSX_LAI_ASSERT( created() );
  return (index >= ilower() && index < iupper()) ? LvArray::integerConversion< localIndex >( index - ilower() ) : -1;
}

globalIndex HypreVector::getGlobalRowID( localIndex const index ) const
{
  GEOSX_LAI_ASSERT( created() );
  return ilower() + index;
}

real64 const * HypreVector::extractLocalVector() const
{
  GEOSX_LAI_ASSERT( ready() );
  return hypre_VectorData( hypre_ParVectorLocalVector ( m_par_vector ) );
}

real64 * HypreVector::extractLocalVector()
{
  GEOSX_LAI_ASSERT( ready() );
  return hypre_VectorData( hypre_ParVectorLocalVector ( m_par_vector ) );
}

void HypreVector::extract( arrayView1d< real64 > const & localVector ) const
{
  GEOSX_LAI_ASSERT_EQ( localSize(), localVector.size() );
  real64 const * const data = extractLocalVector();
  forAll< hypre::execPolicy >( localSize(), [=] GEOSX_HYPRE_HOST_DEVICE ( HYPRE_Int const i )
  {
    localVector[i] = data[i];
  } );
}

globalIndex HypreVector::ilower() const
{
  GEOSX_LAI_ASSERT( created() );
  return LvArray::integerConversion< globalIndex >( hypre_ParVectorFirstIndex( m_par_vector ) );
}

globalIndex HypreVector::iupper() const
{
  GEOSX_LAI_ASSERT( created() );
  return LvArray::integerConversion< globalIndex >( hypre_ParVectorLastIndex( m_par_vector ) ) + 1;
}

MPI_Comm HypreVector::getComm() const
{
  GEOSX_LAI_ASSERT( created() );
  return hypre_IJVectorComm( m_ij_vector );
}

} // end namespace geosx
