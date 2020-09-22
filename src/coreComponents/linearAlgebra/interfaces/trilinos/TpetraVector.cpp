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
 * @file TpetraVector.cpp
 */

#include "TpetraVector.hpp"

#include "codingUtilities/Utilities.hpp"
#include "linearAlgebra/interfaces/trilinos/TpetraUtils.hpp"

#include <Tpetra_Vector.hpp>
#include <KokkosCompat_View.hpp>
#include <MatrixMarket_Tpetra.hpp>

namespace geosx
{

TpetraVector::TpetraVector()
  : VectorBase(),
  m_vector{}
{}

TpetraVector::TpetraVector( TpetraVector const & src )
  : TpetraVector()
{
  *this = src;
}

TpetraVector::TpetraVector( TpetraVector && src ) noexcept
  : TpetraVector()
{
  *this = std::move( src );
}

TpetraVector & TpetraVector::operator=( TpetraVector const & src )
{
  GEOSX_LAI_ASSERT( &src != this );
  GEOSX_LAI_ASSERT( src.ready() );
  if( m_vector )
  {
    Tpetra::deep_copy( *m_vector, *src.m_vector );
  }
  else
  {
    m_vector = std::make_unique< Tpetra_Vector >( *src.m_vector, Teuchos::Copy );
  }
  return *this;
}

TpetraVector & TpetraVector::operator=( TpetraVector && src ) noexcept
{
  GEOSX_LAI_ASSERT( &src != this );
  GEOSX_LAI_ASSERT( src.ready() );
  m_vector = std::move( src.m_vector );
  return *this;
}

TpetraVector::~TpetraVector() = default;

bool TpetraVector::created() const
{
  return bool(m_vector);
}

void TpetraVector::createWithLocalSize( localIndex const localSize,
                                        MPI_Comm const & MPI_PARAM( comm ) )
{
  GEOSX_LAI_ASSERT( closed() );
  GEOSX_LAI_ASSERT_GE( localSize, 0 );

  Teuchos::RCP< Tpetra_Comm > tpetraComm( new Tpetra_Comm( MPI_PARAM( comm ) ) );
  Teuchos::RCP< Tpetra_Map > map( new Tpetra_Map( Teuchos::OrdinalTraits< Tpetra::global_size_t >::invalid(),
                                                  LvArray::integerConversion< size_t >( localSize ),
                                                  0,
                                                  tpetraComm ) );
  m_vector = std::make_unique< Tpetra_Vector >( map );
}

void TpetraVector::createWithGlobalSize( globalIndex const globalSize,
                                         MPI_Comm const & MPI_PARAM( comm ) )
{
  GEOSX_LAI_ASSERT( closed() );
  GEOSX_LAI_ASSERT_GE( globalSize, 0 );

  Teuchos::RCP< Tpetra_Comm > tpetraComm( new Tpetra_Comm( MPI_PARAM( comm ) ) );
  Teuchos::RCP< Tpetra_Map > map( new Tpetra_Map( LvArray::integerConversion< Tpetra::global_size_t >( globalSize ),
                                                  0,
                                                  tpetraComm ) );
  m_vector = std::make_unique< Tpetra_Vector >( map );
}

void TpetraVector::createWithLocalValues( arrayView1d< real64 > const & localValues,
                                          MPI_Comm const & comm )
{
  Teuchos::RCP< Tpetra_Comm > tpetraComm( new Tpetra_Comm( MPI_PARAM( comm ) ) );
  Teuchos::RCP< Tpetra_Map > map( new Tpetra_Map( Teuchos::OrdinalTraits< Tpetra::global_size_t >::invalid(),
                                                  LvArray::integerConversion< size_t >( localValues.size() ),
                                                  0,
                                                  tpetraComm ) );

  // Define what "device" is
#ifdef GEOSX_USE_CUDA
  using KokkosSpaceDevice = Kokkos::CudaSpace; // Not Kokkos::CudaUVMSpace, since we wrap a normal device pointer
  localValues.move( LvArray::MemorySpace::GPU, false ); // Make sure data is current on device
#else
  using KokkosSpaceDevice = Kokkos::HostSpace;
#endif // GEOSX_USE_CUDA

  // Wrap a non-owning Kokkos view around input data
  using InputView = Kokkos::View< real64 * *, Kokkos::LayoutLeft, KokkosSpaceDevice, Kokkos::MemoryUnmanaged >;
  InputView const localValuesView( localValues.data(), localValues.size(), 1 );

  // Sadly, Tpetra currently **requires** the use of UM as default memory space for all its views.
  // The disabled code below would be valid if we were allowed to turn UM off and use regular device allocation.
  // Unfortunately, we are not, and zero-copy construction from a device pointer is impossible to implement.
  // TODO: reevaluate this if Tpetra relaxes the UM requirement.
#if 0

  // Specify a non-owning DualView
  using KokkosDualView = Kokkos::DualView< real64 * *, // this type required by Tpetra::MultiVector, even for a 1D vector
                                           Kokkos::LayoutLeft, // doesn't matter, since the first extent will be 1
                                           Tpetra_Vector::device_type::execution_space,
                                           Kokkos::MemoryUnmanaged >; // tell the view we own the memory

  // DualView constructor requires that host and device views are in sync.
  localValues.move( LvArray::MemorySpace::CPU, false );

  // Directly use device view of input data (shallow copy)
  KokkosDualView::t_dev const values_d = localValuesView;

  // Also create a non-owning view pointing to host data
  KokkosDualView::t_host const values_h( localValues.data( LvArray::MemorySpace::CPU ), localValues.size(), 1 );

#else

  // Use an owning view, since we'll be making a copy
  using KokkosDualView = Tpetra_Vector::dual_view_type;

  // Make a deep copy of input values on device
  KokkosDualView::t_dev const values_d( "", localValuesView.size(), 1 );
  Kokkos::deep_copy( values_d, localValuesView );

  // Also make a host copy as required by DualView
  KokkosDualView::t_host const values_h = Kokkos::create_mirror_view( values_d );

#endif // 0

  // Make a DualView containing host and device views (either shallow or deep) and feed to vector
  KokkosDualView values( values_d, values_h );
  m_vector = std::make_unique< Tpetra_Vector >( map, values );
}

void TpetraVector::open()
{
  GEOSX_LAI_ASSERT( ready() );
  m_closed = false;
  m_vector->sync_host(); // since all element-wise modifications work on host
}

void TpetraVector::close()
{
  GEOSX_LAI_ASSERT( !closed() );
  m_closed = true;
  m_vector->sync_device(); // assuming we've modified host data while vector was open

  // Nothing to do w.r.t. parallel sync since Tpetra::Vector does not support remote assembly.
  // Tpetra::FEMultiVector does, but as of May 2020 it is still lacking proper constructors, documentation, etc.
}

void TpetraVector::reset()
{
  VectorBase::reset();
  m_vector.reset();
}

void TpetraVector::set( globalIndex const globalRowIndex,
                        real64 const value )
{
  GEOSX_LAI_ASSERT( !closed() );
  m_vector->replaceGlobalValue( globalRowIndex, value );
  m_vector->modify_host();
}

void TpetraVector::add( globalIndex const globalRowIndex,
                        real64 const value )
{
  GEOSX_LAI_ASSERT( !closed() );
  m_vector->sumIntoGlobalValue( globalRowIndex, value );
  m_vector->modify_host();
}

void TpetraVector::add( globalIndex const * globalRowIndices,
                        real64 const * values,
                        localIndex const size )
{
  GEOSX_LAI_ASSERT( !closed() );
  for( localIndex i = 0; i < size; ++i )
  {
    m_vector->sumIntoGlobalValue( globalRowIndices[i], values[i] );
  }
  m_vector->modify_host();
}

void TpetraVector::set( globalIndex const * globalRowIndices,
                        real64 const * values,
                        localIndex const size )
{
  GEOSX_LAI_ASSERT( !closed() );
  for( localIndex i = 0; i < size; ++i )
  {
    m_vector->replaceGlobalValue( globalRowIndices[i], values[i] );
  }
  m_vector->modify_host();
}

void TpetraVector::set( arraySlice1d< globalIndex const > const & globalRowIndices,
                        arraySlice1d< real64 const > const & values )
{
  set( globalRowIndices.dataIfContiguous(), values.dataIfContiguous(), globalRowIndices.size() );
}

void TpetraVector::add( arraySlice1d< globalIndex const > const & globalRowIndices,
                        arraySlice1d< real64 const > const & values )
{
  add( globalRowIndices.dataIfContiguous(), values.dataIfContiguous(), globalRowIndices.size() );
}

void TpetraVector::set( real64 const value )
{
  GEOSX_LAI_ASSERT( ready() );
  m_vector->modify_device(); // ensure putScalar runs on device
  m_vector->putScalar( value );
  m_vector->sync_host();
}

void TpetraVector::zero()
{
  set( 0.0 );
}

void TpetraVector::rand( unsigned int const GEOSX_UNUSED_PARAM( seed ) )
{
  GEOSX_LAI_ASSERT( ready() );
  // There's no way currently (May 2020) to set the random seed.
  m_vector->randomize( -1.0, 1.0 );
}

void TpetraVector::scale( real64 const scalingFactor )
{
  GEOSX_LAI_ASSERT( ready() );
  m_vector->scale( scalingFactor );
}

void TpetraVector::reciprocal()
{
  GEOSX_LAI_ASSERT( ready() );
  m_vector->reciprocal( *m_vector );
}

real64 TpetraVector::dot( TpetraVector const & vec ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( vec.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), vec.globalSize() );

  return m_vector->dot( *vec.m_vector );
}

void TpetraVector::copy( TpetraVector const & x )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), x.globalSize() );

  Tpetra::deep_copy( *m_vector, *x.m_vector );
}

void TpetraVector::axpy( real64 alpha,
                         TpetraVector const & x )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), x.globalSize() );

  m_vector->update( alpha, x.unwrapped(), 1. );
}

void TpetraVector::axpby( real64 alpha,
                          TpetraVector const & x,
                          real64 beta )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), x.globalSize() );

  m_vector->update( alpha, x.unwrapped(), beta );
}

real64 TpetraVector::norm1() const
{
  GEOSX_LAI_ASSERT( ready() );
  return m_vector->norm1();
}

real64 TpetraVector::norm2() const
{
  GEOSX_LAI_ASSERT( ready() );
  return m_vector->norm2();
}

real64 TpetraVector::normInf() const
{
  GEOSX_LAI_ASSERT( ready() );
  return m_vector->normInf();
}

globalIndex TpetraVector::globalSize() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_vector->getGlobalLength();
}

localIndex TpetraVector::localSize() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_vector->getLocalLength();
}

globalIndex TpetraVector::ilower() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_vector->getMap()->getMinGlobalIndex();
}

globalIndex TpetraVector::iupper() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_vector->getMap()->getMaxGlobalIndex() + 1;
}

real64 TpetraVector::get( globalIndex const globalRow ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT_GE( globalRow, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), globalRow );

  // TODO: deprecate?

  m_vector->sync_host();
  return m_vector->getLocalViewHost().data()[ getLocalRowID( globalRow ) ];
}

void TpetraVector::get( arraySlice1d< globalIndex const > const & globalIndices,
                        arraySlice1d< real64 > const & values ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT_GE( values.size(), globalIndices.size() );
  GEOSX_LAI_ASSERT_GE( *std::min_element( globalIndices.dataIfContiguous(), globalIndices.dataIfContiguous() + globalIndices.size() ), ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), *std::max_element( globalIndices.dataIfContiguous(),
                                                    globalIndices.dataIfContiguous() + globalIndices.size() ) );

  // TODO: deprecate?

  m_vector->sync_host();
  real64 const * const localVector = m_vector->getLocalViewHost().data();
  for( localIndex i = 0; i < globalIndices.size(); ++i )
  {
    values[i] = localVector[ getLocalRowID( globalIndices[i] ) ];
  }
}

localIndex TpetraVector::getLocalRowID( globalIndex const globalRow ) const
{
  GEOSX_LAI_ASSERT( created() );
  return m_vector->getMap()->getLocalElement( globalRow );
}

globalIndex TpetraVector::getGlobalRowID( localIndex const localRow ) const
{
  GEOSX_LAI_ASSERT( created() );
  return m_vector->getMap()->getGlobalElement( localRow );
}

real64 const * TpetraVector::extractLocalVector() const
{
  GEOSX_LAI_ASSERT( ready() );

  // TODO: deprecate?

  // TODO: This will always return pointer in vector's target memory space (i.e. device with CUDA).
  //       If we want to get data in a different space, we'll need separate methods or a parameter.
  return m_vector->getLocalViewDevice().data();
}

real64 * TpetraVector::extractLocalVector()
{
  GEOSX_LAI_ASSERT( ready() );

  // TODO: deprecate?

  // TODO: This will always return pointer in vector's target memory space (i.e. device with CUDA).
  //       If we want to get data in a different space, we'll need separate methods or a parameter.

  // Mark data dirty since we are returning a pointer to non-const to the user
  m_vector->modify_device();
  return m_vector->getLocalViewDevice().data();
}

void TpetraVector::extract( arrayView1d< real64 > const & localVector ) const
{
  GEOSX_LAI_ASSERT_EQ( localSize(), localVector.size() );
  m_vector->sync_device();
  real64 const * const data = m_vector->getLocalViewDevice().data();
  forAll< parallelDevicePolicy<> >( localSize(), [=] GEOSX_DEVICE ( localIndex const k )
  {
    localVector[k] = data[k];
  } );
}

MPI_Comm TpetraVector::getComm() const
{
  GEOSX_LAI_ASSERT( created() );
#ifdef GEOSX_USE_MPI
  return *dynamic_cast< Tpetra_Comm const & >( *m_vector->getMap()->getComm() ).getRawMpiComm();
#else
  return MPI_COMM_GEOSX;
#endif
}

void TpetraVector::print( std::ostream & os ) const
{
  GEOSX_LAI_ASSERT( ready() );
  Teuchos::RCP< Teuchos::FancyOStream > stream = Teuchos::getFancyOStream( Teuchos::rcp( &os, false ) );
  m_vector->describe( *stream, Teuchos::VERB_EXTREME );
}

void TpetraVector::write( string const & filename,
                          LAIOutputFormat format ) const
{
  GEOSX_LAI_ASSERT( ready() );
  switch( format )
  {
    case LAIOutputFormat::MATRIX_MARKET:
    {
      using Tpetra_CrsMatrix = Tpetra::CrsMatrix< Tpetra_Vector::scalar_type, Tpetra_Vector::local_ordinal_type, Tpetra_Vector::global_ordinal_type >;
      using Writer = Tpetra::MatrixMarket::Writer< Tpetra_CrsMatrix >;
      Writer::writeDenseFile( filename, *m_vector );
      break;
    }
    default:
    {
      GEOSX_ERROR( "Unsupported vector output format" );
    }
  }
}

TpetraVector::Tpetra_Vector const & TpetraVector::unwrapped() const
{
  GEOSX_LAI_ASSERT( created() );
  return *m_vector;
}

TpetraVector::Tpetra_Vector & TpetraVector::unwrapped()
{
  GEOSX_LAI_ASSERT( created() );
  return *m_vector;
}


} // namespace geosx
