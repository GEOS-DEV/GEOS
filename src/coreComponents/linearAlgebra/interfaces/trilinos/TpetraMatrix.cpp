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
 * @file TpetraMatrix.cpp
 */

#include "TpetraMatrix.hpp"

#include "codingUtilities/Utilities.hpp"
#include "linearAlgebra/interfaces/trilinos/TpetraUtils.hpp"

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_RowMatrixTransposer.hpp>
#include <TpetraExt_MatrixMatrix.hpp>
#include <TpetraExt_TripleMatrixMultiply.hpp>
#include <MatrixMarket_Tpetra.hpp>

namespace geosx
{

TpetraMatrix::TpetraMatrix()
  : LinearOperator(),
  MatrixBase()
{ }

TpetraMatrix::TpetraMatrix( TpetraMatrix const & src )
  : TpetraMatrix()
{
  GEOSX_LAI_ASSERT( src.ready() );
  m_matrix = std::make_unique< Tpetra_CrsMatrix >( *src.m_matrix );
  m_src_map = std::make_unique< Tpetra_Map >( *m_matrix->getDomainMap() );
  m_dst_map = std::make_unique< Tpetra_Map >( *m_matrix->getRangeMap() );
  m_assembled = true;
}

TpetraMatrix::~TpetraMatrix() = default;

void TpetraMatrix::createWithGlobalSize( globalIndex const globalRows,
                                         globalIndex const globalCols,
                                         localIndex const maxEntriesPerRow,
                                         MPI_Comm const & comm )
{
  GEOSX_LAI_ASSERT( closed() );
  GEOSX_LAI_ASSERT_GE( globalRows, 0 );
  GEOSX_LAI_ASSERT_GE( globalCols, 0 );
  GEOSX_LAI_ASSERT_GE( maxEntriesPerRow, 0 );

  reset();

  Teuchos::RCP< Tpetra_Comm > tpetraComm( new Tpetra_Comm( MPI_PARAM( comm ) ) );
  m_dst_map = std::make_unique< Tpetra_Map >( LvArray::integerConversion< Tpetra::global_size_t >( globalRows ),
                                              0,
                                              tpetraComm );
  m_src_map = std::make_unique< Tpetra_Map >( LvArray::integerConversion< Tpetra::global_size_t >( globalCols ),
                                              0,
                                              tpetraComm );
  m_matrix = std::make_unique< Tpetra_CrsMatrix >( Teuchos::RCP< Tpetra_Map >( m_dst_map.get(), false ),
                                                   LvArray::integerConversion< size_t >( maxEntriesPerRow ) );
}

void TpetraMatrix::createWithLocalSize( localIndex const localRows,
                                        localIndex const localCols,
                                        localIndex const maxEntriesPerRow,
                                        MPI_Comm const & comm )
{
  GEOSX_LAI_ASSERT( closed() );
  GEOSX_LAI_ASSERT_GE( localRows, 0 );
  GEOSX_LAI_ASSERT_GE( localCols, 0 );
  GEOSX_LAI_ASSERT_GE( maxEntriesPerRow, 0 );

  reset();

  Teuchos::RCP< Tpetra_Comm > tpetraComm( new Tpetra_Comm( MPI_PARAM( comm ) ) );
  m_dst_map = std::make_unique< Tpetra_Map >( Teuchos::OrdinalTraits< Tpetra::global_size_t >::invalid(),
                                              LvArray::integerConversion< size_t >( localRows ),
                                              0,
                                              tpetraComm );
  m_src_map = std::make_unique< Tpetra_Map >( Teuchos::OrdinalTraits< Tpetra::global_size_t >::invalid(),
                                              LvArray::integerConversion< size_t >( localCols ),
                                              0,
                                              tpetraComm );
  m_matrix = std::make_unique< Tpetra_CrsMatrix >( Teuchos::RCP< Tpetra_Map >( m_dst_map.get(), false ),
                                                   LvArray::integerConversion< size_t >( maxEntriesPerRow ) );
}

bool TpetraMatrix::created() const
{
  return bool(m_matrix);
}

void TpetraMatrix::reset()
{
  MatrixBase::reset();
  m_matrix.reset();
  m_dst_map.reset();
  m_src_map.reset();
}

void TpetraMatrix::set( real64 const value )
{
  GEOSX_LAI_ASSERT( ready() );
  open();
  m_matrix->setAllToScalar( value );
  close();
}

void TpetraMatrix::zero()
{
  set( 0 );
}

void TpetraMatrix::open()
{
  GEOSX_LAI_ASSERT( created() && closed() );
  m_matrix->resumeFill();
  m_closed = false;
}

void TpetraMatrix::close()
{
  GEOSX_LAI_ASSERT( !closed() );
  // TODO: pass a Teuchos::ParameterList( { "No Nonlocal Changes", true } )
  m_matrix->fillComplete( Teuchos::RCP< Tpetra_Map >( m_src_map.get(), false ),
                          Teuchos::RCP< Tpetra_Map >( m_dst_map.get(), false ) );
  m_assembled = true;
  m_closed = true;
}

void TpetraMatrix::add( globalIndex const rowIndex,
                        globalIndex const colIndex,
                        real64 const value )
{
  GEOSX_LAI_ASSERT( modifiable() );
  localIndex const numEntries = m_matrix->sumIntoGlobalValues( rowIndex, 1, &value, &colIndex );
  GEOSX_LAI_ASSERT_EQ( numEntries, 1 );
}

void TpetraMatrix::set( globalIndex const rowIndex,
                        globalIndex const colIndex,
                        real64 const value )
{
  GEOSX_LAI_ASSERT( modifiable() );
  localIndex const numEntries = m_matrix->replaceGlobalValues( rowIndex, 1, &value, &colIndex );
  GEOSX_LAI_ASSERT_EQ( numEntries, 1 );
}

void TpetraMatrix::insert( globalIndex const rowIndex,
                           globalIndex const colIndex,
                           real64 const value )
{
  GEOSX_LAI_ASSERT( insertable() );
  m_matrix->insertGlobalValues( rowIndex, 1, &value, &colIndex );
}

void TpetraMatrix::add( globalIndex const rowIndex,
                        globalIndex const * colIndices,
                        real64 const * values,
                        localIndex const size )
{
  GEOSX_LAI_ASSERT( modifiable() );
  localIndex const numEntries = m_matrix->sumIntoGlobalValues( rowIndex, size, values, colIndices );
  GEOSX_LAI_ASSERT_EQ( numEntries, size );
}

void TpetraMatrix::set( globalIndex const rowIndex,
                        globalIndex const * colIndices,
                        real64 const * values,
                        localIndex const size )
{
  GEOSX_LAI_ASSERT( modifiable() );
  localIndex const numEntries = m_matrix->replaceGlobalValues( rowIndex, size, values, colIndices );
  GEOSX_LAI_ASSERT_EQ( numEntries, size );
}

void TpetraMatrix::insert( globalIndex const rowIndex,
                           globalIndex const * colIndices,
                           real64 const * values,
                           localIndex const size )
{
  GEOSX_LAI_ASSERT( insertable() );
  m_matrix->insertGlobalValues( rowIndex, size, values, colIndices );
}

void TpetraMatrix::add( globalIndex const rowIndex,
                        arraySlice1d< globalIndex const > const & colIndices,
                        arraySlice1d< real64 const > const & values )
{
  GEOSX_LAI_ASSERT( modifiable() );
  localIndex const numEntries = m_matrix->sumIntoGlobalValues( rowIndex, colIndices.size(), values, colIndices );
  GEOSX_LAI_ASSERT_EQ( numEntries, colIndices.size() );
}

void TpetraMatrix::set( globalIndex const rowIndex,
                        arraySlice1d< globalIndex const > const & colIndices,
                        arraySlice1d< real64 const > const & values )
{
  GEOSX_LAI_ASSERT( modifiable() );
  localIndex const numEntries = m_matrix->replaceGlobalValues( rowIndex, colIndices.size(), values, colIndices );
  GEOSX_LAI_ASSERT_EQ( numEntries, colIndices.size() );
}

void TpetraMatrix::insert( globalIndex const rowIndex,
                           arraySlice1d< globalIndex const > const & colIndices,
                           arraySlice1d< real64 const > const & values )
{
  GEOSX_LAI_ASSERT( insertable() );
  m_matrix->insertGlobalValues( rowIndex, colIndices.size(), values, colIndices );
}

void TpetraMatrix::add( arraySlice1d< globalIndex const > const & rowIndices,
                        arraySlice1d< globalIndex const > const & colIndices,
                        arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values )
{
  GEOSX_LAI_ASSERT( modifiable() );
  for( localIndex i = 0; i < rowIndices.size(); ++i )
  {
    add( rowIndices[i], colIndices, values[i] );
  }
}

void TpetraMatrix::set( arraySlice1d< globalIndex const > const & rowIndices,
                        arraySlice1d< globalIndex const > const & colIndices,
                        arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values )
{
  GEOSX_LAI_ASSERT( modifiable() );
  for( localIndex i = 0; i < rowIndices.size(); ++i )
  {
    set( rowIndices[i], colIndices, values[i] );
  }
}

void TpetraMatrix::insert( arraySlice1d< globalIndex const > const & rowIndices,
                           arraySlice1d< globalIndex const > const & colIndices,
                           arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values )
{
  for( localIndex i = 0; i < rowIndices.size(); ++i )
  {
    insert( rowIndices[i], colIndices, values[i] );
  }
}

void TpetraMatrix::add( globalIndex const * rowIndices,
                        globalIndex const * colIndices,
                        real64 const * values,
                        localIndex const numRows,
                        localIndex const numCols )
{
  for( localIndex i = 0; i < numRows; ++i )
  {
    add( rowIndices[i], colIndices, values + i * numCols, numCols );
  }
}

void TpetraMatrix::set( globalIndex const * rowIndices,
                        globalIndex const * colIndices,
                        real64 const * values,
                        localIndex const numRows,
                        localIndex const numCols )
{
  for( localIndex i = 0; i < numRows; ++i )
  {
    set( rowIndices[i], colIndices, values + i * numCols, numCols );
  }
}

void TpetraMatrix::insert( globalIndex const * rowIndices,
                           globalIndex const * colIndices,
                           real64 const * values,
                           localIndex const numRows,
                           localIndex const numCols )
{
  for( localIndex i = 0; i < numRows; ++i )
  {
    insert( rowIndices[i], colIndices, values + i * numCols, numCols );
  }
}

void TpetraMatrix::apply( TpetraVector const & src,
                          TpetraVector & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( dst.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalRows(), dst.globalSize() );
  GEOSX_LAI_ASSERT_EQ( numGlobalCols(), src.globalSize() );

  m_matrix->apply( src.unwrapped(),
                   dst.unwrapped(),
                   Teuchos::ETransp::NO_TRANS,
                   1.0,
                   0.0 );
}

void TpetraMatrix::applyTranspose( TpetraVector const & src,
                                   TpetraVector & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( dst.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalCols(), dst.globalSize() );
  GEOSX_LAI_ASSERT_EQ( numGlobalRows(), src.globalSize() );

  m_matrix->apply( src.unwrapped(),
                   dst.unwrapped(),
                   Teuchos::ETransp::TRANS,
                   1.0,
                   0.0 );
}

void TpetraMatrix::multiply( TpetraMatrix const & src,
                             TpetraMatrix & dst ) const
{
  this->multiply( false, src, false, dst );
}

void TpetraMatrix::leftMultiplyTranspose( TpetraMatrix const & src,
                                          TpetraMatrix & dst ) const
{
  this->multiply( true, src, false, dst );
}

void TpetraMatrix::rightMultiplyTranspose( TpetraMatrix const & src,
                                           TpetraMatrix & dst ) const
{
  src.multiply( false, *this, true, dst );
}

void TpetraMatrix::create( Tpetra_CrsMatrix const & src )
{
  GEOSX_LAI_ASSERT( closed() );
  reset();

  m_matrix = std::make_unique< Tpetra_CrsMatrix >( src );
  m_src_map = std::make_unique< Tpetra_Map >( *m_matrix->getDomainMap() );
  m_dst_map = std::make_unique< Tpetra_Map >( *m_matrix->getRangeMap() );
  m_assembled = true;
}

void TpetraMatrix::multiplyRAP( TpetraMatrix const & R,
                                TpetraMatrix const & P,
                                TpetraMatrix & dst ) const
{
  // TODO: Tpetra::TripleMatrixMultiply::MultiplyRAP gives wrong results, investigate
#if 0
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( R.ready() );
  GEOSX_LAI_ASSERT( P.ready() );
  GEOSX_LAI_ASSERT_EQ( R.numLocalCols(), numLocalRows() );
  GEOSX_LAI_ASSERT_EQ( P.numLocalRows(), numLocalCols() );

  // TODO: estimate num nonzeros per row? or it does not matter since matrix is recreated anyway?
  dst.createWithLocalSize( R.numLocalRows(), P.numLocalCols(), 1, getComm() );

  Tpetra::TripleMatrixMultiply::MultiplyRAP( R.unwrapped(),
                                             false,
                                             *m_matrix,
                                             false,
                                             P.unwrapped(),
                                             false,
                                             dst.unwrapped() );
  dst.m_assembled = true;
#else
  MatrixBase::multiplyRAP( R, P, dst );
#endif
}

void TpetraMatrix::multiplyPtAP( TpetraMatrix const & P,
                                 TpetraMatrix & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( P.ready() );
  GEOSX_LAI_ASSERT_EQ( P.numLocalRows(), numLocalRows() );
  GEOSX_LAI_ASSERT_EQ( P.numLocalRows(), numLocalCols() );

  // TODO: estimate num nonzeros per row? or it does not matter since matrix is recreated anyway?
  dst.createWithLocalSize( P.numLocalCols(), P.numLocalCols(), 1, getComm() );

  Tpetra::TripleMatrixMultiply::MultiplyRAP( P.unwrapped(),
                                             true,
                                             *m_matrix,
                                             false,
                                             P.unwrapped(),
                                             false,
                                             dst.unwrapped() );
  dst.m_assembled = true;
}

void TpetraMatrix::gemv( real64 const alpha,
                         TpetraVector const & x,
                         real64 const beta,
                         TpetraVector & y,
                         bool useTranspose ) const
{
  GEOSX_LAI_ASSERT( ready() );
  m_matrix->apply( x.unwrapped(),
                   y.unwrapped(),
                   useTranspose ? Teuchos::ETransp::TRANS : Teuchos::ETransp::NO_TRANS,
                   alpha,
                   beta );
}

void TpetraMatrix::scale( real64 const scalingFactor )
{
  GEOSX_LAI_ASSERT( ready() );

  if( isEqual( scalingFactor, 1.0 ) )
  {
    return;
  }

  open();
  m_matrix->scale( scalingFactor );
  close();
}

void TpetraMatrix::leftScale( TpetraVector const & vec )
{
  GEOSX_LAI_ASSERT( ready() );
  m_matrix->leftScale( vec.unwrapped() );
}

void TpetraMatrix::rightScale( TpetraVector const & vec )
{
  GEOSX_LAI_ASSERT( ready() );
  m_matrix->rightScale( vec.unwrapped() );
}

void TpetraMatrix::leftRightScale( TpetraVector const & vecLeft,
                                   TpetraVector const & vecRight )
{
  leftScale( vecLeft );
  rightScale( vecRight );
}

void TpetraMatrix::transpose( TpetraMatrix & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );

  using TpetraTransposer = Tpetra::RowMatrixTransposer< Tpetra_CrsMatrix::scalar_type,
                                                        Tpetra_CrsMatrix::local_ordinal_type,
                                                        Tpetra_CrsMatrix::global_ordinal_type >;

  TpetraTransposer transposer( Teuchos::RCP< Tpetra_CrsMatrix >( m_matrix.get(), false ) );
  Teuchos::RCP< Tpetra_CrsMatrix > trans = transposer.createTranspose();

  dst.reset();
  dst.m_matrix.reset( trans.release().get() );
  dst.m_src_map = std::make_unique< Tpetra_Map >( *dst.m_matrix->getDomainMap() );
  dst.m_dst_map = std::make_unique< Tpetra_Map >( *dst.m_matrix->getRangeMap() );
  dst.m_assembled = true;
}

real64 TpetraMatrix::clearRow( globalIndex const globalRow,
                               bool const keepDiag,
                               real64 const diagValue )
{
  // TODO: Deprecate/remove this method?

  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_ASSERT_GE( globalRow, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), globalRow );

  Teuchos::ArrayView< Tpetra_CrsMatrix::local_ordinal_type const > indicesView;
  Teuchos::ArrayView< Tpetra_CrsMatrix::scalar_type const > valuesView;
  m_matrix->getLocalRowView( getLocalRowID( globalRow ), indicesView, valuesView );

  array1d< globalIndex > colIndices( indicesView.size() );
  array1d< real64 > values( valuesView.size() );

  bool const square = numGlobalRows() == numGlobalCols();

  real64 oldDiag = 0.0;
  for( localIndex j = 0; j < indicesView.size(); ++j )
  {
    colIndices[j] = m_matrix->getColMap()->getGlobalElement( indicesView[j] );
    if( square && colIndices[j] == globalRow )
    {
      oldDiag = valuesView[j];
      values[j] = keepDiag ? oldDiag : diagValue;
    }
  }

  set( globalRow, colIndices, values );
  return oldDiag;
}

void TpetraMatrix::addEntries( TpetraMatrix const & src, real64 const scale )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( numGlobalRows() == src.numGlobalRows() );
  GEOSX_LAI_ASSERT( numGlobalCols() == src.numGlobalCols() );

  open();
  Tpetra::MatrixMatrix::Add( src.unwrapped(), false, scale, *m_matrix, 1.0 );
  close();
}

void TpetraMatrix::addDiagonal( TpetraVector const & src )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( numGlobalRows() == numGlobalCols() );
  GEOSX_LAI_ASSERT( numLocalRows() == src.localSize() );

  open();

  // TODO: this is only correct if values are on host
  real64 const * const localValues = src.extractLocalVector();
  for( localIndex i = 0; i < numLocalRows(); ++i )
  {
    globalIndex const globalRow = getGlobalRowID( i );
    add( globalRow, globalRow, localValues[i] );
  }

  close();
}

localIndex TpetraMatrix::maxRowLength() const
{
  GEOSX_LAI_ASSERT( assembled() );
  return m_matrix->getGlobalMaxNumRowEntries();
}

localIndex TpetraMatrix::localRowLength( localIndex const localRowIndex ) const
{
  GEOSX_LAI_ASSERT( assembled() );
  return m_matrix->getNumEntriesInLocalRow( localRowIndex );
}

localIndex TpetraMatrix::globalRowLength( globalIndex const globalRowIndex ) const
{
  GEOSX_LAI_ASSERT( assembled() );
  return m_matrix->getNumEntriesInGlobalRow( globalRowIndex );
}

void TpetraMatrix::getRowCopy( globalIndex const globalRow,
                               arraySlice1d< globalIndex > const & colIndices,
                               arraySlice1d< real64 > const & values ) const
{
  // TODO: Deprecate/remove this method?

  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT_GE( globalRow, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), globalRow );

  Teuchos::ArrayView< Tpetra_CrsMatrix::local_ordinal_type const > indicesView;
  Teuchos::ArrayView< Tpetra_CrsMatrix::scalar_type const > valuesView;
  m_matrix->getLocalRowView( getLocalRowID( globalRow ), indicesView, valuesView );

  GEOSX_LAI_ASSERT_GE( colIndices.size(), indicesView.size() );
  GEOSX_LAI_ASSERT_GE( values.size(), indicesView.size() );

  localIndex const numEntries = LvArray::integerConversion< localIndex >( indicesView.size() );
  for( localIndex i = 0; i < numEntries; ++i )
  {
    colIndices[i] = LvArray::integerConversion< globalIndex >( m_src_map->getGlobalElement( indicesView[i] ) );
    values[i] = valuesView[i];
  }
}

real64 TpetraMatrix::getDiagValue( globalIndex const globalRow ) const
{
  // TODO: Deprecate/remove this method?

  GEOSX_LAI_ASSERT( assembled() );
  GEOSX_LAI_ASSERT_GE( globalRow, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), globalRow );

  Teuchos::ArrayView< Tpetra_CrsMatrix::local_ordinal_type const > indicesView;
  Teuchos::ArrayView< Tpetra_CrsMatrix::scalar_type const > valuesView;
  m_matrix->getLocalRowView( getLocalRowID( globalRow ), indicesView, valuesView );

  localIndex const numEntries = LvArray::integerConversion< localIndex >( indicesView.size() );
  for( localIndex j = 0; j < numEntries; ++j )
  {
    if( m_src_map->getGlobalElement( indicesView[j] ) == globalRow )
    {
      return valuesView[j];
    }
  }

  return 0.0;
}

void TpetraMatrix::extractDiagonal( TpetraVector & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( dst.ready() );
  GEOSX_LAI_ASSERT_EQ( dst.localSize(), numLocalRows() );

  m_matrix->getLocalDiagCopy( dst.unwrapped() );
}

TpetraMatrix::Tpetra_CrsMatrix const & TpetraMatrix::unwrapped() const
{
  GEOSX_LAI_ASSERT( created() );
  return *m_matrix;
}

TpetraMatrix::Tpetra_CrsMatrix & TpetraMatrix::unwrapped()
{
  GEOSX_LAI_ASSERT( created() );
  return *m_matrix;
}

globalIndex TpetraMatrix::numGlobalRows() const
{
  GEOSX_LAI_ASSERT( created() );
  return LvArray::integerConversion< globalIndex >( m_dst_map->getGlobalNumElements() );
}

globalIndex TpetraMatrix::numGlobalCols() const
{
  GEOSX_LAI_ASSERT( created() );
  return LvArray::integerConversion< globalIndex >( m_src_map->getGlobalNumElements() );
}

globalIndex TpetraMatrix::ilower() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_dst_map->getMinGlobalIndex();
}

globalIndex TpetraMatrix::iupper() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_dst_map->getMaxGlobalIndex() + 1;
}

globalIndex TpetraMatrix::jlower() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_src_map->getMinGlobalIndex();
}

globalIndex TpetraMatrix::jupper() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_src_map->getMaxGlobalIndex() + 1;
}

localIndex TpetraMatrix::numLocalNonzeros() const
{
  GEOSX_LAI_ASSERT( assembled() );
  return LvArray::integerConversion< localIndex >( m_matrix->getNodeNumEntries() );
}

globalIndex TpetraMatrix::numGlobalNonzeros() const
{
  GEOSX_LAI_ASSERT( assembled() );
  return LvArray::integerConversion< globalIndex >( m_matrix->getGlobalNumEntries() );
}

real64 TpetraMatrix::normInf() const
{
  GEOSX_LAI_ASSERT( ready() );

  using LocalMatrix = Tpetra_CrsMatrix::local_matrix_type;
  LocalMatrix const & localMatrix = m_matrix->getLocalMatrix();

  real64 localNorm = 0.0;
  Kokkos::parallel_reduce( Kokkos::RangePolicy< Tpetra_CrsMatrix::execution_space >( 0, localMatrix.numRows() ),
                           KOKKOS_LAMBDA ( localIndex const i, real64 & maxVal )
    {
      KokkosSparse::SparseRowViewConst< LocalMatrix > const rowView = localMatrix.rowConst( i );
      real64 rowSumAbs = 0.0;
      for( int k = 0; k < rowView.length; ++k )
      {
        rowSumAbs += fabs( rowView.value( k ) );
      }
      maxVal = fmax( maxVal, rowSumAbs );
    }, Kokkos::Max< real64, Kokkos::HostSpace >( localNorm ) );

  return MpiWrapper::Max( localNorm, getComm() );
}

real64 TpetraMatrix::norm1() const
{
  GEOSX_LAI_ASSERT( ready() );

  TpetraMatrix matTrans;
  transpose( matTrans );
  return matTrans.normInf();
}

real64 TpetraMatrix::normFrobenius() const
{
  GEOSX_LAI_ASSERT( ready() );
  return m_matrix->getFrobeniusNorm();
}

localIndex TpetraMatrix::getLocalRowID( globalIndex const globalRow ) const
{
  GEOSX_LAI_ASSERT( created() );
  return m_dst_map->getLocalElement( globalRow );
}

globalIndex TpetraMatrix::getGlobalRowID( localIndex const localRow ) const
{
  GEOSX_LAI_ASSERT( created() );
  GEOSX_LAI_ASSERT_GE( localRow, 0 );
  GEOSX_LAI_ASSERT_GT( numLocalRows(), localRow );
  return m_dst_map->getGlobalElement( localRow );
}

localIndex TpetraMatrix::numLocalCols() const
{
  GEOSX_LAI_ASSERT( created() );
  return LvArray::integerConversion< localIndex >( m_src_map->getNodeNumElements() );
}

localIndex TpetraMatrix::numLocalRows() const
{
  GEOSX_LAI_ASSERT( created() );
  return LvArray::integerConversion< localIndex >( m_dst_map->getNodeNumElements() );
}

MPI_Comm TpetraMatrix::getComm() const
{
  GEOSX_LAI_ASSERT( created() );
#ifdef GEOSX_USE_MPI
  return *dynamic_cast< Tpetra_Comm const & >( *m_matrix->getRowMap()->getComm() ).getRawMpiComm();
#else
  return MPI_COMM_GEOSX;
#endif
}

void TpetraMatrix::print( std::ostream & os ) const
{
  GEOSX_LAI_ASSERT( ready() );
  Teuchos::RCP< Teuchos::FancyOStream > stream = Teuchos::getFancyOStream( Teuchos::rcp( &os, false ) );
  m_matrix->describe( *stream, Teuchos::VERB_EXTREME );
}

void TpetraMatrix::write( string const & filename,
                          LAIOutputFormat const format ) const
{
  GEOSX_LAI_ASSERT( ready() );
  switch( format )
  {
    case LAIOutputFormat::NATIVE_ASCII:
    {
      std::ofstream ofs( filename );
      print( ofs );
      break;
    }
    case LAIOutputFormat::MATRIX_MARKET:
    {
      using Writer = Tpetra::MatrixMarket::Writer< Tpetra_CrsMatrix >;
      Writer::writeSparseFile( filename, *m_matrix );
      break;
    }
    default:
    {
      GEOSX_ERROR( "Unsupported matrix output format" );
    }
  }
}

void TpetraMatrix::multiply( bool const transA,
                             TpetraMatrix const & B,
                             bool const transB,
                             TpetraMatrix & C ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( B.ready() );

  C.createWithLocalSize( transA ? numLocalCols() : numLocalRows(),
                         transB ? B.numLocalRows() : B.numLocalCols(),
                         1, // TODO: estimate entries per row?
                         getComm() );

  Tpetra::MatrixMatrix::Multiply( unwrapped(),
                                  transA,
                                  B.unwrapped(),
                                  transB,
                                  C.unwrapped(),
                                  true );
  C.m_assembled = true;
}

} // end geosx namespace
