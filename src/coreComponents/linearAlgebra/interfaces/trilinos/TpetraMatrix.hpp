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
 * @file TpetraMatrix.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_TPETRAMATRIX_HPP
#define GEOSX_LINEARALGEBRA_INTERFACES_TPETRAMATRIX_HPP

#include "linearAlgebra/interfaces/LinearOperator.hpp"
#include "linearAlgebra/interfaces/MatrixBase.hpp"
#include "linearAlgebra/interfaces/trilinos/TpetraVector.hpp"

#include <Tpetra_CrsMatrix_fwd.hpp>
#include <Tpetra_Map_fwd.hpp>

namespace geosx
{

/**
 * @brief Wrapper class for Trilinos/Tpetra's CrsMatrix class.
 */
class TpetraMatrix : public virtual LinearOperator< TpetraVector >,
  private MatrixBase< TpetraMatrix, TpetraVector >
{
public:

  /// Compatible vector type
  using Vector = TpetraVector;

  /**
   * @name Constructor/Destructor methods
   */
  ///@{

  /**
   * @brief Empty matrix constructor.
   *
   * Create an empty (distributed) matrix.
   */
  TpetraMatrix();

  /**
   * @brief Copy constructor.
   * @param[in] src the matrix to be copied
   */
  TpetraMatrix( TpetraMatrix const & src );

  /**
   * @brief Virtual destructor.
   */
  virtual ~TpetraMatrix() override;

  ///@}

  /**
   * @name MatrixBase interface
   */
  ///@{

  using MatrixBase::createWithLocalSize;
  using MatrixBase::createWithGlobalSize;
  using MatrixBase::create;
  using MatrixBase::closed;
  using MatrixBase::assembled;
  using MatrixBase::insertable;
  using MatrixBase::modifiable;
  using MatrixBase::ready;
  using MatrixBase::residual;

  virtual void createWithLocalSize( localIndex const localRows,
                                    localIndex const localCols,
                                    localIndex const maxEntriesPerRow,
                                    MPI_Comm const & comm ) override;

  virtual void createWithGlobalSize( globalIndex const globalRows,
                                     globalIndex const globalCols,
                                     localIndex const maxEntriesPerRow,
                                     MPI_Comm const & comm ) override;

  virtual void open() override;

  virtual void close() override;

  virtual bool created() const override;

  virtual void reset() override;

  virtual void set( real64 const value ) override;

  virtual void zero() override;

  virtual void add( globalIndex const rowIndex,
                    globalIndex const colIndex,
                    real64 const value ) override;

  virtual void set( globalIndex const rowIndex,
                    globalIndex const colIndex,
                    real64 const value ) override;

  virtual void insert( globalIndex const rowIndex,
                       globalIndex const colIndex,
                       real64 const value ) override;

  virtual void add( globalIndex const rowIndex,
                    globalIndex const * colIndices,
                    real64 const * values,
                    localIndex const size ) override;

  virtual void set( globalIndex const rowIndex,
                    globalIndex const * colIndices,
                    real64 const * values,
                    localIndex const size ) override;

  virtual void insert( globalIndex const rowIndex,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex const size ) override;

  virtual void add( globalIndex const rowIndex,
                    arraySlice1d< globalIndex const > const & colIndices,
                    arraySlice1d< real64 const > const & values ) override;

  virtual void set( globalIndex const rowIndex,
                    arraySlice1d< globalIndex const > const & colIndices,
                    arraySlice1d< real64 const > const & values ) override;

  virtual void insert( globalIndex const rowIndex,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice1d< real64 const > const & values ) override;

  virtual void add( arraySlice1d< globalIndex const > const & rowIndices,
                    arraySlice1d< globalIndex const > const & colIndices,
                    arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values ) override;

  virtual void set( arraySlice1d< globalIndex const > const & rowIndices,
                    arraySlice1d< globalIndex const > const & colIndices,
                    arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values ) override;

  virtual void insert( arraySlice1d< globalIndex const > const & rowIndices,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values ) override;

  virtual void add( globalIndex const * rowIndices,
                    globalIndex const * colIndices,
                    real64 const * values,
                    localIndex const numRows,
                    localIndex const numCols ) override;

  virtual void set( globalIndex const * rowIndices,
                    globalIndex const * colIndices,
                    real64 const * values,
                    localIndex const numRows,
                    localIndex const numCols ) override;

  virtual void insert( globalIndex const * rowIndices,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex const numRows,
                       localIndex const numCols ) override;

  virtual void apply( TpetraVector const & src,
                      TpetraVector & dst ) const override;

  virtual void applyTranspose( TpetraVector const & src,
                               TpetraVector & dst ) const override;

  virtual void multiply( TpetraMatrix const & src,
                         TpetraMatrix & dst ) const override;

  virtual void leftMultiplyTranspose( TpetraMatrix const & src,
                                      TpetraMatrix & dst ) const override;

  virtual void rightMultiplyTranspose( TpetraMatrix const & src,
                                       TpetraMatrix & dst ) const override;

  virtual void multiplyRAP( TpetraMatrix const & R,
                            TpetraMatrix const & P,
                            TpetraMatrix & dst ) const override;

  virtual void multiplyPtAP( TpetraMatrix const & P,
                             TpetraMatrix & dst ) const override;

  virtual void gemv( real64 const alpha,
                     TpetraVector const & x,
                     real64 const beta,
                     TpetraVector & y,
                     bool useTranspose = false ) const override;

  virtual void scale( real64 const scalingFactor ) override;

  virtual void leftScale( TpetraVector const & vec ) override;

  virtual void rightScale( TpetraVector const & vec ) override;

  virtual void leftRightScale( TpetraVector const & vecLeft,
                               TpetraVector const & vecRight ) override;

  virtual void transpose( TpetraMatrix & dst ) const override;

  virtual real64 clearRow( globalIndex const row,
                           bool const keepDiag = false,
                           real64 const diagValue = 0.0 ) override;

  virtual void addEntries( TpetraMatrix const & src, real64 const scale = 1.0 ) override;

  virtual void addDiagonal( TpetraVector const & src ) override;

  virtual localIndex maxRowLength() const override;

  virtual localIndex localRowLength( localIndex localRowIndex ) const override;

  virtual localIndex globalRowLength( globalIndex globalRowIndex ) const override;

  virtual void getRowCopy( globalIndex globalRow,
                           arraySlice1d< globalIndex > const & colIndices,
                           arraySlice1d< real64 > const & values ) const override;

  virtual real64 getDiagValue( globalIndex globalRow ) const override;

  virtual void extractDiagonal( TpetraVector & dst ) const override;

  virtual globalIndex numGlobalRows() const override;

  virtual globalIndex numGlobalCols() const override;

  virtual localIndex numLocalRows() const override;

  virtual localIndex numLocalCols() const override;

  virtual globalIndex ilower() const override;

  virtual globalIndex iupper() const override;

  virtual globalIndex jlower() const override;

  virtual globalIndex jupper() const override;

  virtual localIndex numLocalNonzeros() const override;

  virtual globalIndex numGlobalNonzeros() const override;

  virtual real64 normInf() const override;

  virtual real64 norm1() const override;

  virtual real64 normFrobenius() const override;

  virtual localIndex getLocalRowID( globalIndex const globalRow ) const override;

  virtual globalIndex getGlobalRowID( localIndex const localRow ) const override;

  virtual MPI_Comm getComm() const override;

  virtual void print( std::ostream & os = std::cout ) const override;

  virtual void write( string const & filename,
                      LAIOutputFormat const format = LAIOutputFormat::MATRIX_MARKET ) const override;

  ///@}

  /**
   * @brief Alias for Tpetra map template instantiation used by this class.
   */
  using Tpetra_Map = Tpetra::Map< int, globalIndex >;

  /**
   * @brief Alias for specific Tpetra matrix template instantiation wrapped by this class.
   *
   * @note This uses Tpetra's default execution/memory space. When built with CUDA support,
   * this will be equal to Kokkos::Cuda, so we won't be able to create a host-only vector.
   * If we want both in the same executable, we'll have to make adjustments to our LAI approach.
   */
  using Tpetra_CrsMatrix = Tpetra::CrsMatrix< real64, int, globalIndex >;

  /**
   * @brief Returns a const pointer to the underlying matrix.
   * @return const pointer to the underlying matrix
   */
  Tpetra_CrsMatrix const & unwrapped() const;

  /**
   * @brief Returns a non-const pointer to the underlying matrix.
   * @return non-const pointer to the underlying matrix
   */
  Tpetra_CrsMatrix & unwrapped();

private:

  /**
   * @brief Perform a matrix matrix product with Parallel Matrix
   */
  void multiply( bool const transA,
                 TpetraMatrix const & B,
                 bool const transB,
                 TpetraMatrix & C ) const;

  /**
   * @brief Create the matrix by copying data from an Epetra_CrsMatrix
   * @param src the source matrix
   */
  void create( Tpetra_CrsMatrix const & src );

  /// Pointer to the underlying Epetra_CrsMatrix.
  std::unique_ptr< Tpetra_CrsMatrix > m_matrix;

  /// Map representing the parallel partitioning of a source vector (x in y=Ax)
  std::unique_ptr< Tpetra_Map > m_src_map;

  /// Map representing the parallel partitioning of a destination vector (y in y=Ax)
  std::unique_ptr< Tpetra_Map > m_dst_map;
};

} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_INTERFACES_TPETRAMATRIX_HPP
