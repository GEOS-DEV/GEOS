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
 * @file EpetraMatrix.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_EPETRAMATRIX_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_EPETRAMATRIX_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/EpetraVector.hpp"
#include "linearAlgebra/interfaces/LinearOperator.hpp"
#include "linearAlgebra/interfaces/MatrixBase.hpp"

class Epetra_Map;
class Epetra_FECrsMatrix;

namespace geosx
{

/**
 * \class EpetraMatrix
 * \brief This class creates and provides basic support for the Epetra_CrsMatrix
 *        matrix object type used in Trilinos.
 */
class EpetraMatrix final: public LinearOperator<EpetraVector>,
                          private MatrixBase<EpetraMatrix, EpetraVector>
{
public:

  /**
   * @name Constructor/Destructor methods
   */
  ///@{

  /**
   * @brief Empty matrix constructor.
   *
   * Create an empty (distributed) matrix.
   */
  EpetraMatrix();

  /**
   * @brief Copy constructor.
   *
   * Create new matrix from matrix <tt>src</tt>.
   */
  EpetraMatrix( EpetraMatrix const & src );

  /**
   * @brief Virtual destructor.
   */
  virtual ~EpetraMatrix() override;

  ///@}

  /**
   * @name MatrixBase interface
   */
  ///@{

  using MatrixBase::createWithLocalSize;
  using MatrixBase::createWithGlobalSize;
  using MatrixBase::closed;
  using MatrixBase::assembled;
  using MatrixBase::insertable;
  using MatrixBase::modifiable;
  using MatrixBase::ready;

  void createWithLocalSize( localIndex const localRows,
                            localIndex const localCols,
                            localIndex const maxEntriesPerRow,
                            MPI_Comm const & comm ) override;

  void createWithGlobalSize( globalIndex const globalRows,
                             globalIndex const globalCols,
                             localIndex const maxEntriesPerRow,
                             MPI_Comm const & comm ) override;

  void open() override;

  void close() override;

  bool created() const override;

  void reset() override;

  void set( real64 const value ) override;

  void zero() override;

  void add( globalIndex const rowIndex,
            globalIndex const colIndex,
            real64 const value ) override;

  void set( globalIndex const rowIndex,
            globalIndex const colIndex,
            real64 const value ) override;

  void insert( globalIndex const rowIndex,
               globalIndex const colIndex,
               real64 const value ) override;

  void add( globalIndex const rowIndex,
            globalIndex const * colIndices,
            real64 const * values,
            localIndex const size ) override;

  void set( globalIndex const rowIndex,
            globalIndex const * colIndices,
            real64 const * values,
            localIndex const size ) override;

  void insert( globalIndex const rowIndex,
               globalIndex const * colIndices,
               real64 const * values,
               localIndex const size ) override;

  void add( globalIndex const rowIndex,
            arraySlice1d< globalIndex const > const & colIndices,
            arraySlice1d< real64 const > const & values ) override;

  void set( globalIndex const rowIndex,
            arraySlice1d< globalIndex const > const & colIndices,
            arraySlice1d< real64 const > const & values ) override;

  void insert( globalIndex const rowIndex,
               arraySlice1d< globalIndex const > const & colIndices,
               arraySlice1d< real64 const > const & values ) override;

  void add( arraySlice1d< globalIndex const > const & rowIndices,
            arraySlice1d< globalIndex const > const & colIndices,
            arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values ) override;

  void set( arraySlice1d< globalIndex const > const & rowIndices,
            arraySlice1d< globalIndex const > const & colIndices,
            arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values ) override;

  void insert( arraySlice1d< globalIndex const > const & rowIndices,
               arraySlice1d< globalIndex const > const & colIndices,
               arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values ) override;

  void add( arraySlice1d< globalIndex const > const & rowIndices,
            arraySlice1d< globalIndex const > const & colIndices,
            arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & values ) override;

  void set( arraySlice1d< globalIndex const > const & rowIndices,
            arraySlice1d< globalIndex const > const & colIndices,
            arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & values ) override;

  void insert( arraySlice1d< globalIndex const > const & rowIndices,
               arraySlice1d< globalIndex const > const & colIndices,
               arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & values ) override;

  void add( globalIndex const * rowIndices,
            globalIndex const * colIndices,
            real64 const * values,
            localIndex const numRows,
            localIndex const numCols ) override;

  void set( globalIndex const * rowIndices,
            globalIndex const * colIndices,
            real64 const * values,
            localIndex const numRows,
            localIndex const numCols ) override;

  void insert( globalIndex const * rowIndices,
               globalIndex const * colIndices,
               real64 const * values,
               localIndex const numRows,
               localIndex const numCols ) override;

  void multiply( EpetraVector const & src,
                 EpetraVector & dst ) const override;

  void multiply( EpetraMatrix const & src,
                 EpetraMatrix & dst,
                 bool const closeResult = true ) const override;

  void leftMultiplyTranspose( EpetraMatrix const & src,
                              EpetraMatrix & dst,
                              bool const closeResult = true ) const override;

  void rightMultiplyTranspose( EpetraMatrix const & src,
                               EpetraMatrix & dst,
                               bool const closeResult = true ) const override;

  void gemv( real64 const alpha,
             EpetraVector const & x,
             real64 const beta,
             EpetraVector & y,
             bool useTranspose = false ) const override;

  void scale( real64 const scalingFactor ) override;

  void leftScale( EpetraVector const & vec ) override;

  void rightScale( EpetraVector const & vec ) override;

  void leftRightScale( EpetraVector const & vecLeft,
                       EpetraVector const & vecRight ) override;

  void clearRow( globalIndex const row,
                 real64 const diagValue = 0 ) override;

  localIndex maxRowLength() const override;

  localIndex localRowLength( localIndex localRowIndex ) const override;

  localIndex globalRowLength( globalIndex globalRowIndex ) const override;

  void getRowCopy( globalIndex globalRow,
                   arraySlice1d< globalIndex > const & colIndices,
                   arraySlice1d< real64 > const & values ) const override;

  real64 getDiagValue( globalIndex globalRow ) const override;

  globalIndex numGlobalRows() const override;

  globalIndex numGlobalCols() const override;

  localIndex numLocalRows() const override;

  localIndex numLocalCols() const override;

  globalIndex ilower() const override;

  globalIndex iupper() const override;

  localIndex numLocalNonzeros() const override;

  globalIndex numGlobalNonzeros() const override;

  real64 normInf() const override;

  real64 norm1() const override;

  real64 normFrobenius() const override;

  localIndex getLocalRowID( globalIndex const index ) const override;

  globalIndex getGlobalRowID( localIndex const index ) const override;

  MPI_Comm getComm() const override;

  void print( std::ostream & os = std::cout ) const override;

  void write( string const & filename,
              LAIOutputFormat const format ) const override;

  ///@}

  /**
   * @brief Returns a pointer to the underlying matrix.
   */
  Epetra_FECrsMatrix const & unwrapped() const;

  /**
   * @brief Returns a pointer to the underlying matrix.
   */
  Epetra_FECrsMatrix & unwrapped();

private:

  /**
   * @brief Perform a matrix matrix product with Parallel Matrix
   */
  void multiply( bool const transA,
                 EpetraMatrix const & B,
                 bool const transB,
                 EpetraMatrix & C,
                 bool const closeResult ) const;

  /// Pointer to the underlying Epetra_CrsMatrix.
  std::unique_ptr< Epetra_FECrsMatrix > m_matrix;

  /// Map representing the parallel partitioning of a source vector (x in y=Ax)
  std::unique_ptr< Epetra_Map > m_src_map;

  /// Map representing the parallel partitioning of a destination vector (y in y=Ax)
  std::unique_ptr< Epetra_Map > m_dst_map;
};

} // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_EPETRAMATRIX_HPP_*/
