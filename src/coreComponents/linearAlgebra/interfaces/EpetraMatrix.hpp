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
class EpetraMatrix : public MatrixBase<EpetraMatrix, EpetraVector>
{
public:

  using Base = MatrixBase<EpetraMatrix, EpetraVector>;

  //! @name Constructor/Destructor Methods
  //@{

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

  void createWithLocalSize( localIndex const localSize,
                            localIndex const maxEntriesPerRow,
                            MPI_Comm const & comm ) final;

  void createWithGlobalSize( globalIndex const globalSize,
                             localIndex const maxEntriesPerRow,
                             MPI_Comm const & comm ) final;

  void createWithLocalSize( localIndex const localRows,
                            localIndex const localCols,
                            localIndex const maxEntriesPerRow,
                            MPI_Comm const & comm ) final;

  void createWithGlobalSize( globalIndex const globalRows,
                             globalIndex const globalCols,
                             localIndex const maxEntriesPerRow,
                             MPI_Comm const & comm ) final;

  void open() final;

  void close() final;

  bool isCreated() const final;

  void set( real64 const value ) final;

  void zero() final;

  void add( globalIndex const rowIndex,
            globalIndex const colIndex,
            real64 const value ) final;

  void set( globalIndex const rowIndex,
            globalIndex const colIndex,
            real64 const value ) final;

  void insert( globalIndex const rowIndex,
               globalIndex const colIndex,
               real64 const value ) final;

  void add( globalIndex const rowIndex,
            globalIndex const * colIndices,
            real64 const * values,
            localIndex const size ) final;

  void set( globalIndex const rowIndex,
            globalIndex const * colIndices,
            real64 const * values,
            localIndex const size ) final;

  void insert( globalIndex const rowIndex,
               globalIndex const * colIndices,
               real64 const * values,
               localIndex const size ) final;

  void add( globalIndex const rowIndex,
            arraySlice1d< globalIndex const > const & colIndices,
            arraySlice1d< real64 const > const & values ) final;

  void set( globalIndex const rowIndex,
            arraySlice1d< globalIndex const > const & colIndices,
            arraySlice1d< real64 const > const & values ) final;

  void insert( globalIndex const rowIndex,
               arraySlice1d< globalIndex const > const & colIndices,
               arraySlice1d< real64 const > const & values ) final;

  void add( arraySlice1d< globalIndex const > const & rowIndices,
            arraySlice1d< globalIndex const > const & colIndices,
            arraySlice2d< real64 const, 1 > const & values ) final;

  void set( arraySlice1d< globalIndex const > const & rowIndices,
            arraySlice1d< globalIndex const > const & colIndices,
            arraySlice2d< real64 const, 1 > const & values ) final;

  void insert( arraySlice1d< globalIndex const > const & rowIndices,
               arraySlice1d< globalIndex const > const & colIndices,
               arraySlice2d< real64 const, 1 > const & values ) final;

  void add( arraySlice1d< globalIndex const > const & rowIndices,
            arraySlice1d< globalIndex const > const & colIndices,
            arraySlice2d< real64 const, 0 > const & values ) final;

  void set( arraySlice1d< globalIndex const > const & rowIndices,
            arraySlice1d< globalIndex const > const & colIndices,
            arraySlice2d< real64 const, 0 > const & values ) final;

  void insert( arraySlice1d< globalIndex const > const & rowIndices,
               arraySlice1d< globalIndex const > const & colIndices,
               arraySlice2d< real64 const, 0 > const & values ) final;

  void add( globalIndex const * rowIndices,
            globalIndex const * colIndices,
            real64 const * values,
            localIndex const numRows,
            localIndex const numCols ) final;

  void set( globalIndex const * rowIndices,
            globalIndex const * colIndices,
            real64 const * values,
            localIndex const numRows,
            localIndex const numCols ) final;

  void insert( globalIndex const * rowIndices,
               globalIndex const * colIndices,
               real64 const * values,
               localIndex const numRows,
               localIndex const numCols ) final;

  void multiply( EpetraVector const & src,
                 EpetraVector & dst ) const final;

  void multiply( EpetraMatrix const & src,
                 EpetraMatrix & dst,
                 bool const closeResult = true ) const final;

  void leftMultiplyTranspose( EpetraMatrix const & src,
                              EpetraMatrix & dst,
                              bool const closeResult = true ) const final;

  void rightMultiplyTranspose( EpetraMatrix const & src,
                               EpetraMatrix & dst,
                               bool const closeResult = true ) const final;

  void gemv( real64 const alpha,
             EpetraVector const & x,
             real64 const beta,
             EpetraVector & y,
             bool useTranspose = false ) const final;

  void scale( real64 const scalingFactor ) final;

  void leftScale( EpetraVector const & vec ) final;

  void rightScale( EpetraVector const & vec ) final;

  void leftRightScale( EpetraVector const & vecLeft,
                       EpetraVector const & vecRight ) final;

  void clearRow( globalIndex const row,
                 real64 const diagValue = 0 ) final;

  void getRowCopy( globalIndex globalRow,
                   array1d< globalIndex > & colIndices,
                   array1d< real64 > & values ) const final;

  real64 getDiagValue( globalIndex globalRow ) const final;

  globalIndex globalRows() const final;

  globalIndex globalCols() const final;

  localIndex localRows() const final;

  localIndex localCols() const final;

  globalIndex ilower() const final;

  globalIndex iupper() const final;

  localIndex localNonzeros() const final;

  globalIndex globalNonzeros() const final;

  real64 normInf() const final;

  real64 norm1() const final;

  real64 normFrobenius() const final;

  localIndex getLocalRowID( globalIndex const index ) const final;

  globalIndex getGlobalRowID( localIndex const index ) const final;

  MPI_Comm getComm() const final;

  void print( std::ostream & os = std::cout ) const final;

  void write( string const & filename,
              bool const mtxFormat = true ) const final;

  ///@}

  /**
   * @brief Returns a pointer to the underlying matrix.
   */
  Epetra_FECrsMatrix * unwrappedPointer() const;

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

/**
 * @brief Stream insertion operator for EpetraMatrix
 * @param os the output stream
 * @param matrix the matrix to be printed
 * @return reference to the output stream
 */
std::ostream & operator<<( std::ostream & os,
                           EpetraMatrix const & matrix );

} // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_EPETRAMATRIX_HPP_*/
