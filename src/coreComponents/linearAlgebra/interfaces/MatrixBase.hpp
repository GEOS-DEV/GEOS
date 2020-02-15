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
 * @file MatrixBase.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_MATRIXBASE_HPP
#define GEOSX_LINEARALGEBRA_INTERFACES_MATRIXBASE_HPP

#include "linearAlgebra/interfaces/LinearOperator.hpp"

namespace geosx
{

enum class MatrixOutputFormat
{
  NATIVE_ASCII,
  NATIVE_BINARY,
  MATLAB_ASCII,
  MATLAB_BINARY,
  MATRIX_MARKET
};

/**
 * @brief Common base template for all matrix wrapper types.
 * @tparam MATRIX derived matrix type
 * @tparam VECTOR compatible vector type
 *
 * This class template provides a common interface for all derived matrix
 * wrapper types. Most methods are pure abstract in order to get the compiler
 * to enforce a common interface in all derived classes; however there is no
 * runtime polymorphism or virtual dispatch - derived classes are not related
 * between each other, and pointers/references to base should never be used -
 * the destructor is protected to enforce that.
 *
 * As an added benefit, the documentation for matrix interface also lives here
 * and does not need to be duplicated across matrix wrappers. Derived classes
 * should still document specific functions in case a particular LA package's
 * behavior deviates from expectations or has unexpected performance impacts.
 * In that case, @c \@copydoc tag can be used to copy over the documentation.
 */
template<typename MATRIX, typename VECTOR>
class MatrixBase : public LinearOperator<VECTOR>
{
public:

  using Matrix = MATRIX;
  using Vector = VECTOR;

  /**
   * @name Create Methods
   */
  ///@{

  /**
   * @brief Create a square matrix from local number of rows.
   *
   * @param localSize local number of rows for square matrix.
   * @param maxEntriesPerRow Maximum number of non-zero entries per row.
   * @param comm MPI communicator.
   *
   */
  virtual void
  createWithLocalSize( localIndex const localSize,
                       localIndex const maxEntriesPerRow,
                       MPI_Comm const & comm ) = 0;

  /**
   * @brief Create a square matrix from global number of rows.
   *
   * Create a square matrix with an (approximately) even partitioning of rows.
   *
   * @param globalSize Global dimensions for a square matrix.
   * @param maxEntriesPerRow Maximum number of non-zero entries per row.
   * @param comm MPI communicator.
   *
   */
  virtual void
  createWithGlobalSize( globalIndex const globalSize,
                        localIndex const maxEntriesPerRow,
                        MPI_Comm const & comm ) = 0;

  /**
   * @brief Create a rectangular matrix from number of rows/columns.
   *
   * @param comm MPI communicator.
   * @param localRows Local number of rows.
   * @param localCols Local number of columns.
   * @param maxEntriesPerRow Maximum number of entries per row (hint).
   */
  virtual void
  createWithLocalSize( localIndex const localRows,
                       localIndex const localCols,
                       localIndex const maxEntriesPerRow,
                       MPI_Comm const & comm ) = 0;

  /**
   * @brief Create a rectangular matrix from number of rows/columns.
   *
   * @param comm MPI communicator.
   * @param globalRows Global number of rows.
   * @param globalCols Global number of columns.
   * @param maxEntriesPerRow Maximum number of entries per row (hint).
   */
  virtual void
  createWithGlobalSize( globalIndex const globalRows,
                        globalIndex const globalCols,
                        localIndex const maxEntriesPerRow,
                        MPI_Comm const & comm ) = 0;

  ///@}

  /**
   * @name Open/close methods
   */
  ///@{

  /**
   * @brief Open matrix for adding new entries.
   *
   * @note Adding entries that result in modifications of sparsity pattern may not be allowed
   *       by most implementations. An error will be raised in that case.
   */
  virtual void open() = 0;

  /**
   * @brief Assemble and compress the matrix.
   *
   * Compresses the matrix to CSR format with contiguous memory on each processor. Prevents from
   * adding new entries in the sparsity pattern but allows for modification of existing entries.
   */
  virtual void close() = 0;

  /**
   * @brief Query matrix open status
   * @return @p true if matrix has been opened and has not been closed since; @p false otherwise
   */
  inline bool isOpen() const { return m_open; }

  /**
   * @brief Query matrix assembled status
   * @return @p true if matrix has been opened and closed since creation; @p false otherwise
   */
  inline bool isAssembled() const { return m_assembled; }

  /**
   * @brief Query matrix creation status
   * @return @p true if matrix has been created
   */
  virtual bool isCreated() const = 0;

  /**
   * @brief Reset the matrix to default state
   */
  virtual void reset()
  {
    m_open = false;
    m_assembled = false;
  }

  ///@}

  /**
   * @name Global modification methods
   */
  ///@{

  /**
   * @brief Set all non-zero elements to a value.
   */
  virtual void set( real64 const value ) = 0;

  /**
   * @brief Set all elements to zero.
   */
  virtual void zero() = 0;

  ///@}

  /**
   * @name Add/Set/Insert Methods
   *
   * @p insert method can only be used while the matrix is being filled for the first time,
   * e.g. during sparsity pattern construction. @p add and @p set methods can only be used
   * after the matrix has been assembled, and both only accept entries that already exist
   * in the sparsity pattern.
   *
   * @note Caution: these methods are not thread-safe.
   */
  ///@{

  /**
   * @brief Add to one element.
   *
   * @param rowIndex Global row index.
   * @param colIndex Global column index.
   * @param value Value to add to prescribed location.
   *
   */
  virtual void add( globalIndex const rowIndex,
                    globalIndex const colIndex,
                    real64 const value ) = 0;

  /**
   * @brief Set one element.
   *
   * @param rowIndex Global row index.
   * @param colIndex Global column index.
   * @param value Value to set at prescribed location.
   *
   */
  virtual void set( globalIndex const rowIndex,
                    globalIndex const colIndex,
                    real64 const value ) = 0;

  /**
   * @brief Insert one element.
   *
   * @param rowIndex Global row index.
   * @param colIndex Global column index.
   * @param value Value to insert at prescribed location.
   *
   */
  virtual void insert( globalIndex const rowIndex,
                       globalIndex const colIndex,
                       real64 const value ) = 0;

  /**
   * @brief Add elements to one row using c-style arrays
   *
   * @param rowIndex Global row index.
   * @param colIndices Global column indices
   * @param values Values to add to prescribed locations.
   * @param size Number of elements
   */
  virtual void add( globalIndex const rowIndex,
                    globalIndex const * colIndices,
                    real64 const * values,
                    localIndex const size ) = 0;

  /**
   * @brief Set elements to one row using c-style arrays
   *
   * @param rowIndex Global row index.
   * @param colIndices Global column indices
   * @param values Values to add to prescribed locations.
   * @param size Number of elements
   */
  virtual void set( globalIndex const rowIndex,
                    globalIndex const * colIndices,
                    real64 const * values,
                    localIndex const size ) = 0;

  /**
   * @brief Insert elements to one row using c-style arrays
   *
   * @param rowIndex Global row index.
   * @param colIndices Global column indices
   * @param values Values to add to prescribed locations.
   * @param size Number of elements
   */
  virtual void insert( globalIndex const rowIndex,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex const size ) = 0;

  /**
   * @brief Add elements to one row using array1d
   *
   * @param rowIndex Global row index.
   * @param colIndices Global column indices
   * @param values Values to add to prescribed locations.
   */
  virtual void add( globalIndex const rowIndex,
                    arraySlice1d<globalIndex const> const & colIndices,
                    arraySlice1d<real64 const> const & values ) = 0;

  /**
   * @brief Set elements of one row using array1d
   *
   * @param rowIndex Global row index.
   * @param colIndices Global column indices
   * @param values Values to add to prescribed locations.
   */
  virtual void set( globalIndex const rowIndex,
                    arraySlice1d<globalIndex const> const & colIndices,
                    arraySlice1d<real64 const> const & values ) = 0;

  /**
   * @brief Insert elements of one row using array1d
   *
   * @param rowIndex Global row index.
   * @param colIndices Global column indices
   * @param values Values to add to prescribed locations.
   */
  virtual void insert( globalIndex const rowIndex,
                       arraySlice1d<globalIndex const> const & colIndices,
                       arraySlice1d<real64 const> const & values ) = 0;

  /**
   * @brief Add dense matrix.
   *
   * @param rowIndices Global row indices.
   * @param colIndices Global col indices
   * @param values Dense local matrix of values.
   *
   * @note Row major layout assumed in values
   */
  virtual void add( arraySlice1d<globalIndex const> const & rowIndices,
                    arraySlice1d<globalIndex const> const & colIndices,
                    arraySlice2d<real64 const, 1> const & values ) = 0;

  /**
   * @brief Set dense matrix.
   *
   * @param rowIndices Global row indices.
   * @param colIndices Global col indices
   * @param values Dense local matrix of values.
   *
   * @note Row major layout assumed in values
   */
  virtual void set( arraySlice1d<globalIndex const> const & rowIndices,
                    arraySlice1d<globalIndex const> const & colIndices,
                    arraySlice2d<real64 const, 1> const & values ) = 0;

  /**
   * @brief Insert dense matrix.
   *
   * @param rowIndices Global row indices.
   * @param colIndices Global col indices
   * @param values Dense local matrix of values.
   *
   * @note Row major layout assumed in values
   */
  virtual void insert( arraySlice1d<globalIndex const> const & rowIndices,
                       arraySlice1d<globalIndex const> const & colIndices,
                       arraySlice2d<real64 const, 1> const & values ) = 0;

  /**
   * @brief Add dense matrix.
   *
   * @param rowIndices Global row indices.
   * @param colIndices Global col indices
   * @param values Dense local matrix of values.
   *
   * @note Column major layout assumed in values
   */
  virtual void add( arraySlice1d<globalIndex const> const & rowIndices,
                    arraySlice1d<globalIndex const> const & colIndices,
                    arraySlice2d<real64 const, 0> const & values ) = 0;

  /**
   * @brief Set dense matrix.
   *
   * @param rowIndices Global row indices.
   * @param colIndices Global col indices
   * @param values Dense local matrix of values.
   *
   * @note Column major layout assumed in values
   */
  virtual void set( arraySlice1d<globalIndex const> const & rowIndices,
                    arraySlice1d<globalIndex const> const & colIndices,
                    arraySlice2d<real64 const, 0> const & values ) = 0;

  /**
   * @brief Insert dense matrix.
   *
   * @param rowIndices Global row indices.
   * @param colIndices Global col indices
   * @param values Dense local matrix of values.
   *
   * @note Column major layout assumed in values
   */
  virtual void insert( arraySlice1d<globalIndex const> const & rowIndices,
                       arraySlice1d<globalIndex const> const & colIndices,
                       arraySlice2d<real64 const, 0> const & values ) = 0;

  /**
   * @brief Add dense matrix.
   *
   * @param rowIndices Global row indices.
   * @param colIndices Global col indices
   * @param values Dense local matrix of values.
   * @param numRows Number of row indices.
   * @param numCols Number of column indices.
   *
   * @note Row major layout assumed in values
   */
  virtual void add( globalIndex const * rowIndices,
                    globalIndex const * colIndices,
                    real64 const * values,
                    localIndex const numRows,
                    localIndex const numCols ) = 0;

  /**
   * @brief Set dense matrix.
   *
   * @param rowIndices Global row indices.
   * @param colIndices Global col indices
   * @param values Dense local matrix of values.
   * @param numRows Number of row indices.
   * @param numCols Number of column indices.
   *
   * @note Row major layout assumed in values
   */
  virtual void set( globalIndex const * rowIndices,
                    globalIndex const * colIndices,
                    real64 const * values,
                    localIndex const numRows,
                    localIndex const numCols ) = 0;

  /**
   * @brief Insert dense matrix.
   *
   * @param rowIndices Global row indices.
   * @param colIndices Global col indices
   * @param values Dense local matrix of values.
   * @param numRows Number of row indices.
   * @param numCols Number of column indices.
   *
   * @note Row major layout assumed in values
   */
  virtual void insert( globalIndex const * rowIndices,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex const numRows,
                       localIndex const numCols ) = 0;

  ///@}

  /**
   * @name Linear Algebra Methods
   */
  ///@{

  /**
   * @brief Apply operator to a vector
   * @param src Input vector (x).
   * @param dst Output vector (b).
   */
  virtual void
  multiply( Vector const & src,
            Vector & dst ) const override = 0;

  /**
   * @brief Matrix/Matrix multiplication.
   *
   * Compute <tt>this * B = C<tt>.
   *
   * @param src Input matrix (B).
   * @param dst Output matrix (C).
   * @param closeResult whether to close @p dst for additional entries.
   *
   * Note that the output matrix C should have the same
   * row-map as this.  If close() has already been called
   * on C, then C's sparsity pattern must already contain
   * the nonzero entries produced by the product this*B.
   */
  virtual void multiply( Matrix const & src,
                         Matrix & dst,
                         bool const closeResult = true ) const = 0;

  /**
   * @brief Matrix/Matrix transpose multiplication.
   *
   * Compute <tt>this^T * B = C<tt>.
   *
   * @param src Input matrix (B).
   * @param dst Output matrix (C).
   * @param closeResult whether to close @p dst for additional entries.
   *
   * Note that the output matrix C should have the same
   * row-map as this.  If close() has already been called
   * on C, then C's sparsity pattern must already contain
   * the nonzero entries produced by the product this*B.
   */
  virtual void leftMultiplyTranspose( Matrix const & src,
                                      Matrix & dst,
                                      bool const closeResult = true ) const = 0;

  /**
   * @brief Matrix/Matrix transpose multiplication.
   *
   * Compute <tt>B * this^T = C<tt>.
   *
   * @param src Input matrix (B).
   * @param dst Output matrix (C).
   * @param closeResult whether to close @p dst for additional entries.
   *
   * Note that the output matrix C should have the same
   * row-map as this.  If close() has already been called
   * on C, then C's sparsity pattern must already contain
   * the nonzero entries produced by the product this*B.
   */
  virtual void rightMultiplyTranspose( Matrix const & src,
                                       Matrix & dst,
                                       bool const closeResult = true ) const = 0;

  /**
   * @brief Compute gemv <tt>y = alpha*A*x + beta*y</tt>.
   *
   * @note The naming convention follows the BLAS library.
   *
   * @param alpha Scalar factor for added matvec product.
   * @param x Input vector.
   * @param beta Scalar factor for right hand side.
   * @param y Output vector.
   * @param useTranspose Boolean, set to true to use <tt>A^T</tt>.
   *
   */
  virtual void gemv( real64 const alpha,
                     Vector const & x,
                     real64 const beta,
                     Vector & y,
                     bool useTranspose = false ) const = 0;

  /**
   * @brief Multiply all elements by scalingFactor.
   * @param scalingFactor Scaling factor.
   */
  virtual void scale( real64 const scalingFactor ) = 0;

  /**
   * @brief Pre-multiplies (left) with diagonal matrix consisting of the values in vec.
   * @param vec Vector to pre-multiply with.
   */
  virtual void leftScale( Vector const & vec ) = 0;

  /**
   * @brief Post-multiplies (right) with diagonal matrix consisting of the values in vec.
   * @param vec Vector to post-multiply with.
   */
  virtual void rightScale( Vector const & vec ) = 0;

  /**
   * @brief Post-multiplies (right) with diagonal matrix consisting of the values in vecRight
   * and pre-multiplies (left) with diagonal matrix consisting of the values in vec.
   * @param vec vecLeft to pre-multiply with.
   * @param vec vecRight to post-multiply with.
   *
   */
  virtual void leftRightScale( Vector const & vecLeft,
                               Vector const & vecRight ) = 0;

  /**
   * @brief Clear a row, and optionally set diagonal element to <tt>diagValue</tt>.
   * @param row globalIndex of the row to be cleared.
   * @param diagValue (Optional) set diagonal element to desired value.
   *
   */
  virtual void clearRow( globalIndex const row,
                         real64 const diagValue = 0.0 ) = 0;

  ///@}

  /**
   * @name Accessors Methods
   */
  ///@{

  /**
   * @brief Returns the number of nozero entries in the longest
   * row of the matrix.
   */
  virtual localIndex maxRowLength() const = 0;

  /**
   * @brief Get row length via local row index.
   * @return the number of nonzero entries in the row
   *
   * TODO: Breaks the goal of hiding local row indexing from user.
   *       Revise use cases to use ilower() and iupper().
   */
  virtual localIndex localRowLength( localIndex localRowIndex ) const = 0;

  /**
   * @brief Get row length via global row index.
   * @return the number of nonzero entries in the row
   */
  virtual localIndex globalRowLength( globalIndex globalRowIndex ) const = 0;

  /**
   * @brief Returns a copy of the data in row <tt>globalRow</tt>.
   * Note that the input arrays will be resized internally to fit the number of entries.
   */
  virtual void getRowCopy( globalIndex globalRow,
                           array1d< globalIndex > & colIndices,
                           array1d< real64 > & values ) const = 0;

  /**
   * @brief get diagonal element value on a given row
   * @param globalRow global row index
   * @return value of diagonal element on the row
   */
  virtual real64 getDiagValue( globalIndex globalRow ) const = 0;

  /**
   * @brief Returns the number of global rows.
   */
  virtual globalIndex globalRows() const = 0;

  /**
   * @brief Returns the number of global columns.
   */
  virtual globalIndex globalCols() const = 0;

  /**
   * @brief Return the local number of columns on each processor
   */
  virtual localIndex localRows() const = 0;

  /**
   * @brief Return the local number of columns on each processor
   */
  virtual localIndex localCols() const = 0;

  /**
   * @brief Returns the index of the first global row owned by that processor.
   */
  virtual globalIndex ilower() const = 0;

  /**
   * @brief Returns the next index after last global row owned by that processor.
   *
   * @note The intention is for [ilower; iupper) to be used as a half-open index range
   */
  virtual globalIndex iupper() const = 0;

  /**
   * @brief Returns the number of nonzeros in the local portion of the matrix
   */
  virtual localIndex localNonzeros() const = 0;

  /**
   * @brief Returns the total number of nonzeros in the matrix
   */
  virtual globalIndex globalNonzeros() const = 0;

  /**
   * @brief Returns the infinity norm of the matrix.
   */
  virtual real64 normInf() const = 0;

  /**
   * @brief Returns the one norm of the matrix.
   */
  virtual real64 norm1() const = 0;

  /**
   * @brief Returns the Frobenius norm of the matrix.
   */
  virtual real64 normFrobenius() const = 0;

  /**
   * @brief Map a global row index to local row index
   */
  virtual localIndex getLocalRowID( globalIndex const index ) const = 0;

  /**
   * @brief Map a local row index to global row index
   */
  virtual globalIndex getGlobalRowID( localIndex const index ) const = 0;

  /**
   * @brief Get the MPI communicator the matrix was created with
   * @return MPI communicator passed in @p create...()
   *
   * @note when build without MPI, may return anything
   *       (MPI_Comm will be a mock type defined in MpiWrapper)
   */
  virtual MPI_Comm getComm() const = 0;

  ///@}

  /**
   * @name I/O Methods
   */
  ///@{

  /**
   * @brief Print the matrix in Trilinos format to the terminal.
   */
  virtual void print( std::ostream & os = std::cout ) const = 0;

  /**
   * @brief Write the matrix to filename in a matlab-compatible format.
   * @param mtxFormat if @p true, MatrixMarket format is used, otherwise MATLAB
   *
   * Within octave / matlab:
   * >> load filename
   * >> M = spconvert(filename_root)
   */
  virtual void write( string const & filename,
                      MatrixOutputFormat const format ) const = 0;

  ///@}

protected:

  MatrixBase()
  : m_open(false),
    m_assembled(false)
  {}

  ~MatrixBase() = default;

  /// Flag indicating whether the matrix is currently open for adding new entries
  bool m_open;

  /// Flag indicating whether the matrix (sparsity pattern) has been assembled
  bool m_assembled;

};

}

#endif //GEOSX_LINEARALGEBRA_INTERFACES_MATRIXBASE_HPP
