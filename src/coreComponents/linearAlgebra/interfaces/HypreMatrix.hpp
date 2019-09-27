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
 * @file HypreMatrix.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_HYPREMATRIX_HPP_
#define GEOSX_LINEARALGEBRA_HYPREMATRIX_HPP_

#include "common/DataTypes.hpp"
#include "HypreVector.hpp"

#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_mv.h"


namespace geosx
{

/**
 * \class HypreMatrix
 * \brief This class ...
 */
class HypreMatrix
{
public:

  //! @name Constructor/Destructor Methods
  //@{

  /**
   * @brief Empty matrix constructor.
   *
   * Create an empty (distributed) matrix.
   */

  HypreMatrix();

  /**
   * @brief Copy constructor.
   *
   * Create new matrix from matrix <tt>src</tt>.
   */
  HypreMatrix( HypreMatrix const &src );

  /**
   * @brief Virtual destructor.
   */
  virtual ~HypreMatrix();

  //@}

  //! @name Create Methods
  //@{

  /**
   * @brief Create a matrix from an existing HYPRE_IJMatrix.
   *
   * TODO not implemented yet when the generation of the sparsity pattern will be decided.
   *
   * \param ... .
   */
  void create( );

  /**
   * @brief Create a square matrix from local number of rows.
   *
   * \param localSize local number of rows for square matrix.
   * \param maxEntriesPerRow Maximum number of non-zero entries per row.
   * \param comm MPI communicator.
   *
   */
  void createWithLocalSize( localIndex const localSize,
                            localIndex const maxEntriesPerRow,
                            MPI_Comm const & comm = MPI_COMM_WORLD);

  /**
   * @brief Create a square matrix from global number of rows.
   *
   * Create a square matrix with an (approximately) even partitioning of rows.
   *
   * \param globalSize Global dimensions for a square matrix.
   * \param maxEntriesPerRow Maximum number of non-zero entries per row.
   * \param comm MPI communicator.
   *
   */
  void createWithGlobalSize( globalIndex const globalSize,
                             localIndex const maxEntriesPerRow,
                             MPI_Comm const & comm = MPI_COMM_WORLD);

  /**
   * @brief Create a rectangular matrix from number of rows/columns.
   *
   * \param comm MPI communicator.
   * \param localRows Local number of rows.
   * \param localCols Local number of columns.
   * \param maxEntriesPerRow Maximum number of entries per row (hint).
   */
  void createWithLocalSize( localIndex const localRows,
                            localIndex const localCols,
                            localIndex const maxEntriesPerRow,
                            MPI_Comm const & comm = MPI_COMM_WORLD );

  /**
   * @brief Create a rectangular matrix from number of rows/columns.
   *
   * \param comm MPI communicator.
   * \param globalRows Global number of rows.
   * \param globalCols Global number of columns.
   * \param maxEntriesPerRow Maximum number of entries per row (hint).
   */
  void createWithGlobalSize( globalIndex const globalRows,
                             globalIndex const globalCols,
                             localIndex const maxEntriesPerRow,
                             MPI_Comm const & comm = MPI_COMM_WORLD );

  /**
   * @brief Reinitialize the matrix.
   *
   * Keeps the parallel partitioning and the sparsity pattern but sets all elements to user-defined value.
   *
   */
  void set( real64 const value );


  /**
   * @brief Reinitialize the matrix.
   *
   * Keeps the parallel partitioning and the sparsity pattern but sets all elements to zero.
   *
   */
  void zero();

  /**
   * @brief Empty function for Trilinos implementation. Is required when the HYPRE library is used.
   *
   */
  void open();

  /**
   * @brief Assemble and compress the matrix.
   *
   * Compresses the matrix to CSR format with contiguous memory on each processor. Prevents from
   * adding new entries in the sparsity pattern but allows for modification of existing entries.
   *
   */
  void close();

  //@}

  /** @name Add/Set/Insert Methods
   *
   * TRILINOS logic:
   * The add and set methods assume entries already exist in the sparsity pattern.
   * Insert methods allow for dynamic allocation, but will temporarily use
   * extra memory if one attempts to insert multiple values to the same location.
   *
   * HYPRE logic:
   * The add and set methods can be used also if the sparsity pattern has not been
   * finalized. In Hypre the insert method is an alias for set
   *
   * Caution: In Trilinos these methods are not thread-safe.  //TODO: add thread safety
   */
  //@{

  /**
   * @brief Add to one element.
   *
   * \param rowIndex Global row index.
   * \param colIndex Global column index.
   * \param value Value to add to prescribed location.
   *
   */
  void add( globalIndex const rowIndex,
            globalIndex const colIndex,
            real64 const value );

  /**
   * @brief Set one element.
   *
   * \param rowIndex Global row index.
   * \param colIndex Global column index.
   * \param value Value to set at prescribed location.
   *
   */
  void set( globalIndex const rowIndex,
            globalIndex const colIndex,
            real64 const value );

  /**
   * @brief Insert one element.
   *
   * \param rowIndex Global row index.
   * \param colIndex Global column index.
   * \param value Value to insert at prescribed location.
   *
   */
  void insert( globalIndex const rowIndex,
               globalIndex const colIndex,
               real64 const value );

  /**
   * @brief Add elements to one row using c-style arrays
   *
   * \param rowIndex Global row index.
   * \param colIndices Global column indices
   * \param values Values to add to prescribed locations.
   * \param size Number of elements
   */
  void add( globalIndex const rowIndex,
            globalIndex const * colIndices,
            real64 const * values,
            localIndex const size);

  /**
   * @brief Set elements to one row using c-style arrays
   *
   * \param rowIndex Global row index.
   * \param colIndices Global column indices
   * \param values Values to add to prescribed locations.
   * \param size Number of elements
   */
  void set( globalIndex const rowIndex,
            globalIndex const * colIndices,
            real64 const * values,
            localIndex const size);

  /**
   * @brief Insert elements to one row using c-style arrays
   *
   * \param rowIndex Global row index.
   * \param colIndices Global column indices
   * \param values Values to add to prescribed locations.
   * \param size Number of elements
   */
  void insert( globalIndex const rowIndex,
               globalIndex const * colIndices,
               real64 const * values,
               localIndex const size);

  /**
   * @brief Add elements to one row using array1d
   *
   * \param rowIndex Global row index.
   * \param colIndices Global column indices
   * \param values Values to add to prescribed locations.
   */
  void add( globalIndex const rowIndex,
            array1d<globalIndex> const & colIndices,
            array1d<real64> const & values);

  /**
   * @brief Set elements of one row using array1d
   *
   * \param rowIndex Global row index.
   * \param colIndices Global column indices
   * \param values Values to add to prescribed locations.
   */
  void set( globalIndex const rowIndex,
            array1d<globalIndex> const & colIndices,
            array1d<real64> const & values);

  /**
   * @brief Insert elements of one row using array1d
   *
   * \param rowIndex Global row index.
   * \param colIndices Global column indices
   * \param values Values to add to prescribed locations.
   */
  void insert( globalIndex const rowIndex,
               array1d<globalIndex> const & colIndices,
               array1d<real64> const & values);

  /**
   * @brief Add dense matrix.
   *
   * \param rowIndices Global row indices.
   * \param colIndices Global col indices
   * \param values Dense local matrix of values.
   */
  void add( array1d<globalIndex> const & rowIndices,
            array1d<globalIndex> const & colIndices,
            array2d<real64> const & values);

  /**
   * @brief Set dense matrix.
   *
   * \param rowIndices Global row indices.
   * \param colIndices Global col indices
   * \param values Dense local matrix of values.
   */
  void set( array1d<globalIndex> const & rowIndices,
            array1d<globalIndex> const & colIndices,
            array2d<real64> const & values);

  /**
   * @brief Insert dense matrix.
   *
   * \param rowIndices Global row indices.
   * \param colIndices Global col indices
   * \param values Dense local matrix of values.
   */
  void insert( array1d<globalIndex> const & rowIndices,
               array1d<globalIndex> const & colIndices,
               array2d<real64> const & values);

  //@}

  //! @name Linear Algebra Methods
  //@{
  /**
   * @brief Matrix/Vector multiplication.
   *
   * Compute <tt>Ax = b<tt>.
   *
   * \param src Input vector (x).
   * \param dst Output vector (b).
   *
   */
  void multiply( HypreVector const &src,
                 HypreVector &dst ) const;

  /**
   * @brief Matrix/Matrix multiplication.
   *
   * Compute <tt>this * B = C<tt>.
   *
   * \param src Input matrix (B).
   * \param dst Output matrix (C).
   * \param closeResult whether to close @p dst for additional entries.
   *
   * Note that the output matrix C should have the same
   * row-map as this.  If close() has already been called
   * on C, then C's sparsity pattern must already contain
   * the nonzero entries produced by the product this*B.
   */
  void multiply( HypreMatrix const & src,
		         HypreMatrix & dst,
                 bool const closeResult = true ) const;

  /**
   * @brief Compute residual <tt>r = b - A*x</tt>.
   *
   * \param x Input solution.
   * \param b Input right hand side.
   * \param r Output residual.
   *
   */
  void residual( HypreVector const &x,
                 HypreVector const &b,
                 HypreVector &r ) const;

  /**
   * @brief Compute gemv <tt>y = alpha*A*x + beta*y</tt>.
   *
   * @note The naming convention follows the BLAS library.
   *
   * \param alpha Scalar factor for added matvec product.
   * \param x Input vector.
   * \param beta Scalar factor for right hand side.
   * \param y Output vector.
   * \param useTranspose Boolean, set to true to use <tt>A^T</tt>.
   *
   */
  void gemv( real64 const alpha,
             HypreVector const &x,
             real64 const beta,
             HypreVector &y,
             bool useTranspose=false );

  /**
   * @brief Multiply all elements by scalingFactor.
   *
   * \param scalingFactor Scaling factor.
   *
   */
  void scale( real64 const scalingFactor );

  /**
   * @brief Pre-multiplies (left) with diagonal matrix consisting of the values in vec.
   *
   * \param vec Vector to pre-multiply with.
   *
   */
  void leftScale( HypreVector const &vec );

//  /**
//   * @brief Post-multiplies (right) with diagonal matrix consisting of the values in vec.
//   *
//   * \param vec Vector to post-multiply with.
//   *
//   */
//  void rightScale( EpetraVector const &vec );
//
//  /**
//   * @brief Post-multiplies (right) with diagonal matrix consisting of the values in vecRight
//   * and pre-multiplies (left) with diagonal matrix consisting of the values in vec.
//   *
//   * \param vec vecLeft to pre-multiply with.
//   * \param vec vecRight to post-multiply with.
//   *
//   */
//  void leftRightScale( EpetraVector const &vecLeft,
//                       EpetraVector const &vecRight );
//
//  /**
//   * @brief Clear a row, and optionally set diagonal element to <tt>diagValue</tt>.
//   *
//   * \param row globalIndex of the row to be cleared.
//   * \param diagValue (Optional) set diagonal element to desired value.
//   *
//   */
//  void clearRow( globalIndex const row,
//                 real64 const diagValue = 0 );
//
  //@}

  //! @name Accessors Methods
  //@{

  /**
   * @brief Returns a copy of the data in row <tt>globalRow</tt>.
   * Note that the input arrays will be resized internally to fit the number of entries.
   */
  void getRowCopy( globalIndex globalRow,
                   array1d<globalIndex> & colIndices,
                   array1d<real64> & values) const;

  /**
   * @brief Returns a pointer to the underlying HYPRE_IJMatrix object.
   */
  HYPRE_IJMatrix const * unwrappedPointer() const;

  HYPRE_IJMatrix * unwrappedPointer();

  operator HYPRE_IJMatrix()
  {
    return (HYPRE_IJMatrix) m_ij_mat;
  }

  operator HYPRE_ParCSRMatrix()
  {
    return (HYPRE_ParCSRMatrix) m_parcsr_mat;
  }

  /**
   * @brief Returns the number of global rows.
   */
  globalIndex globalRows() const;

  /**
   * @brief Returns the number of global columns.
   */
  globalIndex globalCols() const;

  /**
   * @brief Returns the index of the first global row owned by that processor.
   */
  globalIndex ilower() const;

  /**
   * @brief Returns the index of the last global row owned by that processor.
   */
  globalIndex iupper() const;

//  /**
//   * @brief Returns the infinity norm of the matrix.
//   */
//  real64 normInf() const;
//
//  /**
//   * @brief Returns the one norm of the matrix.
//   */
//  real64 norm1() const;

  /**
   * @brief Returns the Frobenius norm of the matrix.
   */
  real64 normFrobenius() const;

//  /**
//   * @brief Returns true is the matrix has been assembled, false if not.
//   */
//  bool isAssembled() const;
  //@}

//  //! @name I/O Methods
//  //@{
//  /**
//   * @brief Print the matrix in Trilinos format to the terminal.
//   */
//  void print() const;

  /**
   * @brief Write the matrix to filename in a HYPRE format.
   */
  void write(string const & filename) const;

  //@}

private:

  /**
   * Boolean value, true if the matrix sparsity pattern has been fixed.
   */
  bool m_is_pattern_fixed = false;

  /**
   * Boolean value, true if the matrix had been finalized, false if not.
   */
  bool m_is_ready_to_use = false;

  /**
   * Pointer to underlying HYPRE_IJMatrix type.
   */
  HYPRE_IJMatrix m_ij_mat = nullptr;

  /**
   * Pointer to underlying HYPRE_ParCSRMatrix type.
   */
  HYPRE_ParCSRMatrix m_parcsr_mat = nullptr;

//  /**
//   * Pointer to the underlying Epetra_CrsMatrix.
//   */
//  std::unique_ptr<Epetra_FECrsMatrix> m_matrix = nullptr;
//
//  /*
//   * Map representing the parallel partitioning of a source vector (x in y=Ax)
//   */
//  std::unique_ptr<Epetra_Map> m_src_map = nullptr;
//
//  /*
//   * Map representing the parallel partitioning of a destination vector (y in y=Ax)
//   */
//  std::unique_ptr<Epetra_Map> m_dst_map = nullptr;
};

} // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_HYPREMATRIX_HPP_*/
