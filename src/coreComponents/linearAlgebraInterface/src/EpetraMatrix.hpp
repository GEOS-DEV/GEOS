/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file EpetraMatrix.hpp
 */

#ifndef LAI_EPETRAMATRIX_HPP_
#define LAI_EPETRAMATRIX_HPP_

#include "InterfaceTypes.hpp"

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_FECrsGraph.h>
#include <Epetra_FECrsMatrix.h>
#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_RowMatrixOut.h>

#include "common/DataTypes.hpp"
#include "EpetraVector.hpp"

namespace geosx
{

/**
 * \class EpetraMatrix
 * \brief This class creates and provides basic support for the Epetra_CrsMatrix
 *        matrix object type used in Trilinos.
 */
class EpetraMatrix
{
public:

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
  EpetraMatrix( EpetraMatrix const &src );

  /**
   * @brief Virtual destructor.
   */
  virtual ~EpetraMatrix() = default;

  //@}
  //! @name Create Methods
  //@{

  /**
   * @brief Create a matrix from an existing Epetra_CrsGraph.
   *
   * TODO change this to whatever format the sparsity pattern will be.
   *
   * \param Epetra_FECrsGraph existing graph.
   */
  void create( Epetra_FECrsGraph const &graph );

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
                            MPI_Comm const & comm = MPI_COMM_WORLD );

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
                             MPI_Comm const & comm );

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
                            MPI_Comm const & comm );

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
                             MPI_Comm const & comm );

  /**
   * @brief Create a rectangular matrix from number of local rows/global columns.
   * \param comm MPI communicator.
   * \param localRows Local number of rows.
   * \param globalCols Global number of columns.
   * \param maxEntriesPerRow Maximum number of entries per row (hint).
   */
  void createWithLocalRowGlobalCol( localIndex const localRows,
                                    globalIndex const globalCols,
                                    localIndex const maxEntriesPerRow,
                                    MPI_Comm const & comm );

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
   * The add and set methods assume entries already exist in the sparsity pattern.
   * Insert methods allow for dynamic allocation, but will temporarily use
   * extra memory if one attempts to insert multiple values to the same location.
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
            localIndex const size );

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
            localIndex const size );

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
               localIndex const size );

  /**
   * @brief Add elements to one row using array1d
   *
   * \param rowIndex Global row index.
   * \param colIndices Global column indices
   * \param values Values to add to prescribed locations.
   */
  void add( globalIndex const rowIndex,
            array1d<globalIndex> const & colIndices,
            array1d<real64> const & values );

  /**
   * @brief Set elements of one row using array1d
   *
   * \param rowIndex Global row index.
   * \param colIndices Global column indices
   * \param values Values to add to prescribed locations.
   */
  void set( globalIndex const rowIndex,
            array1d<globalIndex> const & colIndices,
            array1d<real64> const & values );

  /**
   * @brief Insert elements of one row using array1d
   *
   * \param rowIndex Global row index.
   * \param colIndices Global column indices
   * \param values Values to add to prescribed locations.
   */
  void insert( globalIndex const rowIndex,
               array1d<globalIndex> const & colIndices,
               array1d<real64> const & values );

  /**
   * @brief Add dense matrix.
   *
   * \param rowIndices Global row indices.
   * \param colIndices Global col indices
   * \param values Dense local matrix of values.
   */
  void add( array1d<globalIndex> const & rowIndices,
            array1d<globalIndex> const & colIndices,
            array2d<real64> const & values );

  /**
   * @brief Set dense matrix.
   *
   * \param rowIndices Global row indices.
   * \param colIndices Global col indices
   * \param values Dense local matrix of values.
   */
  void set( array1d<globalIndex> const & rowIndices,
            array1d<globalIndex> const & colIndices,
            array2d<real64> const & values );

  /**
   * @brief Insert dense matrix.
   *
   * \param rowIndices Global row indices.
   * \param colIndices Global col indices
   * \param values Dense local matrix of values.
   */
  void insert( array1d<globalIndex> const & rowIndices,
               array1d<globalIndex> const & colIndices,
               array2d<real64> const & values );

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
  void multiply( EpetraVector const &src,
                 EpetraVector &dst ) const;


  /**
   * @brief Matrix/Matrix multiplication.
   *
   * Compute <tt>this*B = C<tt>.
   *
   * \param src Input matrix (B).
   * \param dst Output matrix (C).
   *
   * Note that the output matrix C should have the same
   * row-map as this.  If close() has already been called
   * on C, then C's sparsity pattern must already contain
   * the nonzero entries produced by the product this*B.
   */
  void multiply( EpetraMatrix const &src,
                 EpetraMatrix &dst ) const;


  /**
   * @brief Compute residual <tt>r = Ax - b</tt>.
   *
   * \param x Input solution.
   * \param b Input right hand side.
   * \param r Output residual.
   *
   */
  void residual( EpetraVector const &x,
                 EpetraVector const &b,
                 EpetraVector &r ) const;

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
             EpetraVector const &x,
             real64 const beta,
             EpetraVector &y,
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
  void leftScale( EpetraVector const &vec );

  /**
   * @brief Post-multiplies (right) with diagonal matrix consisting of the values in vec.
   *
   * \param vec Vector to post-multiply with.
   *
   */
  void rightScale( EpetraVector const &vec );

  /**
   * @brief Post-multiplies (right) with diagonal matrix consisting of the values in vecRight
   * and pre-multiplies (left) with diagonal matrix consisting of the values in vec.
   *
   * \param vec vecLeft to pre-multiply with.
   * \param vec vecRight to post-multiply with.
   *
   */
  void leftRightScale( EpetraVector const &vecLeft,
                       EpetraVector const &vecRight );

  /**
   * @brief Clear a row, and optionally set diagonal element to <tt>diagValue</tt>.
   *
   * \param row globalIndex of the row to be cleared.
   * \param diagValue (Optional) set diagonal element to desired value.
   *
   */
  void clearRow( globalIndex const row,
                 real64 const diagValue = 0 );

  //@}

  //! @name Accessors Methods
  //@{

  /**
   * @brief Returns a copy of the data in row <tt>globalRow</tt>.
   * Note that the input arrays will be resized internally to fit the number of entries.
   */
  void getRowCopy( globalIndex globalRow,
                   array1d<globalIndex> & colIndices,
                   array1d<real64> & values ) const;

  /**
   * @brief Returns a pointer to the underlying matrix.
   */
  Epetra_FECrsMatrix* unwrappedPointer() const;

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

  /**
   * @brief Returns the infinity norm of the matrix.
   */
  real64 normInf() const;

  /**
   * @brief Returns the one norm of the matrix.
   */
  real64 norm1() const;

  /**
   * @brief Returns the Frobenius norm of the matrix.
   */
  real64 normFrobenius() const;

  /**
   * @brief Returns true is the matrix has been assembled, false if not.
   */
  bool isAssembled() const;
  //@}

  //! @name I/O Methods
  //@{
  /**
   * @brief Print the matrix in Trilinos format to the terminal.
   */
  void print() const;

  /**
   * @brief Write the matrix to filename in a matlab-compatible format.
   *
   * Within octave / matlab:
   * >> load filename
   * >> M = spconvert(filename_root)
   */
  void write( string const & filename,
              bool const mtxFormat = true ) const;

  /**
   * @brief Performe a matrix matrix product with Parallel Matrix
   */
  void MatrixMatrixMultiply( bool const transA,
                             EpetraMatrix const &B,
                             bool const transB,
                             EpetraMatrix &C,
                             bool const call_FillComplete = true ) const;

  /**
   * @brief Map a global row index to local row index
   */
  localIndex getLocalRowID( globalIndex const index ) const;

  /**
   * @brief Map a local row index to global row index
   */
  localIndex getGlobalRowID( localIndex const index ) const;

  /**
   * @brief Map a local row index to global row index
   */
  localIndex getGlobalRowID( globalIndex const index ) const;

  /**
   * @brief Return the local number of columns on each processor
   */
  localIndex numMyCols() const;

  /**
   * @brief Print the given parallel matrix in Matrix Market format (MTX file)
   */
  void printParallelMatrix( string const & fileName ) const;

  //@}

private:

  /**
   * Boolean value, true if the matrix had been finalized, false if not.
   */
  bool assembled = false;

  /**
   * Pointer to the underlying Epetra_CrsMatrix.
   */
  std::unique_ptr<Epetra_FECrsMatrix> m_matrix = nullptr;

  /*
   * Map representing the parallel partitioning of a source vector (x in y=Ax)
   */
  std::unique_ptr<Epetra_Map> m_src_map = nullptr;

  /*
   * Map representing the parallel partitioning of a destination vector (y in y=Ax)
   */
  std::unique_ptr<Epetra_Map> m_dst_map = nullptr;
};

} // namespace geosx

#endif /* LAI_EPETRAMATRIX_HPP_ */
