/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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
 *
 *  Created on: Jul 24, 2018
 *  Author: Matthias Cremon
 *
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
  void createWithLocalSize( trilinosTypes::lid const localSize,
                            trilinosTypes::lid const maxEntriesPerRow = 1,
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
  void createWithGlobalSize( trilinosTypes::gid const globalSize,
                             trilinosTypes::lid const maxEntriesPerRow = 1,
                             MPI_Comm const & comm = MPI_COMM_WORLD);

  /**
   * @brief Create a rectangular matrix from number of rows/columns.
   *
   * \param comm MPI communicator.
   * \param localRows Local number of rows.
   * \param localCols Local number of columns.
   * \param maxEntriesPerRow Maximum number of entries per row (hint).
   */
  void createWithLocalSize( trilinosTypes::lid const localRows,
                            trilinosTypes::lid const localCols,
                            trilinosTypes::lid const maxEntriesPerRow = 1,
                            MPI_Comm const & comm = MPI_COMM_WORLD );

  /**
   * @brief Create a rectangular matrix from number of rows/columns.
   *
   * \param comm MPI communicator.
   * \param globalRows Global number of rows.
   * \param globalCols Global number of columns.
   * \param maxEntriesPerRow Maximum number of entries per row (hint).
   */
  void createWithGlobalSize( trilinosTypes::gid const globalRows,
                             trilinosTypes::gid const globalCols,
                             trilinosTypes::lid const maxEntriesPerRow = 1,
                             MPI_Comm const & comm = MPI_COMM_WORLD );


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
  void add( trilinosTypes::gid const rowIndex,
            trilinosTypes::gid const colIndex,
            real64 const value );

  /**
   * @brief Set one element.
   *
   * \param rowIndex Global row index.
   * \param colIndex Global column index.
   * \param value Value to set at prescribed location.
   *
   */
  void set( trilinosTypes::gid const rowIndex,
            trilinosTypes::gid const colIndex,
            real64 const value );

  /**
   * @brief Insert one element.
   *
   * \param rowIndex Global row index.
   * \param colIndex Global column index.
   * \param value Value to insert at prescribed location.
   *
   */
  void insert( trilinosTypes::gid const rowIndex,
               trilinosTypes::gid const colIndex,
               real64 const value );

  /**
   * @brief Add elements to one row using c-style arrays
   *
   * \param rowIndex Global row index.
   * \param colIndices Global column indices
   * \param values Values to add to prescribed locations.
   * \param size Number of elements 
   */
  void add( trilinosTypes::gid const rowIndex,
            trilinosTypes::gid const * colIndices,
            real64 const * values,
            trilinosTypes::lid const size);

  /**
   * @brief Set elements to one row using c-style arrays
   *
   * \param rowIndex Global row index.
   * \param colIndices Global column indices
   * \param values Values to add to prescribed locations.
   * \param size Number of elements 
   */
  void set( trilinosTypes::gid const rowIndex,
            trilinosTypes::gid const * colIndices,
            real64 const * values,
            trilinosTypes::lid const size);

  /**
   * @brief Insert elements to one row using c-style arrays
   *
   * \param rowIndex Global row index.
   * \param colIndices Global column indices
   * \param values Values to add to prescribed locations.
   * \param size Number of elements 
   */
  void insert( trilinosTypes::gid const rowIndex,
               trilinosTypes::gid const * colIndices,
               real64 const * values,
               trilinosTypes::lid const size);

  /**
   * @brief Add elements to one row using array1d
   *
   * \param rowIndex Global row index.
   * \param colIndices Global column indices
   * \param values Values to add to prescribed locations.
   */
  void add( trilinosTypes::gid const rowIndex,
            array1d<trilinosTypes::gid> const & colIndices,
            array1d<real64> const & values);

  /**
   * @brief Set elements of one row using array1d
   *
   * \param rowIndex Global row index.
   * \param colIndices Global column indices
   * \param values Values to add to prescribed locations.
   */
  void set( trilinosTypes::gid const rowIndex,
            array1d<trilinosTypes::gid> const & colIndices,
            array1d<real64> const & values);

  /**
   * @brief Insert elements of one row using array1d
   *
   * \param rowIndex Global row index.
   * \param colIndices Global column indices
   * \param values Values to add to prescribed locations.
   */
  void insert( trilinosTypes::gid const rowIndex,
               array1d<trilinosTypes::gid> const & colIndices,
               array1d<real64> const & values);

  /**
   * @brief Add dense matrix.
   *
   * \param rowIndices Global row indices.
   * \param colIndices Global col indices
   * \param values Dense local matrix of values.
   */
  void add( array1d<trilinosTypes::gid> const & rowIndices,
            array1d<trilinosTypes::gid> const & colIndices,
            array2d<real64> const & values);

  /**
   * @brief Set dense matrix.
   *
   * \param rowIndices Global row indices.
   * \param colIndices Global col indices
   * \param values Dense local matrix of values.
   */
  void set( array1d<trilinosTypes::gid> const & rowIndices,
            array1d<trilinosTypes::gid> const & colIndices,
            array2d<real64> const & values);

  /**
   * @brief Insert dense matrix.
   *
   * \param rowIndices Global row indices.
   * \param colIndices Global col indices
   * \param values Dense local matrix of values.
   */
  void insert( array1d<trilinosTypes::gid> const & rowIndices,
               array1d<trilinosTypes::gid> const & colIndices,
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
  void multiply( EpetraVector const &src,
                 EpetraVector &dst ) const;

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
   * @brief Clear a row and multiplies the diagonal term by <tt>factor</tt>.
   *
   * \param row Index of the row to be cleared.
   * \param factor Scaling factor for diagonal element.
   *
   */
  void clearRow( trilinosTypes::gid const row,
                 real64 const factor );

  //@}

  //! @name Accessors Methods
  //@{

  // TODO: replace with array1d versions

  /**
   * @brief Returns the row <tt>GlobalRow</tt>. The number of non zeros in the row is <tt>NumEntries</tt>
   * , the values are sent to <tt>vecValues</tt> and the column indices in <tt>vecIndices</tt>.
   */
  void getRow( trilinosTypes::gid GlobalRow,
               trilinosTypes::lid &NumEntries,
               real64* Values,
               trilinosTypes::gid* Indices ) const;

  /**
   * @brief Returns the row <tt>GlobalRow</tt>. The number of non zeros in the row is <tt>NumEntries</tt>
   * , the values are sent to <tt>vecValues</tt> and the column indices in <tt>vecIndices</tt>.
   */
  void getRow( trilinosTypes::gid GlobalRow,
               trilinosTypes::lid &NumEntries,
               std::vector<real64> &vecValues,
               std::vector<trilinosTypes::gid> &vecIndices ) const;

  /**
   * @brief Returns the row <tt>localRow</tt>. The number of non zeros in the row is <tt>NumEntries</tt>
   * , the values are sent to <tt>vecValues</tt> and the column indices in <tt>vecIndices</tt>.
   */
  void getLocalRow( trilinosTypes::lid myRow,
                    trilinosTypes::lid & NumEntries,
                    real64 * & Values,
                    trilinosTypes::lid * & Indices ) const;

  /**
   * @brief Returns the row <tt>localRow</tt>. The number of non zeros in the row is <tt>NumEntries</tt>
   * , the values are sent to <tt>vecValues</tt> and the column indices in <tt>vecIndices</tt>.
   */
  void getLocalRow( trilinosTypes::lid myRow,
                    trilinosTypes::lid &NumEntries,
                    std::vector<real64> &vecValues,
                    std::vector<trilinosTypes::lid> &vecIndices ) const;

  /**
   * @brief Returns a pointer to the underlying matrix.
   */
  Epetra_FECrsMatrix* unwrappedPointer() const;

  /**
   * @brief Returns the number of global rows.
   */
  trilinosTypes::gid globalRows() const;

  /**
   * @brief Returns the number of global columns.
   */
  trilinosTypes::gid globalCols() const;

  /**
   * @brief Returns the number of local rows.
   */
  trilinosTypes::lid localRows() const;

  /**
   * @brief Returns the number of local columns.
   */
  trilinosTypes::lid localCols() const;

  /**
   * @brief Returns the index of the first global row owned by that processor.
   */
  trilinosTypes::gid ilower() const;

  /**
   * @brief Returns the index of the last global row owned by that processor.
   */
  trilinosTypes::gid iupper() const;

  /**
   * @brief Wrapper for LRID function. Returns the local index of the corresponding global index.
   * Returns -1 if the row is not owned by the processor.
   */
  trilinosTypes::lid rowLID( trilinosTypes::gid const GID ) const;

  /**
   * @brief Wrapper for GRID64 function. Returns the global index of the corresponding local index.
   * Returns -1 if the row is not owned by the processor.
   */
  trilinosTypes::gid rowGID( trilinosTypes::lid const LID ) const;

  //TODO: Do we ever need LCID and GCID64 functions as well?

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
