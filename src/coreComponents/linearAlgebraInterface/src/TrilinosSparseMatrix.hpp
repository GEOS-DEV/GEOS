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
 * @file EpetraSparseMatrix.hpp
 *
 *  Created on: Jul 24, 2018
 *  Author: Matthias Cremon
 *
 */

#ifndef EPETRASPARSEMATRIX_HPP_
#define EPETRASPARSEMATRIX_HPP_

#include "InterfaceTypes.hpp"
#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <EpetraExt_MatrixMatrix.h>
#include "common/DataTypes.hpp"
#include "TrilinosVector.hpp"

namespace geosx
{

/**
 * \class EpetraSparseMatrix
 * \brief This class creates and provides basic support for the Epetra_CrsMatrix
 *        matrix object type used in Trilinos.
 */
class EpetraSparseMatrix
{
public:

  //! @name Constructor/Destructor Methods
  //@{
  /**
   * @brief Empty matrix constructor.
   *
   * Create an empty (distributed) matrix.
   */
  EpetraSparseMatrix();

  /**
   * @brief Copy constructor.
   *
   * Create new matrix from matrix <tt>in_mat</tt>.
   */
  EpetraSparseMatrix( EpetraSparseMatrix const &in_mat );

  /**
   * @brief Virtual destructor.
   */
  virtual ~EpetraSparseMatrix() = default;
  //@}

  //! @name Create/Finalize/Reinitialize Methods
  //@{
  /**
   * @brief Create a square matrix from number of rows.
   *
   * \param comm MPI communicator.
   * \param m_nRowGlobal Maximum number of entries per row.
   * \param nMaxEntriesPerRow Maximum number of entries per row.
   *
   */
  void create( MPI_Comm const comm,
               trilinosTypes::gid const m_nRowGlobal,
               trilinosTypes::lid const nMaxEntriesPerRow );

  /**
   * @brief Create a rectangular matrix from number of rows/columns.
   *
   * \param comm MPI communicator.
   * \param m_nRowGlobal Global number of rows.
   * \param m_nColGlobal Global number of columns.
   * \param nMaxEntriesPerRow Maximum number of entries per row.
   */
  void create( MPI_Comm const comm,
               trilinosTypes::gid const m_nRowGlobal,
               trilinosTypes::gid const m_nColGlobal,
               trilinosTypes::lid const nMaxEntriesPerRow = 0 );

  /**
   * @brief Create a square matrix from Epetra_Map.
   *
   * Prepare the matrix to be filled. Takes as inputs an existing Epetra_Map and
   * a hint on the maximum number of entries per row. This method is meant to be called
   * after the matrix has been declared from the empty constructor.
   *
   * \param input_map Epetra_Map.
   * \param nMaxEntriesPerRow Maximum number of entries per row.
   *
   */
  void create( Epetra_Map const &input_map,
               trilinosTypes::lid const nMaxEntriesPerRow );

  /**
   * @brief Create a rectangular matrix from two existing Epetra_Map, row and column maps.
   *
   * \param row_map Epetra_Map for rows.
   * \param col_map Epetra_Map for columns.
   * \param nMaxEntriesPerRow Maximum number of entries per row.
   */
  void create( Epetra_Map const &row_map,
               Epetra_Map const &col_map,
               trilinosTypes::lid const nMaxEntriesPerRow = 0 );

  /**
   * @brief Create a matrix from an existing Epetra_CrsGraph.
   *
   * TODO change that to whatever format the sparsity pattern will be.
   *
   * \param Epetra_CrsGraph existing graph.
   */
  void create( Epetra_CrsGraph const &graph );

  /**
   * @brief Create a matrix from an existing Epetra_CrsMatrix.
   *
   * \param Epetra_CrsMatrix existing matrix.
   */
  void create( Epetra_CrsMatrix &matrix );

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

  //! @name Insertion/Replace/SumInto Methods
  //@{
  /**
   * @brief Add to row of elements.
   *
   * Adds the values <tt>values</tt> to row <tt>iRow</tt>, at locations specified
   * by <tt>cols</tt>. <tt>nCols</tt> is the number of entries (columns) in the row, and the
   * size of both <tt>values</tt> and <tt>cols</tt>.
   *
   * \param iRow Global row index.
   * \param nCols Number of columns to modify.
   * \param values Values to add to prescribed locations.
   * \param cols Global column indices in which to add the values.
   */
  void add( trilinosTypes::gid const iRow,
            trilinosTypes::lid const nCols,
            real64 const *values,
            trilinosTypes::gid const *cols );

  /**
   * @brief Add dense matrix.
   *
   * Adds the matrix <tt>values</tt> to the sparse matrix, at locations specified
   * by <tt>indices</tt>.
   *
   * \param indices Vector of indices.
   * \param values Values.
   */
  void add( array1d<trilinosTypes::gid> const rowIndices,
            array1d<trilinosTypes::gid> const colIndices,
            array2d<real64> const values);

  /**
   * @brief Add to one element.
   *
   * Adds the value <tt>value</tt> to location (<tt>iRow</tt>,<tt>iCol</tt>).
   *
   * \param iRow Global row index.
   * \param iCol Global column index.
   * \param value Value to add to prescribed locations.
   *
   */
  void add( trilinosTypes::gid const iRow,
            trilinosTypes::gid const iCol,
            real64 const value );

  /**
   * @brief Set row of elements.
   *
   * Sets the values <tt>values</tt> of row <tt>iRow</tt>, at locations specified
   * by <tt>cols</tt>. <tt>nCols</tt> is the number of entries (columns) in the row, and the
   * size of both <tt>values</tt> and <tt>cols</tt>.
   *
   * \param iRow Global row index.
   * \param nCols Number of columns to modify.
   * \param values Values to set in prescribed locations.
   * \param cols Global column indices in which to set the values.
   *
   */
  void set( trilinosTypes::gid const iRow,
            trilinosTypes::lid const nCols,
            real64 const *values,
            trilinosTypes::gid const *cols );

  /**
   * @brief Set one element.
   *
   * Sets the value of the location (<tt>iRow</tt>,<tt>iCol</tt>) to <tt>value</tt>.
   *
   * \param iRow Global row index.
   * \param iCol Global column index.
   * \param value Value to set at prescribed locations.
   *
   */
  void set( trilinosTypes::gid const iRow,
            trilinosTypes::gid const iCol,
            real64 const value );

  /**
   * @brief Insert to row of elements.
   *
   * Inserts the values <tt>values</tt> to row <tt>iRow</tt>, at locations specified
   * by <tt>cols</tt>. <tt>nCols</tt> is the number of entries (columns) in the row, and the
   * size of both <tt>values</tt> and <tt>cols</tt>.
   *
   * TODO remove the possibility to dynamically construct the sparsity pattern.
   *
   * \param iRow Global row index.
   * \param nCols Number of columns to modify.
   * \param values Values to add to prescribed locations.
   * \param cols Global column indices in which to add the values.
   */
  void insert( trilinosTypes::gid const iRow,
               trilinosTypes::lid const nCols,
               real64 const *values,
               trilinosTypes::gid const *cols );

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
  Epetra_CrsMatrix* getPointer() const;

  /**
   * @brief Returns the number of global rows.
   */
  trilinosTypes::gid globalRows() const;

  /**
   * @brief Returns the number of global columns.
   */
  trilinosTypes::gid globalCols() const;

  /**
   * @brief Returns the number of unique columns (can be used to check if matrix is square).
   */
  trilinosTypes::gid uniqueCols() const;

  /**
   * @brief Returns the index of the first global row owned by that processor.
   */
  trilinosTypes::gid ilower() const;

  /**
   * @brief Returns the index of the last global row owned by that processor.
   */
  trilinosTypes::gid iupper() const;

  /**
   * @brief Returns the number of local rows.
   */
  int myRows() const;

  /**
   * @brief Returns the number of local columns.
   */
  int myCols() const;

  /**
   * @brief Returns the row map.
   */
  Epetra_Map const & RowMap() const;

  /**
   * @brief Returns the (usually overlapping) column map.
   */
  Epetra_Map const & ColMap() const;

  /**
   * @brief Returns the (1-to-1) domain column map.
   */
  Epetra_Map const & DomainMap() const;

  /**
   * @brief Wrapper for LID function. Returns the local map of the corresponding global index.
   * Returns -1 if the global row is not owned by the processor.
   */
  trilinosTypes::lid rowMapLID( trilinosTypes::gid const GID ) const;

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
   * @brief Boolean value, true if the matrix had been finalized, false if not.
   */
  bool assembled = false;

  /**
   * @brief Pointer to the underlying Epetra_CrsMatrix.
   */
  std::unique_ptr<Epetra_CrsMatrix> m_matrix = nullptr;

};

} // namespace geosx

#endif /* EpetraSPARSEMATRIX_HPP_ */
