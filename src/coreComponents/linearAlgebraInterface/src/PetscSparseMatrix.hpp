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
 * @file PetscSparseMatrix.hpp
 */

#ifndef LAI_PETSCSPARSEMATRIX_HPP_
#define LAI_PETSCSPARSEMATRIX_HPP_

#include "InterfaceTypes.hpp"

#include <petscvec.h>
#include <petscmat.h>

#include "common/DataTypes.hpp"
#include "PetscVector.hpp"

namespace geosx
{

/**
 * \class PetscSparseMatrix
 * \brief This class creates and provides basic support for the Mat
 *        matrix object type used in PETSc.
 */
class PetscSparseMatrix
{

public:

	//! @name Constructor/Destructor Methods
  //@{

  /**
   * @brief Empty matrix constructor.
   *
   * Create an empty (distributed) matrix.
   */
	PetscSparseMatrix();

	/**
   * @brief Copy constructor.
   *
   * Create new matrix from matrix <tt>in_matrix</tt>.
   */
	PetscSparseMatrix(PetscSparseMatrix const &in_matrix);

	/**
   * @brief Virtual destructor.
   */
	virtual ~PetscSparseMatrix() = default;
  //@}

  //! @name Create Methods
  //@{

  // Hannah: skip?
  //void create( Epetra_FECrsGraph const &graph );

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
   * @brief Reinitialize the matrix.
   *
   * Keeps the parallel partitioning and the sparsity pattern but sets all elements to zero.
   *
   */
  void zero();

  /**
   * @brief Empty function for Petsc implementation. Is required when the HYPRE library is used.
   *
   */
  void open();

  /**
   * @brief Assemble and compress the matrix.
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
   * Hannah: is this thread safe in PETSc?
   */
  //@{

   // change single element
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
   // Hannah: what is the difference between set and insert?

  // change multiple elements with arrays
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

  // change multiple elements with array1d
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

   // Hannah: dense matrix?
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
  void multiply( PetscVector const &src,
                 PetscVector &dst ) const;

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
  void multiply( PetscMatrix const &src,
                 PetscMatrix &dst ) const;

  /**
   * @brief Compute residual <tt>r = Ax - b</tt>.
   *
   * \param x Input solution.
   * \param b Input right hand side.
   * \param r Output residual.
   *
   */
  void residual( PetscVector  const &x,
                 PetscVector  const &b,
                 PetscVector &res ) const;

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
            PetscVector  const &x,
            real64 const beta,
            PetscVector  &y,
            bool useTranspose=false);

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
  void leftScale( PetscVector const &vec );

  /**
   * @brief Post-multiplies (right) with diagonal matrix consisting of the values in vec.
   *
   * \param vec Vector to post-multiply with.
   *
   */
  void rightScale( PetscVector const &vec );

  /**
   * @brief Post-multiplies (right) with diagonal matrix consisting of the values in vecRight
   * and pre-multiplies (left) with diagonal matrix consisting of the values in vec.
   *
   * \param vec vecLeft to pre-multiply with.
   * \param vec vecRight to post-multiply with.
   *
   */
  void leftRightScale( PetscVector const &vecLeft,
                       PetscVector const &vecRight );

  /**
   * @brief Clear a row and multiplies the diagonal term by <tt>factor</tt>.
   *
   * \param row Index of the row to be cleared.
   * \param factor Scaling factor for diagonal element.
   *
   */
  void clearRow( globalIndex const row,
                 real64 const diagValue=0 );

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
   // Hannah: to do

  /**
   * @brief Returns the row <tt>GlobalRow</tt>. The number of non zeros in the row is <tt>NumEntries</tt>
   * , the values are sent to <tt>Values</tt> and the column indices in <tt>Indices</tt>.
   */
  void getRow( globalIndex GlobalRow,
               globalIndex &NumEntries,
               real64* Values,
               globalIndex* Indices ) const;

   /**
   * @brief Returns a pointer to the underlying matrix.
   */
  const Mat* unwrappedPointer() const;

  /**
   * @brief Returns underlying PETSc matrix.
   */
  Mat getMat()

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
   * @brief Print the matrix in PETSc format to the terminal.
   */
  void print() const;

  /**
   * @brief Write the matrix to filename in a matlab-compatible format.
   *
   * Within octave / matlab:
   * >> load filename
   * >> M = spconvert(filename_root)
   */
  void write( string const & filename ) const;

  //@}

// Hannah: protected?
private:

   /**
   * Boolean value, true if the matrix had been finalized, false if not.
   */
	bool assembled = false;

  /**
   * Underlying Petsc object.
   */
	Mat _mat;

};

} // namespace geosx

#endif /* LAI_PETSCSPARSEMATRIX_HPP_ */