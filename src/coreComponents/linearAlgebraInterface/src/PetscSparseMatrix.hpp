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
 *
 *  Created on: Jan 28, 2019
 *  Author: Hannah Morgan
 *
 */

#ifndef PETSCSPARSEMATRIX_HPP_
#define PETSCSPARSEMATRIX_HPP_

#include "InterfaceTypes.hpp"
#include "common/DataTypes.hpp"
#include <petscvec.h>
#include <petscmat.h>
#include "PetscVector.hpp"

namespace geosx
{

/**
 * \class PetscSparseMatrix
 * \brief This class creates and provides basic support for the Mat
 *        matrix object type used in Petsc.
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
				       int const m_nRowGlobal,
				       int const nMaxEntriesPerRow );

  /**
   * @brief Create a rectangular matrix from number of rows/columns.
   *
   * \param comm MPI communicator.
   * \param m_nRowGlobal Global number of rows.
   * \param m_nColGlobal Global number of columns.
   * \param nMaxEntriesPerRow Maximum number of entries per row.
   */
  void create( MPI_Comm const comm,
               int const m_nRowGlobal,
               int const m_nColGlobal,
               int const nMaxEntriesPerRow);

   /**
   * @brief Create a square matrix from number of rows.
   *
   * \param comm MPI communicator.
   * \param m_nRowGlobal Maximum number of entries per row.
   * \param nMaxEntriesPerRow vector of maximum entries per row.
   *
   */
   void create( MPI_Comm const comm,
                int const m_nRowGlobal,
                std::vector<int> const nMaxEntriesPerRow );
  
   /**
   * @brief Create a rectangular matrix from number of rows/columns.
   *
   * \param comm MPI communicator.
   * \param m_nRowGlobal Global number of rows.
   * \param m_nColGlobal Global number of columns.
   * \param nMaxEntriesPerRow vector of maximum entries per row.
   */
   void create( MPI_Comm const comm,
                int const m_nRowGlobal,
                int const m_nColGlobal,
                std::vector<int> const nMaxEntriesPerRow );

  /**
   * @brief Create a matrix from an existing PetscSparseMatrix.
   *
   * \param PetscSparseMatrix existing matrix.
   */
  void create( PetscSparseMatrix &matrix );

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
  void add( int const iRow,
            int const nCols,
            real64 const *values,
            int const *cols );

  /* add values to sparse matrix located at indices */
  // void add( array1d<int> const rowIndices,
  //           array1d<int> const colIndices,
  //           array2d<real64> const values);

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
  void add( int const iRow,
            int const iCol,
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
  void set( int const iRow,
            int const nCols,
            real64 const *values,
            int const *cols );

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
  void set( int const iRow,
            int const iCol,
            real64 const value );

  /**
   * @brief Insert to row of elements.
   *
   * Inserts the values <tt>values</tt> to row <tt>iRow</tt>, at locations specified
   * by <tt>cols</tt>. <tt>nCols</tt> is the number of entries (columns) in the row, and the
   * size of both <tt>values</tt> and <tt>cols</tt>.
   *
   * \param iRow Global row index.
   * \param nCols Number of columns to modify.
   * \param values Values to add to prescribed locations.
   * \param cols Global column indices in which to add the values.
   */
  void insert( int const iRow,
               int const nCols,
               real64 const *values,
               int const *cols );

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
            bool useTranspose);

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
  void clearRow( int const row,
                 real64 const factor );

  //@}

  //! @name Accessors Methods
  //@{

  /**
   * @brief Returns the row <tt>GlobalRow</tt>. The number of non zeros in the row is <tt>NumEntries</tt>
   * , the values are sent to <tt>vecValues</tt> and the column indices in <tt>vecIndices</tt>.
   */
  void getRow( int GlobalRow,
               int &NumEntries,
               const real64* Values,
               const int* Indices ) const;

  /*
   * Get global row myRow
   * - numEntries: number of nonzeros 
   * - vecValues: vector of values
   * - vecIndices: vector of column indices */
  // void getRow( int GlobalRow,
  //              int &NumEntries,
  //              std::vector<real64> &vecValues,
  //              std::vector<int> &vecIndices ) const;

  /**
   * @brief Returns the row <tt>localRow</tt>. The number of non zeros in the row is <tt>NumEntries</tt>
   * , the values are sent to <tt>vecValues</tt> and the column indices in <tt>vecIndices</tt>.
   */
  void getLocalRow( int myRow,
                    int & NumEntries,
                    const real64 * & Values,
                    const int * & Indices ) const;

  /*
   * Get local row myRow
   * - numEntries: number of nonzeros 
   * - vecValues: vector of values
   * - vecIndices: vector of column indices */
  // void getLocalRow( int myRow,
  //                   int &NumEntries,
  //                   std::vector<real64> &vecValues,
  //                   std::vector<int> &vecIndices ) const;

  /**
   * @brief Returns the number of global rows.
   */
  int globalRows() const;

  /**
   * @brief Returns the number of global columns.
   */
  int globalCols() const;

  // /**
  //  * @brief Returns the number of unique columns (can be used to check if matrix is square).
  //  */
  // int uniqueCols() const;

  /**
   * @brief Returns the index of the first global row owned by that processor.
   */
  int ilower() const;

  /**
   * @brief Returns the index of the last global row owned by that processor.
   */
  int iupper() const;

  /**
   * @brief Returns the number of local rows.
   */
  int myRows() const;

  /**
   * @brief Returns the number of local columns.
   */
  int myCols() const;

  /* @brief Return local row number from global row number if owned
     by process, return -1 if the global row is not owned by the processor. */
  int rowMapLID( int const GID ) const;

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
   * @brief Print the matrix in Petsc format to the terminal.
   */
  void print() const;

  /**
   * @brief Returns a pointer to the underlying matrix.
   */
  const Mat* getPointer() const;

  /* Get underlying Petsc object */
  Mat getMat();

protected:

	bool assembled = false;

  // underlying Petsc object
	Mat _mat;

};

} // namespace geosx

#endif /* PetscSPARSEMATRIX_HPP_ */