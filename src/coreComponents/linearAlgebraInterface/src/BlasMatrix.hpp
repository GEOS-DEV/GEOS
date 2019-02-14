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
 * @file BlasMatrix.hpp
 */

#ifndef CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLASMATRIX_HPP_
#define CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLASMATRIX_HPP_

#include "common/DataTypes.hpp"
#include "Logger.hpp"

#include "cblas.h"
#include "lapacke.h"

namespace geosx
{
/**
 * \class BlasMatrix
 * \brief This class creates and provides basic support for for manipulating
 *        a BLAS-style column-major matrix.
 */

class BlasMatrix
{
public:

  //----------------------------------------------------------------------------
  //! @name Constructor/Destructor Methods
  //@{

  /**
   * @brief Default constructor.
   */
  BlasMatrix();

  /**
   * @brief Shaped constructor; defines a variable-sized matrix.
   *
   * \param IN
   * nRows - Number of matrix rows.
   * \param IN
   * nCols - Number of matrix cols.
   *
   * @note All values are initialized to 0.
   */
  BlasMatrix( localIndex numRows,
              localIndex numCols );

  /**
   * @brief Shape constructor; defines a square matrix.
   *
   * \param IN
   * order - Matrix order
   *
   * @note All values are initialized to 0.
   */
  BlasMatrix( localIndex order );

  /**
   * @brief Copy constructor.
   *
   * \param IN
   * src - BlasMatrix
   *
   */
  BlasMatrix( BlasMatrix const & src );

  /**
   * @brief Matrix destructor.
   */
  ~BlasMatrix();

  //@}

  //----------------------------------------------------------------------------
  //! @name Shaping/sizing methods
  //@{

  /**
   * @brief Resize matrix. All entries set to zero
   *
   * \param IN
   * nRows - New number of matrix rows.
   * \param IN
   * nCols - New number of matrix cols.
   *
   */
  void resize( localIndex numRows,
               localIndex numCols );

  /**
   * @brief Resize matrix. All entries set to zero
   *
   * \param IN
   * order - Matrix order
   *
   */
  void resize( localIndex order );

  //@}

  //----------------------------------------------------------------------------
  //! @name Mathematical methods
  //@{

  /**
   * @brief Computes determinant.
   *
   * The function is implemented for a matrix up to order three. The matrix must
   * be square.
   */
  real64 determinant();

  /**
   * @brief Compute inverse; \a this = <tt>M</tt><sup>-1</sup>.
   *
   * Assign the inverse of the given matrix to \a this.
   *
   * \param IN
   * <tt>M</tt> - Dense matrix.
   *
   * @warning
   * Assumes \a this already has the same size as <tt>M</tt>.
   *
   * @note This function is hardcoded for square matrices up to order four.
   * For dimensions larger than four, the function calls LAPACK functions DGETRF
   * and DGETRI using Trilinos/Epetra LAPACK Wrapper Class.
   */
  void computeInverse( BlasMatrix& dst );

  /**
   * @brief Compute inverse; \a this = <tt>M</tt><sup>-1</sup>.
   *
   * Assign the inverse of the given matrix to \a dst and return the
   * determinant of \a this.
   *
   * \param IN
   * \a dst - Dense matrix containg the inverse.
   *
   * \param IN
   * \a determinant - Determinant of \a this.
   *
   * @warning
   * Assumes \a dst already has the same size as \a this.
   *
   * @note This function is hardcoded for square matrices up to order four.
   * For dimensions larger than four, the function calls LAPACK functions DGETRF
   * and DGETRI using lapacke interface.
   */
  void computeInverse( BlasMatrix& dst, real64& det );

  /**
   * @brief Matrix-Matrix sum;
   * \a this = scalarA*<tt>A</tt> + \a this.
   *
   * Computes (scalarA*<tt>A</tt> + \a this) and overwrites the result on the
   * current \a this, with optional scaling.
   *
   * \param IN
   * <tt>A</tt> - BlasMatrix matrix.
   * \param [IN]
   * scalarA - Optional scalar to multiply with <tt>A</tt>.
   *
   * @warning
   * Assumes that <tt>A</tt> and \a this have the same size.
   */
  void MatAdd( BlasMatrix const & A,
               const real64 scalarA = 1. );

  /**
   * @brief Matrix-Matrix product;
   * \a dst = scalarDst * \a dst + scalarThisSrc * \a this * \a src.
   *
   * Computes matrix-matrix product with optional scaling and accumulation.
   *
   * \param IN
   * \a src - Source BlasMatrix.
   * \param IN
   * \a dst - Destination BlasMatrix.
   * \param [IN]
   * scalarThisSrc - Optional scalar to multiply with \a this * \a src.
   * \param [IN]
   * scalarDst - Optional parameter to control the accumulation.
   *
   * @warning
   * Assumes that \a this and \a src have compatible sizes and that \a dst
   * already has the right size.
   *
   */
  void matrixMultiply(BlasMatrix const &src,
                      BlasMatrix &dst,
                      real64 const scalarThisSrc=1.,
                      real64 const scalarDst=0.);

  /**
   * @brief transpose(Matrix)-Matrix product;
   * \a dst = scalarDst * \a dst + scalarThisSrc * \a this<sup>T</sup> * \a src.
   *
   * Computes transpose(matrix)-matrix product with optional scaling and accumulation.
   *
   * \param IN
   * \a src - Source BlasMatrix.
   * \param IN
   * \a dst - Destination BlasMatrix.
   * \param [IN]
   * scalarThisSrc - Optional scalar to multiply with \a this<sup>T</sup> * \a src.
   * \param [IN]
   * scalarDst - Optional parameter to control the accumulation.
   *
   * @warning
   * Assumes that \a this and \a src have compatible sizes and that \a dst
   * already has the right size.
   *
   */
  void TmatrixMultiply(BlasMatrix const &src,
                       BlasMatrix &dst,
                       real64 const scalarThisSrc=1.,
                       real64 const scalarDst=0.);

  /**
   * @brief Matrix-transpose(Matrix) product;
   * \a dst = scalarDst * \a dst + scalarThisSrc * \a this * \a src<sup>T</sup>.
   *
   * Computes matrix-transpose(matrix) product with optional scaling and accumulation.
   *
   * \param IN
   * \a src - Source BlasMatrix.
   * \param IN
   * \a dst - Destination BlasMatrix.
   * \param [IN]
   * scalarThisSrc - Optional scalar to multiply with \a this * \a src<sup>T</sup>.
   * \param [IN]
   * scalarDst - Optional parameter to control the accumulation.
   *
   * @warning
   * Assumes that \a this and \a src have compatible sizes and that \a dst
   * already has the right size.
   *
   */
  void matrixTMultiply(BlasMatrix const &src,
                       BlasMatrix &dst,
                       real64 const scalarThisSrc=1.,
                       real64 const scalarDst=0.);

  /**
   * @brief transpose(Matrix)-transpose(Matrix) product;
   * \a dst = scalarDst * \a dst +
   *          scalarThisSrc * \a this<sup>T</sup> * \a src<sup>T</sup>.
   *
   * Computes transpose(matrix)matrix-transpose(matrix) product with optional
   * scaling and accumulation.
   *
   * \param IN
   * \a src - Source BlasMatrix.
   * \param IN
   * \a dst - Destination BlasMatrix.
   * \param [IN]
   * scalarThisSrc - Optional scalar to multiply with \a this<sup>T</sup> * \a src<sup>T</sup>.
   * \param [IN]
   * scalarDst - Optional parameter to control the accumulation.
   *
   * @warning
   * Assumes that \a this and \a src have compatible sizes and that \a dst
   * already has the right size.
   *
   */
  void TmatrixTMultiply(BlasMatrix const &src,
                        BlasMatrix &dst,
                        real64 const scalarThisSrc=1.,
                        real64 const scalarDst=0.);

//
//  /**
//   * @brief Matrix-Vector product;
//   * \a <tt>y</tt> = \a this * <tt>x</tt>.
//   *
//   * Computes matrix-vector product with optional scaling and accumulation.
//   *
//   * \param IN
//   * <tt>x</tt> - Dense vector.
//   * \param OUT
//   * <tt>y</tt> - Dense vector.
//   *
//   * @warning
//   * Assumes that <tt>x</tt> and <tt>y</tt> have compatible sizes.
//   */
//  void MatVecMult(SerialDenseMatrix& x, SerialDenseMatrix& y);
//
//  /**
//   * @brief Matrix-Vector product;
//   * \a <tt>y</tt> = \a this <sup>T</sup> * <tt>x</tt>.
//   *
//   * Computes matrix-vector product with optional scaling and accumulation using
//   * the transpose of \a this.
//   *
//   * \param IN
//   * <tt>x</tt> - Dense vector.
//   * \param OUT
//   * <tt>y</tt> - Dense vector.
//   *
//   * @warning
//   * Assumes that <tt>x</tt> and <tt>y</tt> have compatible sizes.
//   */
//  void MatTVecMult(SerialDenseMatrix& x, SerialDenseMatrix& y);
//
//
//  /**
//   * @brief In-place scalar-matrix product;
//   * \a this = scalarThis* \a this
//   *
//   * \param IN
//   * scalarThis - Scalar to multiply with \a this.
//   *
//   * @note This function calls BLAS function DSCAL using Trilinos/Epetra BLAS
//   *       Wrapper Class if built with Trilinos.
//   */
//  void Scale(double scalarThis);
//
//  //@}
//
  //----------------------------------------------------------------------------
  //! @name Data Accessor methods
  //@{

  /**
   * @brief Coefficient access function.
   *
   * The parentheses operator returns the reference to the coefficient in
   * position (row, column).
   */
  inline real64 & operator ()( localIndex row,
                               localIndex column );

  /**
   * @brief Constant coefficient access function.
   *
   * The parentheses operator returns the reference to the coefficient in
   * position (row, column).
   */
  inline real64 const & operator ()( localIndex iRow,
                                     localIndex jCol ) const;

  /**
   * @brief Returns number of matrix rows.
   */
  localIndex getNumRows() const
  {
    return m_numRows;
  }
  ;

  /**
   * @brief Returns number of matrix columns.
   */
  localIndex getNumCols() const
  {
    return m_numCols;
  }
  ;

//  /**
//   * @brief Assign scalar \a value to all matrix entries.
//   *
//   * \param IN
//   * <tt>value</tt> - Scalar to be assigned to all matrix entries.
//   */
//  void assign_value(double value);

  //@}

  //----------------------------------------------------------------------------
  //! @name I/O methods
  //@{

  /**
   * @brief Print service method; defines behavior of ostream << operator.
   */
  void print();

  //@}

protected:

  localIndex m_numRows = 0; ///< Number of rows of the matrix.
  localIndex m_numCols = 0; ///< Number of columns of the matrix.
  array1d<real64> m_values; ///< array1d storing matrix entries (column-major)
};

// Inline methods
inline real64 & BlasMatrix::operator ()( localIndex iRow,
                                         localIndex jCol )
{
  GEOS_ASSERT_MSG( 0 <= iRow &&
                       iRow <= m_numRows &&
                       0 <= jCol &&
                       jCol <= m_numCols,
                   "Requested value out of bounds" );
  return m_values[jCol * m_numRows + iRow];
}

inline real64 const & BlasMatrix::operator ()( localIndex iRow,
                                               localIndex jCol ) const
                                               {
  GEOS_ASSERT_MSG( 0 <= iRow &&
                       iRow <= m_numRows &&
                       0 <= jCol &&
                       jCol <= m_numCols,
                   "Requested value out of bounds" );
  return m_values[jCol * m_numRows + iRow];
}

} // namespace geosx

#endif /* CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLASMATRIX_HPP_ */
