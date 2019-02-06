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
 * @file LapackMatrix.hpp
 */

#ifndef CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_LAPACKMATRIX_HPP_
#define CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_LAPACKMATRIX_HPP_

#include "common/DataTypes.hpp"
#include "Logger.hpp"

#include "cblas.h"
#include "lapacke.h"

namespace geosx
{
/**
 * \class LapackMatrix
 * \brief This class creates and provides basic support for for real-valued,
 *        double-precision dense rectangular. Column-major storage is used.
 */
class LapackMatrix
{
public:

  //----------------------------------------------------------------------------
  //! @name Constructor/Destructor Methods
  //@{

  /**
   * @brief Empty constructor.
   */
  LapackMatrix();

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
  LapackMatrix( localIndex nRows,
                localIndex nCols );

  /**
   * @brief Matrix destructor.
   */
  virtual ~LapackMatrix()
  {

  }

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
  void resize( localIndex nRows,
               localIndex nCols );
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
  double determinant();

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
  void invert(LapackMatrix& src);

//  /**
//   * @brief Matrix-Matrix sum;
//   * \a this = scalarA*<tt>A</tt> + \a this.
//   *
//   * Computes (scalarA*<tt>A</tt> + \a this) and overwrites the result on the
//   * current \a this, with optional scaling.
//   *
//   * \param IN
//   * <tt>A</tt> - Dense matrix.
//   * \param [IN]
//   * scalarA - Optional scalar to multiply with <tt>A</tt>.
//   *
//   * @warning
//   * Assumes that <tt>A</tt> and \a this have the same size.
//   *
//   * @note This function calls BLAS function DAXPY using Trilinos/Epetra BLAS
//   *       Wrapper Class if built with Trilinos.
//   */
//  void MatAdd(SerialDenseMatrix& A, const double scalarA=1.);
//
//  /**
//   * @brief Matrix-Matrix product;
//   * \a this = scalarThis* \a this + scalarAB*<tt>A</tt>*<tt>B</tt>.
//   *
//   * Computes matrix-matrix product with optional scaling and accumulation.
//   *
//   * \param IN
//   * <tt>A</tt> - Dense matrix.
//   * \param IN
//   * <tt>B</tt> - Dense matrix.
//   * \param [IN]
//   * scalarAB - Optional scalar to multiply with <tt>A</tt>*<tt>B</tt>.
//   * \param [IN]
//   * scalarThis - Optional parameter to control the accumulation.
//   *
//   * @warning
//   * Assumes that <tt>A</tt> and <tt>B</tt> have compatible sizes and that
//   * \a this already has the right size.
//   *
//   * @note This function calls BLAS function DGEMM using Trilinos/Epetra BLAS
//   *       Wrapper Class if built with Trilinos.
//   */
//  void MatMatMult(SerialDenseMatrix& A, SerialDenseMatrix& B,
//                  const double scalarAB=1., const double scalarThis=0.);
//
//  /**
//   * @brief Matrix-Matrix product;
//   * \a this = scalarThis* \a this +
//   *           scalarAB*<tt>A</tt>*<tt>B</tt><sup>T</sup>.
//   *
//   * Computes matrix-matrix product with optional scaling and accumulation
//   * using the transpose of <tt>B</tt>.
//   *
//   * \param IN
//   * <tt>A</tt> - Dense matrix.
//   * \param IN
//   * <tt>B</tt> - Dense matrix.
//   * \param [IN]
//   * scalarAB - Optional scalar to multiply with
//   * <tt>A</tt>*<tt>B</tt><sup>T</sup>.
//   * \param [IN]
//   * scalarThis - Optional parameter to control the accumulation.
//   *
//   * @warning
//   * Assumes that <tt>A</tt> and <tt>B</tt><sup>T</sup> have compatible sizes
//   * and that \a this already has the right size.
//   *
//   * @note This function calls BLAS function DGEMM using Trilinos/Epetra BLAS
//   *       Wrapper Class if built with Trilinos.
//   */
//  void MatMatTMult(SerialDenseMatrix& A, SerialDenseMatrix& B,
//                   const double scalarAB=1., const double scalarThis=0.);
//
//  /**
//   * @brief Matrix-Matrix product;
//   * \a this = scalarThis* \a this +
//   *           scalarAB*<tt>A</tt><sup>T</sup>*<tt>B</tt>.
//   *
//   * Computes matrix-matrix product with optional scaling and accumulation
//   * using the transpose of <tt>A</tt>.
//   *
//   * \param IN
//   * <tt>A</tt> - Dense matrix.
//   * \param IN
//   * <tt>B</tt> - Dense matrix.
//   * \param [IN]
//   * scalarAB - Optional scalar to multiply with
//   * <tt>A</tt><sup>T</sup>*<tt>B</tt>.
//   * \param [IN]
//   * scalarThis - Optional parameter to control the accumulation.
//   *
//   * @warning
//   * Assumes that <tt>A</tt><sup>T</sup> and <tt>B</tt> have compatible sizes
//   * and that \a this already has the right size.
//   *
//   * @note This function calls BLAS function DGEMM using Trilinos/Epetra BLAS
//   *       Wrapper Class if built with Trilinos.
//   */
//  void MatTMatMult(SerialDenseMatrix& A, SerialDenseMatrix& B,
//                   const double scalarAB=1., const double scalarThis=0.);
//
//  /**
//   * @brief Matrix-Matrix product;
//   * \a this = scalarThis* \a this +
//   *           scalarAB*<tt>A</tt><sup>T</sup>*<tt>B</tt><sup>T</sup>.
//   *
//   * Computes matrix-matrix product with optional scaling and accumulation
//   * using the transpose of <tt>A</tt> and <tt>B</tt>.
//   *
//   * \param IN
//   * <tt>A</tt> - Dense matrix.
//   * \param IN
//   * <tt>B</tt> - Dense matrix.
//   * \param [IN]
//   * scalarAB - Optional scalar to multiply with
//   * A</tt><sup>T</sup>*<tt>B</tt><sup>T</sup>.
//   * \param [IN]
//   * scalarThis - Optional parameter to control the accumulation.
//   *
//   * @warning
//   * Assumes that <tt>A</tt><sup>T</sup> and <tt>B</tt><sup>T</sup> have
//   * compatible sizes and that \a this already has the right size.
//   *
//   * @note This function calls BLAS function DGEMM using Trilinos/Epetra BLAS
//   *       Wrapper Class if built with Trilinos.
//   */
//  void MatTMatTMult(SerialDenseMatrix& A, SerialDenseMatrix& B,
//                    const double scalarAB=1., const double scalarThis=0.);
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
   * @brief Element access function.
   *
   * The parentheses operator returns the element in the ith row (RowIndex) and
   * jth column (ColIndex).
   */
  real64& operator()(localIndex RowIndex, localIndex ColIndex);

  /**
   * @brief Returns number of matrix rows.
   */
  localIndex get_nRows() const
  {
    return m_nRows;
  };

  /**
   * @brief Returns number of matrix columns.
   */
  localIndex get_coefficientPtr() const
  {
    return m_nCols;
  };

  /**
   * @brief Returns number of matrix columns.
   */
  localIndex get_nCols() const
  {
    return m_nCols;
  };

//  /**
//   * @brief Assign scalar \a value to all matrix entries.
//   *
//   * \param IN
//   * <tt>value</tt> - Scalar to be assigned to all matrix entries.
//   */
//  void assign_value(double value);

  //@}
//
//  //----------------------------------------------------------------------------
//  //! @name I/O methods
//  //@{
//
//  /**
//    * @brief Print service method; defines behavior of ostream << operator.
//    */
//  void print();
//
//  //@}

protected:

  localIndex m_nRows = 0;
  localIndex m_nCols = 0;
  array1d<real64> m_Elements;
};

// inlined definitions of op()
inline real64& LapackMatrix::operator () (localIndex RowIndex,
                                          localIndex ColIndex)
{
  return m_Elements[ColIndex*m_nRows + RowIndex];
}

}// namespace geosx

#endif /* CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_LAPACKMATRIX_HPP_ */
