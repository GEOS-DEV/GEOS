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
 * @file BlasLapackLA.hpp
 */
#ifndef CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLASLAPACKLA_HPP_
#define CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLASLAPACKLA_HPP_

#include "common/DataTypes.hpp"
#include "Logger.hpp"


namespace geosx
{

class BlasLapackLA
{
  /**
   * \class BlasLapackLA
   * \brief This class contains a collection of BLAS and LAPACK linear
   *        algebra operations (dense) for GEOSX array1d and array2d
   */

public:

  //----------------------------------------------------------------------------
  //! @name Mathematical methods
  //@{

  /**
   * @brief Returns the 1-norm of the vector.
   */
  real64 vectorNorm1( array1d<real64> const & X ) const;

  /**
   * @brief Returns the two norm of the vector.
   */
  real64 vectorNorm2( array1d<real64> const & X ) const;

  /**
   * @brief Infinity-norm of the vector.
   */
  real64 vectorNormInf( array1d<real64> const & X ) const;

  /**
   * @brief Computes matrix determinant.
   *
   * The matrix must be square.
   *
   * @note
   * This function is hardcoded for square matrices up to order four.
   * For dimensions larger than four, the determinant is computed using
   * LAPACK's function DGETRF. To avoid matrix transposition/copy
   * due to the row major ordering used in GEOSX for array2d, the determinant
   * is computed for the transpose matrix, i.e. assuming column major
   * ordering, for best performance.
   */
  real64 determinant( array2d<real64> const & A ) const;

  /**
   * @brief Returns the infinity norm of the matrix.
   *
   * @note
   * Row major ordering is used for GEOSX array2d. Since LAPACK native
   * routines are using a column major ordering (Fortran), the infinity
   * norm is computed as the one norm of the transpose matrix, i.e. assuming
   * column major ordering, for best performance.
   */
  real64 matrixNormInf(array2d<real64> const & A) const;

  /**
   * @brief Returns the one norm of the matrix.
   *
   * @note
   * Row major ordering is used for GEOSX array2d. Since LAPACK native
   * routines are using a column major ordering (Fortran), the one norm
   * is computed as the infinity norm of the transpose matrix, i.e. assuming
   * column major ordering, for best performance.
   */
  real64 matrixNorm1(array2d<real64> const & A) const;

  /**
   * @brief Returns the Frobenius norm of the matrix.
   *
   * @note
   * Row major ordering is used for GEOSX array2d. Since LAPACK native
   * routines are using a column major ordering (Fortran), the one norm
   * is computed for the transpose matrix, i.e. assuming column major
   * ordering, for best performance.
   */
  real64 matrixNormFrobenius(array2d<real64> const & A) const;

  /**
   * @brief Vector-Vector sum;
   * <tt>y</tt> = alpha*<tt>x</tt> + <tt>y</tt>.
   *
   * Computes (alpha*<tt>x</tt> + <tt>y</tt>) and overwrites the result on
   * <tt>y</tt>, with optional scaling.
   *
   * \param IN
   * <tt>x</tt> - GEOSX array1d.
   * \param [IN]
   * alpha - Optional scalar to multiply with <tt>Vec</tt>.
   *
   * \param INOUT
   * <tt>y</tt> - GEOSX array1d.
   *
   * @warning
   * Assumes that <tt>x</tt> and <tt>y</tt> have the same size.
   */
  void vectorVectorAdd( array1d<real64> const & X,
                        array1d<real64> & Y,
                        real64 const alpha = 1. ) const;

  /**
   * @brief Matrix-Matrix sum;
   * <tt>B</tt> = alpha*<tt>A</tt> + <tt>B</tt>.
   *
   * Computes (alpha*<tt>A</tt> + <tt>B</tt>) and overwrites the result on
   * <tt>B</tt>, with optional scaling.
   *
   * \param IN
   * <tt>A</tt> -  GEOSX array2d.
   * \param [IN]
   * alpha - Optional scalar to multiply with <tt>A</tt>.
   *
   * \param INout
   * <tt>B</tt> -  GEOSX array2d.
   *
   * @warning
   * Assumes that <tt>A</tt> and <tt>B</tt> have the same size.
   */
  void matrixMatrixAdd( array2d<real64> const & A,
                        array2d<real64> & B,
                        real64 const alpha = 1. ) const;

  /**
   * @brief In-place scalar-vector product;
   * <tt>x</tt> = alpha*<tt>x<tt>
   *
   * \param IN
   * scalarThis - Scalar to multiply with \a this.
   */
  void vectorScale( array1d<real64> & X,
                    real64 alpha ) const;

  /**
   * @brief In-place scalar-matrix product;
   * <tt>A</tt> = alpha*<tt>A<tt>
   *
   * \param IN
   * alpha - Scalar to multiply with <tt>A</tt>.
   */
  void matrixScale( array2d<real64> & A,
                    real64 alpha) const;

  /**
   * @brief Dot product of two vectors.
   *
   * \param IN
   * <tt>x</tt> - GEOSX array1d.
   * \param IN
   * <tt>y</tt> - GEOSX array1d.
   *
   */
  real64 vectorDot( array1d<real64> const & X,
                    array1d<real64> const & Y) const;

  /**
   * @brief Matrix-Vector product;
   * <tt>Y</tt> = alpha*<tt>A</tt>*<tt>X</tt> + beta*<tt>Y<tt>.
   *
   * Computes matrix-vector product with optional scaling and accumulation.
   *
   * \param IN
   * <tt>X</tt> - GEOSX array1d.
   * \param [IN]
   * alpha - Optional scalar to multiply with <tt>A</tt>*<tt>X</tt>.
   * \param [IN]
   * beta - Optional parameter to control the accumulation.
   *
   * \param INOUT
   * <tt>Y</tt> - GEOSX array1d.
   *
   * @warning
   * Assumes that <tt>X</tt> and <tt>X</tt> have compatible sizes
   * with <tt>A<tt>.
   */
  void matrixVectorMultiply(array2d<real64> const & A,
                            array1d<real64> const & X,
                            array1d<real64>  & Y,
                            real64 const alpha=1.,
                            real64 const beta=0.) const;

  /**
   * @brief transpose(Matrix)-Vector product;
   * <tt>Y</tt> = alpha*<tt>A</tt><sup>T</sup>*<tt>X</tt> + beta*<tt>Y<tt>.
   *
   * Computes transpose(matrix)-vector product with optional scaling and accumulation.
   *
   * \param IN
   * <tt>X</tt> - GEOSX array1d.
   * \param [IN]
   * alpha - Optional scalar to multiply with <tt>A</tt><sup>T</sup>*<tt>X</tt>.
   * \param [IN]
   * beta - Optional parameter to control the accumulation.
   *
   * \param INOUT
   * <tt>Y</tt> - GEOSX array1d.
   *
   * @warning
   * Assumes that <tt>X</tt> and <tt>X</tt> have compatible sizes
   * with <tt>A<tt><sup>T</sup>.
   */
  void matrixTVectorMultiply(array2d<real64> const & A,
                             array1d<real64> const & X,
                             array1d<real64>  & Y,
                             real64 const alpha=1.,
                             real64 const beta=0.) const;

  /**
   * @brief Matrix-Matrix product;
   * * <tt>C</tt> = alpha*<tt>A<tt>*<tt>B<tt> + beta**<tt>C<tt>.
   *
   * Computes matrix-matrix product with optional scaling and accumulation.
   *
   * \param IN
   * <tt>A<tt> - Source GEOSX array2d.
   * \param IN
   * <tt>B<tt> - Source GEOSX array2d.
   * \param [IN]
   * alpha - Optional scalar to multiply with <tt>A<tt>*<tt>B<tt>.
   * \param [IN]
   * beta - Optional parameter to control the accumulation.
   *
   * \param INOUT
   * <tt>C</tt> - Destination GEOSX array2d.
   *
   * @warning
   * Assumes that <tt>A<tt> and <tt>B<tt> have compatible sizes and that
   * <tt>C<tt> already has the right size.
   *
   */
  void matrixMatrixMultiply( array2d<real64> const & A,
                             array2d<real64> const & B,
                             array2d<real64> & C,
                             real64 const alpha=1.,
                             real64 const beta=0.) const;

  /**
   * @brief transpose(Matrix)-Matrix product;
   * * <tt>C</tt> = alpha*<tt>A<tt><sup>T</sup>*<tt>B<tt> + beta**<tt>C<tt>.
   *
   * Computes transpose(matrix)-matrix product with optional scaling and accumulation.
   *
   * \param IN
   * <tt>A<tt> - Source GEOSX array2d.
   * \param IN
   * <tt>B<tt> - Source GEOSX array2d.
   * \param [IN]
   * alpha - Optional scalar to multiply with <tt>A<tt><sup>T</sup>*<tt>B<tt>.
   * \param [IN]
   * beta - Optional parameter to control the accumulation.
   *
   * \param INOUT
   * <tt>C</tt> - Destination GEOSX array2d.
   *
   * @warning
   * Assumes that <tt>A<tt><sup>T</sup> and <tt>B<tt> have compatible sizes and that
   * <tt>C<tt> already has the right size.
   *
   */
  void matrixTMatrixMultiply( array2d<real64> const & A,
                              array2d<real64> const & B,
                              array2d<real64> & C,
                              real64 const alpha=1.,
                              real64 const beta=0.) const;

  /**
   * @brief Matrix-transpose(Matrix) product;
   * * <tt>C</tt> = alpha*<tt>A<tt>*<tt>B<tt><sup>T</sup> + beta**<tt>C<tt>.
   *
   * Computes matrix-transpose(matrix) product with optional scaling and accumulation.
   *
   * \param IN
   * <tt>A<tt> - Source GEOSX array2d.
   * \param IN
   * <tt>B<tt> - Source GEOSX array2d.
   * \param [IN]
   * alpha - Optional scalar to multiply with <tt>A<tt>*<tt>B<tt><sup>T</sup>.
   * \param [IN]
   * beta - Optional parameter to control the accumulation.
   *
   * \param INOUT
   * <tt>C</tt> - Destination GEOSX array2d.
   *
   * @warning
   * Assumes that <tt>A<tt> and <tt>B<tt><sup>T</sup> have compatible sizes and that
   * <tt>C<tt> already has the right size.
   *
   */
  void matrixMatrixTMultiply( array2d<real64> const & A,
                              array2d<real64> const & B,
                              array2d<real64> & C,
                              real64 const alpha=1.,
                              real64 const beta=0.) const;

  /**
   * @brief transpose(Matrix)-transpose(Matrix) product;
   * * <tt>C</tt> = alpha*<tt>A<tt><sup>T</sup>*<tt>B<tt><sup>T</sup>
   *                + beta**<tt>C<tt>.
   *
   * Computes transpose(matrix)-transpose(matrix) product with optional
   * scaling and accumulation.
   *
   * \param IN
   * <tt>A<tt> - Source GEOSX array2d.
   * \param IN
   * <tt>B<tt> - Source GEOSX array2d.
   * \param [IN]
   * alpha - Optional scalar to multiply with <tt>A<tt><sup>T</sup>*<tt>B<tt><sup>T</sup>.
   * \param [IN]
   * beta - Optional parameter to control the accumulation.
   *
   * \param INOUT
   * <tt>C</tt> - Destination GEOSX array2d.
   *
   * @warning
   * Assumes that <tt>A<tt><sup>T</sup> and <tt>B<tt><sup>T</sup> have compatible sizes and that
   * <tt>C<tt> already has the right size.
   *
   */
  void matrixTMatrixTMultiply( array2d<real64> const & A,
                               array2d<real64> const & B,
                               array2d<real64> & C,
                               real64 const alpha=1.,
                               real64 const beta=0.) const;

  /**
   * @brief Compute inverse; <tt>Ainv<tt> = <tt>A</tt><sup>-1</sup>.
   *
   * Assign the inverse of the given square matrix <tt>A<tt> to <tt>Ainv<tt>.
   *
   * \param IN
   * <tt>A</tt> - GEOSX array2d.
   *
   * \param INOUT
   * <tt>Ainv</tt> - GEOSX array2d.
   *
   * @warning
   * Assumes <tt>Ainv<tt> already has the same size as <tt>A</tt>.
   *
   * @note This function is hardcoded for matrices up to order three.
   * For dimensions larger than three, the function calls LAPACK
   * functions DGETRF and DGETRI. Because of the row major
   * ordering used by GEOSX array2d, the inverse of <tt>A<tt> is practically
   * computed as the transpose matrix of the transpose matrix inverse using
   * lapack operation based on column-major ordering. This removes the need
   * for any copy/transposition that would be required operating with the
   * row-major layout.
   */
  void matrixInverse( array2d<real64> const & A,
                       array2d<real64> & Ainv ) const;

  /**
   * @brief Compute inverse; <tt>Ainv<tt> = <tt>A</tt><sup>-1</sup>.
   *
   * Assign the inverse of the given matrix <tt>A<tt> to <tt>Ainv<tt> and
   * return also the determinant of <tt>A<tt>.
   *
   * \param IN
   * <tt>A</tt> - GEOSX array2d.
   *
   * \param INOUT
   * <tt>Ainv</tt> - GEOSX array2d.
   *
   * \param INOUT
   * <tt>detA</tt> - Determinant of input matrix <tt>A</tt>
   *
   * @warning
   * Assumes <tt>Ainv<tt> already has the same size as <tt>A</tt>.
   *
   * @note This function is hardcoded for matrices up to order three.
   * For dimensions larger than three, the function calls LAPACK
   * functions DGETRF and DGETRI. Because of the row major
   * ordering used by GEOSX array2d, the inverse of <tt>A<tt> is practically
   * computed as the transpose matrix of the transpose matrix inverse using
   * lapack operation based on column-major ordering. This removes the need
   * for any copy/transposition that would be required operating with the
   * row-major layout.
   */
  void matrixInverse( array2d<real64> const & A,
                      array2d<real64> & Ainv,
                      real64 & detA) const;

  /**
   * @brief Vector copy;
   * <tt>y</tt> = <tt>x<tt>
   *
   * \param IN
   * <tt>x</tt> - GEOSX array1d.
   *
   * \param INOUT
   * <tt>y</tt> - GEOSX array1d.
   *
   * @warning
   * Assumes that <tt>x</tt> and <tt>y</tt> have the same size.
   *
   */
  void vectorCopy( array1d<real64> const & X,
                   array1d<real64> & Y ) const;

  /**
   * @brief Vector copy;
   * <tt>B</tt> = <tt>A<tt>
   *
   * \param IN
   * <tt>A</tt> - GEOSX array2d.
   *
   * \param INOUT
   * <tt>B</tt> - GEOSX array2d.
   *
   * @warning
   * Assumes that <tt>A</tt> and <tt>B</tt> have the same size.
   *
   */
  void matrixCopy( array2d<real64> const & A,
                   array2d<real64> & B ) const;
  //@}

  //----------------------------------------------------------------------------
  //! @name I/O methods
  //@{

  /**
   * @brief Print service method for GEOSX array1d.
   */
  void printVector(array1d<real64> const & X) const;

  /**
   * @brief Print service method for GEOSX array2d.
   */
  void printMatrix(array2d<real64> const & X) const;

  //@}

};

}

#endif /* CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLASLAPACKLA_HPP_ */
