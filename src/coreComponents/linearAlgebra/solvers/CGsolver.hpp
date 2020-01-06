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
 * @file CGsolver.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_SOLVERS_CGSOLVER_HPP_
#define GEOSX_LINEARALGEBRA_SOLVERS_CGSOLVER_HPP_

namespace geosx
{

template< typename LAI > class BlockMatrixView;
template< typename LAI > class BlockVectorView;

/**
 * \class CGsolver
 * \brief This class creates and provides basic support for block
 *        CG (templated on the LA interface).
 * \note  The notation is consistent with "Iterative Methods for
 *        Linear and Non-Linear Equations" from C.T. Kelley (1995)
 *        and "Iterative Methods for Sparse Linear Systems"
 *        from Y. Saad (2003).
 */

template< typename LAI >
class CGsolver
{

  using ParallelMatrix = typename LAI::ParallelMatrix;
  using ParallelVector = typename LAI::ParallelVector;

public:

  //! @name Constructor/Destructor Methods
  //@{
  /**
   * @brief Empty constructor.
   *
   * Creates a solver object.
   */
  CGsolver();

  /**
   * @brief Virtual destructor.
   */
  ~CGsolver() = default;
  //@}

  /**
   * @brief Solve the system <tt>M^{-1}(Ax - b) = 0</tt> with CG
   * using monolithic GEOSX matrices.
   *
   * \param A system matrix.
   * \param x system solution (input = initial guess, output = solution).
   * \param b system right hand side.
   * \param M preconditioner.
   */
  void solve( ParallelMatrix const &A,
              ParallelVector &x,
              ParallelVector const &b,
              ParallelMatrix const &M );

  /**
   * @brief Solve the system <tt>M^{-1}(Ax - b) = 0</tt> with CG
   * using block GEOSX matrices.
   *
   * \param A system block matrix.
   * \param x system block solution (input = initial guess, output = solution).
   * \param b system block right hand side.
   * \param M block preconditioner.
   */
  void solve( BlockMatrixView<LAI> const &A,
              BlockVectorView<LAI> &x,
              BlockVectorView<LAI> const &b,
              BlockMatrixView<LAI> const &M );

private:

};

} // namespace GEOSX

#endif /*GEOSX_LINEARALGEBRA_SOLVERS_CGSOLVER_HPP_*/
