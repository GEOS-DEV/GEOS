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

#ifndef GEOSX_LINEARALGEBRA_SOLVERS_KRYLOVSOLVER_HPP_
#define GEOSX_LINEARALGEBRA_SOLVERS_KRYLOVSOLVER_HPP_

#include "linearAlgebra/utilities/BlockVectorView.hpp"
#include "linearAlgebra/utilities/BlockVector.hpp"
#include "linearAlgebra/utilities/BlockOperatorView.hpp"

namespace geosx
{

template< typename Vector > class LinearOperator;
template< typename Vector > class BlockVectorView;

//// StatsStruct
//struct KrylovConvegernceStats
//{
//  localIndex numIterations;
//  array1d<real64> relativeResidual;
//  bool convergenceFlag;
//};

namespace internal
{

// Helper struct to get temporary stored vector type from vector view types
template< typename VECTOR >
struct VectorStorage
{
  using type = VECTOR;
};

template< typename VECTOR >
struct VectorStorage< BlockVectorView< VECTOR > >
{
  using type = BlockVector< VECTOR >;
};

}

/**
 * @brief Base class for Krylov solvers
 * @tparam VECTOR type of vector handled by this solver
 */
template< typename VECTOR >
class KrylovSolver : public LinearOperator< VECTOR >
{
public:

  /// Alias for template parameter
  using Vector = typename LinearOperator< VECTOR >::Vector;

  /**
   * @brief Constructor
   */
  KrylovSolver( LinearOperator< Vector > const & A,
                LinearOperator< Vector > const & M,
                real64 const tolerance,
                localIndex const maxIterations,
                integer const verbosity );

  /**
   * @brief Virtual destructor (does nothing)
   */
  virtual ~KrylovSolver() override;

  /**
   * @brief Solve the system <tt>M^{-1}(Ax - b) = 0</tt> with CG
   * using monolithic GEOSX matrices.
   *
   * @param A system matrix.
   * @param x system solution (input = initial guess, output = solution).
   * @param b system right hand side.
   * @param M preconditioner.
   */
  virtual void
  solve( Vector const & b,
         Vector & x ) const = 0;

  virtual void
  multiply( Vector const & src,
            Vector & dst ) const override final;

protected:

  /// Alias for vector type that can be used for temporaries
  using VectorTemp = typename internal::VectorStorage< VECTOR >::type;

  /// reference to the operator to be solved
  LinearOperator<Vector> const & m_operator;

  /// reference to the preconditioning operator
  LinearOperator<Vector> const & m_precond;

  /// relative residual norm reduction tolerance
  real64 m_tolerance;

  /// maximum number of Krylov iterations
  localIndex m_maxIterations;

  /// solver verbosity level
  integer m_verbosity;

};

} //namespace geosx

#endif //GEOSX_LINEARALGEBRA_SOLVERS_KRYLOVSOLVER_HPP_
