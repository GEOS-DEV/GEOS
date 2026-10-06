/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file EigenSolverBase.hpp
 */

#ifndef GEOS_LINEARALGEBRA_SOLVERS_EIGENSOLVERBASE_HPP_
#define GEOS_LINEARALGEBRA_SOLVERS_EIGENSOLVERBASE_HPP_

#include "linearAlgebra/common/LinearOperator.hpp"
#include "linearAlgebra/utilities/EigenSolverParameters.hpp"
#include "linearAlgebra/utilities/MultiVectorOperations.hpp"

#include <memory>
#include <vector>

namespace geos
{

/**
 * @brief Outcome of a generalized eigenvalue solve.
 */
struct EigenSolverResult
{
  /// True if all the requested eigenpairs satisfied the convergence criterion
  bool converged = false;

  /// Number of eigenpairs that satisfied the convergence criterion
  integer numConverged = 0;

  /// Number of restarts (Arnoldi) or iterations (LOBPCG)
  integer numIterations = 0;

  /// Number of applications of the shift-and-invert operator or of the preconditioner
  integer numOperatorApplications = 0;

  /// Wall-clock time of the solve (seconds)
  real64 solveTime = 0.0;

  /// Eigenvalues lambda, ascending. These are Rayleigh quotients of the returned vectors.
  array1d< real64 > eigenvalues;

  /// Relative residuals ||K x - lambda M x||_2 / ( |lambda - shift| ||M x||_2 ) of the returned eigenpairs
  array1d< real64 > residuals;
};

/**
 * @brief Operators defining the generalized symmetric problem K x = lambda M x.
 * @tparam VECTOR type of the vectors
 *
 * K is symmetric positive semi-definite and M is symmetric positive semi-definite (typically diagonal). Rows of
 * K and M that are zero (constrained degrees of freedom) are tolerated: they produce zero entries in the
 * eigenvectors.
 */
template< typename VECTOR >
struct GeneralizedEigenProblem
{
  /// Stiffness-like operator K
  LinearOperator< VECTOR > const & stiffness;

  /// Mass-like operator M
  LinearOperator< VECTOR > const & mass;

  /// Applies (K - shift M)^{-1} (needed by spectral-transformation methods such as Arnoldi)
  LinearOperator< VECTOR > const * shiftedInverse = nullptr;

  /// Applies an approximation of (K - shift M)^{-1} (needed by preconditioned methods such as LOBPCG)
  LinearOperator< VECTOR > const * preconditioner = nullptr;

  /// Dimension of the space on which the problem lives, if smaller than the vector size (e.g. when constrained
  /// unknowns are kept in the vectors with zero rows). Zero means the vector size.
  globalIndex numUnknowns = 0;

  /// Optional vector with 1 on the unknowns of the problem and 0 on the constrained ones, used to generate
  /// starting vectors that vanish on constrained unknowns. Null if there is none.
  VECTOR const * freeMask = nullptr;

  /// Known eigenvectors (e.g. the rigid-body modes of a free structure). They do not need to be orthonormal,
  /// but they must span an invariant subspace of the pencil. The solver deflates them: it computes the
  /// remaining eigenpairs in the M-orthogonal complement, and returns the constraints as the first
  /// eigenvectors, so that numEigenvalues counts them.
  std::vector< VECTOR const * > constraints = {};
};

/**
 * @brief Interface of the generalized symmetric eigensolvers.
 * @tparam VECTOR type of the vectors
 *
 * Eigenvectors are M-orthonormal. Algorithms only access the problem through operator applications, so new
 * methods are added by deriving from this class and extending the factory.
 */
template< typename VECTOR >
class GeneralizedEigenSolver
{
public:

  /// Type of the vectors
  using Vector = VECTOR;

  /// Type of the problem description
  using Problem = GeneralizedEigenProblem< VECTOR >;

  /**
   * @brief Factory method.
   * @param parameters solver parameters
   * @return an owning pointer to the requested solver
   */
  static std::unique_ptr< GeneralizedEigenSolver< VECTOR > >
  create( EigenSolverParameters const & parameters );

  /**
   * @brief Constructor.
   * @param parameters solver parameters
   */
  explicit GeneralizedEigenSolver( EigenSolverParameters const & parameters ):
    m_params( parameters )
  {}

  /// Destructor
  virtual ~GeneralizedEigenSolver() = default;

  /**
   * @brief Compute the eigenpairs of the pencil (K, M) closest to the shift.
   * @param[in] problem the operators of the pencil
   * @param[in] prototype a vector used to size and distribute the returned eigenvectors
   * @param[out] modes the M-orthonormal eigenvectors, ordered as the eigenvalues in the result
   * @return the eigenvalues and convergence information
   */
  virtual EigenSolverResult solve( Problem const & problem,
                                   Vector const & prototype,
                                   std::vector< Vector > & modes ) const = 0;

  /// @return the solver parameters
  EigenSolverParameters const & parameters() const { return m_params; }

protected:

  /**
   * @brief Fill eigenvalues and residuals of the result from the vectors in @p modes and sort by eigenvalue.
   * @param[in] problem the operators of the pencil
   * @param[in,out] modes the eigenvectors (reordered)
   * @param[in,out] result the result, whose eigenvalues and residuals are (over)written
   *
   * Eigenvalues are the Rayleigh quotients of the M-normalized vectors. Residuals are measured on the
   * original pencil, independently of the algorithm used to build the vectors.
   */
  void finalizeResult( Problem const & problem,
                       std::vector< Vector > & modes,
                       EigenSolverResult & result ) const;

public:

  /**
   * @brief Orthonormal basis of the constraint space of a problem.
   */
  struct ConstraintSpace
  {
    /// M-orthonormal vectors
    std::vector< VECTOR > vectors;

    /// Their images by M
    std::vector< VECTOR > massVectors;

    /// @return the number of independent constraints
    integer size() const { return LvArray::integerConversion< integer >( vectors.size() ); }

    /**
     * @brief Remove from a vector its M-orthogonal projection on the constraint space.
     * @param[in,out] z the vector to project
     */
    void project( VECTOR & z ) const
    {
      if( vectors.empty() )
      {
        return;
      }
      std::vector< VECTOR const * > basis;
      std::vector< VECTOR const * > massBasis;
      for( size_t i = 0; i < vectors.size(); ++i )
      {
        basis.push_back( &vectors[i] );
        massBasis.push_back( &massVectors[i] );
      }
      array2d< real64 > products;
      std::vector< real64 > coefficients( vectors.size() );
      for( int pass = 0; pass < 2; ++pass )
      {
        multiVectorOperations::dots( massBasis, std::vector< VECTOR const * >{ & z }, products );
        for( size_t i = 0; i < vectors.size(); ++i )
        {
          coefficients[i] = -products( i, 0 );
        }
        multiVectorOperations::combine( basis, coefficients, z, true );
      }
    }
  };

  /**
   * @brief M-orthonormalize the constraints of a problem.
   * @param[in] problem the operators of the pencil
   * @param[in] prototype a vector used to size and distribute the new vectors
   * @return the constraint space (dependent constraints are dropped)
   */
  static ConstraintSpace makeConstraintSpace( Problem const & problem, Vector const & prototype );

  /**
   * @brief Create a vector with the same distribution as a prototype.
   * @param[in] prototype the vector to mimic
   * @return a new (zero) vector
   */
  static Vector makeVector( Vector const & prototype );

protected:

  /// Solver parameters
  EigenSolverParameters m_params;
};

} // namespace geos

#endif /* GEOS_LINEARALGEBRA_SOLVERS_EIGENSOLVERBASE_HPP_ */
