/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_LINEARALGEBRA_SOLVERS_PRECONDITIONERTWOLEVEL_HPP_
#define GEOSX_LINEARALGEBRA_SOLVERS_PRECONDITIONERTWOLEVEL_HPP_

#include "linearAlgebra/common/LinearOperator.hpp"
#include "linearAlgebra/common/PreconditionerBase.hpp"

namespace geosx
{

struct twoLevelStructuredMesh;

/**
 * @brief Common interface for two-level preconditioning operator
 * @tparam LAI linear algebra interface providing vectors, matrices and solvers
 */
template< typename LAI >
class PreconditionerTwoLevel : public PreconditionerBase< LAI >
{
public:

  /// Alias for base type
  using Base = PreconditionerBase< LAI >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename Base::Matrix;

  /**
   * @brief Constructor.
   * @param numComp number of components in the field
   * @param precond the actual preconditioner to apply to filtered matrix (ownership transferred)
   */
  PreconditionerTwoLevel ( Matrix const & mat,
                           Matrix const & matCoarse );

  virtual ~PreconditionerTwoLevel() override;

  /**
   * @brief Clean up the preconditioner setup.
   *
   * Releases memory used and allows the matrix to be deleted cleanly.
   * This method should be called before the matrix used to compute the preconditioner
   * goes out of scope or is re-created. Some implementations require the matrix
   * to outlive the preconditioner (for example, Trilinos/ML may crash the program if
   * deleted after the matrix).
   *
   * @note Should be properly overridden in derived classes, which may call this method.
   */
   virtual void clear() override
   {
     m_coarseSolver.reset();
     m_coarseMat = nullptr;
     m_diagInv.reset();
   }

  /**
   * @brief Apply operator to a vector.
   *
   * @param src Input vector (src).
   * @param dst Output vector (dst).
   */
  virtual void apply( Vector const & src,
                      Vector & dst ) const override;

private:

  // Nested mesh data structure
  std::unique_ptr< twoLevelStructuredMesh > m_mesh;

  // Coarse operator
  Matrix const * m_coarseMat{};

  // Coarse solver
  std::unique_ptr< PreconditionerBase< LAI > > m_coarseSolver;

  // The inverse of the diagonal of the fine operator (Jacobi smoother)
  Vector m_diagInv;
};

}

#endif //GEOSX_LINEARALGEBRA_SOLVERS_PRECONDITIONERTWOLEVEL_HPP_
