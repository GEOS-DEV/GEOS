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

/**
 * @file BlockPreconditioner.hpp
 */

#ifndef GEOS_LINEARALGEBRA_SOLVERS_BLOCKPRECONDITIONER_HPP_
#define GEOS_LINEARALGEBRA_SOLVERS_BLOCKPRECONDITIONER_HPP_

#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/common/PreconditionerBase.hpp"
#include "linearAlgebra/utilities/BlockOperator.hpp"
#include "linearAlgebra/utilities/BlockVector.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

namespace geos
{

/*
 * Since formulas in Doxygen are broken with 1.8.13 and certain versions of ghostscript,
 * keeping this documentation in a separate comment block for now. Should be moved into
 * documentation of BlockPreconditioner class.
 *
 * This class implements a 2x2 block preconditioner of the form:
 * @f$
 * M^{-1} =
 * \begin{bmatrix}
 * I_{1} & -M_{1}^{-1}A_{12} \\
 *       &  I_{2}
 * \end{bmatrix}
 * \begin{bmatrix}
 * M_{1}^{-1} &            \\
 *            & M_{2}^{-1}
 * \end{bmatrix}
 * \begin{bmatrix}
 *  I_{1}            &       \\
 * -A_{21}M_{1}^{-1} & I_{2}
 * \end{bmatrix}
 * @f$
 * with
 * @f$ M_1^{-1} ~= A_{11}^{-1} @f$, @f$ M_2^{-1} ~= S_{22}^{-1} @f$, and
 * @f$ S_{22} ~= A_{22} - A_{21} A_{11}^{-1} A_{12} @f$ being the (optional) approximate Schur complement.
 *
 * The first step (predictor solve with @f$ M_{1} @f$) is optional. Enabling it results in a more powerful
 * preconditioner, but involves two applications of @f$ M_{1}^{-1} @f$ instead of just one on each solve.
 *
 * The user provides individual block solvers @f$ M_{1} @f$ and @f$ M_{2} @f$ as well as
 * a description of the block split of monolithic matrix in terms of DOF components.
 */

/**
 * @brief General 2x2 block preconditioner.
 * @tparam LAI type of linear algebra interface providing matrix/vector types
 */
template< typename LAI >
class BlockPreconditioner : public PreconditionerBase< LAI >
{
public:

  /// Alias for the base type
  using Base = PreconditionerBase< LAI >;

  /// Alias for the vector type
  using Vector = typename Base::Vector;

  /// Alias for the matrix type
  using Matrix = typename Base::Matrix;

  /**
   * @brief Constructor.
   * @param params block preconditioner parameters
   */
  explicit BlockPreconditioner( LinearSolverParameters::Block params );

  /**
   * @brief Setup data for one of the two blocks.
   * @param blockIndex index of the block to set up
   * @param blockDofs choice of DoF components (from a monolithic system)
   * @param solver instance of the inner preconditioner for the block
   * @param scaling user-provided row scaling coefficient for this block
   *
   * @note While not strictly required, it is generally expected that @p blockDofs
   *       of the two blocks are non-overlapping and their union includes all
   *       DoF components in the monolithic system.
   */
  void setupBlock( localIndex const blockIndex,
                   std::vector< DofManager::SubComponent > blockDofs,
                   std::unique_ptr< PreconditionerBase< LAI > > solver,
                   real64 const scaling = 1.0 );

  /**
   * @brief Setup data for one of the two blocks.
   * @param blockIndex index of the block to set up
   * @param blockDofs choice of DoF components (from a monolithic system)
   * @param solver instance of the inner preconditioner for the block
   * @param scaling user-provided row scaling coefficient for this block
   *
   * @note This version does not take ownership of @p solver that must have its lifetime guaranteed elsewhere.
   *
   * @note While not strictly required, it is generally expected that @p blockDofs
   *       of the two blocks are non-overlapping and their union includes all
   *       DoF components in the monolithic system.
   */
  void setupBlock( localIndex const blockIndex,
                   std::vector< DofManager::SubComponent > blockDofs,
                   PreconditionerBase< LAI > * const solver,
                   real64 const scaling = 1.0 );

  /**
   * @brief Set user-provided prolongation operator for a block
   * @param blockIndex index of the block to set up
   * @param P the operator
   */
  void setProlongation( localIndex const blockIndex,
                        Matrix const & P );

  /**
   * @name PreconditionerBase interface methods
   */
  ///@{

  /**
   * @brief Compute the preconditioner from a matrix
   * @param mat the matrix to precondition
   */
  virtual void setup( Matrix const & mat ) override;

  /**
   * @brief Apply operator to a vector
   * @param src Input vector (x).
   * @param dst Output vector (b).
   *
   * @warning @p src and @p dst cannot alias the same vector (some implementations may allow this).
   */
  virtual void apply( Vector const & src, Vector & dst ) const override;

  virtual void clear() override;

  ///@}

  /**
   * @brief @return the block operator consisting of extracted matrix blocks
   */
  BlockOperator< Vector, Matrix > const & blocks() const
  {
    return m_matBlocks;
  }

private:

  /**
   * @brief Initialize/resize internal data structures for a new linear system.
   * @param mat the new system matrix
   */
  void reinitialize( Matrix const & mat );

  /**
   * @brief Apply block scaling to system blocks (which must be already extracted).
   */
  void applyBlockScaling();

  /**
   * @brief Compute and apply the Schur complement to (1,1)-block.
   */
  void computeSchurComplement();

  /// Block preconditioner parameters
  LinearSolverParameters::Block m_params;

  /// Description of dof components making up each of the two main blocks
  std::array< std::vector< DofManager::SubComponent >, 2 > m_blockDofs{};

  /// Pointers to prolongation operators for each sub-block
  std::array< Matrix const *, 2 > m_prolongators{};

  /// Prolongation operators for each sub-block
  std::array< Matrix, 2 > m_prolongatorsOwned{};

  /// Matrix blocks
  BlockOperator< Vector, Matrix > m_matBlocks;

  /// Individual block preconditioner pointers
  std::array< PreconditionerBase< LAI > *, 2 > m_solvers{};

  /// Individual block preconditioner operators (when owned)
  std::array< std::unique_ptr< PreconditionerBase< LAI > >, 2 > m_solversOwned{};

  /// Scaling of each block
  std::array< real64, 2 > m_scaling{};

  /// Internal vector of block residuals
  mutable BlockVector< Vector > m_rhs;

  /// Internal vector of block solutions
  mutable BlockVector< Vector > m_sol;
};

} //namespace geos

#endif //GEOS_LINEARALGEBRA_SOLVERS_BLOCKPRECONDITIONER_HPP_
