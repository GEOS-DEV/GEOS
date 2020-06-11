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
 * @file HyprePreconditioner.hpp
 */

#ifndef GEOSX_HYPREPRECONDITIONER_HPP
#define GEOSX_HYPREPRECONDITIONER_HPP

#include "linearAlgebra/solvers/PreconditionerBase.hpp"
#include "linearAlgebra/interfaces/hypre/HypreInterface.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

#include <memory>

/**
 * @name Hypre forward declarations.
 *
 * Forward declare hypre's solver structs and pointer aliases in order
 * to avoid including hypre headers and leaking into the rest of GEOSX.
 */
///@{

/// Hypre solver struct forward declaration
extern "C" struct hypre_Solver_struct;

/// Solver pointer alias
using HYPRE_Solver = hypre_Solver_struct *;

///@}

namespace geosx
{

/// Forward-declared struct that hosts pointers to preconditioner functions
struct HyprePrecFuncs;

/**
 * @brief Wrapper around hypre-based preconditioners.
 */
class HyprePreconditioner final : public PreconditionerBase< HypreInterface >
{
public:

  /// Alias for base type
  using Base = PreconditionerBase< HypreInterface >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename Base::Matrix;

  /// Allow for partial overload of Base::compute()
  using Base::compute;

  /**
   * @brief Constructor.
   * @param params preconditioner parameters
   */
  explicit HyprePreconditioner( LinearSolverParameters params );

  /**
   * @brief Destructor.
   */
  virtual ~HyprePreconditioner() override;

  /**
   * @brief Compute the preconditioner from a matrix.
   * @param mat the matrix to precondition.
   */
  virtual void compute( Matrix const & mat ) override;

  /**
   * @brief Apply operator to a vector
   * @param src Input vector (x).
   * @param dst Output vector (b).
   *
   * @warning @p src and @p dst cannot alias the same vector.
   */
  virtual void apply( Vector const & src, Vector & dst ) const override;

  virtual void clear() override;

  /**
   * @brief Access the underlying implementation.
   * @return reference to hypre preconditioner
   */
  HYPRE_Solver const & unwrapped() const;

  /**
   * @brief Access the underlying implementation.
   * @return reference to container of preconditioner functions.
   *
   * Intended for use by HypreSolver.
   */
  HyprePrecFuncs const & unwrappedFuncs() const;

private:

  void createAMG();

  void createILU();

  void createILUT();

  /// Parameters for all preconditioners
  LinearSolverParameters m_parameters;

  /// Pointer to the Hypre implementation
  HYPRE_Solver m_precond;

  /// Pointers to hypre functions to setup/solve/destroy preconditioner
  std::unique_ptr< HyprePrecFuncs > m_functions;
};

}

#endif //GEOSX_HYPREPRECONDITIONER_HPP
