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
 * @file TrilinosTpetraPreconditioner.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSTPETRAPRECONDITIONER_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSTPETRAPRECONDITIONER_HPP_

#include "linearAlgebra/solvers/PreconditionerBase.hpp"
#include "linearAlgebra/interfaces/trilinos/TrilinosTpetraInterface.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

#include <Tpetra_Operator_fwd.hpp>
#include <Tpetra_RowMatrix_fwd.hpp>

#include <memory>

namespace geosx
{

/**
 * @brief Wrapper around Trilinos Tpetra-based preconditioners.
 */
class TrilinosTpetraPreconditioner final : public PreconditionerBase< TrilinosTpetraInterface >
{
public:

  /// Alias for base type
  using Base = PreconditionerBase< TrilinosTpetraInterface >;

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
  explicit TrilinosTpetraPreconditioner( LinearSolverParameters params );

  /**
   * @brief Destructor.
   */
  virtual ~TrilinosTpetraPreconditioner() override;

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

  /// Instantiation of Tpetra::Operator template held by this class
  using Tpetra_Operator = Tpetra::Operator< real64, int, globalIndex >;

  /// Instantiation of Tpetra::RowMatrix template used by this class
  using Tpetra_RowMatrix = Tpetra::RowMatrix< real64, int, globalIndex >;

  /**
   * @brief Access the underlying Epetra preconditioning operator.
   * @return reference to the Epetra operator
   */
  Tpetra_Operator const & unwrapped() const;

  /**
   * @copydoc unwrapped() const
   */
  Tpetra_Operator & unwrapped();

private:

  /// Parameters for all preconditioners
  LinearSolverParameters m_parameters;

  /// Pointer to the Trilinos implementation
  std::unique_ptr< Tpetra_Operator > m_precond;
};

} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSTPETRAPRECONDITIONER_HPP_
