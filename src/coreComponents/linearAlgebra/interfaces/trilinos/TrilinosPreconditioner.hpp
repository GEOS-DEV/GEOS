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
 * @file PreconditionerAMG.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSAMG_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSAMG_HPP_

#include "linearAlgebra/solvers/PreconditionerBase.hpp"
#include "linearAlgebra/interfaces/trilinos/TrilinosInterface.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

#include <memory>

class Epetra_Operator;

namespace Teuchos
{
class ParameterList;
}

namespace geosx
{

class TrilinosPreconditioner final : public PreconditionerBase< TrilinosInterface >
{
public:

  using Base = PreconditionerBase< TrilinosInterface >;
  using Vector = typename Base::Vector;
  using Matrix = typename Base::Matrix;

  explicit TrilinosPreconditioner( LinearSolverParameters params );

  virtual ~TrilinosPreconditioner() override;

  virtual void compute( Matrix const & mat ) override;

  virtual void apply( Vector const & src, Vector & dst ) const override;

  virtual void clear() override;

  Epetra_Operator const & unwrapped() const;

  Epetra_Operator & unwrapped();

private:

  /// Parameters for all preconditioners
  LinearSolverParameters m_parameters;

  /// Pointer to the Trilinos implementation
  std::unique_ptr< Epetra_Operator > m_precond;
};

}

#endif //GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSAMG_HPP_
