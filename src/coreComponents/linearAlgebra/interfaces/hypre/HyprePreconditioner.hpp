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

class HyprePreconditioner final : public PreconditionerBase< HypreInterface >
{
public:

  using Base = PreconditionerBase< HypreInterface >;
  using Vector = typename Base::Vector;
  using Matrix = typename Base::Matrix;

  explicit HyprePreconditioner( LinearSolverParameters params );

  virtual ~HyprePreconditioner() override;

  virtual void compute( Matrix const & mat ) override;

  virtual void apply( Vector const & src, Vector & dst ) const override;

  virtual void clear() override;

  HypreMatrix const & getPreconditioningMatrix() const;

  HYPRE_Solver const & unwrapped() const;

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
