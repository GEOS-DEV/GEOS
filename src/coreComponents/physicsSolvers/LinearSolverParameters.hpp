/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file LinearSolverParameters.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_LINEARSOLVERPARAMETERS_HPP_
#define GEOSX_PHYSICSSOLVERS_LINEARSOLVERPARAMETERS_HPP_

#include "dataRepository/Group.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

namespace geosx
{

/**
 * @brief Linear solver parameters with Group capabilities
 *
 * This class is a derived version of LinearSolverParameters with
 * dataRepository::Group capabilities to allow for XML input.
 *
 * Developer note: This class exposes LAI solver settings to external users.
 * As a general philosophy, only a subset of frequently tuned parameters should
 * be exposed.  Many advanced parameters can be set by the PhysicsSolver itself
 * since it has knowledge of the underlying problem (e.g. isSymmetric = true,
 * dofsPerNode = 3, etc.).  While we want to enable power users to tune the
 * solvers, most users prefer a short list of options with good default settings
 */
class LinearSolverParametersInput : public dataRepository::Group
{
public:

  LinearSolverParametersInput() = delete;

  /// Constructor
  LinearSolverParametersInput( std::string const & name, Group * const parent );

  /// Copy constructor
  LinearSolverParametersInput( LinearSolverParametersInput && ) = default;

  /// Destructor
  virtual ~LinearSolverParametersInput() override = default;

  /// Catalog name
  static string CatalogName() { return "LinearSolverParameters"; }

  /// Postprocessing of input
  virtual void PostProcessInput() override;

  LinearSolverParameters const & get() const
  { return m_parameters; }

  LinearSolverParameters & get()
  { return m_parameters; }

  /// Keys appearing in XML
  struct viewKeyStruct
  {
    static constexpr auto solverTypeString         = "solverType";         ///< Solver type key
    static constexpr auto preconditionerTypeString = "preconditionerType"; ///< Preconditioner type key
    static constexpr auto stopIfErrorString        = "stopIfError";        ///< stop if error key

    static constexpr auto directCheckResidualString = "directCheckResidual";  ///< direct solver check residual key
    static constexpr auto directEquilString         = "directEquil";          ///< direct solver equilibrate key
    static constexpr auto directColPermString       = "directColPerm";        ///< direct solver columns permutation key
    static constexpr auto directRowPermString       = "directRowPerm";        ///< direct solver rows permutation key
    static constexpr auto directReplTinyPivotString = "directReplTinyPivot";  ///< direct solver replace tiny pivot key
    static constexpr auto directIterRefString       = "directIterRef";        ///< direct solver iterative refinement key
    static constexpr auto directParallelString      = "directParallel";       ///< direct solver parallelism key

    static constexpr auto krylovMaxIterString     = "krylovMaxIter";     ///< Krylov max iterations key
    static constexpr auto krylovMaxRestartString  = "krylovMaxRestart";  ///< Krylov max iterations key
    static constexpr auto krylovTolString         = "krylovTol";         ///< Krylov tolerance key
    static constexpr auto krylovAdaptiveTolString = "krylovAdaptiveTol"; ///< Krylov adaptive tolerance key
    static constexpr auto krylovWeakTolString     = "krylovWeakestTol";  ///< Krylov weakest tolerance key

    static constexpr auto amgNumSweepsString     = "amgNumSweeps";       ///< AMG number of sweeps key
    static constexpr auto amgSmootherString      = "amgSmootherType";    ///< AMG smoother type key
    static constexpr auto amgCoarseString        = "amgCoarseSolver";    ///< AMG coarse solver key
    static constexpr auto amgThresholdString     = "amgThreshold";       ///< AMG threshold key
    static constexpr auto amgNullSpaceTypeString = "amgNullSpaceType";   ///< AMG near null space type key

    static constexpr auto iluFillString      = "iluFill";       ///< ILU fill key
    static constexpr auto iluThresholdString = "iluThreshold";  ///< ILU threshold key
  } viewKeys;

private:

  LinearSolverParameters m_parameters;

};

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_LINEARSOLVERPARAMETERS_HPP_
