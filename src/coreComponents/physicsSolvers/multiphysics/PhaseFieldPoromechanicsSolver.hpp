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
 * @file PhaseFieldPoromechanicsSolver.hpp
 *
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_PHASEFIELDPOROMECHANICSSOLVER_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_PHASEFIELDPOROMECHANICSSOLVER_HPP_

#include "physicsSolvers/multiphysics/CoupledSolver.hpp"
#include "physicsSolvers/multiphysics/PhaseFieldFractureSolver.hpp"
#include "physicsSolvers/multiphysics/SinglePhasePoromechanics.hpp"
#include "physicsSolvers/simplePDE/PhaseFieldDamageFEM.hpp"

namespace geos
{

class PhaseFieldPoromechanicsSolver : public CoupledSolver< SinglePhasePoromechanics< SinglePhaseBase >, PhaseFieldDamageFEM >
{
public:

  using Base = CoupledSolver< SinglePhasePoromechanics< SinglePhaseBase >, PhaseFieldDamageFEM >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;

  PhaseFieldPoromechanicsSolver( const string & name,
                                 Group * const parent );

  ~PhaseFieldPoromechanicsSolver() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName()
  {
    return "PhaseFieldPoromechanics";
  }

  string getCatalogName() const override { return catalogName(); }

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "PhaseFieldPoromechanics"; }

  enum class SolverType : integer
  {
    Poromechanics = 0,
    Damage = 1
  };

  virtual void resetStateToBeginningOfStep( DomainPartition & domain ) override final;

  virtual void postInputInitialization() override final;

  /**
   * @brief accessor for the pointer to the poromechanics solver
   * @return a pointer to the poromechanics solver
   */
  SinglePhasePoromechanics< SinglePhaseBase > * poromechanicsSolver() const
  {
    return std::get< toUnderlying( SolverType::Poromechanics ) >( m_solvers );
  }

  /**
   * @brief accessor for the pointer to the flow solver
   * @return a pointer to the flow solver
   */
  PhaseFieldDamageFEM * damageSolver() const
  {
    return std::get< toUnderlying( SolverType::Damage ) >( m_solvers );
  }

  virtual void mapSolutionBetweenSolvers( DomainPartition & Domain, integer const idx ) override final;

  void applyDamageOnTractionBC( DomainPartition & domain );

protected:

  virtual void initializePostInitialConditionsPreSubGroups() override final {}

};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_PHASEFIELDPOROMECHANICSSOLVER_HPP_ */
