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
 * @file PhaseFieldFractureSolver.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_PHASEFIELDFRACTURESOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_PHASEFIELDFRACTURESOLVER_HPP_

#include "physicsSolvers/multiphysics/CoupledSolver.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/simplePDE/PhaseFieldDamageFEM.hpp"

namespace geosx
{

class PhaseFieldFractureSolver : public CoupledSolver< SolidMechanicsLagrangianFEM, PhaseFieldDamageFEM >
{
public:

  using Base = CoupledSolver< SolidMechanicsLagrangianFEM, PhaseFieldDamageFEM >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;
  using Base::m_couplingType;

  PhaseFieldFractureSolver( const string & name,
                            Group * const parent );

  ~PhaseFieldFractureSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName()
  {
    return "PhaseFieldFracture";
  }

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "PhaseFieldFracture"; }

  enum class SolverType : integer
  {
    SolidMechanics = 0,
    Damage = 1
  };

  /**
   * @brief accessor for the pointer to the solid mechanics solver
   * @return a pointer to the solid mechanics solver
   */
  SolidMechanicsLagrangianFEM * solidMechanicsSolver() const
  {
    return std::get< toUnderlying( SolverType::SolidMechanics ) >( m_solvers );
  }

  /**
   * @brief accessor for the pointer to the flow solver
   * @return a pointer to the flow solver
   */
  PhaseFieldDamageFEM * damageSolver() const
  {
    return std::get< toUnderlying( SolverType::Damage ) >( m_solvers );
  }

  virtual void registerDataOnMesh( Group & MeshBodies ) override final;

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override final;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual void mapSolutionBetweenSolvers( DomainPartition & Domain, int const idx ) override final;

  struct viewKeyStruct : Base::viewKeyStruct
  {
    constexpr static char const * totalMeanStressString() { return "totalMeanStress"; }
    constexpr static char const * oldTotalMeanStressString() { return "oldTotalMeanStress"; }
  };

protected:
  virtual void postProcessInput() override final {}

  virtual void initializePostInitialConditionsPreSubGroups() override final {}

};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_PHASEFIELDFRACTURESOLVER_HPP_ */
