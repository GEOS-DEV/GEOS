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
 * @file ThermalSinglePhasePoromechanicsFixedStressSolver.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_THERMALSINGLEPHASEPOROMECHANICSFIXEDSTRESSSOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_THERMALSINGLEPHASEPOROMECHANICSFIXEDSTRESSSOLVER_HPP_

#include "physicsSolvers/multiphysics/CoupledSolver.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"

namespace geosx
{

class ThermalSinglePhasePoromechanicsFixedStressSolver : public CoupledSolver< SolidMechanicsLagrangianFEM, SinglePhaseBase >
{
public:

  using Base = CoupledSolver< SolidMechanicsLagrangianFEM, SinglePhaseBase >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;

  ThermalSinglePhasePoromechanicsFixedStressSolver( const string & name,
                                                    Group * const parent );

  ~ThermalSinglePhasePoromechanicsFixedStressSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName()
  {
    return "ThermalSinglePhasePoromechanicsFixedStress";
  }

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "ThermalSinglePhasePoromechanicsFixedStress"; }

  enum class SolverType : integer
  {
    SolidMechanics = 0,
    Flow = 1
  };

  virtual void postProcessInput() override final;

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
  SinglePhaseBase * flowSolver() const
  {
    return std::get< toUnderlying( SolverType::Flow ) >( m_solvers );
  }

protected:

  virtual void initializePreSubGroups() override;

};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_THERMALSINGLEPHASEPOROMECHANICSFIXEDSTRESSSOLVER_HPP_ */
