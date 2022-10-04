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

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_PhaseFieldFractureSOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_PhaseFieldFractureSOLVER_HPP_

#include "codingUtilities/EnumStrings.hpp"
#include "physicsSolvers/SolverBase.hpp"

namespace geosx
{

/**
 * @brief class that implements a FEM-based solver for phase-field fracture problems
 *
 * This class implements a FEM-based solver for the phase-field description of brittle fracture.
 * In a nutshell, a regularized fracture surface is represented by a continuous variable, that attains
 * values from 0 to 1. This variable is viewed as a damaged parameter. It's evolution is governed
 * by the minimization of an energy functional that generalizes the theory of Griffith. This minimization
 * principle gives rise to a pair of PDEs describing the balance of linear momentum and the evolution
 * of the damage.
 *
 * The balance of linear momentum is implemented by modifying the stress calculators on the solid
 * consitutive objects (to include the effects of damage). The standard SolidMechanicsLagrangianFEM
 * solver can then be used to solve the stress equilibrium equation.
 *
 * The equation for the evolution of damage is implemented in the PhaseFieldDamageFEM solver.
 *
 * This class, PhaseFieldFractureSolver is responsible for implementing the solution scheme
 * that effects the coupling between these two solvers. This solution scheme is an alternate
 * minimization approach that solves each of the problems in an alternate fashion until convergence
 * is obtain. The subcycling option can be turned off if one wants a weaker coupling, without the
 * convergence check (single-pass solution). This can be useful for debugging purposes.
 *
 * In the future, a monolithic scheme can be consider, but those tend to require special techiniques
 * to avoid issues associated with the lack of convexity of the energy functional.
 *
 * Some references:
 *
 * Miehe, Christian; Hofacker, Martina; Welschinger, Fabian. A phase field model for rate-independent crack
 * propagation: Robust algorithmic implementation based on operator splits.
 * Computer Methods in Applied Mechianics and Engineering, v. 199, n. 45-48, p. 2765-2778, 2010
 *
 * Borden, Micheal J., et al. A phase-field description of dynamic brittle fracture.
 * Computer Methods in Applied Mechanics and Engineering, v. 217, p. 77-95, 2012
 *
 * Bourdin, Blaise; Francfort, Gille A.; Marigo, Jean-Jacques. The variational approach to fracture.
 * Journal of Elasticity, v. 91, n. 1-3, p. 5-148, 2008.
 *
 *
 */
class PhaseFieldFractureSolver : public SolverBase
{
public:
  PhaseFieldFractureSolver( const string & name,
                            Group * const parent );
  ~PhaseFieldFractureSolver() override;

  /**
   * @brief define a name for the PhaseFieldFractureSolver in the object catalog (interface with XML input file)
   * @return string used as XML tag in the input file
   */
  static string catalogName()
  {
    return "PhaseFieldFracture";
  }

  //documented on base class
  virtual void registerDataOnMesh( Group & MeshBodies ) override final;

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override final;

  virtual void
  implicitStepComplete( real64 const & time_n,
                        real64 const & dt,
                        DomainPartition & domain ) override final;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual real64
  solverStep( real64 const & time_n,
              real64 const & dt,
              int const cycleNumber,
              DomainPartition & domain ) override;

  /**
   * @brief implements staggered solution scheme for phase-field fracture, alternating between damage and mechanics solvers
   *
   * @param time_n time at the beginning of the step
   * @param dt the perscribed timestep
   * @param cycleNumber the current cycle number
   * @param domain the domain object
   * @return return the timestep that was achieved during the step
   *
   * The splitOperatorStep solves the coupled u-d problem by separating the coupled system in a mechanics problem and a damage problem.
   * It then solves these problems in an alternate fashion, transfering the solutions back and forth until convergence is obtained.
   * This approach circunvents the issue associated with the lack of convexity of the energy functional. However, it may display very
   * slow convergence sometimes.
   */
  real64 splitOperatorStep( real64 const & time_n,
                            real64 const & dt,
                            integer const cycleNumber,
                            DomainPartition & domain );

  /**
   * @brief this function updates the damage values at quadrature points in the constitutive object
   * @param domain the domain object
   *
   * The value of the damage at the quadrature points is needed to perform some constitutive updates. In the current implementation,
   * a m_damage( k,q ) exists in the SolidBase constitutive object. This field, defined at quadrature points must then be updated after
   * each solution step of the damage solver to maintain the coupling.
   */
  void mapDamageToQuadrature( DomainPartition & domain );

  /**
   * @brief enumeration with coupling options for phase-field fracture
   *
   * In this enumeration, to keep consitency with other solver, FixedStress is the staggered approach, and TightlyCoupled is the
   * monolithic one.
   */
  enum class CouplingTypeOption : integer
  {
    FixedStress,
    TightlyCoupled
  };

  /**
   * @brief struct with keys used in the xml file
   *
   */
  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static char const * couplingTypeOptionString() { return "couplingTypeOption"; }

    constexpr static char const * totalMeanStressString() { return "totalMeanStress"; }
    constexpr static char const * oldTotalMeanStressString() { return "oldTotalMeanStress"; }

    constexpr static char const * solidSolverNameString() { return "solidSolverName"; }
    constexpr static char const * damageSolverNameString() { return "damageSolverName"; }
    constexpr static char const * subcyclingOptionString() { return "subcycling"; }
  };

protected:

  //documentation in SolverBase
  virtual void postProcessInput() override final;

  virtual void initializePostInitialConditionsPreSubGroups() override final;

private:

  ///string with the name of the solid mechanics solver
  string m_solidSolverName;
  ///string with the name of the damage solver
  string m_damageSolverName;
  ///type of coupling between damage and solid mechancis
  CouplingTypeOption m_couplingTypeOption;
  ///option of turning off the convergence check in the staggered solution, useful for debugging
  integer m_subcyclingOption;

};
///register enum type in the data repository
ENUM_STRINGS( PhaseFieldFractureSolver::CouplingTypeOption,
              "FixedStress",
              "TightlyCoupled" );

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_PhaseFieldFractureSOLVER_HPP_ */
