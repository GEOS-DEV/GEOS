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
 * @file SinglePhasePoromechanicsLagrangianContactSolver.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSLAGRANGIANCONTACTSOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSLAGRANGIANCONTACTSOLVER_HPP_

#include "codingUtilities/EnumStrings.hpp"
#include "physicsSolvers/SolverBase.hpp"
#include "physicsSolvers/multiphysics/SinglePhasePoromechanicsSolver.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/contact/ContactExtrinsicData.hpp"

namespace geosx
{


class LagrangianContactSolver;
class LagrangianContactFlowSolver;
class SinglePhaseBase;

class SinglePhasePoromechanicsLagrangianContactSolver : public SinglePhasePoromechanicsSolver
{
public:
  SinglePhasePoromechanicsLagrangianContactSolver( const string & name,
                                                   Group * const parent );
  ~SinglePhasePoromechanicsLagrangianContactSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName() { return "SinglePhasePoromechanicsLagrangianContact"; }

  virtual void setupSystem( DomainPartition & domain,
                            DofManager & dofManager,
                            CRSMatrix< real64, globalIndex > & localMatrix,
                            ParallelVector & rhs,
                            ParallelVector & solution,
                            bool const setSparsity = true ) override;

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  virtual void
  assembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  virtual void
  applyBoundaryConditions( real64 const time_n,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override;

  virtual real64
  calculateResidualNorm( DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override;

  virtual void
  solveLinearSystem( DofManager const & dofManager,
                     ParallelMatrix & matrix,
                     ParallelVector & rhs,
                     ParallelVector & solution ) override;

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override;

  virtual void updateState( DomainPartition & domain ) override;


  virtual void
  implicitStepComplete( real64 const & time_n,
                        real64 const & dt,
                        DomainPartition & domain ) override;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual real64
  solverStep( real64 const & time_n,
              real64 const & dt,
              int const cycleNumber,
              DomainPartition & domain ) override;

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static char const * dPerm_dDisplacementString() { return "dPerm_dDisplacement"; }
    constexpr static char const * contactSolverNameString() { return "contactSolverName"; }
    constexpr static char const * contactFlowSolverNameString() { return "contactFlowSolverName"; }
    constexpr static char const * flowSolverNameString() { return "fluidSolverName"; }
    constexpr static char const * porousMaterialNamesString() { return "porousMaterialNames"; }
  };

  virtual bool
  lineSearch( real64 const & time_n,
              real64 const & dt,
              integer const cycleNumber,
              DomainPartition & domain,
              DofManager const & dofManager,
              CRSMatrixView< real64, globalIndex const > const & localMatrix,
              ParallelVector & rhs,
              ParallelVector & solution,
              real64 const scaleFactor,
              real64 & lastResidual ) override;

protected:

  virtual void postProcessInput() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

  string m_contactSolverName;
  string m_contactFlowSolverName;

  // pointer to the flow sub-solver
  SinglePhaseBase * m_flowSolver;

  // pointer to the contact mechanics sub-solver
  LagrangianContactSolver * m_contactSolver;
  LagrangianContactFlowSolver * m_contactFlowSolver;

private:

  integer m_activeSetIter = 0;
  void createPreconditioner( DomainPartition const & domain );

};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSLAGRANGIANCONTACTSOLVER_HPP_ */
