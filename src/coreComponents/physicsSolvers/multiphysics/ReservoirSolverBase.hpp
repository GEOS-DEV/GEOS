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
 * @file ReservoirSolverBase.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_RESERVOIRSOLVERBASE_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_RESERVOIRSOLVERBASE_HPP_

#include "physicsSolvers/SolverBase.hpp"

namespace geosx
{

class FlowSolverBase;
class WellSolverBase;

class ReservoirSolverBase : public SolverBase
{
public:

  /**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  ReservoirSolverBase( const string & name,
                       Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~ReservoirSolverBase() override;

  /// deleted copy constructor
  ReservoirSolverBase( ReservoirSolverBase const & ) = delete;

  /// default move constructor
  ReservoirSolverBase( ReservoirSolverBase && ) = default;

  /// deleted assignment operator
  ReservoirSolverBase & operator=( ReservoirSolverBase const & ) = delete;

  /// deleted move operator
  ReservoirSolverBase & operator=( ReservoirSolverBase && ) = delete;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName() { return "Reservoir"; }

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override final;

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  setupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               ParallelVector & rhs,
               ParallelVector & solution,
               bool const setSparsity = true ) override;

  virtual void
  assembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  virtual void
  applyBoundaryConditions( real64 const time,
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
  solveSystem( DofManager const & dofManager,
               ParallelMatrix & matrix,
               ParallelVector & rhs,
               ParallelVector & solution ) override;

  virtual real64
  scalingForSystemSolution( DomainPartition const & domain,
                            DofManager const & dofManager,
                            arrayView1d< real64 const > const & localSolution ) override;

  virtual bool
  checkSystemSolution( DomainPartition const & domain,
                       DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor ) override;

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override;

  virtual void updateState( DomainPartition & domain ) override;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;


  virtual void
  implicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain ) override;

  virtual real64
  solverStep( real64 const & time_n,
              real64 const & dt,
              int const cycleNumber,
              DomainPartition & domain ) override;

  /**@}*/

  /**
   * @Brief assembles the perforation rate terms
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void assembleCouplingTerms( real64 const time_n,
                                      real64 const dt,
                                      DomainPartition const & domain,
                                      DofManager const & dofManager,
                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                      arrayView1d< real64 > const & localRhs ) = 0;

  FlowSolverBase * getFlowSolver() const { return m_flowSolver; }

  WellSolverBase * getWellSolver() const { return m_wellSolver; }

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    // solver that assembles the reservoir equations
    constexpr static char const * flowSolverNameString() { return "flowSolverName"; }

    // solver that assembles the well
    constexpr static char const * wellSolverNameString() { return "wellSolverName"; }
  };


protected:

  virtual void initializePostInitialConditionsPreSubGroups() override;

  virtual void postProcessInput() override;

  void addCouplingNumNonzeros( DomainPartition & domain,
                               DofManager & dofManager,
                               arrayView1d< localIndex > const & rowLengths ) const;

  /**
   * @Brief add the sparsity pattern induced by the perforations
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param pattern the sparsity pattern
   */
  virtual void addCouplingSparsityPattern( DomainPartition const & domain,
                                           DofManager const & dofManager,
                                           SparsityPatternView< globalIndex > const & pattern ) const = 0;

  /**
   * @brief Setup stored views into domain data for the current step
   */
  virtual void resetViews( DomainPartition & domain );

  /// solver that assembles the reservoir equations
  string m_flowSolverName;

  /// solver that assembles the well equations and compute perforation rates
  string m_wellSolverName;

  /// pointer to the flow sub-solver
  FlowSolverBase * m_flowSolver;

  /// pointer to the well sub-solver
  WellSolverBase * m_wellSolver;

};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_RESERVOIRSOLVERBASE_HPP_ */
