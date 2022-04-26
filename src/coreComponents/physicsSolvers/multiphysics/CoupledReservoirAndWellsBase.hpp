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
 * @file CoupledReservoirAndWellsBase.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_COUPLEDRESERVOIRANDWELLSBASE_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_COUPLEDRESERVOIRANDWELLSBASE_HPP_

#include "physicsSolvers/multiphysics/CoupledSolver.hpp"

namespace geosx
{

template< typename RESERVOIR_SOLVER, typename WELL_SOLVER >
class CoupledReservoirAndWellsBase : public CoupledSolver< RESERVOIR_SOLVER, WELL_SOLVER >
{
public:

  using AbstractBase = CoupledSolver< RESERVOIR_SOLVER, WELL_SOLVER >;
  using AbstractBase::m_solvers;
  using AbstractBase::m_dofManager;
  using AbstractBase::m_localMatrix;
  using AbstractBase::m_rhs;
  using AbstractBase::m_solution;

  enum class SolverType : integer
  {
    Reservoir = 0,
    Well = 1
  };

  /**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  CoupledReservoirAndWellsBase ( const string & name,
                                 dataRepository::Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~CoupledReservoirAndWellsBase () override;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

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
  solveLinearSystem( DofManager const & dofManager,
                     ParallelMatrix & matrix,
                     ParallelVector & rhs,
                     ParallelVector & solution ) override;

  /**@}*/

  /**
   * @brief assembles the perforation rate terms
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void
  assembleCouplingTerms( real64 const time_n,
                         real64 const dt,
                         DomainPartition const & domain,
                         DofManager const & dofManager,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs ) = 0;

  /**
   * @brief accessor for the pointer to the reservoir solver
   * @return a pointer to the reservoir solver
   */
  RESERVOIR_SOLVER * getReservoirSolver() const;

  /**
   * @brief accessor for the pointer to the well solver
   * @return a pointer to the well solver
   */
  WELL_SOLVER * getWellSolver() const;

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    // solver that assembles the reservoir equations
    constexpr static char const * reservoirSolverNameString() { return "flowSolverName"; }

    // solver that assembles the well
    constexpr static char const * wellSolverNameString() { return "wellSolverName"; }
  };


protected:

  virtual void
  initializePostInitialConditionsPreSubGroups() override;

  /**
   * @Brief loop over perforations and increase Jacobian matrix row lengths for reservoir and well elements accordingly
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param rowLengths the row-by-row length
   */
  void
  addCouplingNumNonzeros( DomainPartition & domain,
                          DofManager & dofManager,
                          arrayView1d< localIndex > const & rowLengths ) const;

  /**
   * @Brief add the sparsity pattern induced by the perforations
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param pattern the sparsity pattern
   */
  virtual void
  addCouplingSparsityPattern( DomainPartition const & domain,
                              DofManager const & dofManager,
                              SparsityPatternView< globalIndex > const & pattern ) const = 0;

  /// solver that assembles the reservoir equations
  string m_reservoirSolverName;

  /// solver that assembles the well equations and compute perforation rates
  string m_wellSolverName;

};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_COUPLEDRESERVOIRANDWELLSBASE_HPP_ */
