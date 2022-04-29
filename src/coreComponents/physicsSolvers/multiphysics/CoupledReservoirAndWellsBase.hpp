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

#include "common/TimingMacros.hpp"

namespace geosx
{

namespace coupledReservoirAndWellsInternal
{

/**
 * @brief Utility function for the implementation details of initializePostInitialConditionsPreSubGroups
 */
void
initializePostInitialConditionsPreSubGroups( SolverBase * const solver );

/**
 * @brief Utility function for the implementation details of addCouplingNumZeros
 * @param solver the coupled solver
 * @param domain the physical domain object
 * @param dofManager degree-of-freedom manager associated with the linear system
 * @param rowLengths the row-by-row length
 * @param resNumDof number of reservoir element dofs
 * @param wellNumDof number of well element dofs
 * @param resElemDofName name of the reservoir element dofs
 * @param wellElemDofName name of the well element dofs
 */
void
addCouplingNumNonzeros( SolverBase const * const solver,
                        DomainPartition & domain,
                        DofManager & dofManager,
                        arrayView1d< localIndex > const & rowLengths,
                        integer const resNumDof,
                        integer const wellNumDof,
                        string const & resElemDofName,
                        string const & wellElemDofName );


/**
 * @brief Utility function for the implementation details of assembleSystem
 * @param reservoirSolver the reservoir solver
 * @param wellSolver the well solver
 * @param time_n the time at the beginning of the time step
 * @param dt the time step size
 * @param domain the domain partition
 * @param dofManager the dof manager
 * @param localMatrix the local Jacobian matrix
 * @param localRhs the local residual
 */
void
assembleSinglePhysicsSystems( SolverBase * const reservoirSolver,
                              SolverBase * const wellSolver,
                              real64 const time_n,
                              real64 const dt,
                              DomainPartition & domain,
                              DofManager const & dofManager,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs );

/**
 * @brief Utility function for the implementation details of solveLinearSystem
 * @param solver the coupled solver
 * @param dofManager the degree of freedom manager
 * @param matrix the Jacobian matrix
 * @param rhs the residual
 * @param solution the solution of the Jacobian system
 */
void
solveLinearSystem( SolverBase * const solver,
                   DofManager const & dofManager,
                   ParallelMatrix & matrix,
                   ParallelVector & rhs,
                   ParallelVector & solution );

}

template< typename RESERVOIR_SOLVER, typename WELL_SOLVER >
class CoupledReservoirAndWellsBase : public CoupledSolver< RESERVOIR_SOLVER, WELL_SOLVER >
{
public:

  using Base = CoupledSolver< RESERVOIR_SOLVER, WELL_SOLVER >;
  using Base::m_solvers;
  using Base::m_names;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;

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
                                 dataRepository::Group * const parent )
    : Base( name, parent )
  {
    this->template getWrapper< string >( Base::viewKeyStruct::discretizationString() ).
      setInputFlag( dataRepository::InputFlags::FALSE );
  }

  /**
   * @brief default destructor
   */
  virtual ~CoupledReservoirAndWellsBase () override {}

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
               bool const setSparsity = true ) override
  {
    GEOSX_MARK_FUNCTION;

    GEOSX_UNUSED_VAR( setSparsity );

    dofManager.setDomain( domain );

    Base::setupDofs( domain, dofManager );
    dofManager.reorderByRank();

    // Set the sparsity pattern without reservoir-well coupling
    SparsityPattern< globalIndex > patternDiag;
    dofManager.setSparsityPattern( patternDiag );

    // Get the original row lengths (diagonal blocks only)
    array1d< localIndex > rowLengths( patternDiag.numRows() );
    for( localIndex localRow = 0; localRow < patternDiag.numRows(); ++localRow )
    {
      rowLengths[localRow] = patternDiag.numNonZeros( localRow );
    }

    // Add the number of nonzeros induced by coupling on perforations
    addCouplingNumNonzeros( domain, dofManager, rowLengths.toView() );

    // Create a new pattern with enough capacity for coupled matrix
    SparsityPattern< globalIndex > pattern;
    pattern.resizeFromRowCapacities< parallelHostPolicy >( patternDiag.numRows(), patternDiag.numColumns(), rowLengths.data() );

    // Copy the original nonzeros
    for( localIndex localRow = 0; localRow < patternDiag.numRows(); ++localRow )
    {
      globalIndex const * cols = patternDiag.getColumns( localRow ).dataIfContiguous();
      pattern.insertNonZeros( localRow, cols, cols + patternDiag.numNonZeros( localRow ) );
    }

    // Add the nonzeros from coupling
    addCouplingSparsityPattern( domain, dofManager, pattern.toView() );

    // Finally, steal the pattern into a CRS matrix
    localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );
    localMatrix.setName( this->getName() + "/localMatrix" );

    rhs.setName( this->getName() + "/rhs" );
    rhs.create( dofManager.numLocalDofs(), MPI_COMM_GEOSX );

    solution.setName( this->getName() + "/solution" );
    solution.create( dofManager.numLocalDofs(), MPI_COMM_GEOSX );
  }


  virtual void
  assembleSystem( real64 const time_n,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override
  {
    coupledReservoirAndWellsInternal::
      assembleSinglePhysicsSystems( getReservoirSolver(),
                                    getWellSolver(),
                                    time_n, dt,
                                    domain,
                                    dofManager,
                                    localMatrix,
                                    localRhs );

    // assemble perforation rates in J_WR, J_RW, J_RR and J_WW
    assembleCouplingTerms( time_n, dt,
                           domain,
                           dofManager,
                           localMatrix,
                           localRhs );
  }

  virtual void
  solveLinearSystem( DofManager const & dofManager,
                     ParallelMatrix & matrix,
                     ParallelVector & rhs,
                     ParallelVector & solution ) override
  {
    GEOSX_MARK_FUNCTION;

    coupledReservoirAndWellsInternal::
      solveLinearSystem( this,
                         dofManager,
                         matrix,
                         rhs,
                         solution );
  }

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
  RESERVOIR_SOLVER * getReservoirSolver() const { return std::get< toUnderlying( SolverType::Reservoir ) >( m_solvers ); }

  /**
   * @brief accessor for the pointer to the well solver
   * @return a pointer to the well solver
   */
  WELL_SOLVER * getWellSolver() const { return std::get< toUnderlying( SolverType::Well ) >( m_solvers ); }

protected:

  virtual void
  postProcessInput() override
  {
    Base::postProcessInput();

    getWellSolver()->setFlowSolverName( m_names[toUnderlying( SolverType::Reservoir )] );
  }

  virtual void
  initializePostInitialConditionsPreSubGroups() override
  {
    Base::initializePostInitialConditionsPreSubGroups( );

    coupledReservoirAndWellsInternal::initializePostInitialConditionsPreSubGroups( this );
  }

  /**
   * @Brief loop over perforations and increase Jacobian matrix row lengths for reservoir and well elements accordingly
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param rowLengths the row-by-row length
   */
  void
  addCouplingNumNonzeros( DomainPartition & domain,
                          DofManager & dofManager,
                          arrayView1d< localIndex > const & rowLengths ) const
  {
    coupledReservoirAndWellsInternal::
      addCouplingNumNonzeros( this,
                              domain,
                              dofManager,
                              rowLengths,
                              getWellSolver()->numDofPerResElement(),
                              getWellSolver()->numDofPerWellElement(),
                              getWellSolver()->resElementDofName(),
                              getWellSolver()->wellElementDofName() );
  }

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

};


} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_COUPLEDRESERVOIRANDWELLSBASE_HPP_ */
