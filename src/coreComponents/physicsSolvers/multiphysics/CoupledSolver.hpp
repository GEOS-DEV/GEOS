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
 * @file CoupledSolver.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_COUPLEDSOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_COUPLEDSOLVER_HPP_

#include "physicsSolvers/SolverBase.hpp"

#include <tuple>

namespace geosx
{

template< typename ... SOLVERS >
class CoupledSolver : public SolverBase
{

public:

  /**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  CoupledSolver( const string & name,
                 Group * const parent )
    : SolverBase( name, parent )
  {
    forEachArgInTuple( m_solvers, [&]( auto solver, auto idx )
    {
      using SolverType = TYPEOFPTR( solver );
      string const key = SolverType::coupledSolverAttributePrefix() + "SolverName";
      registerWrapper( key, &m_names[idx()] ).
        setInputFlag( dataRepository::InputFlags::REQUIRED ).
        setDescription( "Name of the " + SolverType::coupledSolverAttributePrefix() + " solver used by the coupled solver" );
    } );
  }

  /// deleted copy constructor
  CoupledSolver( CoupledSolver const & ) = delete;

  /// default move constructor
  CoupledSolver( CoupledSolver && ) = default;

  /// deleted assignment operator
  CoupledSolver & operator=( CoupledSolver const & ) = delete;

  /// deleted move operator
  CoupledSolver & operator=( CoupledSolver && ) = delete;


  /**
   * @brief Utility function to set the subsolvers pointers using the names provided by the user
   */
  void
  setSubSolvers()
  {
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto idx )
    {
      using SolverPtr = TYPEOFREF( solver );
      using SolverType = TYPEOFPTR( SolverPtr {} );
      solver = this->getParent().template getGroupPointer< SolverType >( m_names[idx()] );
      GEOSX_THROW_IF( solver == nullptr,
                      GEOSX_FMT( "Could not find solver '{}' of type {}",
                                 m_names[idx()], LvArray::system::demangleType< SolverType >() ),
                      InputError );
    } );
  }


  /**
   * @brief Utility function to set the coupling between degrees of freedom
   * @param[in] domain the domain partition
   * @param[in] dofManager the dof manager
   */
  virtual void
  setupCoupling( DomainPartition const & domain,
                 DofManager & dofManager ) const
  { GEOSX_UNUSED_VAR( domain, dofManager ); }

  /**
   * @brief Utility function to compute coupling terms
   * @param[in] time_n the time at the beginning of the time step
   * @param[in] dt the time step size
   * @param[in] domain the domain partition
   * @paran[in] dofManager the degree of freedom manager
   * @param[in] localMatrix the local matrix
   * @param[in] localRhs the local rhs
   */
  virtual void
  assembleCouplingTerms( real64 const time_n,
                         real64 const dt,
                         DomainPartition const & domain,
                         DofManager const & dofManager,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs )
  { GEOSX_UNUSED_VAR( time_n, dt, domain, dofManager, localMatrix, localRhs ); }

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override
  {
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->setupDofs( domain, dofManager );
    } );

    setupCoupling( domain, dofManager );
  }

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override
  {
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->implicitStepSetup( time_n, dt, domain );
    } );
  }

  virtual void
  implicitStepComplete( real64 const & time_n,
                        real64 const & dt,
                        DomainPartition & domain ) override
  {
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->implicitStepComplete( time_n, dt, domain );
    } );
  }

  virtual void
  assembleSystem( real64 const time_n,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override
  {
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->assembleSystem( time_n, dt, domain, dofManager, localMatrix, localRhs );
    } );

    assembleCouplingTerms( time_n, dt, domain, dofManager, localMatrix, localRhs );
  }

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override
  {
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->applySystemSolution( dofManager, localSolution, scalingFactor, domain );
    } );
  }

  virtual void
  updateState( DomainPartition & domain ) override
  {
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->updateState( domain );
    } );
  }

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override
  {
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->resetStateToBeginningOfStep( domain );
    } );
  }


  virtual real64
  solverStep( real64 const & time_n,
              real64 const & dt,
              int const cycleNumber,
              DomainPartition & domain ) override
  {
    GEOSX_MARK_FUNCTION;

    real64 dt_return = dt;

    // setup the coupled linear system
    setupSystem( domain, m_dofManager, m_localMatrix, m_rhs, m_solution );

    // setup reservoir and well systems
    implicitStepSetup( time_n, dt, domain );

    // currently the only method is implicit time integration
    dt_return = nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

    // complete time step in reservoir and well systems
    implicitStepComplete( time_n, dt_return, domain );

    return dt_return;
  }


  virtual real64
  calculateResidualNorm( real64 const & time_n,
                         real64 const & dt,
                         DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override
  {
    real64 norm = 0.0;
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      real64 const singlePhysicsNorm = solver->calculateResidualNorm( time_n, dt, domain, dofManager, localRhs );
      norm += singlePhysicsNorm * singlePhysicsNorm;
    } );

    return sqrt( norm );
  }

  virtual void
  applyBoundaryConditions( real64 const time_n,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override
  {
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );
    } );
  }

  virtual bool
  checkSystemSolution( DomainPartition const & domain,
                       DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor ) override
  {
    bool validSolution = true;
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      bool const validSinglePhysicsSolution = solver->checkSystemSolution( domain, dofManager, localSolution, scalingFactor );
      validSolution = validSolution && validSinglePhysicsSolution;
    } );
    return validSolution;
  }

  virtual real64
  scalingForSystemSolution( DomainPartition const & domain,
                            DofManager const & dofManager,
                            arrayView1d< real64 const > const & localSolution ) override
  {
    real64 scalingFactor = 1e9;
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      real64 const singlePhysicsScalingFactor = solver->scalingForSystemSolution( domain, dofManager, localSolution );
      scalingFactor = LvArray::math::min( scalingFactor, singlePhysicsScalingFactor );
    } );
    return scalingFactor;
  }

  /**@}*/

protected:

  virtual void
  postProcessInput() override
  {
    setSubSolvers();
  }

  /// Pointers of the single-physics solvers
  std::tuple< SOLVERS *... > m_solvers;

  /// Names of the single-physics solvers
  std::array< string, sizeof...( SOLVERS ) > m_names;

};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_COUPLEDSOLVER_HPP_ */
