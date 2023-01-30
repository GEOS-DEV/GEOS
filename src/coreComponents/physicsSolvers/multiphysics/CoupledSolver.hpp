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

    this->getWrapper< string >( SolverBase::viewKeyStruct::discretizationString() ).
      setInputFlag( dataRepository::InputFlags::FALSE );
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
    /// Fully-coupled assembly.

    // 1. we sync the nonlinear convergence history. The coupled solver parameters are the one being
    // used. We want to propagate the info to subsolvers. It can be important for solvers that
    // have special treatment for specific iterations.
    synchronizeNonLinearParameters();

    // 2. Assemble matrix blocks of each individual solver
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->assembleSystem( time_n, dt, domain, dofManager, localMatrix, localRhs );
    } );

    // 3. Assemble coupling blocks
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

  /// This method is meant to be kept final. Derived CoupledSolvers are expected, if needed,
  /// to override fullyCoupledSolverStep and/or sequentiallyCoupledSolverStep.
  real64
  solverStep( real64 const & time_n,
              real64 const & dt,
              int const cycleNumber,
              DomainPartition & domain ) override final
  {
    GEOSX_MARK_FUNCTION;

    if( getNonlinearSolverParameters().m_couplingType == NonlinearSolverParameters::CouplingType::FullyImplicit )
    {
      return fullyCoupledSolverStep( time_n, dt, cycleNumber, domain );
    }
    else if( getNonlinearSolverParameters().m_couplingType == NonlinearSolverParameters::CouplingType::Sequential )
    {
      return sequentiallyCoupledSolverStep( time_n, dt, cycleNumber, domain );
    }
    else
    {
      GEOSX_ERROR( "Invalid coupling type option." );
      return 0;
    }

  }


  virtual real64
  calculateResidualNorm( DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override
  {
    real64 norm = 0.0;
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      real64 const singlePhysicsNorm = solver->calculateResidualNorm( domain, dofManager, localRhs );
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

  virtual real64
  setNextDtBasedOnStateChange( real64 const & currentDt,
                               DomainPartition & domain ) override
  {
    real64 nextDt = LvArray::NumericLimits< real64 >::max;
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      real64 const singlePhysicsNextDt =
        solver->setNextDtBasedOnStateChange( currentDt, domain );
      nextDt = LvArray::math::min( singlePhysicsNextDt, nextDt );
    } );
    return nextDt;
  }


  /**@}*/

protected:

  /**
   * @brief Fully coupled solution approach solution step.
   *
   * @param time_n the current time
   * @param dt timestep size
   * @param cycleNumber
   * @param domain the domain partition
   * @return real64 size of the accepted timestep
   */
  virtual real64 fullyCoupledSolverStep( real64 const & time_n,
                                         real64 const & dt,
                                         int const cycleNumber,
                                         DomainPartition & domain )
  {
    GEOSX_MARK_FUNCTION;

    real64 dtReturn = dt;

    // setup the coupled linear system
    if( !m_systemSetupDone )
    {
      setupSystem( domain, m_dofManager, m_localMatrix, m_rhs, m_solution );
      m_systemSetupDone = true;
    }

    // setup reservoir and well systems
    implicitStepSetup( time_n, dt, domain );

    // currently the only method is implicit time integration
    dtReturn = nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

    // complete time step
    implicitStepComplete( time_n, dtReturn, domain );

    return dtReturn;
  }
  /**
   * @brief Sequentially coupled solver step. It solves a nonlinear system of
   * equations using a sequential approach.
   *
   * @param time_n the current time
   * @param dt timestep size
   * @param cycleNumber
   * @param domain the domain partition
   * @return real64 size of the accepted timestep
   */
  virtual real64 sequentiallyCoupledSolverStep( real64 const & time_n,
                                                real64 const & dt,
                                                int const cycleNumber,
                                                DomainPartition & domain )
  {
    GEOSX_MARK_FUNCTION;

    real64 dtReturn = dt;

    real64 dtReturnTemporary;

    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      if( !solver->systemSetupDone() )
      {
        solver->setupSystem( domain,
                             solver->getDofManager(),
                             solver->getLocalMatrix(),
                             solver->getSystemRhs(),
                             solver->getSystemSolution() );
        m_systemSetupDone = true;
      }

      solver->implicitStepSetup( time_n, dt, domain );

    } );

    NonlinearSolverParameters & solverParams = getNonlinearSolverParameters();
    integer & iter = solverParams.m_numNewtonIterations;
    iter = 0;
    bool isConverged = false;
    /// Sequential coupling loop
    while( iter < solverParams.m_maxIterNewton )
    {
      if( iter == 0 )
      {
        // Reset the states of all solvers if any of them had to restart
        forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
        {
          solver->resetStateToBeginningOfStep( domain );
        } );
        resetStateToBeginningOfStep( domain );
      }

      forEachArgInTuple( m_solvers, [&]( auto & solver, auto idx )
      {
        GEOSX_LOG_LEVEL_RANK_0( 1, GEOSX_FMT( "  Iteration {:2}: {}", iter+1, solver->getName() ) );
        dtReturnTemporary = solver->nonlinearImplicitStep( time_n,
                                                           dtReturn,
                                                           cycleNumber,
                                                           domain );

        mapSolutionBetweenSolvers( domain, idx() );

        if( dtReturnTemporary < dtReturn )
        {
          iter = 0;
          dtReturn = dtReturnTemporary;
        }
      } );

      isConverged = checkSequentialConvergence( iter );

      if( isConverged )
      {
        break;
      }
      // Add convergence check:
      ++iter;
    }

    GEOSX_ERROR_IF( !isConverged, getName() << "::sequentiallyCoupledSolverStep did not converge!" );

    implicitStepComplete( time_n, dt, domain );

    return dtReturn;
  }

  /**
   * @brief Maps the solution obtained from one solver to the fields used by the other solver(s)
   *
   * @param Domain the domain parition
   * @param solverType the index of the solver withing this coupled solver.
   */
  virtual void mapSolutionBetweenSolvers( DomainPartition & Domain, integer const solverType )
  {
    GEOSX_UNUSED_VAR( Domain, solverType );
  }

  bool checkSequentialConvergence( int const & iter ) const
  {
    bool isConverged = true;
    if( getNonlinearSolverParameters().m_subcyclingOption == 0 )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "***** Single Pass solver, no subcycling *****\n" );
    }
    else
    {
      // TODO: a better convergence check could/should be found.
      forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
      {
        if( solver->getNonlinearSolverParameters().m_numNewtonIterations > solver->getNonlinearSolverParameters().m_minIterNewton )
        {
          isConverged = false;
        }
      } );
      if( isConverged )
      {
        GEOSX_LOG_LEVEL_RANK_0( 1, "***** The iterative coupling has converged in " << iter + 1 << " iterations! *****\n" );
      }
    }
    return isConverged;
  }

  virtual void
  postProcessInput() override
  {
    setSubSolvers();
  }

  struct viewKeyStruct : SolverBase::viewKeyStruct {};

  void synchronizeNonLinearParameters()
  {
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->getNonlinearSolverParameters().m_numNewtonIterations =
        m_nonlinearSolverParameters.m_numNewtonIterations;
    } );
  }

  /// Pointers of the single-physics solvers
  std::tuple< SOLVERS *... > m_solvers;

  /// Names of the single-physics solvers
  std::array< string, sizeof...( SOLVERS ) > m_names;
};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_COUPLEDSOLVER_HPP_ */
