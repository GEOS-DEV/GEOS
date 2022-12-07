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

/// Declare strings associated with enumeration values.

/**
 * @brief Coupling type.
  */
enum class CouplingType : integer
{
  FIM,        ///< Fully-implicit coupling
  Sequential  ///< Sequential coupling
};

ENUM_STRINGS( CouplingType,
              "FIM",
              "Sequential" );

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
    : SolverBase( name, parent ),
    m_subcyclingOption( 0 )
  {
    /// GEOS mainly uses FIM coupling so let's define FIM as the default.
    registerWrapper( viewKeyStruct::couplingTypeString(), &m_couplingType ).
       setInputFlag( dataRepository::InputFlags::OPTIONAL ).
       setApplyDefaultValue( CouplingType::FIM ).
       setDescription( "Type of coupling. Options are: Sequential and FIM" );

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


  real64
  solverStep( real64 const & time_n,
              real64 const & dt,
              int const cycleNumber,
              DomainPartition & domain ) override
  {
    GEOSX_MARK_FUNCTION;

    if( m_couplingType == CouplingType::FIM )
    {
      return fullyCoupledSolverStep( time_n, dt, cycleNumber, domain );
    }
    else if( m_couplingType == CouplingType::Sequential )
    {
      return sequentiallyCoupledSolverStep( time_n, dt, cycleNumber, domain );
    }else
    {
      GEOSX_ERROR("Invalid coupling type option.");
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

  real64 fullyCoupledSolverStep( real64 const & time_n,
                                 real64 const & dt,
                                 int const cycleNumber,
                                 DomainPartition & domain )
  {
    GEOSX_MARK_FUNCTION;

    real64 dtReturn = dt;

    // setup the coupled linear system
    setupSystem( domain, m_dofManager, m_localMatrix, m_rhs, m_solution );

    // setup reservoir and well systems
    implicitStepSetup( time_n, dt, domain );

    // currently the only method is implicit time integration
    dtReturn = nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

    // complete time step
    implicitStepComplete( time_n, dtReturn, domain );

    return dtReturn;
  }

  real64 sequentiallyCoupledSolverStep( real64 const & time_n,
                                        real64 const & dt,
                                        int const cycleNumber,
                                        DomainPartition & domain )
  {
    GEOSX_MARK_FUNCTION;

    real64 dtReturn = dt;

    real64 dtReturnTemporary;

    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->setupSystem( domain,
                           solver->getDofManager(),
                           solver->getLocalMatrix(),
                           solver->getSystemRhs(),
                           solver->getSystemSolution() );

      solver->implicitStepSetup( time_n, dt, domain );

    } );

    // need to check what to do with this. this->implicitStepSetup( time_n, dt, domain );
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
        GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1 << solver->getName() );
        dtReturnTemporary = solver->nonlinearImplicitStep( time_n,
                                                           dtReturn,
                                                           cycleNumber,
                                                           domain );

        if( solver->getNonlinearSolverParameters().m_numNewtonIterations == 0 && iter > 0 )
        {
          GEOSX_LOG_LEVEL_RANK_0( 1, "***** The iterative coupling has converged in " << iter << " iterations! *****\n" );
          isConverged = true;
        }
        else if( m_subcyclingOption == 0 && iter > 0 )
        {
          GEOSX_LOG_LEVEL_RANK_0( 1, "***** Single Pass solver, no subcycling *****\n" );
          isConverged = true;
        }

        mapSolutionBetweenSolvers( domain, idx() );

        if( dtReturnTemporary < dtReturn )
        {
          iter = 0;
          dtReturn = dtReturnTemporary;
        }
      } );

      if ( isConverged )
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

  virtual void mapSolutionBetweenSolvers( DomainPartition & Domain, integer const idx)
  {}

  virtual void
  postProcessInput() override
  {
    setSubSolvers();

    // We need to set the minimum number of newton's iterations to 0 for the sequentially
    // coupled approach to converge. 
    if( m_couplingType == CouplingType::Sequential )
    {
      forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
      {
        solver->getNonlinearSolverParameters().m_minIterNewton = 0;
      } );
    }
  }


  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static char const * couplingTypeString() { return "couplingType"; }
  };

  /// Pointers of the single-physics solvers
  std::tuple< SOLVERS *... > m_solvers;

  /// Names of the single-physics solvers
  std::array< string, sizeof...( SOLVERS ) > m_names;

  /// Type of coupling
  CouplingType m_couplingType;

  int m_subcyclingOption;

};



} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_COUPLEDSOLVER_HPP_ */
