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

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_COUPLEDSOLVER_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_COUPLEDSOLVER_HPP_

#include "physicsSolvers/SolverBase.hpp"
#include "constitutive/solid/PorousSolid.hpp"

#include <tuple>

namespace geos
{

template< typename ... SOLVERS >
class CoupledSolver : public SolverBase
{

public:

  /**
   * @brief main constructor for CoupledSolver Objects
   * @param name the name of this instantiation of CoupledSolver in the repository
   * @param parent the parent group of this instantiation of CoupledSolver
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

    // rm laterr
    // m_disp( nodeManager.getField< fields::solidMechanics::totalDisplacement >() );
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
      GEOS_THROW_IF( solver == nullptr,
                     GEOS_FMT( "Could not find solver '{}' of type {}",
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
  { GEOS_UNUSED_VAR( domain, dofManager ); }

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
  { GEOS_UNUSED_VAR( time_n, dt, domain, dofManager, localMatrix, localRhs ); }

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
    GEOS_MARK_FUNCTION;

    if( getNonlinearSolverParameters().couplingType() == NonlinearSolverParameters::CouplingType::FullyImplicit )
    {
      return fullyCoupledSolverStep( time_n, dt, cycleNumber, domain );
    }
    else if( getNonlinearSolverParameters().couplingType() == NonlinearSolverParameters::CouplingType::Sequential )
    {
      return sequentiallyCoupledSolverStep( time_n, dt, cycleNumber, domain );
    }
    else
    {
      GEOS_ERROR( "Invalid coupling type option." );
      return 0;
    }

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

  virtual void cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override
  {
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->cleanup( time_n, cycleNumber, eventCounter, eventProgress, domain );
    } );
    SolverBase::cleanup( time_n, cycleNumber, eventCounter, eventProgress, domain );
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
    return SolverBase::solverStep( time_n, dt, cycleNumber, domain );
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
    GEOS_MARK_FUNCTION;

    real64 dtReturn = dt;

    real64 dtReturnTemporary;

    Timestamp const meshModificationTimestamp = getMeshModificationTimestamp( domain );

    // First call Coupled Solver setup  (important for poromechanics initialization for sequentially coupled)
    implicitStepSetup( time_n, dt, domain );

    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {

      // Only build the sparsity pattern if the mesh has changed
      if( meshModificationTimestamp > solver->getSystemSetupTimestamp() )
      {
        solver->setupSystem( domain,
                             solver->getDofManager(),
                             solver->getLocalMatrix(),
                             solver->getSystemRhs(),
                             solver->getSystemSolution() );
        solver->setSystemSetupTimestamp( meshModificationTimestamp );
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
          solver->getSolverStatistics().initializeTimeStepStatistics(); // initialize counters for subsolvers
        } );
        resetStateToBeginningOfStep( domain );
      }

      // Increment the solver statistics for reporting purposes
      // Pass a "0" as argument (0 linear iteration) to skip the output of linear iteration stats at the end
      m_solverStatistics.logNonlinearIteration( 0 );

      // Nonlinear Acceleration (Aitken)
      beforeOuterIter( iter, domain );

      // Solve the subproblems nonlinearly
      forEachArgInTuple( m_solvers, [&]( auto & solver, auto idx )
      {
        GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "  Iteration {:2}: {}", iter+1, solver->getName() ) );
        dtReturnTemporary = solver->nonlinearImplicitStep( time_n,
                                                           dtReturn,
                                                           cycleNumber,
                                                           domain );

        mapSolutionBetweenSolvers( domain, idx() );

        if( idx() == 1 )
        {
          // Nonlinear Acceleration (Aitken): record unaccelerated averageMeanTotalStressIncrement
          afterGeomechanicsInnerLoop( domain );
        }

        if( dtReturnTemporary < dtReturn )
        {
          iter = 0;
          dtReturn = dtReturnTemporary;
        }
      } );

      // Check convergence of the outer loop
      isConverged = checkSequentialConvergence( iter,
                                                time_n,
                                                dtReturn,
                                                domain );
      if( isConverged )
      {
        // Save Time step statistics for the subsolvers
        forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
        {
          solver->getSolverStatistics().saveTimeStepStatistics();
        } );
        break;
      }
      else
      {
        // Nonlinear Acceleration (Aitken)
        afterOuterIter( iter, domain );
      }
      // Add convergence check:
      ++iter;
    }

    GEOS_ERROR_IF( !isConverged, getName() << "::sequentiallyCoupledSolverStep did not converge!" );

    implicitStepComplete( time_n, dt, domain );

    return dtReturn;
  }

  /**
   * @brief Maps the solution obtained from one solver to the fields used by the other solver(s)
   *
   * @param domain the domain partition
   * @param solverType the index of the solver withing this coupled solver.
   */
  virtual void mapSolutionBetweenSolvers( DomainPartition & domain, integer const solverType )
  {
    GEOS_UNUSED_VAR( domain, solverType );
  }

  bool checkSequentialConvergence( int const & iter,
                                   real64 const & time_n,
                                   real64 const & dt,
                                   DomainPartition & domain ) const
  {
    NonlinearSolverParameters const & params = getNonlinearSolverParameters();
    bool isConverged = true;

    if( params.m_subcyclingOption == 0 )
    {
      GEOS_LOG_LEVEL_RANK_0( 1, "***** Single Pass solver, no subcycling *****\n" );
    }
    else
    {
      if( params.sequentialConvergenceCriterion() == NonlinearSolverParameters::SequentialConvergenceCriterion::ResidualNorm )
      {
        GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "  Iteration {:2}: outer-loop convergence check", iter+1 ) );
        real64 residualNorm = 0;

        // loop over all the single-physics solvers
        forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
        {

          solver->getLocalMatrix().toViewConstSizes().zero();
          solver->getSystemRhs().zero();
          arrayView1d< real64 > const localRhs = solver->getSystemRhs().open();

          // for each solver, we have to recompute the residual (and Jacobian, although not necessary)
          solver->assembleSystem( time_n,
                                  dt,
                                  domain,
                                  solver->getDofManager(),
                                  solver->getLocalMatrix().toViewConstSizes(),
                                  localRhs );
          solver->applyBoundaryConditions( time_n,
                                           dt,
                                           domain,
                                           solver->getDofManager(),
                                           solver->getLocalMatrix().toViewConstSizes(),
                                           localRhs );
          solver->getSystemRhs().close();

          // once this is done, we recompute the single-physics residual
          real64 const singlePhysicsNorm =
            solver->calculateResidualNorm( time_n,
                                           dt,
                                           domain,
                                           solver->getDofManager(),
                                           solver->getSystemRhs().values() );
          residualNorm += singlePhysicsNorm * singlePhysicsNorm;
        } );

        // finally, we perform the convergence check on the multiphysics residual
        residualNorm = sqrt( residualNorm );
        GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "    ( R ) = ( {:4.2e} ) ; ", residualNorm ) );

        isConverged = ( residualNorm < params.m_newtonTol );

      }
      else if( params.sequentialConvergenceCriterion() == NonlinearSolverParameters::SequentialConvergenceCriterion::NumberOfNonlinearIterations )
      {
        forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
        {
          NonlinearSolverParameters const & singlePhysicsParams = solver->getNonlinearSolverParameters();
          if( singlePhysicsParams.m_numNewtonIterations > singlePhysicsParams.m_minIterNewton )
          {
            isConverged = false;
          }
        } );
      }
      else
      {
        GEOS_ERROR( "Invalid sequential convergence criterion." );
      }

      if( isConverged )
      {
        GEOS_LOG_LEVEL_RANK_0( 1, "***** The iterative coupling has converged in " << iter + 1 << " iteration(s)! *****\n" );
      }
    }
    return isConverged;
  }

  virtual void
  postProcessInput() override
  {
    setSubSolvers();

    bool const isSequential = getNonlinearSolverParameters().couplingType() == NonlinearSolverParameters::CouplingType::Sequential;
    bool const usesLineSearch = getNonlinearSolverParameters().m_lineSearchAction != NonlinearSolverParameters::LineSearchAction::None;
    GEOS_THROW_IF( isSequential && usesLineSearch,
                   GEOS_FMT( "`{}`: line search is not supported by the coupled solver when {} is set to `{}`. Please set {} to `{}` to remove this error",
                             getName(),
                             NonlinearSolverParameters::viewKeysStruct::couplingTypeString(),
                             EnumStrings< NonlinearSolverParameters::CouplingType >::toString( NonlinearSolverParameters::CouplingType::Sequential ),
                             NonlinearSolverParameters::viewKeysStruct::lineSearchActionString(),
                             EnumStrings< NonlinearSolverParameters::LineSearchAction >::toString( NonlinearSolverParameters::LineSearchAction::None ) ),
                   InputError );
  }

  void
  synchronizeNonLinearParameters()
  {
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->getNonlinearSolverParameters().m_numNewtonIterations =
        m_nonlinearSolverParameters.m_numNewtonIterations;
    } );
  }

  /* Implementation of Nonlinear Acceleration (Aitken) of averageMeanTotalStressIncrement */

  void recordAverageMeanTotalStressIncrement( DomainPartition & domain,
                                              std::vector< real64 > & s )
  {
    // s denotes averageMeanTotalStressIncrement
    s.resize( 0 );
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 arrayView1d< string const > const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                            auto & subRegion )
      {
        // get the solid model (to access stress increment)
        string const solidName = subRegion.template getReference< string >( "porousMaterialNames" );
        geos::constitutive::CoupledSolidBase & solid = getConstitutiveModel< geos::constitutive::CoupledSolidBase >( subRegion, solidName );

        arrayView1d< const real64 > const & averageMeanTotalStressIncrement_k = solid.getAverageMeanTotalStressIncrement_k();
        for( localIndex k = 0; k < localIndex( averageMeanTotalStressIncrement_k.size() ); k++ )
        {
          s.push_back( averageMeanTotalStressIncrement_k[ k ] );
        }
      } );
    } );
  }

  void applyAcceleratedAverageMeanTotalStressIncrement( DomainPartition & domain,
                                                        std::vector< real64 > & s )
  {
    // s denotes averageMeanTotalStressIncrement
    integer i = 0;
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 arrayView1d< string const > const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                            auto & subRegion )
      {
        // get the solid model (to access stress increment)
        string const solidName = subRegion.template getReference< string >( "porousMaterialNames" );
        geos::constitutive::CoupledSolidBase & solid = getConstitutiveModel< geos::constitutive::CoupledSolidBase >( subRegion, solidName );
        auto & porosityModel = dynamic_cast< geos::constitutive::BiotPorosity const & >( solid.getBasePorosityModel() );
        arrayView1d< real64 > const & averageMeanTotalStressIncrement_k = solid.getAverageMeanTotalStressIncrement_k();
        for( localIndex k = 0; k < localIndex( averageMeanTotalStressIncrement_k.size() ); k++ )
        {
          porosityModel.updateAverageMeanTotalStressIncrement( k, s[ i ] );
          i++;
        }
      } );
    } );
  }

  std::vector< real64 > addTwoVecs( const std::vector< real64 > & vec1,
                                    const std::vector< real64 > & vec2,
                                    const real64 sign )
  {
    assert( vec1.size() == vec2.size() );
    std::vector< real64 > result;
    for( size_t i = 0; i < vec1.size(); i++ )
    {
      result.push_back( vec1[ i ] + sign * vec2[ i ] );
    }
    return result;
  }

  std::vector< real64 > scalarMultiplyAVec( const std::vector< real64 > & vec,
                                            const real64 scalarMult )
  {
    std::vector< real64 > result;
    for( size_t i = 0; i < vec.size(); i++ )
    {
      result.push_back( scalarMult * vec[ i ] );
    }
    return result;
  }

  real64 dotTwoVecs( const std::vector< real64 > & vec1,
                     const std::vector< real64 > & vec2 )
  {
    assert( vec1.size() == vec2.size );
    real64 result = 0;
    for( size_t i = 0; i < vec1.size(); i++ )
    {
      result += vec1[ i ] * vec2[ i ];
    }
    return result;
  }

  real64 computeAitkenRelaxationFactor( const std::vector< real64 > & s0,
                                        const std::vector< real64 > & s1,
                                        const std::vector< real64 > & s1_tilde,
                                        const std::vector< real64 > & s2_tilde,
                                        const real64 omega0 )
  {
    std::vector< real64 > r1 = addTwoVecs( s1_tilde, s0, -1.0 );
    std::vector< real64 > r2 = addTwoVecs( s2_tilde, s1, -1.0 );

    // diff = r2 - r1
    std::vector< real64 > diff = addTwoVecs( r2, r1, -1.0 );

    real64 denom = dotTwoVecs( diff, diff );
    real64 numer = dotTwoVecs( r1, diff );

    real64 omega1 = 1;
    if( !isZero( denom ) )
    {
      omega1 = -1.0 * omega0 * numer / denom;
    }
    return omega1;
  }

  std::vector< real64 > computeUpdate( const std::vector< real64 > & s1,
                                       const std::vector< real64 > & s2_tilde,
                                       const real64 omega1 )
  {
    return addTwoVecs( scalarMultiplyAVec( s1, 1.0 - omega1 ),
                       scalarMultiplyAVec( s2_tilde, omega1 ),
                       1.0 );
  }

  void beforeOuterIter( integer const & iter,
                        DomainPartition & domain )
  {
    if( iter == 0 )
    {
      recordAverageMeanTotalStressIncrement( domain, m_s1 );
    }
    else
    {
      m_s0 = m_s1;
      m_s1 = m_s2;
      m_s1_tilde = m_s2_tilde;
      m_omega0 = m_omega1;
    }
  }

  void afterGeomechanicsInnerLoop( DomainPartition & domain )
  {
    recordAverageMeanTotalStressIncrement( domain, m_s2_tilde );
  }

  void afterOuterIter( integer const & iter,
                       DomainPartition & domain )
  {
    if( iter == 0 )
    {
      m_s2 = m_s2_tilde;
      m_omega1 = 1.0;
    }
    else
    {
      m_omega1 = computeAitkenRelaxationFactor( m_s0, m_s1, m_s1_tilde, m_s2_tilde, m_omega0 );
      m_s2 = computeUpdate( m_s1, m_s2_tilde, m_omega1 );
      applyAcceleratedAverageMeanTotalStressIncrement( domain, m_s2 );
    }
  }

protected:

  /// Pointers of the single-physics solvers
  std::tuple< SOLVERS *... > m_solvers;

  /// Names of the single-physics solvers
  std::array< string, sizeof...( SOLVERS ) > m_names;

  /// member variables needed for Nonlinear Acceleration (Aitken)
  std::vector< real64 > m_s0;
  std::vector< real64 > m_s1;
  std::vector< real64 > m_s1_tilde;
  std::vector< real64 > m_s2;
  std::vector< real64 > m_s2_tilde;
  real64 m_omega0;
  real64 m_omega1;
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_COUPLEDSOLVER_HPP_ */
