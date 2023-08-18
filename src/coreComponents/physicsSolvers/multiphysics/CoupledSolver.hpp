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

// #include "constitutive/ConstitutivePassThru.hpp"
// #include "common/DataTypes.hpp"
// #include "constitutive/solid/SolidBase.hpp"
// #include "common/MpiWrapper.hpp"
// #include "mesh/ElementRegionManager.hpp"
// #include "mesh/CellElementSubRegion.hpp"
// #include "mesh/FaceElementSubRegion.hpp"
// #include "mesh/InterObjectRelation.hpp"
// #include "codingUtilities/Utilities.hpp"
// #include "constitutive/ConstitutiveManager.hpp"
// #include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
// #include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
// #include "constitutive/contact/ContactBase.hpp"
// #include "constitutive/NullModel.hpp"
// #include "mesh/DomainPartition.hpp"
// #include "mesh/MeshBody.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"
// #include "physicsSolvers/multiphysics/PoromechanicsFields.hpp"

// #include "finiteElement/BilinearFormUtilities.hpp"
// #include "finiteElement/LinearFormUtilities.hpp"
// #include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
// #include "linearAlgebra/interfaces/hypre/HypreVector.hpp"
// #include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"

#include "constitutive/solid/PorousSolid.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
// #include "constitutive/solid/porosity/PorosityFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"

// #include "linearAlgebra/interfaces/VectorBase.hpp"

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

    std::cout << "Beginning of sequentially Coupled Solver Step \n";

    real64 dtReturn = dt;

    real64 dtReturnTemporary;

    Timestamp const meshModificationTimestamp = getMeshModificationTimestamp( domain );

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

    // remove laterr
    std::vector<real64> res_v;
    std::vector<real64> omega_v;
    std::vector<real64> tot_str_v;
    
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

      // Increment the solver statistics for reporting purposes
      // Pass a "0" as argument (0 linear iteration) to skip the output of linear iteration stats at the end
      m_solverStatistics.logNonlinearIteration( 0 );

      // Nonlinear Acceleration
      aitkenRelaxation( iter, "BEFORE_OUTER_ITER", domain );

      // if (iter == 0)
      // {
      //   int start_ind = 0;
      //   int end_ind = m_subRegion_sizes[ 0 ];
      //   std::cout << "size of subRegion_sizes = " << int( m_subRegion_sizes.size() );
      //   std::cout << "start_ind = " << start_ind << ", end_ind = " << end_ind << std::endl;
      //   for (int i = 0; i < int( m_subRegion_sizes.size() ); i++)
      //   {
      //     std::vector<real64> v( m_s1.begin() + start_ind, m_s1.begin() + end_ind );
      //     saveVecToTextFile( v, "initialStress" + std::to_string( i ) + ".txt", " copying initial tot stress ");
      //     start_ind += m_subRegion_sizes[ i ]; 
      //     end_ind += m_subRegion_sizes[ i + 1 ];
      //     std::cout << "start_ind = " << start_ind << ", end_ind = " << end_ind << std::endl;
      //   }
      // }

      // std::fstream f;
      // f.open( "iter.txt", std::ios_base::out );
      // if ( !f )
      // {
      //   std::cout << "file not created" << std::endl;
      // }
      // else
      // {
      //   std::cout << "creating iter.txt" << std::endl;
      //     f << std::to_string(iter) << std::endl;
      //   f.close();
      // }

      // Solve the subproblems nonlinearly
      forEachArgInTuple( m_solvers, [&]( auto & solver, auto idx )
      {
        GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "  Iteration {:2}: {}", iter+1, solver->getName() ) );
        dtReturnTemporary = solver->nonlinearImplicitStep( time_n,
                                                           dtReturn,
                                                           cycleNumber,
                                                           domain );

        mapSolutionBetweenSolvers( domain, idx() );

        if (idx() == 1)
        {
          // Nonlinear Acceleration
          aitkenRelaxation( iter, "AFTER_OUTER_ITER", domain );
          tot_str_v.push_back( m_s2[ 33 ] );
          saveVecToTextFile( m_s2, "accTtlStrs.txt", " !! ");
        }
        mapSolutionBetweenSolvers( domain, idx() );

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
                                                domain,
                                                res_v );                                          

      // Nonlinear Acceleration
      // aitkenRelaxation( iter, "AFTER_OUTER_ITER", domain );
      omega_v.push_back( m_omega1 );
 
      // std::vector<real64> tot_str;
      // tot_str.push_back( recordAverageMeanTotalStressIncrement( domain )[33] );
      // saveVecToTextFile( tot_str, "tot_str"+ std::to_string(iter)+".txt", " !! ");

      // int start_index = 0;
      // int end_index = m_subRegion_sizes[ 0 ];
      // std::cout << "start_index = " << start_index << ", end_index = " << end_index << std::endl;
      // for (int i = 0; i < int( m_subRegion_sizes.size() ); i++)
      // {
      //   std::vector<real64> v( m_s2.begin() + start_index, m_s2.begin() + end_index );
      //   saveVecToTextFile( v, "accTtlStrs" + std::to_string( i ) + ".txt", " copying tot stress ");
      //   start_index += m_subRegion_sizes[ i ]; 
      //   end_index += m_subRegion_sizes[ i + 1 ];
      //   std::cout << "start_index = " << start_index << ", end_index = " << end_index << std::endl;
      // }

      // std::cout << "size of subRegion_sizes = " << m_subRegion_sizes.size() << std::endl;
      // for (int i = 0; i < int( m_subRegion_sizes.size() ); i++)
      // {
      //   std::cout << "size of subRegion_size " << i << " = " << m_subRegion_sizes[i] << std::endl;
      // }

      // saveVecToTextFile( m_s2, "accTtlStrs.txt", " !! ");
      if( isConverged )
      {
        saveVecToTextFile( res_v, "res.txt", "Copying res to res.txt");
        saveVecToTextFile( omega_v, "relaxationFactor.txt", "Copying relaxation factor to relaxationFactor.txt");
        saveVecToTextFile( tot_str_v, "TotalStressIncrement.txt", "Copying totalStressIncrement to TotalStressIncrement.txt");

        break;
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
                                   DomainPartition & domain,
                                   std::vector<real64> & res_v ) const
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

        // remove laterr
        res_v.push_back(residualNorm);     
        // std::cout << "Hello!" << std::endl;
        // std::cout << res_v.size() << std::endl;

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

  /* Helper funtions for Aitken's acceleration */
  void saveVecToTextFile( const std::vector<real64> & vec, 
                          const std::string filename,
                          const std::string consoleMessage)
  {
    std::fstream f;
    f.open( filename, std::ios_base::out );
    if ( !f )
    {
      std::cout << "file not created" << std::endl;
    }
    else
    {
      std::cout << consoleMessage << std::endl;
      f << vec.size() << std::endl;
      for ( size_t i = 0; i < vec.size(); i++ )
      {
        f << vec[ i ] << std::endl;
      }
      f.close();
    }
  }
  void recordSubRegionSizes( DomainPartition & domain )
  { m_subRegion_sizes.resize( 0 );
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

        // arrayView2d< real64 const > const meanStressIncrement_k = solid.getMeanEffectiveStressIncrement_k();
        arrayView1d< real64 > const & averageMeanStressIncrement_k = solid.getAverageMeanEffectiveStressIncrement_k();
        m_subRegion_sizes.push_back( int( averageMeanStressIncrement_k.size() ) );
      } );
    } );
  }

  std::vector<real64> recordAverageMeanTotalStressIncrement( DomainPartition & domain )
  {
    // s denotes averageMeanTotalStressIncrement
    std::vector<real64> s;
    // m_subRegion_sizes.resize(0);
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

        // arrayView2d< real64 const > const meanStressIncrement_k = solid.getMeanEffectiveStressIncrement_k();
        arrayView1d< real64 > const & averageMeanStressIncrement_k = solid.getAverageMeanEffectiveStressIncrement_k();
        arrayView1d< real64 const > const & biotCoefficient = solid.getBiotCoefficient();
        // arrayView1d< real64 const > const bulkModulus = solid.getBulkModulus();
        arrayView1d< real64 const > const & pres_k = subRegion.template getField< fields::flow::pressure_k >();
        arrayView1d< real64 > const & pres_n = subRegion.template getField< fields::flow::pressure_n >();
        for (int j = 0; j < averageMeanStressIncrement_k.size(); j++)
        {
          s.push_back( averageMeanStressIncrement_k[ j ] - biotCoefficient[ j ] * ( pres_k[ j ] - pres_n[ j ]) );
        }
        // m_subRegion_sizes.push_back( int( averageMeanStressIncrement_k.size() ) );
      } );
    } );
    return s;
  }

  // // std::vector<real64> recordVolumetricStrain( DomainPartition & domain )
  // // {
  // //   NodeManager const & nodeManager = domain.getMeshBody( 0 ).getBaseDiscretization().getNodeManager();
  // //   // NodeManager const & nodeManager = meshLevel.getNodeManager();
  // //   std::vector<real64> totalDisp;
  // //   if( nodeManager.hasField< fields::solidMechanics::incrementalDisplacement >() )
  // //   {
  // //     fields::solidMechanics::arrayViewConst2dLayoutincrementalDisplacement const & uhat =
  // //       nodeManager.getField< fields::solidMechanics::incrementalDisplacement >();
  // //     size_t nNodes = uhat.size();
  // //     std::cout << "uhat vector size = " << nNodes << std::endl;
  // //   }
  // // }

  std::vector<real64> recordTotalDisp( DomainPartition & domain )
  {
    NodeManager const & nodeManager = domain.getMeshBody( 0 ).getBaseDiscretization().getNodeManager();
    std::vector<real64> totalDisp;
    if( nodeManager.hasField< fields::solidMechanics::totalDisplacement >() )
    {
      fields::solidMechanics::arrayViewConst2dLayoutTotalDisplacement const & totalDisplacement =
        nodeManager.getField< fields::solidMechanics::totalDisplacement >();
      size_t nNodes = size_t( totalDisplacement.size() / 3 );
      std::cout << "totalDisplacement vector size = " << nNodes << std::endl;
      // totalDisp.resize( nNodes ); 
      for ( size_t i = 0; i < nNodes; i++ )
      {
        for ( size_t j = 0; j < 3; j++ )
        {
          totalDisp.push_back( totalDisplacement[ i ][ j ] );
        }
      }    
    }
    std::cout << "vector size = " << totalDisp.size() << std::endl;
    return totalDisp;
      // saveVecToTextFile( m_totalDisp, nNodes, "totalDisp.txt" );
  }

  std::vector<real64> addTwoVecs( const std::vector<real64> & vec1,  
                                  const std::vector<real64> & vec2,
                                  const real64 sign)
  {
    assert( vec1.size() == vec2.size() );
    std::vector<real64> result;
    for ( size_t i = 0; i < vec1.size(); i++ )
    {
      result.push_back( vec1[ i ] + sign * vec2[ i ] );
    }
    return result;
  }

  std::vector<real64> scalarMultiplyAVec( const std::vector<real64> & vec,
                                          const real64 scalarMult)
  {
    std::vector<real64> result;
    for ( size_t i = 0; i < vec.size(); i++ )
    {
      result.push_back( scalarMult * vec[ i ]);
    }
    return result;
  }

  real64 dotTwoVecs( const std::vector<real64> & vec1,
                     const std::vector<real64> & vec2 )
  {
    assert( vec1.size() == vec2.size );
    real64 result = 0;
    for ( size_t i = 0; i < vec1.size(); i++ )
    {
      result += vec1[ i ] * vec2[ i ];
    }
    return result;
  }

  real64 computeAitkenRelaxationFactor( const std::vector<real64> & s0,
                                        const std::vector<real64> & s1, 
                                        const std::vector<real64> & s1_tilde,
                                        const std::vector<real64> & s2_tilde, 
                                        const real64 omega0 )
  {
    // saveVecToTextFile( m_s0, "s0.txt", "Copying s0");
    // saveVecToTextFile( m_s1_tilde, "s1_tilde.txt", "Copying s1_tilde");
    // saveVecToTextFile( m_s1, "s1.txt", "Copying s1");
    // saveVecToTextFile( m_s2_tilde, "s2_tilde.txt", "Copying s2_tilde");

    std::cout << "Computing relaxation factor, omega:" << std::endl; 

    std::vector<real64> r1, r2;
    if ( m_scaling_flag )
    {
      std::vector<real64> s0_scaled, s1_tilde_scaled, s1_scaled, s2_tilde_scaled;
      performScaling( s0_scaled, s1_tilde_scaled, s1_scaled, s2_tilde_scaled );
      r1 = addTwoVecs( s1_tilde_scaled, s0_scaled, -1.0 );
      r2 = addTwoVecs( s2_tilde_scaled, s1_scaled, -1.0 );

      // saveVecToTextFile( s0_scaled, "s0_scaled.txt", "Copying s0_scaled");
      // saveVecToTextFile( s1_tilde_scaled, "s1_tilde_scaled.txt", "Copying s1_tilde_scaled");
      // saveVecToTextFile( s1_scaled, "s1_scaled.txt", "Copying s1_scaled");
      // saveVecToTextFile( s2_tilde_scaled, "s2_tilde_scaled.txt", "Copying s2_tilde_scaled");
      std::cout << "xDisp min = " << m_min_xDisp << ", xDisp max = " << m_max_xDisp << std::endl;
      std::cout << "yDisp min = " << m_min_yDisp << ", yDisp max = " << m_max_yDisp << std::endl;
      std::cout << "zDisp min = " << m_min_zDisp << ", zDisp max = " << m_max_zDisp << std::endl;
    }
    else
    {
      r1 = addTwoVecs( s1_tilde, s0, -1.0 );
      r2 = addTwoVecs( s2_tilde, s1, -1.0 );
    }

    // diff = r2 - r1
    std::vector<real64> diff = addTwoVecs( r2, r1, -1.0 );

    real64 denom = dotTwoVecs( diff, diff );
    real64 numer = dotTwoVecs( r1, diff );

    real64 omega1 = 1;
    if ( !isZero( denom ) )
    {
      omega1 = -1.0 * omega0 * numer / denom;
      // if ( omega1 > 5.0 )
      // {
      //   omega1 = 5.0;
      // }
      // if ( omega1 < -5 )
      // {
      //   omega1 = -5;
      // }
    }
    // omega1 = 1.0;
    std::cout << "omega1 = " << omega1 << std::endl;
    return omega1;
  }

  std::vector<real64> computeUpdate( const std::vector<real64> & s1,
                                     const std::vector<real64> & s2_tilde, 
                                     const real64 omega1 )
  {
    return addTwoVecs( scalarMultiplyAVec( s1, 1.0 - omega1 ), 
                       scalarMultiplyAVec( s2_tilde, omega1 ), 
                       1.0 );
  }

  void applyUpdateToTheSystemSolution( DomainPartition & domain )
  {
    std::cout << "Applying update to the solution:" << std::endl;
    assert( m_s2_tilde.size() == m_s2.size() );
    size_t size_ = m_s2.size();
    array1d< real64 > localSolution_array( size_ );

      for (size_t i = 0; i < size_; i++)
      {
        localSolution_array[ i ] = m_s2[ i ] - m_s2_tilde[ i ];
      }
      arrayView1d< real64 > const & localSolution = localSolution_array.toView();

    forEachArgInTuple( m_solvers, [&]( auto & solver, auto idx )
    {
      // TODO: use solver name instead of indexes, idx() = 1 should be geomechanics solver
      if( idx() == 1 )
      {
        solver->applySystemSolution( solver->getDofManager(), localSolution, 1.0, domain );
        solver->updateState( domain );
        // NodeManager const & nodeManager = domain.getMeshBody( 0 ).getBaseDiscretization().getNodeManager();
        // std::cout << "size of local Solution = " << solver->getDofManager().numGlobalDofs() << std::endl;
      }
      mapSolutionBetweenSolvers( domain, idx() );
    } );
    
  }

  // std::vector<real64> updatePorosityWithAccTotalStressIncrement ( DomainPartition & domain )
  // {
  //   std::vector<real64> poro;
  //   int i = 0;
  //   forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
  //                                                               MeshLevel & mesh,
  //                                                               arrayView1d< string const > const & regionNames )
  //   {
  //     mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
  //                                                                                           auto & subRegion )
  //     {
  //       // get the solid model (to access stress increment)
  //       string const solidName = subRegion.template getReference< string >( "porousMaterialNames" );
  //       geos::constitutive::CoupledSolidBase & solid = getConstitutiveModel< geos::constitutive::CoupledSolidBase >( subRegion, solidName );
  //       geos::constitutive::PorosityBase const & poroModel = solid.getBasePorosityModel();

  //       arrayView1d< real64 const > const & biotCoefficient = solid.getBiotCoefficient();
  //       arrayView1d< real64 const > const & bulkModulus = solid.getBulkModulus();
  //       arrayView2d< real64 const > const & poro_n = solid.getPorosity_n();
  //       // arrayView1d< real64 const > const & poro_ref = solid.getReferencePorosity();
  //       arrayView1d< real64 const > const & pres = subRegion.template getField< fields::flow::pressure >();
  //       // arrayView1d< real64 const > const & pres_k = subRegion.template getField< fields::flow::pressure_k >();
  //       arrayView1d< real64 > const & pres_n = subRegion.template getField< fields::flow::pressure_n >();
  //       // std::cout << "pres size 0 = " << pres.size(0) << std::endl;
  //       // std::cout << "pres size 1 = " << pres.size(1) << std::endl;
  //       std::cout << "pres_n size 0 = " << pres_n.size(0) << std::endl;
  //       std::cout << "pres_n size 1 = " << pres_n.size(1) << std::endl;
  //       std::cout << "poro_n size 0  = " << poro_n.size(0) << std::endl;
  //       std::cout << "poro_n size 1  = " << poro_n.size(1) << std::endl;
  //       for (int j = 0; j < pres_n.size(); j++)
  //       {
  //         for (int q = 0; q < poro_n.size( 1 ); q++)
  //         {
  //           real64 porosity = poro_n[j][q] + biotCoefficient[j] / bulkModulus[j] * m_s2[i] + 
  //                         biotCoefficient[j] * biotCoefficient[j] / bulkModulus[j] * ( pres[j] - pres_n[j] );
  //           poroModel.savePorosity( j, q, porosity, biotCoefficient[j] * biotCoefficient[j] / bulkModulus[j] );
  //           poro.push_back( porosity );
  //         }
  //         ++i;
  //       }
  //     } );
  //   } );
  //   std::cout << "i = " << i << std::endl;
  //   return poro;
  // }

  // template< typename POROUSWRAPPER_TYPE >
  // void updatePorosityAndPermeabilityFromPressureAndTemperatureAitken( POROUSWRAPPER_TYPE porousWrapper,
  //                                                                     CellElementSubRegion & subRegion,
  //                                                                     arrayView1d< real64 const > const & pressure,
  //                                                                     arrayView1d< real64 const > const & pressure_k,
  //                                                                     arrayView1d< real64 const > const & pressure_n,
  //                                                                     arrayView1d< real64 const > const & temperature,
  //                                                                     arrayView1d< real64 const > const & temperature_k,
  //                                                                     arrayView1d< real64 const > const & temperature_n,
  //                                                                     real64 const omega )
  // {
  //   forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_DEVICE ( localIndex const k )
  //   {

  //     for( localIndex q = 0; q < porousWrapper.numGauss(); ++q )
  //     {
  //       porousWrapper.updateStateFromPressureAndTemperature( k, q,
  //                                                           pressure[k],
  //                                                           pressure_k[k],
  //                                                           pressure_n[k],
  //                                                           temperature[k],
  //                                                           temperature_k[k],
  //                                                           temperature_n[k],
  //                                                           omega );
  //     }
  //   } );
  // }

  // void updatePorosityAndPermeabilityAitken( CellElementSubRegion & subRegion,
  //                                           real64 const omega ) const
  // {
  //   GEOS_MARK_FUNCTION;
  //   arrayView1d< real64 const > const & pressure = subRegion.getField< fields::flow::pressure >();
  //   arrayView1d< real64 const > const & pressure_n = subRegion.getField< fields::flow::pressure_n >();

  //   arrayView1d< real64 const > const & temperature = subRegion.getField< fields::flow::temperature >();
  //   arrayView1d< real64 const > const & temperature_n = subRegion.getField< fields::flow::temperature_n >();

  //   string const & solidName = subRegion.getReference< string >( "solidNames" );
  //   geos::constitutive::CoupledSolidBase & porousSolid = subRegion.template getConstitutiveModel< geos::constitutive::CoupledSolidBase >( solidName );

  //   GEOS_THROW_IF( subRegion.hasField< fields::flow::pressure_k >() !=
  //                 subRegion.hasField< fields::flow::temperature_k >(),
  //                 GEOS_FMT( "`{}` and `{}` must be either both existing or both non-existing on subregion {}",
  //                           fields::flow::pressure_k::key(), fields::flow::temperature_k::key(), subRegion.getName() ),
  //                 std::runtime_error );

  //   constitutive::ConstitutivePassThru< geos::constitutive::CoupledSolidBase >::execute( porousSolid, [=, &subRegion] ( auto & castedPorousSolid )
  //   {
  //     typename TYPEOFREF( castedPorousSolid ) ::KernelWrapper porousWrapper = castedPorousSolid.createKernelUpdates();

  //     if( subRegion.hasField< fields::flow::pressure_k >() && // for sequential simulations
  //         subRegion.hasField< fields::flow::temperature_k >() )
  //     {
  //       arrayView1d< real64 const > const & pressure_k = subRegion.getField< fields::flow::pressure_k >();
  //       arrayView1d< real64 const > const & temperature_k = subRegion.getField< fields::flow::temperature_k >();

  //       updatePorosityAndPermeabilityFromPressureAndTemperatureAitken( porousWrapper, subRegion,
  //                                                                     pressure, pressure_k, pressure_n,
  //                                                                     temperature, temperature_k, temperature_n, 
  //                                                                     omega );
  //     }
  //   } );
  // }

    // TODO: instead of copying, use pointers
  void beforeOuterIter( integer const & iter,
                        DomainPartition & domain )
  {
    std::cout << "Aitken Relaxation (Before outer iteration):" << std::endl;
    if ( iter == 0)
    {
      std::cout << "Iter = 0 so m_s1 is recorded" << std::endl;
      // m_s1 = recordTotalDisp( domain );
      m_s1 = recordAverageMeanTotalStressIncrement( domain );
      saveVecToTextFile( m_s1, "accTtlStrs.txt", " !! ");
      // recordSubRegionSizes( domain );
    }
    else 
    {
      std::cout << "iter = " << iter << std::endl;
      std::cout << "omega_n = " << m_omega0 << std::endl;
      std::cout << "omega_n+1 = " << m_omega1 << std::endl;
      m_s0 = m_s1;
      m_s1 = m_s2;
      m_s1_tilde = m_s2_tilde;
      m_omega0 = m_omega1;
      std::cout << "omega^n after update " << m_omega0 << std::endl;
    }
    
  }

  void afterOuterIter( integer const & iter,
                       DomainPartition & domain )
  {
    std::cout << "Aitken Relaxation (after outer iteration):" << std::endl;
    // m_s2_tilde = recordTotalDisp( domain );
    m_s2_tilde = recordAverageMeanTotalStressIncrement( domain );
    if ( iter >= 0)
    {
      m_s2 = m_s2_tilde;
      m_omega1 = 1.0;
    }
    else 
    {
      m_omega1 = computeAitkenRelaxationFactor( m_s0, m_s1, m_s1_tilde, m_s2_tilde, m_omega0 );
      m_s2 = computeUpdate( m_s1, m_s2_tilde, m_omega1 );
      
      
      if ( m_scaling_flag )
      {
        // scale back
        scaleBack( m_s2);
      }

      // std::vector<real64> beforeUpdate = recordTotalDisp( domain );
      // saveVecToTextFile( beforeUpdate, "beforeUpdate.txt", "Copying beforeUpdate");

      // applyUpdateToTheSystemSolution( domain );
      // std::vector<real64> poro = updatePorosityWithAccTotalStressIncrement ( domain );
      // saveVecToTextFile( poro, "Acc_poro.txt", "Copying Acc poro!");

    // forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
    //                                                              MeshLevel & mesh,
    //                                                              arrayView1d< string const > const & regionNames )
    // {

    //   mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
    //                                                                                         auto & subRegion )
    //   {
    //     updatePorosityAndPermeabilityAitken( subRegion, m_omega1 );
    //   } );
    // } );
      

      // std::vector<real64> afterUpdate = recordTotalDisp( domain );
      // saveVecToTextFile( afterUpdate, "afterUpdate.txt", "Copying afterUpdate");

      // saveVecToTextFile( m_s2_tilde, "s2_tilde.txt", "Copying s2_tilde");
      // saveVecToTextFile( m_s2, "s2.txt", "Copying s2");
      // std::cout << "min m_s2 = " << *std::min_element( m_s2.begin(), m_s2.end() ) << std::endl;
      // std::cout << "max m_s2 = " << *std::max_element( m_s2.begin(), m_s2.end() ) << std::endl;
    }
  }

  void aitkenRelaxation( integer const & iter, 
                         std::string const when,
                         DomainPartition & domain )
  {
    if ( when.compare( "BEFORE_OUTER_ITER" ) == 0 )
    {
      beforeOuterIter( iter, domain );
    }
    else if ( when.compare( "AFTER_OUTER_ITER" ) == 0 )
    {
      afterOuterIter( iter, domain );
    }
  }

  /* Min-max normalization */
  void computeMinMax( const int start_index,
                      real64 & min_,
                      real64 & max_ )
  {
    min_ = 0;
    max_ = 0;
    for (int i = start_index; i < int ( m_s0.size() ); i += 3)
    {
      min_ = std::min( min_, m_s0[ i ] );
      min_ = std::min( min_, m_s1_tilde[ i ] );
      min_ = std::min( min_, m_s1[ i ] );
      min_ = std::min( min_, m_s2_tilde[ i ] );
      max_ = std::max( max_, m_s0[ i ] );
      max_ = std::max( max_, m_s1_tilde[ i ] );
      max_ = std::max( max_, m_s1[ i ] );
      max_ = std::max( max_, m_s2_tilde[ i ] );
    }
  }

  void normalizationMinMax( const std::vector<real64> & v,
                            std::vector<real64> & v_scaled )
  {
    v_scaled.resize( v.size() );
    for (int i = 0; i < int ( v.size() ); i += 3)
    {
      v_scaled[ i ] = ( v[ i ] - m_min_xDisp ) / ( m_max_xDisp - m_min_xDisp );
      v_scaled[ i + 1 ] = ( v[ i + 1 ] - m_min_yDisp ) / ( m_max_yDisp - m_min_yDisp );
      v_scaled[ i + 2 ] = ( v[ i + 2 ] - m_min_zDisp ) / ( m_max_zDisp - m_min_zDisp );
    }
  }

  void performScaling( std::vector<real64> & s0_scaled,
                       std::vector<real64> & s1_tilde_scaled,
                       std::vector<real64> & s1_scaled,
                       std::vector<real64> & s2_tilde_scaled )
  {
    computeMinMax( 0, m_min_xDisp, m_max_xDisp ); // global min and max of total disp in the x-direction is computed
    computeMinMax( 1, m_min_yDisp, m_max_yDisp );
    computeMinMax( 2, m_min_zDisp, m_max_zDisp );

    normalizationMinMax( m_s0, s0_scaled ); // s0 is scaled to [0,1]
    normalizationMinMax( m_s1_tilde, s1_tilde_scaled );
    normalizationMinMax( m_s1, s1_scaled );
    normalizationMinMax( m_s2_tilde, s2_tilde_scaled );
  }

  void scaleBack( std::vector<real64> & v)
  {
    for (int i = 0; i < int( v.size() ); i += 3 )
    {
      v[ i ] = ( m_max_xDisp - m_min_xDisp ) * v[ i ] + m_min_xDisp;
      v[ i + 1 ] = ( m_max_yDisp - m_min_yDisp ) * v[ i + 1 ] + m_min_yDisp;
      v[ i + 2 ] = ( m_max_zDisp - m_min_zDisp ) * v[ i + 2 ] + m_min_zDisp;
    }
  }

protected:

  /// Pointers of the single-physics solvers
  std::tuple< SOLVERS *... > m_solvers;

  /// Names of the single-physics solvers
  std::array< string, sizeof...( SOLVERS ) > m_names;

  // // arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_disp;
  std::vector<real64> m_s0;
  std::vector<real64> m_s1;
  std::vector<real64> m_s1_tilde;
  std::vector<real64> m_s2;
  std::vector<real64> m_s2_tilde;
  real64 m_omega0;
  real64 m_omega1;
  real64 m_max_xDisp;
  real64 m_min_xDisp;
  real64 m_max_yDisp;
  real64 m_min_yDisp;
  real64 m_max_zDisp;
  real64 m_min_zDisp;
  const bool m_scaling_flag = false;
  std::vector<int> m_subRegion_sizes;

  // arrayView2d< real64 const > const m_meanStressIncrement_k;
  // arrayView1d< real64 > const m_averageMeanStressIncrement_k;
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_COUPLEDSOLVER_HPP_ */
