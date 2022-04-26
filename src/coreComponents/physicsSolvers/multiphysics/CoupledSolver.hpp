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
  {}


  /// deleted copy constructor
  CoupledSolver( CoupledSolver const & ) = delete;

  /// default move constructor
  CoupledSolver( CoupledSolver && ) = default;

  /// deleted assignment operator
  CoupledSolver & operator=( CoupledSolver const & ) = delete;

  /// deleted move operator
  CoupledSolver & operator=( CoupledSolver && ) = delete;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override
  {
    forEachArgInTuple( std::tuple< SOLVERS *... >{}, [&]( auto t, auto idx )
    {
      GEOSX_UNUSED_VAR( t );
      auto & solver = std::get< idx() >( m_solvers );
      solver->setupDofs( domain, dofManager );
    } );
  }

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override
  {
    forEachArgInTuple( std::tuple< SOLVERS *... >{}, [&]( auto t, auto idx )
    {
      GEOSX_UNUSED_VAR( t );
      auto & solver = std::get< idx() >( m_solvers );
      solver->implicitStepSetup( time_n, dt, domain );
    } );
  }

  virtual void
  implicitStepComplete( real64 const & time_n,
                        real64 const & dt,
                        DomainPartition & domain ) override
  {
    forEachArgInTuple( std::tuple< SOLVERS *... >{}, [&]( auto t, auto idx )
    {
      GEOSX_UNUSED_VAR( t );
      auto & solver = std::get< idx() >( m_solvers );
      solver->implicitStepComplete( time_n, dt, domain );
    } );
  }

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override
  {
    forEachArgInTuple( std::tuple< SOLVERS *... >{}, [&]( auto t, auto idx )
    {
      GEOSX_UNUSED_VAR( t );
      auto & solver = std::get< idx() >( m_solvers );
      solver->applySystemSolution( dofManager, localSolution, scalingFactor, domain );
    } );
  }

  virtual void
  updateState( DomainPartition & domain ) override
  {
    forEachArgInTuple( std::tuple< SOLVERS *... >{}, [&]( auto t, auto idx )
    {
      GEOSX_UNUSED_VAR( t );
      auto & solver = std::get< idx() >( m_solvers );
      solver->updateState( domain );
    } );
  }

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override
  {
    forEachArgInTuple( std::tuple< SOLVERS *... >{}, [&]( auto t, auto idx )
    {
      GEOSX_UNUSED_VAR( t );
      auto & solver = std::get< idx() >( m_solvers );
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
  calculateResidualNorm( DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override
  {
    real64 norm = 0.0;
    forEachArgInTuple( std::tuple< SOLVERS *... >{}, [&]( auto t, auto idx )
    {
      GEOSX_UNUSED_VAR( t );
      auto & solver = std::get< idx() >( m_solvers );
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
    forEachArgInTuple( std::tuple< SOLVERS *... >{}, [&]( auto t, auto idx )
    {
      GEOSX_UNUSED_VAR( t );
      auto & solver = std::get< idx() >( m_solvers );
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
    forEachArgInTuple( std::tuple< SOLVERS *... >{}, [&]( auto t, auto idx )
    {
      GEOSX_UNUSED_VAR( t );
      auto & solver = std::get< idx() >( m_solvers );
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
    forEachArgInTuple( std::tuple< SOLVERS *... >{}, [&]( auto t, auto idx )
    {
      GEOSX_UNUSED_VAR( t );
      auto & solver = std::get< idx() >( m_solvers );
      real64 const singlePhysicsScalingFactor = solver->scalingForSystemSolution( domain, dofManager, localSolution );
      scalingFactor = LvArray::math::min( scalingFactor, singlePhysicsScalingFactor );
    } );
    return scalingFactor;
  }

  /**@}*/

protected:

  std::tuple< SOLVERS *... > m_solvers;
};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_COUPLEDSOLVER_HPP_ */
