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
 * @file HypreMGR.cpp
 */

#include "HypreMGR.hpp"

#include "linearAlgebra/common/common.hpp"

#include "linearAlgebra/interfaces/hypre/mgrStrategies/CompositionalMultiphaseFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/CompositionalMultiphaseHybridFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/CompositionalMultiphaseReservoirFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/CompositionalMultiphaseReservoirHybridFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/HybridSinglePhasePoromechanics.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/Hydrofracture.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/LagrangianContactMechanics.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/MultiphasePoromechanics.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/MultiphasePoromechanicsReservoirFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhaseHybridFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhasePoromechanics.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhasePoromechanicsReservoirFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhaseReservoirFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhaseReservoirHybridFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/ThermalCompositionalMultiphaseFVM.hpp"

#include "LvArray/src/output.hpp"

namespace geosx
{

void hypre::mgr::createMGR( LinearSolverParameters const & params,
                            DofManager const * const dofManager,
                            HyprePrecWrapper & precond,
                            HypreMGRData & mgrData )
{
  GEOSX_ERROR_IF( dofManager == nullptr, "MGR preconditioner requires a DofManager instance" );

  GEOSX_LAI_CHECK_ERROR( HYPRE_MGRCreate( &precond.ptr ) );

  // Hypre's parameters to use MGR as a preconditioner
  GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetTol( precond.ptr, 0.0 ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetMaxIter( precond.ptr, 1 ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetPrintLevel( precond.ptr, LvArray::math::min( LvArray::integerConversion< HYPRE_Int >( params.logLevel - 1), 0 ) ) );

  array1d< int > const numComponentsPerField = dofManager->numComponentsPerField();
  dofManager->getLocalDofComponentLabels( mgrData.pointMarkers );

  if( params.logLevel >= 1 )
  {
    GEOSX_LOG_RANK_0( numComponentsPerField );
  }
  if( params.logLevel >= 2 )
  {
    GEOSX_LOG_RANK_VAR( mgrData.pointMarkers );
  }

//  using namespace hypre::mgr;
  switch( params.mgr.strategy )
  {
    case LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseFVM:
    {
      setStrategy< CompositionalMultiphaseFVM >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseHybridFVM:
    {
      setStrategy< CompositionalMultiphaseHybridFVM >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseReservoirFVM:
    {
      setStrategy< CompositionalMultiphaseReservoirFVM >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseReservoirHybridFVM:
    {
      setStrategy< CompositionalMultiphaseReservoirHybridFVM >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::thermalCompositionalMultiphaseFVM:
    {
      setStrategy< ThermalCompositionalMultiphaseFVM >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::hybridSinglePhasePoromechanics:
    {
      setStrategy< HybridSinglePhasePoromechanics >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::hydrofracture:
    {
      setStrategy< Hydrofracture >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::lagrangianContactMechanics:
    {
      setStrategy< LagrangianContactMechanics >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::multiphasePoromechanics:
    {
      setStrategy< MultiphasePoromechanics >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::multiphasePoromechanicsReservoirFVM:
    {
      setStrategy< MultiphasePoromechanicsReservoirFVM >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::singlePhaseHybridFVM:
    {
      setStrategy< SinglePhaseHybridFVM >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanics:
    {
      setStrategy< SinglePhasePoromechanics >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanicsReservoirFVM:
    {
      setStrategy< SinglePhasePoromechanicsReservoirFVM >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::singlePhaseReservoirFVM:
    {
      setStrategy< SinglePhaseReservoirFVM >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    case LinearSolverParameters::MGR::StrategyType::singlePhaseReservoirHybridFVM:
    {
      setStrategy< SinglePhaseReservoirHybridFVM >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    default:
    {
      GEOSX_ERROR( "Unsupported MGR strategy: " << params.mgr.strategy );
    }
  }

  GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseSolver( precond.ptr,
                                                   mgrData.coarseSolver.solve,
                                                   mgrData.coarseSolver.setup,
                                                   mgrData.coarseSolver.ptr ) );
  precond.setup = HYPRE_MGRSetup;
  precond.solve = HYPRE_MGRSolve;
  precond.destroy = HYPRE_MGRDestroy;

  // Set custom F-solver based on SDC for mechanics case
  // Requirement: displacement degrees of freedom are the first being eliminated,
  //              i.e. they are F-points for the first MGR level
  //if( params.preconditionerType == LinearSolverParameters::PreconditionerType::mgr && params.mgr.separateComponents )
  if( params.preconditionerType == LinearSolverParameters::PreconditionerType::mgr )
  {
    HYPRE_Int logLevel = (params.logLevel == 2 || params.logLevel >= 4) ? 1 : 0;

    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &mgrData.mechSolver.ptr ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( mgrData.mechSolver.ptr, 0.0 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( mgrData.mechSolver.ptr, 1 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxRowSum( mgrData.mechSolver.ptr, 1.0 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetStrongThreshold( mgrData.mechSolver.ptr, 0.8 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( mgrData.mechSolver.ptr, logLevel ) );
#ifdef GEOSX_USE_HYPRE_CUDA
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCoarsenType( mgrData.mechSolver.ptr, hypre::getAMGCoarseningType( "PMIS" ) ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( mgrData.mechSolver.ptr, hypre::getAMGRelaxationType( LinearSolverParameters::AMG::SmootherType::l1jacobi ) ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumSweeps( mgrData.mechSolver.ptr, 2 ) );
#else
    // Note: commenting out the line below, and uncommenting what follows seems to consistently deteriorate
    // the behavior of the linear solver (number of iterations and wall-clock time). This requires further investigation.
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxOrder( mgrData.mechSolver.ptr, 1 ) );
    // GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( mgrData.mechSolver.ptr, hypre::getAMGRelaxationType(
    // LinearSolverParameters::AMG::SmootherType::l1sgs ) ) );
    // GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumSweeps( mgrData.mechSolver.ptr, 1 ) );
    // GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxWt( mgrData.mechSolver.ptr, 0.8 ) );
#endif
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumFunctions( mgrData.mechSolver.ptr, 3 ) );

    mgrData.mechSolver.setup = HYPRE_BoomerAMGSetup;
    mgrData.mechSolver.solve = HYPRE_BoomerAMGSolve;
    mgrData.mechSolver.destroy = HYPRE_BoomerAMGDestroy;

    // Ignore the setup function here, since we'll be performing it manually in setupSeparateComponent()
    HYPRE_MGRSetFSolver( precond.ptr, mgrData.mechSolver.solve, mgrData.mechSolver.setup, mgrData.mechSolver.ptr );
  }
}

} // namespace geosx
