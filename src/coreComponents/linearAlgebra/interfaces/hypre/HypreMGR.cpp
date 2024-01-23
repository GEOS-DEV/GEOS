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


#include "linearAlgebra/interfaces/hypre/mgrStrategies/CompositionalMultiphaseFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/CompositionalMultiphaseHybridFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/CompositionalMultiphaseReservoirFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/CompositionalMultiphaseReservoirHybridFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/HybridSinglePhasePoromechanics.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/Hydrofracture.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/LagrangianContactMechanics.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/MultiphasePoromechanics.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/MultiphasePoromechanicsReservoirFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/ReactiveCompositionalMultiphaseOBL.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhaseHybridFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhasePoromechanics.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhasePoromechanicsEmbeddedFractures.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhasePoromechanicsReservoirFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhaseReservoirFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhaseReservoirHybridFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/ThermalCompositionalMultiphaseFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/ThermalSinglePhasePoromechanics.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/ThermalMultiphasePoromechanics.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SolidMechanicsEmbeddedFractures.hpp"

#include "LvArray/src/output.hpp"

namespace geos
{

void hypre::mgr::createMGR( LinearSolverParameters const & params,
                            DofManager const * const dofManager,
                            HyprePrecWrapper & precond,
                            HypreMGRData & mgrData )
{
  GEOS_ERROR_IF( dofManager == nullptr, "MGR preconditioner requires a DofManager instance" );

  GEOS_LAI_CHECK_ERROR( HYPRE_MGRCreate( &precond.ptr ) );

  // Set MGR's print level
  HYPRE_Int logLevel = LvArray::integerConversion< HYPRE_Int >( params.logLevel );

  // Update logLevel (1 is reserved in GEOS to log execution times)
  logLevel = (logLevel > 0) ? ((logLevel - 1) & ~0x2) : 0;

  // Hypre's parameters to use MGR as a preconditioner
  GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetTol( precond.ptr, 0.0 ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetMaxIter( precond.ptr, 1 ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetPrintLevel( precond.ptr, logLevel ) );

  array1d< int > const numComponentsPerField = dofManager->numComponentsPerField();
  dofManager->getLocalDofComponentLabels( mgrData.pointMarkers );

  if( params.logLevel >= 1 )
  {
    GEOS_LOG_RANK_0( GEOS_FMT( "        MGR preconditioner: numComponentsPerField = {}", numComponentsPerField ) );
  }
  if( params.logLevel >= 1024 )
  {
    GEOS_LOG_RANK( GEOS_FMT( "        MGR preconditioner: pointMarkers = {}", mgrData.pointMarkers ) );
  }

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
    case LinearSolverParameters::MGR::StrategyType::reactiveCompositionalMultiphaseOBL:
    {
      setStrategy< ReactiveCompositionalMultiphaseOBL >( params.mgr, numComponentsPerField, precond, mgrData );
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
    case LinearSolverParameters::MGR::StrategyType::thermalSinglePhasePoromechanics:
    {
      setStrategy< ThermalSinglePhasePoromechanics >( params.mgr, numComponentsPerField, precond, mgrData );
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
    case LinearSolverParameters::MGR::StrategyType::thermalMultiphasePoromechanics:
    {
      setStrategy< ThermalMultiphasePoromechanics >( params.mgr, numComponentsPerField, precond, mgrData );
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
    case LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanicsEmbeddedFractures:
    {
      setStrategy< SinglePhasePoromechanicsEmbeddedFractures >( params.mgr, numComponentsPerField, precond, mgrData );
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
    case LinearSolverParameters::MGR::StrategyType::solidMechanicsEmbeddedFractures:
    {
      setStrategy< SolidMechanicsEmbeddedFractures >( params.mgr, numComponentsPerField, precond, mgrData );
      break;
    }
    default:
    {
      GEOS_ERROR( "Unsupported MGR strategy: " << params.mgr.strategy );
    }
  }

  GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseSolver( precond.ptr,
                                                  mgrData.coarseSolver.solve,
                                                  mgrData.coarseSolver.setup,
                                                  mgrData.coarseSolver.ptr ) );
  precond.setup = HYPRE_MGRSetup;
  precond.solve = HYPRE_MGRSolve;
  precond.destroy = HYPRE_MGRDestroy;
}

} // namespace geos
