/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file LinearSolverParameters.cpp
 */

#include "LinearSolverParameters.hpp"

namespace geosx
{
using namespace dataRepository;

LinearSolverParametersGroup::LinearSolverParametersGroup( std::string const & name,
                                                          Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
  enableLogLevelInput();

  // note: default parameter values preset by base class

  registerWrapper( viewKeysStruct::solverTypeString, &solverType )->
    setApplyDefaultValue( solverType )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Linear solver type" );

  registerWrapper( viewKeysStruct::preconditionerTypeString, &preconditionerType )->
    setApplyDefaultValue( preconditionerType )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Preconditioner type" );

  registerWrapper( viewKeysStruct::krylovMaxIterString, &krylov.maxIterations )->
    setApplyDefaultValue( krylov.maxIterations )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Maximum iterations allowed" );

  registerWrapper( viewKeysStruct::krylovTolString, &krylov.relTolerance )->
    setApplyDefaultValue( krylov.relTolerance )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Relative convergence tolerance" );

  registerWrapper( viewKeysStruct::krylovAdaptiveTolString, &krylov.useAdaptiveTol )->
    setApplyDefaultValue( krylov.useAdaptiveTol )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Use Eisenstat-Walker adaptive linear tolerance" );

  registerWrapper( viewKeysStruct::krylovWeakTolString, &krylov.weakestTol )->
    setApplyDefaultValue( krylov.weakestTol )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Weakest-allowed tolerance for adaptive method" );

  registerWrapper( viewKeysStruct::amgNumSweepsString, &amg.numSweeps )->
    setApplyDefaultValue( amg.numSweeps )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "AMG smoother sweeps" );

  registerWrapper( viewKeysStruct::amgSmootherString, &amg.smootherType )->
    setApplyDefaultValue( amg.smootherType )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "AMG smoother type" );

  registerWrapper( viewKeysStruct::amgCoarseString, &amg.coarseType )->
    setApplyDefaultValue( amg.coarseType )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "AMG coarsest level solver/smoother type" );

  registerWrapper( viewKeysStruct::amgThresholdString, &amg.threshold )->
    setApplyDefaultValue( amg.threshold )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "AMG strength-of-connection threshold" );

  registerWrapper( viewKeysStruct::iluFillString, &ilu.fill )->
    setApplyDefaultValue( ilu.fill )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "ILU(K) fill factor" );

  registerWrapper( viewKeysStruct::iluThresholdString, &ilu.threshold )->
    setApplyDefaultValue( ilu.threshold )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "ILU(T) threshold factor" );
}

void LinearSolverParametersGroup::PostProcessInput()
{
  logLevel = getLogLevel();
}

REGISTER_CATALOG_ENTRY( Group, LinearSolverParametersGroup, std::string const &, Group * const )

} /* namespace geosx */
