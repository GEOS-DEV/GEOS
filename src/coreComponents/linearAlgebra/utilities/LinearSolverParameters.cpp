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

  // note: default values currently set by base class

  registerWrapper( viewKeysStruct::solverTypeString, &solverType )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Linear solver type" );

  registerWrapper( viewKeysStruct::preconditionerTypeString, &preconditionerType )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Preconditioner type" );

  registerWrapper( viewKeysStruct::krylovTolString, &krylov.tolerance )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Relative convergence tolerance" );

  registerWrapper( viewKeysStruct::krylovAdaptiveTolString, &krylov.useAdaptiveTol )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Use Eisenstat-Walker adaptive linear tolerance" );

  registerWrapper( viewKeysStruct::krylovMaxIterString, &krylov.maxIterations )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Maximum iterations allowed" );

  registerWrapper( viewKeysStruct::amgNumSweepsString, &amg.numSweeps )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "AMG smoother sweeps" );

  registerWrapper( viewKeysStruct::amgSmootherString, &amg.smootherType )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "AMG smoother type" );

  registerWrapper( viewKeysStruct::amgCoarseString, &amg.coarseType )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "AMG coarsest level solver/smoother type" );

  registerWrapper( viewKeysStruct::amgAggregationString, &amg.aggregationThreshold )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "AMG aggregation threshold" );

  registerWrapper( viewKeysStruct::iluFillString, &ilu.fill )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "ILU(K) fill factor" );

  registerWrapper( viewKeysStruct::iluThresholdString, &ilu.threshold )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "ILU(T) threshold factor" );
}

void LinearSolverParametersGroup::PostProcessInput()
{}

REGISTER_CATALOG_ENTRY( Group, LinearSolverParametersGroup, std::string const &, Group * const )

} /* namespace geosx */
