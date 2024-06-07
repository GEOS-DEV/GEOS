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
 * @file LinearSolverParameters.cpp
 */

#include "LinearSolverParameters.hpp"

namespace geos
{
using namespace dataRepository;

LinearSolverParametersInput::LinearSolverParametersInput( string const & name,
                                                          Group * const parent )
  :
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
  enableLogLevelInput();

  registerWrapper( viewKeyStruct::solverTypeString(), &m_parameters.solverType ).
    setApplyDefaultValue( m_parameters.solverType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Linear solver type. Available options are: "
                    "``" + EnumStrings< LinearSolverParameters::SolverType >::concat( "|" ) + "``" );

  registerWrapper( viewKeyStruct::preconditionerTypeString(), &m_parameters.preconditionerType ).
    setApplyDefaultValue( m_parameters.preconditionerType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Preconditioner type. Available options are: "
                    "``" + EnumStrings< LinearSolverParameters::PreconditionerType >::concat( "|" ) + "``" );

  registerWrapper( viewKeyStruct::stopIfErrorString(), &m_parameters.stopIfError ).
    setApplyDefaultValue( m_parameters.stopIfError ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Whether to stop the simulation if the linear solver reports an error" );

  registerWrapper( viewKeyStruct::directCheckResidualString(), &m_parameters.direct.checkResidual ).
    setApplyDefaultValue( m_parameters.direct.checkResidual ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Whether to check the linear system solution residual" );

  registerWrapper( viewKeyStruct::directEquilString(), &m_parameters.direct.equilibrate ).
    setApplyDefaultValue( m_parameters.direct.equilibrate ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Whether to scale the rows and columns of the matrix" );

  registerWrapper( viewKeyStruct::directColPermString(), &m_parameters.direct.colPerm ).
    setApplyDefaultValue( m_parameters.direct.colPerm ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "How to permute the columns. Available options are: "
                    "``" + EnumStrings< LinearSolverParameters::Direct::ColPerm >::concat( "|" ) + "``" );

  registerWrapper( viewKeyStruct::directRowPermString(), &m_parameters.direct.rowPerm ).
    setApplyDefaultValue( m_parameters.direct.rowPerm ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "How to permute the rows. Available options are: "
                    "``" + EnumStrings< LinearSolverParameters::Direct::RowPerm >::concat( "|" ) + "``" );

  registerWrapper( viewKeyStruct::directReplTinyPivotString(), &m_parameters.direct.replaceTinyPivot ).
    setApplyDefaultValue( m_parameters.direct.replaceTinyPivot ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Whether to replace tiny pivots by sqrt(epsilon)*norm(A)" );

  registerWrapper( viewKeyStruct::directIterRefString(), &m_parameters.direct.iterativeRefine ).
    setApplyDefaultValue( m_parameters.direct.iterativeRefine ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Whether to perform iterative refinement" );

  registerWrapper( viewKeyStruct::directParallelString(), &m_parameters.direct.parallel ).
    setApplyDefaultValue( m_parameters.direct.parallel ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Whether to use a parallel solver (instead of a serial one)" );

  registerWrapper( viewKeyStruct::krylovMaxIterString(), &m_parameters.krylov.maxIterations ).
    setApplyDefaultValue( m_parameters.krylov.maxIterations ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum iterations allowed for an iterative solver" );

  registerWrapper( viewKeyStruct::krylovMaxRestartString(), &m_parameters.krylov.maxRestart ).
    setApplyDefaultValue( m_parameters.krylov.maxRestart ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum iterations before restart (GMRES only)" );

  registerWrapper( viewKeyStruct::krylovTolString(), &m_parameters.krylov.relTolerance ).
    setApplyDefaultValue( m_parameters.krylov.relTolerance ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Relative convergence tolerance of the iterative method\n"
                    "If the method converges, the iterative solution :math:`\\mathsf{x}_k` is such that\n"
                    "the relative residual norm satisfies:\n"
                    ":math:`\\left\\lVert \\mathsf{b} - \\mathsf{A} \\mathsf{x}_k \\right\\rVert_2` < ``" +
                    string( viewKeyStruct::krylovTolString() ) + "`` * :math:`\\left\\lVert\\mathsf{b}\\right\\rVert_2`" );

  registerWrapper( viewKeyStruct::krylovAdaptiveTolString(), &m_parameters.krylov.useAdaptiveTol ).
    setApplyDefaultValue( m_parameters.krylov.useAdaptiveTol ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Use Eisenstat-Walker adaptive linear tolerance" );

  registerWrapper( viewKeyStruct::krylovWeakTolString(), &m_parameters.krylov.weakestTol ).
    setApplyDefaultValue( m_parameters.krylov.weakestTol ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Weakest-allowed tolerance for adaptive method" );

  registerWrapper( viewKeyStruct::krylovStrongTolString(), &m_parameters.krylov.strongestTol ).
    setApplyDefaultValue( m_parameters.krylov.strongestTol ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Strongest-allowed tolerance for adaptive method" );

  registerWrapper( viewKeyStruct::adaptiveGammaString(), &m_parameters.krylov.adaptiveGamma ).
    setApplyDefaultValue( m_parameters.krylov.adaptiveGamma ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Gamma parameter for adaptive method" );
    
  registerWrapper( viewKeyStruct::adaptiveExponentString(), &m_parameters.krylov.adaptiveExponent ).
    setApplyDefaultValue( m_parameters.krylov.adaptiveExponent ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Exponent parameter for adaptive method" );

  registerWrapper( viewKeyStruct::amgNumSweepsString(), &m_parameters.amg.numSweeps ).
    setApplyDefaultValue( m_parameters.amg.numSweeps ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "AMG smoother sweeps" );

  registerWrapper( viewKeyStruct::amgSmootherString(), &m_parameters.amg.smootherType ).
    setApplyDefaultValue( m_parameters.amg.smootherType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "AMG smoother type. Available options are: "
                    "``" + EnumStrings< LinearSolverParameters::AMG::SmootherType >::concat( "|" ) + "``" );

  registerWrapper( viewKeyStruct::amgRelaxWeight(), &m_parameters.amg.relaxWeight ).
    setApplyDefaultValue( m_parameters.amg.relaxWeight ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "AMG relaxation factor for the smoother" );

  registerWrapper( viewKeyStruct::amgCoarseString(), &m_parameters.amg.coarseType ).
    setApplyDefaultValue( m_parameters.amg.coarseType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "AMG coarsest level solver/smoother type. Available options are: "
                    "``" + EnumStrings< LinearSolverParameters::AMG::CoarseType >::concat( "|" ) + "``" );

  registerWrapper( viewKeyStruct::amgCoarseningString(), &m_parameters.amg.coarseningType ).
    setApplyDefaultValue( m_parameters.amg.coarseningType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "AMG coarsening algorithm. Available options are: "
                    "``" + EnumStrings< LinearSolverParameters::AMG::CoarseningType >::concat( "|" ) + "``" );

  registerWrapper( viewKeyStruct::amgInterpolationString(), &m_parameters.amg.interpolationType ).
    setApplyDefaultValue( m_parameters.amg.interpolationType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "AMG interpolation algorithm. Available options are: "
                    "``" + EnumStrings< LinearSolverParameters::AMG::InterpType >::concat( "|" ) + "``" );

  registerWrapper( viewKeyStruct::amgInterpMaxNonZerosString(), &m_parameters.amg.interpolationMaxNonZeros ).
    setApplyDefaultValue( m_parameters.amg.interpolationMaxNonZeros ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "AMG interpolation maximum number of nonzeros per row" );

  registerWrapper( viewKeyStruct::amgNumFunctionsString(), &m_parameters.amg.numFunctions ).
    setApplyDefaultValue( m_parameters.amg.numFunctions ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "AMG number of functions" );

  registerWrapper( viewKeyStruct::amgAggressiveNumPathsString(), &m_parameters.amg.aggressiveNumPaths ).
    setApplyDefaultValue( m_parameters.amg.aggressiveNumPaths ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "AMG number of paths for aggressive coarsening" );

  registerWrapper( viewKeyStruct::amgAggressiveNumLevelsString(), &m_parameters.amg.aggressiveNumLevels ).
    setApplyDefaultValue( m_parameters.amg.aggressiveNumLevels ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "AMG number of levels for aggressive coarsening" );

  registerWrapper( viewKeyStruct::amgAggressiveInterpTypeString(), &m_parameters.amg.aggressiveInterpType ).
    setApplyDefaultValue( m_parameters.amg.aggressiveInterpType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "AMG aggressive interpolation algorithm. Available options are: "
                    "``" + EnumStrings< LinearSolverParameters::AMG::AggInterpType >::concat( "|" ) + "``" );

  registerWrapper( viewKeyStruct::amgThresholdString(), &m_parameters.amg.threshold ).
    setApplyDefaultValue( m_parameters.amg.threshold ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "AMG strength-of-connection threshold" );

  registerWrapper( viewKeyStruct::amgSeparateComponentsString(), &m_parameters.amg.separateComponents ).
    setApplyDefaultValue( m_parameters.amg.separateComponents ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "AMG apply separate component filter for multi-variable problems" );

  registerWrapper( viewKeyStruct::amgNullSpaceTypeString(), &m_parameters.amg.nullSpaceType ).
    setApplyDefaultValue( m_parameters.amg.nullSpaceType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "AMG near null space approximation. Available options are:"
                    "``" + EnumStrings< LinearSolverParameters::AMG::NullSpaceType >::concat( "|" ) + "``" );

  registerWrapper( viewKeyStruct::iluFillString(), &m_parameters.ifact.fill ).
    setApplyDefaultValue( m_parameters.ifact.fill ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "ILU(K) fill factor" );

  registerWrapper( viewKeyStruct::iluThresholdString(), &m_parameters.ifact.threshold ).
    setApplyDefaultValue( m_parameters.ifact.threshold ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "ILU(T) threshold factor" );
}

void LinearSolverParametersInput::postProcessInput()
{
  m_parameters.logLevel = getLogLevel();

  static const std::set< integer > binaryOptions = { 0, 1 };

  GEOS_ERROR_IF( binaryOptions.count( m_parameters.stopIfError ) == 0,
                 getWrapperDataContext( viewKeyStruct::stopIfErrorString() ) <<
                 ": option can be either 0 (false) or 1 (true)" );
  GEOS_ERROR_IF( binaryOptions.count( m_parameters.direct.checkResidual ) == 0,
                 getWrapperDataContext( viewKeyStruct::directCheckResidualString() ) <<
                 ": option can be either 0 (false) or 1 (true)" );
  GEOS_ERROR_IF( binaryOptions.count( m_parameters.direct.equilibrate ) == 0,
                 getWrapperDataContext( viewKeyStruct::directEquilString() ) <<
                 ": option can be either 0 (false) or 1 (true)" );
  GEOS_ERROR_IF( binaryOptions.count( m_parameters.direct.replaceTinyPivot ) == 0,
                 getWrapperDataContext( viewKeyStruct::directReplTinyPivotString() ) <<
                 ": option can be either 0 (false) or 1 (true)" );
  GEOS_ERROR_IF( binaryOptions.count( m_parameters.direct.iterativeRefine ) == 0,
                 getWrapperDataContext( viewKeyStruct::directIterRefString() ) <<
                 ": option can be either 0 (false) or 1 (true)" );
  GEOS_ERROR_IF( binaryOptions.count( m_parameters.direct.parallel ) == 0,
                 getWrapperDataContext( viewKeyStruct::directParallelString() ) <<
                 ": option can be either 0 (false) or 1 (true)" );

  GEOS_ERROR_IF_LT_MSG( m_parameters.krylov.maxIterations, 0,
                        getWrapperDataContext( viewKeyStruct::krylovMaxIterString() ) <<
                        ": Invalid value." );
  GEOS_ERROR_IF_LT_MSG( m_parameters.krylov.maxRestart, 0,
                        getWrapperDataContext( viewKeyStruct::krylovMaxRestartString() ) <<
                        ": Invalid value." );

  GEOS_ERROR_IF_LT_MSG( m_parameters.krylov.relTolerance, 0.0,
                        getWrapperDataContext( viewKeyStruct::krylovTolString() ) <<
                        ": Invalid value." );
  GEOS_ERROR_IF_GT_MSG( m_parameters.krylov.relTolerance, 1.0,
                        getWrapperDataContext( viewKeyStruct::krylovTolString() ) <<
                        ": Invalid value." );

  GEOS_ERROR_IF_LT_MSG( m_parameters.ifact.fill, 0,
                        getWrapperDataContext( viewKeyStruct::iluFillString() ) <<
                        ": Invalid value." );
  GEOS_ERROR_IF_LT_MSG( m_parameters.ifact.threshold, 0.0,
                        getWrapperDataContext( viewKeyStruct::iluThresholdString() ) <<
                        ": Invalid value." );

  GEOS_ERROR_IF_LT_MSG( m_parameters.amg.numSweeps, 0,
                        getWrapperDataContext( viewKeyStruct::amgNumSweepsString() ) <<
                        ": Invalid value." );
  GEOS_ERROR_IF_LT_MSG( m_parameters.amg.threshold, 0.0,
                        getWrapperDataContext( viewKeyStruct::amgThresholdString() ) <<
                        ": Invalid value." );
  GEOS_ERROR_IF_GT_MSG( m_parameters.amg.threshold, 1.0,
                        getWrapperDataContext( viewKeyStruct::amgThresholdString() ) <<
                        ": Invalid value." );

  // TODO input validation for other AMG parameters ?
}

REGISTER_CATALOG_ENTRY( Group, LinearSolverParametersInput, string const &, Group * const )

} // namespace geos
