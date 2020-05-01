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
 * @file SystemSolverParameters.cpp
 *
 */

#include "SystemSolverParameters.hpp"

namespace geosx
{
using namespace dataRepository;
SystemSolverParameters::SystemSolverParameters( std::string const & name,
                                                Group * const parent ):
  Group( name, parent ),
  m_solverType( "Klu" ),
  m_krylovTol(),
  m_numKrylovIter(),
  m_kspace(),
  m_ilut_fill( 3.0 ),
  m_ilut_drop(),
  m_useMLPrecond(),
  m_useInnerSolver(),
  m_scalingOption(),
  m_useBicgstab(),
  m_useDirectSolver()
{
  setInputFlags( InputFlags::OPTIONAL );

  // This enables logLevel filtering
  enableLogLevelInput();

  registerWrapper( viewKeysStruct::solverTypeString, &m_solverType )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "" );

  registerWrapper( viewKeysStruct::krylovTolString, &m_krylovTol )->
    setApplyDefaultValue( 1.0e-6 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Desired tolerance for Krylov solve" );

  registerWrapper( viewKeysStruct::useAdaptiveKrylovString, &m_useAdaptiveKrylovTol )->
    setApplyDefaultValue( 0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Enable Eisenstat-Walker adaptive Krylov tolerance" );

  registerWrapper( viewKeysStruct::numKrylovIterString, &m_numKrylovIter )->
    setApplyDefaultValue( 100 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Maximum number of Krylov Iterations" );

  registerWrapper( viewKeysStruct::kspaceString, &m_kspace )->
    setApplyDefaultValue( 300 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "" );

  registerWrapper( viewKeysStruct::ilut_fillString, &m_ilut_fill )->
    setApplyDefaultValue( 3.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "" );

  registerWrapper( viewKeysStruct::ilut_dropString, &m_ilut_drop )->
    setApplyDefaultValue( 0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "" );

  registerWrapper( viewKeysStruct::useMLPrecondString, &m_useMLPrecond )->
    setApplyDefaultValue( 0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "" );

  registerWrapper( viewKeysStruct::useInnerSolverString, &m_useInnerSolver )->
    setApplyDefaultValue( 0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "" );

  registerWrapper( viewKeysStruct::scalingOptionString, &m_scalingOption )->
    setApplyDefaultValue( 0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "" );

  registerWrapper( viewKeysStruct::useBicgstabString, &m_useBicgstab )->
    setApplyDefaultValue( 0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "" );

  registerWrapper( viewKeysStruct::useDirectSolverString, &m_useDirectSolver )->
    setApplyDefaultValue( 0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "" );

  registerWrapper( viewKeysStruct::krylovResidualInitString, &m_krylovResidualInit )->
    setApplyDefaultValue( 0 )->
    setDescription( "Initial Krylov solver residual." );

  registerWrapper( viewKeysStruct::krylovResidualFinalString, &m_krylovResidualFinal )->
    setApplyDefaultValue( 0 )->
    setDescription( "Final Krylov solver residual." );

}

void SystemSolverParameters::PostProcessInput()
{}

REGISTER_CATALOG_ENTRY( Group, SystemSolverParameters, std::string const &, Group * const )

} /* namespace geosx */
