/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TwoPhaseImmiscibleFluid.cpp
 */

#include "constitutive/fluid/multifluid/CO2Brine/functions/PVTFunctionHelpers.hpp" // for readTable
#include "TwoPhaseImmiscibleFluid.hpp"
#include "TwoPhaseImmiscibleFluidFields.hpp"

#include "functions/FunctionManager.hpp"


namespace geos
{

using namespace dataRepository;

namespace constitutive
{


TwoPhaseImmiscibleFluid::TwoPhaseImmiscibleFluid( string const & name, Group * const parent )
  : ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeyStruct::phaseNamesString(), &m_phaseNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "List of fluid phases" );

  // 1) First option: specify PVT tables from one file per phase, read the files line by line, and populate the internal TableFunctions
  registerWrapper( viewKeyStruct::tableFilesString(), &m_tableFiles ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "List of filenames with input PVT tables (one per phase)" );

  // 2) Second option: specify TableFunction names for each phase,
  registerWrapper( viewKeyStruct::densityTableNamesString(), &m_densityTableNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of density TableFuncion names from the Function block. \n"
                    "The user must provide one TableFunction per phase, respecting the order provided in \"phaseNames\"." );

  registerWrapper( viewKeyStruct::viscosityTableNamesString(), &m_viscosityTableNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of viscosity TableFuncion names from the Function block. \n"
                    "The user must provide one TableFunction per phase, respecting the order provided in \"phaseNames\"." );

  registerField< fields::twophaseimmisciblefluid::phaseDensity >( &m_phaseDensity.value );
  registerField< fields::twophaseimmisciblefluid::dPhaseDensity >( &m_phaseDensity.derivs );
  registerField< fields::twophaseimmisciblefluid::phaseDensity_n >( &m_phaseDensity_n );

  registerField< fields::twophaseimmisciblefluid::phaseViscosity >( &m_phaseViscosity.value );
  registerField< fields::twophaseimmisciblefluid::dPhaseViscosity >( &m_phaseViscosity.derivs );
}


void TwoPhaseImmiscibleFluid::allocateConstitutiveData( Group & parent,
                                                        localIndex const numPts )
{
  // Assume sole dependency on pressure, i.e. one derivative
  m_phaseDensity.value.resize( 0, numPts, 2 );
  m_phaseDensity.derivs.resize( 0, numPts, 2, 1 );

  m_phaseDensity_n.resize( 0, numPts, 2 );

  m_phaseViscosity.value.resize( 0, numPts, 2 );
  m_phaseViscosity.derivs.resize( 0, numPts, 2, 1 );

  ConstitutiveBase::allocateConstitutiveData( parent, numPts );
}


void TwoPhaseImmiscibleFluid::postInputInitialization()
{
  ConstitutiveBase::postInputInitialization();

  // Ensure that we have two phases defined
  GEOS_THROW_IF_NE_MSG( m_phaseNames.size(), 2,
                        GEOS_FMT( "invalid number of phases. There should be 2 phases defined by {}",
                                  viewKeyStruct::phaseNamesString() ),
                        InputError, getDataContext() );

  // Input relationships can be provided either as text files or TableFunctions.
  m_tableFiles.empty() ? readInputDataFromTableFunctions() : readInputDataFromFileTableFunctions();

  checkTableConsistency();
}


void TwoPhaseImmiscibleFluid::fillData( integer const ip,
                                        array1d< array1d< real64 > > const & tableValues )
{
  array1d< array1d< real64 > > pressureCoords( 1 );
  pressureCoords[0].resize( tableValues.size() );
  array1d< real64 > density( tableValues.size() );
  array1d< real64 > viscosity( tableValues.size() );

  for( localIndex i = 0; i < tableValues.size(); ++i )
  {
    GEOS_THROW_IF_NE_MSG( tableValues[i].size(), 3,
                          "three columns (pressure, density, and viscosity) are expected",
                          InputError, getDataContext() );

    pressureCoords[0][i] = tableValues[i][0];
    density[i] = tableValues[i][1];
    viscosity[i] = tableValues[i][2];
  }

  string const densityTableName = getName() + "DensityPhase" + GEOS_FMT( "{}", ip );
  string const viscosityTableName = getName() +  "ViscosityPhase" + GEOS_FMT( "{}", ip );
  m_densityTableNames.emplace_back( densityTableName );
  m_viscosityTableNames.emplace_back( viscosityTableName );

  FunctionManager & functionManager = FunctionManager::getInstance();

  TableFunction & tableDensity =
    dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", densityTableName ) );
  tableDensity.setTableCoordinates( pressureCoords, { units::Pressure } );
  tableDensity.setTableValues( density, units::Density );
  tableDensity.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  TableFunction & tableViscosity =
    dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", viscosityTableName ) );
  tableViscosity.setTableCoordinates( pressureCoords, { units::Pressure } );
  tableViscosity.setTableValues( viscosity, units::Viscosity );
  tableViscosity.setInterpolationMethod( TableFunction::InterpolationType::Linear );
}


void TwoPhaseImmiscibleFluid::readInputDataFromFileTableFunctions()
{
  // Check for ambiguous definition
  GEOS_THROW_IF( !(m_densityTableNames.empty() && m_viscosityTableNames.empty()),
                 "input is redundant (both TableFunction names and text files)",
                 InputError, getDataContext() );


  // Check that we have exactly two table files (one per phase)
  GEOS_THROW_IF_NE_MSG( m_tableFiles.size(), 2,
                        "expecting two table files (one per phase)",
                        InputError, getDataContext() );

  array1d< array1d< real64 > > tableValues;
  for( integer ip = 0; ip < 2; ++ip )
  {
    tableValues.clear();
    geos::constitutive::PVTProps::BlackOilTables::readTable( m_tableFiles[ip], 3, tableValues );
    fillData( ip, tableValues );
  }
}


void TwoPhaseImmiscibleFluid::readInputDataFromTableFunctions()
{
  // Check for ambiguous definition
  GEOS_THROW_IF( !m_tableFiles.empty(),
                 "input is redundant (both TableFunction names and text files)",
                 InputError, getDataContext() );

  // Since we are considering a two phase fluid, we should have exactly 2 tables per property
  GEOS_THROW_IF_NE_MSG( m_densityTableNames.size(), 2,
                        "one density table must be provided for each phase",
                        InputError, getDataContext() );

  GEOS_THROW_IF_NE_MSG( m_viscosityTableNames.size(), 2,
                        "one viscosity table must be provided for each phase",
                        InputError, getDataContext() );


  FunctionManager const & functionManager = FunctionManager::getInstance();

  for( integer iph = 0; iph < 2; ++iph )
  {
    GEOS_THROW_IF( !functionManager.hasGroup( m_densityTableNames[iph] ),
                   GEOS_FMT( "density table '{}' not found", m_densityTableNames[iph] ),
                   InputError, getDataContext() );

    GEOS_THROW_IF( !functionManager.hasGroup( m_viscosityTableNames[iph] ),
                   GEOS_FMT( "viscosity table '{}' not found", m_viscosityTableNames[iph] ),
                   InputError, getDataContext() );
  }
}


void TwoPhaseImmiscibleFluid::initializePostSubGroups()
{
  ConstitutiveBase::initializePostSubGroups();

  FunctionManager const & functionManager = FunctionManager::getInstance();
  for( integer iph = 0; iph < 2; ++iph )
  {
    // Grab the tables by name from the function manager,
    // then add them in a list to create their table wrappers when needed
    TableFunction const & densityTable = functionManager.getGroup< TableFunction const >( m_densityTableNames[iph] );
    m_densityTables.emplace_back( &densityTable );
    m_densityTableKernels.emplace_back( m_densityTables[iph]->createKernelWrapper() );

    TableFunction const & viscosityTable = functionManager.getGroup< TableFunction const >( m_viscosityTableNames[iph] );
    m_viscosityTables.emplace_back( &viscosityTable );
    m_viscosityTableKernels.emplace_back( m_viscosityTables[iph]->createKernelWrapper() );
  }
}


void TwoPhaseImmiscibleFluid::checkTableConsistency() const
{
  FunctionManager const & functionManager = FunctionManager::getInstance();
  for( integer iph = 0; iph < 2; ++iph )
  {
    TableFunction const & densityTable = functionManager.getGroup< TableFunction const >( m_densityTableNames[iph] );
    arrayView1d< real64 const > const density = densityTable.getValues();

    for( localIndex i = 1; i < density.size(); ++i )
    {
      GEOS_THROW_IF( density[i] - density[i-1] < 0,
                     GEOS_FMT( "in table '{}' density values must be increasing", densityTable.getName() ),
                     InputError, getDataContext() );
    }
  }
}


TwoPhaseImmiscibleFluid::KernelWrapper
TwoPhaseImmiscibleFluid::createKernelWrapper()
{
  return KernelWrapper( m_densityTableKernels,
                        m_viscosityTableKernels,
                        m_phaseDensity.toView(),
                        m_phaseViscosity.toView());
}


TwoPhaseImmiscibleFluid::KernelWrapper::KernelWrapper(
  arrayView1d< TableFunction::KernelWrapper const > densityTables,
  arrayView1d< TableFunction::KernelWrapper const > viscosityTables,
  PhaseProp::ViewType phaseDensity,
  PhaseProp::ViewType phaseViscosity )
  : m_densityTables( std::move( densityTables )),
  m_viscosityTables( std::move( viscosityTables )),
  m_phaseDensity( std::move( phaseDensity )),
  m_phaseViscosity( std::move( phaseViscosity )) {}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, TwoPhaseImmiscibleFluid, string const &, Group * const )

}  // namespace constitutive
}  // namespace geos
