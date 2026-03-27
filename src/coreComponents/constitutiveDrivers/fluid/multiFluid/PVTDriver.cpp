/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PVTDriver.cpp
 */

#include "PVTDriver.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutiveDrivers/LogLevelsInfo.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidSelector.hpp"
#include "functions/FunctionManager.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

PVTDriver::PVTDriver( const string & name,
                      Group * const parent ):
  ConstitutiveDriver( name, parent )
{

  registerWrapper( viewKeyStruct::fluidNameString(), &m_fluidName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Fluid to test" );

  registerWrapper( viewKeyStruct::feedString(), &m_feed ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Feed composition array [mol fraction]" );

  registerWrapper( viewKeyStruct::pressureFunctionString(), &m_pressureFunctionName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Function controlling pressure time history" );

  registerWrapper( viewKeyStruct::temperatureFunctionString(), &m_temperatureFunctionName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Function controlling temperature time history" );

  registerWrapper( viewKeyStruct::outputMassDensityString(), &m_outputMassDensity ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag to indicate that the mass density of each phase should be output" );

  registerWrapper( viewKeyStruct::outputCompressibilityString(), &m_outputCompressibility ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag to indicate that the total compressibility should be output" );

  registerWrapper( viewKeyStruct::outputPhaseCompositionString(), &m_outputPhaseComposition ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag to indicate that phase compositions should be output" );
}

void PVTDriver::postInputInitialization()
{
  ConstitutiveDriver::postInputInitialization();

  // Validate some inputs
  GEOS_ERROR_IF( m_outputMassDensity != 0 && m_outputMassDensity != 1,
                 "Option can be either 0 (false) or 1 (true)",
                 getWrapperDataContext( viewKeyStruct::outputMassDensityString() ) );

  GEOS_ERROR_IF( m_outputCompressibility != 0 && m_outputCompressibility != 1,
                 "Option can be either 0 (false) or 1 (true)",
                 getWrapperDataContext( viewKeyStruct::outputCompressibilityString() ) );

  GEOS_ERROR_IF( m_outputPhaseComposition != 0 && m_outputPhaseComposition != 1,
                 "Option can be either 0 (false) or 1 (true)",
                 getWrapperDataContext( viewKeyStruct::outputPhaseCompositionString() ) );

  // Check that the functions exist
  FunctionManager & functionManager = FunctionManager::getInstance();
  GEOS_ERROR_IF( !functionManager.hasGroup< TableFunction >( m_pressureFunctionName ),
                 getWrapperDataContext( viewKeyStruct::pressureFunctionString() ) <<
                 "Pressure function with name '" << m_pressureFunctionName << "' not found" );

  GEOS_ERROR_IF( !functionManager.hasGroup< TableFunction >( m_temperatureFunctionName ),
                 getWrapperDataContext( viewKeyStruct::temperatureFunctionString() ) <<
                 "Temperature function with name '" << m_temperatureFunctionName << "' not found" );

  // get number of phases and components
  MultiFluidBase & baseFluid = getFluid();

  m_numPhases = baseFluid.numFluidPhases();
  m_numComponents = baseFluid.numFluidComponents();

  GEOS_ERROR_IF( m_feed.size() != m_numComponents,
                 getWrapperDataContext( viewKeyStruct::feedString() ) <<
                 "Feed must have the same number of components as the fluid. "
                 "Feed has " << m_feed.size() << " and fluid has " << m_numComponents );

  string_array columnNames;
  getColumnNames( columnNames );
  integer const numCols = static_cast< integer >(columnNames.size());

  // initialize functions
  TableFunction & pressureFunction = functionManager.getGroup< TableFunction >( m_pressureFunctionName );
  TableFunction & temperatureFunction = functionManager.getGroup< TableFunction >( m_temperatureFunctionName );

  pressureFunction.initializeFunction();
  temperatureFunction.initializeFunction();

  ArrayOfArraysView< real64 > coordinates = pressureFunction.getCoordinates();
  real64 const minTime = coordinates[0][0];
  real64 const maxTime = coordinates[0][coordinates.sizeOfArray( 0 )-1];

  // Allocate the data
  allocateTable( numCols, minTime, maxTime );

  // set input columns
  integer const numRows = m_table.size( 0 );
  for( integer step = 0; step < numRows; ++step )
  {
    real64 const time = m_table( step, TIME );
    m_table( step, PRES ) = pressureFunction.evaluate( &time );
    m_table( step, TEMP ) = temperatureFunction.evaluate( &time );
  }
}

bool PVTDriver::execute()
{
  // get the fluid out of the constitutive manager.
  // for the moment it is of type MultiFluidBase.

  MultiFluidBase & baseFluid = getFluid();

  // depending on logLevel, print some useful info

  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "Launching PVT Driver" );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Fluid .................. " << m_fluidName );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Type ................... " << baseFluid.getCatalogName() );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  No. of Phases .......... " << m_numPhases );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  No. of Components ...... " << m_numComponents );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Pressure Control ....... " << m_pressureFunctionName );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Temperature Control .... " << m_temperatureFunctionName );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Feed ................... " << m_feed );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Steps .................. " << m_numSteps );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Output ................. " << (m_outputFile.empty() ? "<none>" : m_outputFile) );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Output Mass Density .... " << m_outputMassDensity );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Output Compressibility . " << m_outputCompressibility );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Output Phase Comp. ..... " << m_outputPhaseComposition );

  // storing constitutive data
  conduit::Node node;
  dataRepository::Group rootGroup( "root", node );
  dataRepository::Group discretization( "discretization", &rootGroup );

  // Allocate as many elements as the number of rows
  integer const numRows = m_table.size( 0 );
  discretization.resize( numRows );       // numRows elements
  baseFluid.allocateConstitutiveData( discretization, 1 );   // one quadrature point

  baseFluid.initializeState();

  // pass the fluid through the ConstitutivePassThru to downcast from the
  // base type to a known model type.  the lambda here then executes the
  // appropriate test driver. note that these calls will move data to device if available.

  constitutiveUpdatePassThru( baseFluid, [&] ( auto & selectedFluid )
  {
    using FLUID_TYPE = TYPEOFREF( selectedFluid );
    runTest< FLUID_TYPE >( selectedFluid, m_table );
  } );

  return false;
}

void PVTDriver::getColumnNames( string_array & columnNames ) const
{
  MultiFluidBase const & baseFluid = getFluid();
  auto const phaseNames = baseFluid.phaseNames();
  auto const componentNames = baseFluid.componentNames();

  // Output columns depend on options
  columnNames.emplace_back( "time" );
  columnNames.emplace_back( "pressure" );
  columnNames.emplace_back( "temperature" );
  columnNames.emplace_back( "density" );
  if( m_outputCompressibility != 0 )
  {
    columnNames.emplace_back( "total compressibility" );
  }

  for( auto const & phaseName : phaseNames )
  {
    columnNames.emplace_back( GEOS_FMT( "phase fraction,{}", phaseName ));
  }
  for( auto const & phaseName : phaseNames )
  {
    columnNames.emplace_back( GEOS_FMT( "phase density,{}", phaseName ));
  }
  if( m_outputMassDensity != 0 )
  {
    for( auto const & phaseName : phaseNames )
    {
      columnNames.emplace_back( GEOS_FMT( "phase mass density,{}", phaseName ));
    }
  }
  for( auto const & phaseName : phaseNames )
  {
    columnNames.emplace_back( GEOS_FMT( "phase viscosity,{}", phaseName ));
  }

  if( baseFluid.isThermal())
  {
    for( auto const & phaseName : phaseNames )
    {
      columnNames.emplace_back( GEOS_FMT( "phase enthalpy,{}", phaseName ));
    }
  }

  if( m_outputPhaseComposition != 0 )
  {
    for( auto const & phaseName : phaseNames )
    {
      for( auto const & compName : componentNames )
      {
        columnNames.emplace_back( GEOS_FMT( "phase composition,{},{}", phaseName, compName ));
      }
    }
  }
}

MultiFluidBase & PVTDriver::getFluid()
{
  return getConstitutiveManager().getGroup< MultiFluidBase >( m_fluidName );
}

MultiFluidBase const & PVTDriver::getFluid() const
{
  return getConstitutiveManager().getGroup< MultiFluidBase >( m_fluidName );
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        PVTDriver,
                        string const &, dataRepository::Group * const )

} /* namespace geos */
