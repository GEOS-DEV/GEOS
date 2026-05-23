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

#include "RelpermDriver.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "constitutiveDrivers/LogLevelsInfo.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilitySelector.hpp"
#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"
#include "common/format/StringUtilities.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

RelpermDriver::RelpermDriver( const string & name,
                              Group * const parent )
  : ConstitutiveDriver( name, parent )
{
  registerWrapper( viewKeyStruct::relpermNameString(), &m_relpermName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Relative permeability model to test" );

  registerWrapper( viewKeyStruct::phaseNamesString(), &m_phaseNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "The names of the phases for which saturations are defined" );

  registerWrapper( viewKeyStruct::saturationFunctionsString(), &m_saturationFunctionNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Functions controlling saturation time history of the selected phases" );
}

void RelpermDriver::postInputInitialization()
{
  ConstitutiveDriver::postInputInitialization();

  RelativePermeabilityBase const & baseRelperm = getRelperm();

  integer const numPhases = baseRelperm.numFluidPhases();

  // Must be 2-phase or 3-phase
  GEOS_ERROR_IF( numPhases < 2 || 3 < numPhases,
                 "Number of phases for relative permeability model must be 2 or 3",
                 getWrapperDataContext( viewKeyStruct::relpermNameString() ) );

  // Table functions should be provided for np-1 phases
  integer const numSelectedPhases = m_phaseNames.size();
  GEOS_ERROR_IF( numSelectedPhases != numPhases - 1,
                 GEOS_FMT( "Number of selected phases should be {} not {}", numPhases - 1, numSelectedPhases ),
                 getWrapperDataContext( viewKeyStruct::phaseNamesString() ) );

  integer const numFunctions = m_saturationFunctionNames.size();
  GEOS_ERROR_IF( numFunctions != numSelectedPhases,
                 "Number of saturations functions should match the number of selected phases",
                 getWrapperDataContext( viewKeyStruct::saturationFunctionsString() ) );

  // Check that the phase names are valid
  std::set< string > const relpermPhases( baseRelperm.phaseNames().begin(), baseRelperm.phaseNames().end());
  std::set< string > seenPhases;
  for( auto const & phaseName : m_phaseNames )
  {
    GEOS_ERROR_IF ( relpermPhases.find( phaseName ) == relpermPhases.end(),
                    GEOS_FMT( "Phase {} is not in the list of allowed phases for the relative permeability model", phaseName ),
                    getWrapperDataContext( viewKeyStruct::phaseNamesString() ) );

    GEOS_ERROR_IF ( seenPhases.find( phaseName ) != seenPhases.end(),
                    GEOS_FMT( "Phase {} is repeated in the list of selected phases", phaseName ),
                    getWrapperDataContext( viewKeyStruct::phaseNamesString() ) );

    seenPhases.insert( phaseName );
  }

  // Check that the functions exist
  FunctionManager & functionManager = FunctionManager::getInstance();
  for( auto const & functionName : m_saturationFunctionNames )
  {
    GEOS_ERROR_IF( !functionManager.hasGroup< TableFunction >( functionName ),
                   GEOS_FMT( "Saturation function with name '{}' not found", functionName ),
                   getWrapperDataContext( viewKeyStruct::saturationFunctionsString() ) );
  }

  string_array columnNames;
  getColumnNames( columnNames );
  integer const numCols = static_cast< integer >(columnNames.size());

  // Initialize functions and extract limits
  real64 minTime =  LvArray::NumericLimits< real64 >::max;
  real64 maxTime = -LvArray::NumericLimits< real64 >::max;
  for( auto const & functionName : m_saturationFunctionNames )
  {
    TableFunction & function = functionManager.getGroup< TableFunction >( functionName );
    function.initializeFunction();
    ArrayOfArraysView< real64 > coordinates = function.getCoordinates();
    minTime = LvArray::math::min( minTime, coordinates[0][0] );
    maxTime = LvArray::math::max( maxTime, coordinates[0][coordinates.sizeOfArray( 0 )-1] );
  }

  // Allocate the data
  allocateTable( numCols, minTime, maxTime );

  // Populate the data
  initializeTable( baseRelperm );
}

bool RelpermDriver::execute()
{
  RelativePermeabilityBase & baseRelperm = getRelperm();

  integer const numPhases = baseRelperm.numFluidPhases();

  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "Launching Relperm Driver" );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Relperm ................ " << m_relpermName );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Type ................... " << baseRelperm.getCatalogName() );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  No. of Phases .......... " << numPhases );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Steps .................. " << m_numSteps );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Selected phases ........ " << stringutilities::join( m_phaseNames, ", " ) );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Saturation functions ... " << stringutilities::join( m_saturationFunctionNames, ", " ) );
  if( !m_outputFile.empty())
  {
    GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Output ................. " << m_outputFile );
  }

  // create a dummy discretization with one quadrature point for
  // storing constitutive data
  conduit::Node node;
  dataRepository::Group rootGroup( "root", node );
  dataRepository::Group discretization( "discretization", &rootGroup );

  // Allocate as many elements as the number of rows
  integer const numRows = m_table.size( 0 );
  discretization.resize( numRows );   // numRows elements
  baseRelperm.allocateConstitutiveData( discretization, 1 );   // one quadrature point

  constitutiveUpdatePassThru( baseRelperm, [&]( auto & selectedRelpermModel )
  {
    using RELPERM_TYPE = TYPEOFREF( selectedRelpermModel );
    runTest< RELPERM_TYPE >( selectedRelpermModel, m_table );
  } );

  // move table back to host for output
  m_table.move( LvArray::MemorySpace::host );

  return false;
}

void RelpermDriver::getColumnNames( string_array & columnNames ) const
{
  RelativePermeabilityBase const & baseRelperm = getRelperm();
  bool const has_hysteresis = (dynamic_cast< constitutive::TableRelativePermeabilityHysteresis const * >(&baseRelperm) != nullptr);

  integer const numPhases = baseRelperm.numFluidPhases();
  string_array const & phaseNames = baseRelperm.phaseNames();

  columnNames.emplace_back( "index" );
  for( integer ip = 0; ip < numPhases; ip++ )
  {
    columnNames.emplace_back( GEOS_FMT( "saturation,{}", phaseNames[ip] ));
  }
  if( has_hysteresis )
  {
    for( integer ip = 0; ip < numPhases; ip++ )
    {
      columnNames.emplace_back( GEOS_FMT( "historical saturation,{}", phaseNames[ip] ));
    }
  }
  for( integer ip = 0; ip < numPhases; ip++ )
  {
    columnNames.emplace_back( GEOS_FMT( "relperm,{}", phaseNames[ip] ));
  }
}

void RelpermDriver::initializeTable( RelativePermeabilityBase const & baseRelperm )
{
  integer const numRows = m_table.size( 0 );
  integer const numPhases = baseRelperm.numFluidPhases();

  string_array const & phaseNames = baseRelperm.phaseNames();
  array1d< integer > phaseOrder( numPhases );
  phaseOrder[0] = numPhases * (numPhases - 1) / 2;
  for( integer ip = 0; ip < numPhases-1; ++ip )
  {
    integer const index = static_cast< integer >(std::distance( phaseNames.begin(), std::find( phaseNames.begin(), phaseNames.end(), m_phaseNames[ip] )));
    phaseOrder[ip+1] = index;
    phaseOrder[0] -= index;
  }

  FunctionManager const & functionManager = FunctionManager::getInstance();
  stdVector< TableFunction const * > tableFunctions;
  for( auto const & functionName : m_saturationFunctionNames )
  {
    TableFunction const * function = functionManager.getGroupPointer< TableFunction >( functionName );
    tableFunctions.emplace_back( function );
  }

  // Offset for saturations in table
  constexpr integer SATURATION = 1;

  for( integer step = 0; step < numRows; ++step )
  {
    real64 const time = m_table( step, TIME );
    real64 sumSaturation = 0.0;
    for( integer ip = 1; ip < numPhases; ip++ )
    {
      real64 const saturation = LvArray::math::max( tableFunctions[ip-1]->evaluate( &time ), 0.0 );
      m_table( step, phaseOrder[ip] + SATURATION ) = saturation;
      sumSaturation += saturation;
    }
    if( 1.0 - sumSaturation < -LvArray::NumericLimits< real64 >::epsilon )
    {
      real64 const scale = 1.0 / sumSaturation;
      for( integer ip = 1; ip < numPhases; ip++ )
      {
        m_table( step, phaseOrder[ip] + SATURATION ) *= scale;
      }
      sumSaturation = 1.0;
    }
    m_table( step, phaseOrder[0] + SATURATION ) = 1.0 - sumSaturation;
  }
}

RelativePermeabilityBase & RelpermDriver::getRelperm()
{
  return getConstitutiveManager().getGroup< RelativePermeabilityBase >( m_relpermName );
}

RelativePermeabilityBase const & RelpermDriver::getRelperm() const
{
  return getConstitutiveManager().getGroup< RelativePermeabilityBase >( m_relpermName );
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        RelpermDriver,
                        string const &, dataRepository::Group * const )

}
