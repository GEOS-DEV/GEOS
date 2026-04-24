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

  registerWrapper( viewKeyStruct::historicalSaturationsString(), &m_historicalSaturations ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Historical saturations for each phase." );
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

  // Historical saturations must be the same number as the phases
  if( !m_historicalSaturations.empty())
  {
    GEOS_ERROR_IF( m_historicalSaturations.size() != numPhases,
                   "Number of historical saturations must be the same as the number of phases",
                   getWrapperDataContext( viewKeyStruct::historicalSaturationsString() ) );
  }

  string_array columnNames;
  getColumnNames( columnNames );
  integer const numCols = static_cast< integer >(columnNames.size());

  allocateTable( numCols, numPhases );
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
  if( !m_outputFile.empty())
  {
    GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Output ................. " << m_outputFile );
  }

  initializeTable();

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

void RelpermDriver::allocateTable( integer numColumns, integer numPhases )
{
  // For 3-phase we have m_numSteps+1 points for each of the other two phases
  integer const numRows = (numPhases == 3) ? (m_numSteps+1)*(m_numSteps+1) : (m_numSteps+1);
  m_table.resize( numRows, numColumns );
  for( integer index = 0; index < numRows; ++index )
  {
    m_table( index, TIME ) = index;
  }
}

void RelpermDriver::initializeTable()
{
  RelativePermeabilityBase & baseRelperm = getRelperm();

  using PT = RelativePermeabilityBase::PhaseType;
  integer const ipWater = baseRelperm.getPhaseOrder()[PT::WATER];
  integer const ipOil = baseRelperm.getPhaseOrder()[PT::OIL];
  integer const ipGas = baseRelperm.getPhaseOrder()[PT::GAS];

  integer const numPhases = baseRelperm.numFluidPhases();

  auto const [ipWetting, ipNonWetting] = baseRelperm.wettingAndNonWettingPhaseIndices();
  real64 const min_wetting_saturation = baseRelperm.getPhaseMinVolumeFraction()[ipWetting];
  real64 const min_non_wetting_saturation = baseRelperm.getPhaseMinVolumeFraction()[ipNonWetting];

  real64 const dSw = ( 1.0 - min_wetting_saturation - min_non_wetting_saturation ) / m_numSteps;

  // Offset for saturations in table
  constexpr integer SATURATION = 1;

  // 3-phase branch
  if( numPhases == 3 )
  {
    real64 swat = 0.0;
    real64 sgas = 0.0;
    for( integer ni = 0; ni < m_numSteps + 1; ++ni )
    {
      swat = min_wetting_saturation + ni*dSw;
      for( integer nj = 0; nj < m_numSteps + 1; ++nj )
      {
        sgas = min_non_wetting_saturation + nj*dSw;

        integer index = ni * ( m_numSteps + 1 ) + nj;
        m_table( index, ipWater + SATURATION ) = swat;
        m_table( index, ipGas + SATURATION ) = sgas;
        m_table( index, ipOil + SATURATION ) = 1.0 - swat - sgas;
      }
    }
  }
  else // 2-phase branch
  {
    for( integer ni = 0; ni < m_numSteps + 1; ++ni )
    {
      real64 const s_nw = min_non_wetting_saturation + ni * dSw;
      m_table( ni, ipNonWetting + SATURATION ) = s_nw;
      m_table( ni, ipWetting + SATURATION ) = 1.0 - s_nw;
    }
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
