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

#include "common/MpiWrapper.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutiveDrivers/contact/LogLevelsInfo.hpp"
#include "constitutive/contact/FrictionBase.hpp"
#include "constitutive/contact/FrictionSelector.hpp"

#include "FrictionDriver.hpp"

namespace geos {

using namespace dataRepository;
using namespace constitutive;

FrictionDriver::FrictionDriver(const string& name, Group * const parent)
:
TaskBase(name, parent)
{
  registerWrapper( viewKeyStruct::frictionNameString(), &m_frictionName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Relperm model to test" );

  registerWrapper( viewKeyStruct::numStepsString(), &m_numSteps ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Number of saturation steps to take" );

  registerWrapper( viewKeyStruct::outputString(), &m_outputFile ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( "none" ).
    setDescription( "Output file" );

  registerWrapper( viewKeyStruct::baselineString(), &m_baselineFile ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( "none" ).
    setDescription( "Baseline file" );

  addLogLevel< logInfo::LogOutput >();
}

void FrictionDriver::outputResults()
{
  // TODO: improve file path output to grab command line -o directory
  //       for the moment, we just use the specified m_outputFile directly

  FILE * fp = fopen( m_outputFile.c_str(), "w" );

  fprintf( fp, "# column 1 = time\n" );
  fprintf( fp, "# columns %d-%d = phase vol fractions\n", 2, 1 + m_numPhases );
  fprintf( fp, "# columns %d-%d = phase relperm\n", 2 + m_numPhases, 1 + 2 * m_numPhases );

  if( ( m_numPhases == 2 && m_table.size( 1 ) > 5 ) || m_table.size( 1 ) > 7 )
  {
    fprintf( fp, "# columns %d-%d = phase relperm (hyst)\n", 1 + 2 * m_numPhases, 1 + 3 * m_numPhases );
  }


  for( integer n = 0; n < m_table.size( 0 ); ++n )
  {
    for( integer col = 0; col < m_table.size( 1 ); ++col )
    {
      fprintf( fp, "%.4e ", m_table( n, col ) );
    }
    fprintf( fp, "\n" );
  }
  fclose( fp );


}


void FrictionDriver::postInputInitialization()
{
  ConstitutiveManager
  & constitutiveManager = this->getGroupByPath< ConstitutiveManager >( "/Problem/domain/Constitutive" );
  FrictionBase& baseFriction = constitutiveManager.getGroup< FrictionBase >( m_frictionName );

//   m_numPhases = baseFriction.numSubGroups();

}


bool FrictionDriver::execute( const geos::real64 GEOS_UNUSED_PARAM( time_n ),
                             const geos::real64 GEOS_UNUSED_PARAM( dt ),
                             const geos::integer GEOS_UNUSED_PARAM( cycleNumber ),
                             const geos::integer GEOS_UNUSED_PARAM( eventCounter ),
                             const geos::real64 GEOS_UNUSED_PARAM( eventProgress ),
                             geos::DomainPartition &
                             GEOS_UNUSED_PARAM( domain ) )
{
  // this code only makes sense in serial

  GEOS_THROW_IF( MpiWrapper::commRank() > 0, "RelpermDriver should only be run in serial", std::runtime_error );


  ConstitutiveManager
  & constitutiveManager = this->getGroupByPath< ConstitutiveManager >( "/Problem/domain/Constitutive" );
  FrictionBase
  & baseFriction = constitutiveManager.getGroup< FrictionBase >( m_frictionName );

  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "Launching Relperm Driver" );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Relperm .................. " << m_frictionName );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Type ................... " << baseFriction.getCatalogName() );
//   GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  No. of Phases .......... " << m_numPhases );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Steps .................. " << m_numSteps );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Output ................. " << m_outputFile );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Baseline ............... " << m_baselineFile );

  // create a dummy discretization with one quadrature point for
  // storing constitutive data

  conduit::Node node;
  dataRepository::Group rootGroup( "root", node );
  dataRepository::Group discretization( "discretization", &rootGroup );

  discretization.resize( 1 );   // one element
  baseRelperm.allocateConstitutiveData( discretization, 1 );   // one quadrature point

  constitutiveUpdatePassThru( baseFriction, [&]( auto & selectedFrictionModel )
  {
    using RELPERM_TYPE = TYPEOFREF( selectedFrictionModel );
    resizeTables< RELPERM_TYPE >();
    runTest< RELPERM_TYPE >( selectedFrictionModel, m_table );
  } );

  // move table back to host for output
  m_table.move( LvArray::MemorySpace::host );

  if( m_outputFile != "none" )
  {
    outputResults();
  }

  if( m_baselineFile != "none" )
  {
    compareWithBaseline();
  }

  return false;
}


template< typename RELPERM_TYPE >
void FrictionDriver::resizeTables()
{
  ConstitutiveManager
  & constitutiveManager = this->getGroupByPath< ConstitutiveManager >( "/Problem/domain/Constitutive" );
  FrictionBase
  & baseRelperm = constitutiveManager.getGroup< FrictionBaseBase >( m_frictionName );

  using PT = RelativePermeabilityBase::PhaseType;
  integer const ipWater = baseRelperm.getPhaseOrder()[PT::WATER];
  integer const ipOil = baseRelperm.getPhaseOrder()[PT::OIL];
  integer const ipGas = baseRelperm.getPhaseOrder()[PT::GAS];

  real64 minSw = 0., minSnw = 0.;
  if( baseRelperm.numFluidPhases() > 2 )
  {
    minSw = baseRelperm.getWettingPhaseMinVolumeFraction();
    minSnw = baseRelperm.getNonWettingMinVolumeFraction();
  }
  else
  {
    if( ipWater < 0 )// a.k.a o/g
    {
      minSw = 0;
      minSnw = baseRelperm.getNonWettingMinVolumeFraction();
    }
    else if( ipGas < 0 || ipOil < 0 )// a.k.a w/o or w/g
    {
      minSnw = 0;
      minSw = baseRelperm.getWettingPhaseMinVolumeFraction();
    }
  }

  real64 const dSw = ( 1 - minSw - minSnw ) / m_numSteps;
  // set input columns

  resizeTable< RELPERM_TYPE >();
  // 3-phase branch
  if( m_numPhases > 2 )
  {
    for( integer ni = 0; ni < m_numSteps + 1; ++ni )
    {
      for( integer nj = 0; nj < m_numSteps + 1; ++nj )
      {

        integer index = ni * ( m_numSteps + 1 ) + nj;
        m_table( index, TIME ) = minSw + index * dSw;
        m_table( index, ipWater + 1 ) = minSw + nj * dSw;
        m_table( index, ipGas + 1 ) = minSnw + ni * dSw;
        m_table( index, ipOil + 1 ) =
          1. - m_table( index, ipWater + 1 ) - m_table( index, ipOil + 1 );
      }
    }
  }
  else // 2-phase branch
  {
    for( integer ni = 0; ni < m_numSteps + 1; ++ni )
    {
      integer index = ni;
      m_table( index, TIME ) = minSw + index * dSw;
      if( ipWater < 0 )
      {
        m_table( index, ipGas + 1 ) = minSnw + ni * dSw;
        m_table( index, ipOil + 1 ) = 1. - m_table( index, ipGas + 1 );
      }
      else if( ipGas < 0 )
      {
        m_table( index, ipWater + 1 ) = minSw + ni * dSw;
        m_table( index, ipOil + 1 ) = 1. - m_table( index, ipWater + 1 );
      }
      else if( ipOil < 0 )
      {
        m_table( index, ipWater + 1 ) = minSw + ni * dSw;
        m_table( index, ipGas + 1 ) = 1. - m_table( index, ipWater + 1 );
      }
    }

  }


}


// template< typename RELPERM_TYPE >
// std::enable_if_t< std::is_same< TableRelativePermeabilityHysteresis, RELPERM_TYPE >::value, void >
// RelpermDriver::resizeTable()
// {
//   if( m_numPhases > 2 )
//   {
//     m_table.resize( ( m_numSteps + 1 ) * ( m_numSteps + 1 ), 1 + 3 * m_numPhases );
//   }
//   else
//   {
//     m_table.resize( m_numSteps + 1, 1 + 3 * m_numPhases );
//   }

// }


}