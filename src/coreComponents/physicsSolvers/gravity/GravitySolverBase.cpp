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
 * @file GravitySolverBase.cpp
 */

#include "GravitySolverBase.hpp"
#include <filesystem>


namespace geos
{

using namespace dataRepository;

GravitySolverBase::GravitySolverBase( const std::string & name, Group * const parent )
  : PhysicsSolverBase( name, parent )
{
  registerWrapper( viewKeyStruct::modeString(), &m_modeString )
    .setInputFlag( InputFlags::OPTIONAL )
    .setApplyDefaultValue( "modeling" )
    .setDescription( "Mode: 'modeling' (default) or 'adjoint'" );

  registerWrapper( viewKeyStruct::stationCoordinatesString(), &m_stationCoordinates )
    .setInputFlag( InputFlags::REQUIRED )
    .setSizedFromParent( 0 )
    .setDescription( "Coordinates (x,y,z) of the gravimeter stations" );

  registerWrapper( viewKeyStruct::outputGzString(), &m_outputGz )
    .setInputFlag( InputFlags::OPTIONAL )
    .setApplyDefaultValue( 1 )
    .setDescription( "Flag to dump to disk field recorded at gravimeters" );

  registerWrapper( viewKeyStruct::outputGzBasenameString(), &m_outputGzBasename )
    .setInputFlag( InputFlags::OPTIONAL )
    .setApplyDefaultValue( "gz" )
    .setDescription( "Basename to dump to disk field recorded at gravimeters" );

  registerWrapper( viewKeyStruct::residueString(), &m_residue )
    .setInputFlag( InputFlags::OPTIONAL )
    .setSizedFromParent( 0 )
    .setDescription( "Residue at each station" );

  registerWrapper( viewKeyStruct::gzAtStationsString(), &m_gzAtStations )
    .setInputFlag( InputFlags::FALSE )
    .setSizedFromParent( 0 )
    .setDescription( "Gz value at each station" );
}

GravitySolverBase::~GravitySolverBase() = default;


GravitySolverBase::GravityMode GravitySolverBase::parseMode( const std::string & modeStr )
{
  std::string lowerMode = modeStr;
  std::transform( lowerMode.begin(), lowerMode.end(), lowerMode.begin(),
                  [] ( unsigned char c ) -> unsigned char { return std::tolower( c ); } );

  auto it = modeMap.find( lowerMode );
  if( it != modeMap.end())
    return it->second;

  std::ostringstream oss;
  oss << "Invalid gravity mode string: '" << modeStr << "'. Valid options are: ";
  bool first = true;
  for( const auto & pair : modeMap )
  {
    if( !first )
      oss << ", ";
    oss << pair.first;
    first = false;
  }
  throw std::invalid_argument( oss.str());
}


void GravitySolverBase::reinit()
{
  initializePostInitialConditionsPreSubGroups();
}

void GravitySolverBase::initializePreSubGroups()
{
  PhysicsSolverBase::initializePreSubGroups();
}

void GravitySolverBase::initializePostInitialConditionsPreSubGroups()
{
  PhysicsSolverBase::initializePostInitialConditionsPreSubGroups();

  try
  {
    m_mode = parseMode( m_modeString );
  }
  catch( const std::invalid_argument & e )
  {
    GEOS_THROW( e.what(), InputError );
  }

  if( m_mode == GravityMode::Adjoint )
  {


    const auto stationCount = m_stationCoordinates.size( 0 );
    const auto residueCount = m_residue.size( 0 );

    GEOS_THROW_IF(
      stationCount != residueCount,
      "GravitySolverBase: Residue size (" + std::to_string( residueCount ) +
      ") does not match the number of stations (" + std::to_string( stationCount ) + ")",
      InputError
      );
  }


  if( this->getLogLevel()>0 )
  {
    constexpr auto yesno = [] (int flag) noexcept->const char * { return flag ? "yes" : "no"; };

    LogPart gravitySolverLog( "Gravity Solver: ", MpiWrapper::commRank() == 0 );
    gravitySolverLog.begin();
    GEOS_LOG_RANK_0( "Name:                        " << getName());
    GEOS_LOG_RANK_0( "Mode:                        " << m_modeString );
    GEOS_LOG_RANK_0( "Number of stations:          " << m_stationCoordinates.size( 0 ));
    GEOS_LOG_RANK_0( "Output Gz to file:           " << yesno( m_outputGz ));
    GEOS_LOG_RANK_0( "Output Gz basename:          " << m_outputGzBasename );
    GEOS_LOG_RANK_0( "Log level:                   " << getLogLevel());
    GEOS_LOG_RANK_0( "  Output Gz to logs:         " << yesno( getLogLevel() > 1 ));
    GEOS_LOG_RANK_0( "  Output Adjoint to logs:    " << yesno( getLogLevel() > 2 ));
    GEOS_LOG_RANK_0( "  Output Properties to logs: " << yesno( getLogLevel() > 3 ));
    gravitySolverLog.end();
  }
}

void GravitySolverBase::postInputInitialization()
{
  PhysicsSolverBase::postInputInitialization();

  GEOS_THROW_IF( m_stationCoordinates.size( 1 ) != 3,
                 "GravitySolverBase: Invalid number of physical coordinates for the gravimeter stations",
                 InputError );

  m_gzAtStations.resize( m_stationCoordinates.size( 0 ));
}

real64 GravitySolverBase::solverStep( real64 const & time_n,
                                      real64 const & dt,
                                      integer const cycleNumber,
                                      DomainPartition & domain )
{
  GEOS_LOG_RANK_0( "GravitySolverBase: time_n: " << time_n << " dt: " << dt );
  return explicitStep( time_n, dt, cycleNumber, domain );
}

real64 GravitySolverBase::explicitStep( real64 const & time_n,
                                        real64 const & dt,
                                        integer const cycleNumber,
                                        DomainPartition & domain )
{
  switch( m_mode )
  {
    case GravityMode::Modeling:
      return explicitStepModeling( time_n, dt, cycleNumber, domain );
    case GravityMode::Adjoint:
      return explicitStepAdjoint( time_n, dt, cycleNumber, domain );
    default:
      GEOS_THROW( "GravitySolverBase::explicitStep: Mode not supported", InputError );
  }
}

void GravitySolverBase::saveGz( real64 const & time_n,
                                integer const cycleNumber,
                                string const basename,
                                arrayView1d< real64 > const & gzAtStations )
{
  // Convert to float32
  std::vector< float > tmp32( gzAtStations.size());
  std::transform( gzAtStations.begin(), gzAtStations.end(), tmp32.begin(),
                  [] ( real64 val )-> float { return static_cast< float >(val); } );

  // Binary data
  std::filesystem::path dataFilename = basename + GEOS_FMT( "_{:015}.H@", time_n );
  std::ofstream fout( dataFilename, std::ios::binary );
  if( !fout )
    throw std::ios_base::failure( "Failed to open data file" );
  fout.write( reinterpret_cast< const char * >(tmp32.data()), tmp32.size() * sizeof(float));

  // Header file
  std::filesystem::path headerFilename = basename + GEOS_FMT( "_{:015}.H", time_n );
  std::ofstream fheader( headerFilename );
  if( !fheader )
    throw std::ios_base::failure( "Failed to open header file" );

  fheader << "n1=" << gzAtStations.size() << '\n'
          << "d1=1.\n"
          << "o1=0.\n"
          << "label1=STATION\n"
          << "esize=4\n"
          << "data_format=native_float\n"
          << "data_label=GZ\n"
          << "in=" << dataFilename << '\n'
          << GEOS_FMT( "time={:015}", time_n ) << '\n'
          << "cycle=" << cycleNumber << '\n';
}

} // namespace geos
