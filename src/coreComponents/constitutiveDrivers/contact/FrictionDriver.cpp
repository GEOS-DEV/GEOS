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

#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"
#include <cmath>
#include <type_traits>

#include "FrictionDriver.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

FrictionDriver::FrictionDriver( const string & name, Group * const parent )
  :
  TaskBase( name, parent )
{
  registerWrapper( viewKeyStruct::frictionNameString(), &m_frictionName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Friction model to test" );

  registerWrapper( viewKeyStruct::numStepsString(), &m_numSteps ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Number of sample step to take in both jumps and traction increments" );

  registerWrapper( viewKeyStruct::jumpFunctionString(), &m_jumpFunctionName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the input function representing jump function along world x-axis" );

  registerWrapper( viewKeyStruct::tractionFunctionString(), &m_tractionFunctionName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the input function representing traction function along world x-axis" );

  registerWrapper( viewKeyStruct::thetaString(), &m_theta ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Number of increment step to take in both jumps and traction increments" );

  registerWrapper( viewKeyStruct::phiString(), &m_phi ).
    setInputFlag( InputFlags::INVALID ).
    setDescription( "Number of increment step to take in both jumps and traction increments" );

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

void FrictionDriver::compareWithBaseline()
{
  // open baseline file

  std::ifstream file( m_baselineFile.c_str() );
  GEOS_THROW_IF( !file.is_open(), GEOS_FMT( "Can't seem to open the baseline file ", m_baselineFile ), InputError );

  // discard file header

  string line;
  for( integer row=0; row < m_numColumns; ++row )
  {
    getline( file, line );
  }

  // read data block.  we assume the file size is consistent with m_table,
  // but check for a premature end-of-file. we then compare results value by value.
  // we ignore the newton iteration and residual columns, as those may be platform
  // specific.

  real64 value;
  real64 error;

  for( integer row=0; row < m_table.size( 0 ); ++row )
  {
    for( integer col=0; col < m_table.size( 1 ); ++col )
    {
      GEOS_THROW_IF( file.eof(), "Baseline file appears shorter than internal results", InputError );
      file >> value;

      error = fabs( m_table[row][col]-value ) / ( fabs( value )+1 );
      GEOS_THROW_IF( error > m_baselineTol, GEOS_FMT( "Results do not match baseline at data row {} (row {} with header) and column {}",
                                                      row+1,
                                                      row+10,
                                                      col+1 ),
                     InputError );

    }
  }

  // check we actually reached the end of the baseline file

  file >> value;
  GEOS_THROW_IF( !file.eof(), "Baseline file appears longer than internal results", InputError );

  // success

  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Comparison ........ Internal results consistent with baseline." );

  file.close();
}

void FrictionDriver::outputResults()
{
  // TODO: improve file path output to grab command line -o directory
  //       for the moment, we just use the specified m_outputFile directly

  FILE * fp = fopen( m_outputFile.c_str(), "w" );

  fprintf( fp, "# column 1 = index \n" );
  fprintf( fp, "# column 2-3 = normal and in-plane displacement jump\n" );
  fprintf( fp, "# columns 5-7 = normal and in-place tractions\n" );
  fprintf( fp, "# columns 8 = fracture state (0:Stick,1-2:[new]Slip,3:Open)\n" );
  fprintf( fp, "# columns 9 = tau lim\n" );

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
//   ConstitutiveManager
//   & constitutiveManager = this->getGroupByPath< ConstitutiveManager >( "/Problem/domain/Constitutive" );
//   FrictionBase& baseFriction = constitutiveManager.getGroup< FrictionBase >( m_frictionName );

// //   m_numPhases = baseFriction.numSubGroups();

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

  GEOS_THROW_IF( MpiWrapper::commRank() > 0, "FrictionDriver should only be run in serial", geos::RuntimeError );


  ConstitutiveManager
  & constitutiveManager = this->getGroupByPath< ConstitutiveManager >( "/Problem/domain/Constitutive" );
  FrictionBase
  & baseFriction = constitutiveManager.getGroup< FrictionBase >( m_frictionName );

  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "Launching Friction Driver" );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Friction .................. " << m_frictionName );
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
  baseFriction.allocateConstitutiveData( discretization, 1 );   // one quadrature point

  constitutiveUpdatePassThru( baseFriction, [&]( auto & selectedFrictionModel )
  {
    using FRICTION_TYPE = TYPEOFREF( selectedFrictionModel );
    resizeTables< FRICTION_TYPE >();
    runTest< FRICTION_TYPE >( selectedFrictionModel, m_table );
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


template< typename FRICTION_TYPE >
void FrictionDriver::resizeTables()
{
  ConstitutiveManager
  & constitutiveManager = this->getGroupByPath< ConstitutiveManager >( "/Problem/domain/Constitutive" );
  FrictionBase
  & baseFriction = constitutiveManager.getGroup< FrictionBase >( m_frictionName );

  // initialize table functions
  FunctionManager & functionManager = FunctionManager::getInstance();

  TableFunction & jumpFunction = functionManager.getGroup< TableFunction >( m_jumpFunctionName );
  TableFunction & tractionFunction = functionManager.getGroup< TableFunction >( m_tractionFunctionName );

  jumpFunction.initializeFunction();
  tractionFunction.initializeFunction();

  ArrayOfArraysView< real64 > coordinates = jumpFunction.getCoordinates();
  real64 const minTime = coordinates[0][0];
  real64 const maxTime = coordinates[0][coordinates.sizeOfArray( 0 )-1];
  real64 const dt = (maxTime-minTime) / m_numSteps;

  // set input columns
  resizeTable< FRICTION_TYPE >();

  // set time column
  for( integer k=0; k<m_numSteps+1; ++k )
    for( integer n=0; n<m_numSteps+1; ++n )
    {
      m_table( n + k*(m_numSteps+1), TIME ) = minTime + n*dt;
    }


  //TODO Somewhere ALM.updateConfiguration(domain, dummy);

  real64 cohesion = 0.;
  real64 frictionCoeff = 0.;
  if constexpr ( std::is_same_v< CoulombFriction, FRICTION_TYPE > ) {

    cohesion = dynamic_cast< CoulombFriction * >(&baseFriction)->getCohesion();
    frictionCoeff = dynamic_cast< CoulombFriction * >(&baseFriction)->getFrictionCoeff();
  }

  //All variation
  for( integer nt = 0; nt < m_numSteps+1; ++nt )
  {
    for( integer nj = 0; nj < m_numSteps+1; ++nj )
    {

      integer index = nt * (m_numSteps+1) + nj;
      m_table( index, NTRAC ) =  tractionFunction.evaluate( &m_table( nt, TIME ))*sin( m_theta * M_PI/180 );
      m_table( index, STRAC0 ) = tractionFunction.evaluate( &m_table( nt, TIME ))*cos( m_theta * M_PI/180 );
      m_table( index, STRAC1 ) = tractionFunction.evaluate( &m_table( nt, TIME ))*cos( m_theta * M_PI/180 );

      m_table( index, NJUMP ) = jumpFunction.evaluate( &m_table( nj, TIME ))*sin( m_theta * M_PI/180 );
      m_table( index, SLIP0 ) = jumpFunction.evaluate( &m_table( nj, TIME ))*cos( m_theta * M_PI/180 );
      m_table( index, SLIP1 ) = jumpFunction.evaluate( &m_table( nj, TIME ))*cos( m_theta * M_PI/180 );

      m_table( index, FS ) = fields::contact::FractureState::Stick;

      //Only for Coulomb
      m_table( index, TLIM ) = cohesion - m_table( index, NTRAC ) * frictionCoeff;

    }
  }

}

//TODO updateConfig proxy
// simultaneaous branch
// unsim
// checkconstaint
// *keep track of iter


template< typename FRICTION_TYPE >
void
FrictionDriver::resizeTable()
{
  m_table.resize((m_numSteps + 1)*(m_numSteps + 1), m_numColumns );
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        FrictionDriver,
                        string const &, dataRepository::Group * const )

}
