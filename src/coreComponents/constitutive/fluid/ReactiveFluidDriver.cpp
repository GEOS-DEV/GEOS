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
 * @file ReactiveFluidDriver.cpp
 */

#include "ReactiveFluidDriver.hpp"
#include "fileIO/Outputs/OutputBase.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;


ReactiveFluidDriver::ReactiveFluidDriver( const string & name,
                                          Group * const parent ):
  TaskBase( name, parent )
{
  enableLogLevelInput();

  registerWrapper( viewKeyStruct::fluidNameString(), &m_fluidName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Fluid to test" );

  registerWrapper( viewKeyStruct::feedString(), &m_feed ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Feed composition array: total concentration of the primary species " );

  registerWrapper( viewKeyStruct::pressureFunctionString(), &m_pressureFunctionName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Function controlling pressure time history" );

  registerWrapper( viewKeyStruct::temperatureFunctionString(), &m_temperatureFunctionName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Function controlling temperature time history" );

  registerWrapper( viewKeyStruct::numStepsString(), &m_numSteps ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Number of load steps to take" );

  registerWrapper( viewKeyStruct::outputString(), &m_outputFile ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( "none" ).
    setDescription( "Output file" );

  registerWrapper( viewKeyStruct::baselineString(), &m_baselineFile ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( "none" ).
    setDescription( "Baseline file" );
}


ReactiveFluidDriver::~ReactiveFluidDriver()
{}


void ReactiveFluidDriver::postProcessInput()
{
  // get number of phases and components

  ConstitutiveManager & constitutiveManager = this->getGroupByPath< ConstitutiveManager >( "/Problem/domain/Constitutive" );
  ReactiveMultiFluid & fluid = constitutiveManager.getGroup< ReactiveMultiFluid >( m_fluidName );

  m_numPhases = fluid.numFluidPhases();
  m_numPrimarySpecies = fluid.numPrimarySpecies();
  m_numSecondarySpecies = fluid.numSecondarySpecies();
  m_numComponents = fluid.numFluidComponents();

  // resize data table to fit number of timesteps and concentrations
  // (numRows,numCols) = (numSteps+1,3+numPrimarySpecies + numSecSpecies)
  // column order = time, pressure, temp, primarySpeciesConcentration, secondarySpeciesConcentration
  m_table.resize( m_numSteps+1, m_numPrimarySpecies + m_numSecondarySpecies + 3 );

  // initialize functions

  FunctionManager & functionManager = FunctionManager::getInstance();

  TableFunction & pressureFunction = functionManager.getGroup< TableFunction >( m_pressureFunctionName );
  TableFunction & temperatureFunction = functionManager.getGroup< TableFunction >( m_temperatureFunctionName );

  pressureFunction.initializeFunction();
  temperatureFunction.initializeFunction();

  // determine time increment

  ArrayOfArraysView< real64 > coordinates = pressureFunction.getCoordinates();
  real64 const minTime = coordinates[0][0];
  real64 const maxTime = coordinates[0][coordinates.sizeOfArray( 0 )-1];
  real64 const dt = (maxTime-minTime) / m_numSteps;

  // set input columns
  for( integer n=0; n<m_numSteps+1; ++n )
  {
    m_table( n, TIME ) = minTime + n*dt;
    m_table( n, PRES ) = pressureFunction.evaluate( &m_table( n, TIME ) );
    m_table( n, TEMP ) = temperatureFunction.evaluate( &m_table( n, TIME ) );
  }
}


bool ReactiveFluidDriver::execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                   real64 const GEOSX_UNUSED_PARAM( dt ),
                                   integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                   integer const GEOSX_UNUSED_PARAM( eventCounter ),
                                   real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                                   DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  // this code only makes sense in serial

  GEOSX_THROW_IF( MpiWrapper::commRank() > 0, "ReactiveFluidDriver should only be run in serial", std::runtime_error );

  // get the fluid out of the constitutive manager.
  // for the moment it is of type MultiFluidBase.

  ConstitutiveManager & constitutiveManager = this->getGroupByPath< ConstitutiveManager >( "/Problem/domain/Constitutive" );
  ReactiveMultiFluid & fluid = constitutiveManager.getGroup< ReactiveMultiFluid >( m_fluidName );

  // depending on logLevel, print some useful info

  if( getLogLevel() > 0 )
  {
    GEOSX_LOG_RANK_0( "Launching ReactiveFluid Driver" );
    GEOSX_LOG_RANK_0( "  Fluid .................. " << m_fluidName );
    GEOSX_LOG_RANK_0( "  Type ................... " << fluid.getCatalogName() );
    GEOSX_LOG_RANK_0( "  No. of Phases .......... " << m_numPhases );
    GEOSX_LOG_RANK_0( "  No. of Primary Species ...... " << m_numPrimarySpecies );
    GEOSX_LOG_RANK_0( "  No. of Secondary Species ...... " << m_numSecondarySpecies );
    GEOSX_LOG_RANK_0( "  Pressure Control ....... " << m_pressureFunctionName );
    GEOSX_LOG_RANK_0( "  Temperature Control .... " << m_temperatureFunctionName );
    GEOSX_LOG_RANK_0( "  Steps .................. " << m_numSteps );
    GEOSX_LOG_RANK_0( "  Output ................. " << m_outputFile );
    GEOSX_LOG_RANK_0( "  Baseline ............... " << m_baselineFile );
  }

  // create a dummy discretization with one quadrature point for
  // storing constitutive data

  conduit::Node node;
  dataRepository::Group rootGroup( "root", node );
  dataRepository::Group discretization( "discretization", &rootGroup );

  discretization.resize( 1 );   // one element
  fluid.allocateConstitutiveData( discretization, 1 );   // one quadrature point

  // Eventually we will dispatch but for now there is only one model.
  runTest< ReactiveMultiFluid >( selectedFluid, m_table );

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


template< typename FLUID_TYPE >
void ReactiveFluidDriver::runTest( FLUID_TYPE & fluid, arrayView2d< real64 > const & table )
{
  // get number of phases and components
  localIndex const numPhases = fluid.numFluidPhases();
  localIndex const numPrimarySpecies = fluid.numPrimarySpecies();
  localIndex const numSecondarySpecies = fluid.numSecondarySpecies();

  // get output data views
  arrayView2d< real64 const, multifluid::USD_FLUID > primarySpeciesConcentration   = fluid.primarySpeciesConcentration();
  arrayView2d< real64 const, multifluid::USD_FLUID > secondarySpeciesConcentration   = fluid.secondarySpeciesConcentration();

  // create kernel wrapper

  typename FLUID_TYPE::KernelWrapper kernelWrapper = fluid.createKernelWrapper();

  // set composition to user specified feed
  // it is more convenient to provide input in molar, so perform molar to mass conversion here

  GEOSX_ASSERT_EQ( numPrimarySpecies, m_feed.size() );
  array2d< real64, compflow::LAYOUT_COMP > primarySpeciesTotalConcentrationValues( 1, numPrimarySpecies );

  for( localIndex i = 0; i < numPrimarySpecies; ++i )
  {
    primarySpeciesTotalConcentrationValues[0][i] = m_feed[i];
  }

  arrayView2d< real64 const, compflow::USD_COMP > const primarySpeciesTotalConcentration = compositionValues;

  // perform fluid update using table (P,T) and save resulting compositions, etc.
  // note: column indexing should be kept consistent with output file header below.

  integer numSteps = m_numSteps;
  using ExecPolicy = typename FLUID_TYPE::exec_policy;
  forAll< ExecPolicy >( 1, [=]  GEOSX_HOST_DEVICE ( localIndex const ei )
  {
    for( integer n = 0; n <= numSteps; ++n )
    {
      kernelWrapper.update( ei, 0, table( n, PRES ), table( n, TEMP ), primarySpeciesTotalConcentration[0] );
      for( integer p=0; p<numPrimarySpecies; ++p )
      {
        table( n, TEMP+1+p ) = primarySpeciesConcentration( ei, 0, p );
      }
      for( integer s=0; s<numSecondarySpecies; ++s )
      {
        table( n, TEMP+1+numPrimarySpecies+s ) = secondarySpeciesConcentration( ei, 0, s );
      }
    }
  } );

}


void ReactiveFluidDriver::outputResults()
{
  // TODO: improve file path output to grab command line -o directory
  //       for the moment, we just use the specified m_outputFile directly

  FILE * fp = fopen( m_outputFile.c_str(), "w" );

  fprintf( fp, "# column 1 = time\n" );
  fprintf( fp, "# column 2 = pressure\n" );
  fprintf( fp, "# column 3 = temperature\n" );
  fprintf( fp, "# columns %d-%d = concentrations\n", 4, 3+m_numPrimarySpecies+m_numSecondarySpecies );

  for( integer n=0; n<m_table.size( 0 ); ++n )
  {
    for( integer col=0; col<m_table.size( 1 ); ++col )
    {
      fprintf( fp, "%.4e ", m_table( n, col ) );
    }
    fprintf( fp, "\n" );
  }
  fclose( fp );
}


void ReactiveFluidDriver::compareWithBaseline()
{
  // open baseline file

  std::ifstream file( m_baselineFile.c_str() );
  GEOSX_THROW_IF( !file.is_open(), "Can't seem to open the baseline file " << m_baselineFile, InputError );

  // discard file header

  string line;
  for( integer row=0; row < 7; ++row )
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
      GEOSX_THROW_IF( file.eof(), "Baseline file appears shorter than internal results", std::runtime_error );
      file >> value;

      error = fabs( m_table[row][col]-value ) / ( fabs( value )+1 );
      GEOSX_THROW_IF( error > m_baselineTol, "Results do not match baseline at data row " << row+1
                                                                                          << " (row " << row+m_numColumns << " with header)"
                                                                                          << " and column " << col+1, std::runtime_error );
    }
  }

  // check we actually reached the end of the baseline file

  file >> value;
  GEOSX_THROW_IF( !file.eof(), "Baseline file appears longer than internal results", std::runtime_error );

  // success

  if( getLogLevel() > 0 )
  {
    GEOSX_LOG_RANK_0( "  Comparison ............. Internal results consistent with baseline." );
  }

  file.close();
}


REGISTER_CATALOG_ENTRY( TaskBase,
                        ReactiveFluidDriver,
                        string const &, dataRepository::Group * const )

} /* namespace geosx */
