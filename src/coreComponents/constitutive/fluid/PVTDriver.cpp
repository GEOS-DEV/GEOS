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
 * @file PVTDriver.cpp
 */

#include "PVTDriver.hpp"
#include "fileIO/Outputs/OutputBase.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;


PVTDriver::PVTDriver( const string & name,
                      Group * const parent ):
  TaskBase( name, parent )
{
  enableLogLevelInput();

  registerWrapper( viewKeyStruct::fluidNameString(), &m_fluidName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Fluid to test" );

  registerWrapper( viewKeyStruct::feedString(), &m_feed ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Feed composition array [mol fraction]" );

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


PVTDriver::~PVTDriver()
{}


void PVTDriver::postProcessInput()
{
  // get number of phases and components

  ConstitutiveManager & constitutiveManager = this->getGroupByPath< ConstitutiveManager >( "/Problem/domain/Constitutive" );
  MultiFluidBase & baseFluid = constitutiveManager.getGroup< MultiFluidBase >( m_fluidName );

  m_numPhases = baseFluid.numFluidPhases();
  m_numComponents = baseFluid.numFluidComponents();

  // resize data table to fit number of timesteps and fluid phases:
  // (numRows,numCols) = (numSteps+1,4+3*numPhases)
  // column order = time, pressure, temp, totalDensity, phaseFraction_{1:NP}, phaseDensity_{1:NP}, phaseViscosity_{1:NP}

  m_table.resize( m_numSteps+1, 3*m_numPhases+4 );

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


bool PVTDriver::execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                         real64 const GEOSX_UNUSED_PARAM( dt ),
                         integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                         integer const GEOSX_UNUSED_PARAM( eventCounter ),
                         real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                         DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  // this code only makes sense in serial

  GEOSX_THROW_IF( MpiWrapper::commRank() > 0, "PVTDriver should only be run in serial", std::runtime_error );

  // get the fluid out of the constitutive manager.
  // for the moment it is of type MultiFluidBase.

  ConstitutiveManager & constitutiveManager = this->getGroupByPath< ConstitutiveManager >( "/Problem/domain/Constitutive" );
  MultiFluidBase & baseFluid = constitutiveManager.getGroup< MultiFluidBase >( m_fluidName );

  // depending on logLevel, print some useful info

  if( getLogLevel() > 0 )
  {
    GEOSX_LOG_RANK_0( "Launching PVT Driver" );
    GEOSX_LOG_RANK_0( "  Fluid .................. " << m_fluidName );
    GEOSX_LOG_RANK_0( "  Type ................... " << baseFluid.getCatalogName() );
    GEOSX_LOG_RANK_0( "  No. of Phases .......... " << m_numPhases );
    GEOSX_LOG_RANK_0( "  No. of Components ...... " << m_numComponents );
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
  baseFluid.allocateConstitutiveData( discretization, 1 );   // one quadrature point

  // pass the solid through the ConstitutivePassThru to downcast from the
  // base type to a known model type.  the lambda here then executes the
  // appropriate test driver. note that these calls will move data to device if available.

  constitutiveUpdatePassThru( baseFluid, [&] ( auto & selectedFluid )
  {
    using FLUID_TYPE = TYPEOFREF( selectedFluid );
    runTest< FLUID_TYPE >( selectedFluid, m_table );
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


template< typename FLUID_TYPE >
void PVTDriver::runTest( FLUID_TYPE & fluid, arrayView2d< real64 > const & table )
{
  // get number of phases and components

  localIndex const NC = fluid.numFluidComponents();
  localIndex const NP = fluid.numFluidPhases();

  // prefer output in mass

  fluid.setMassFlag( true );

  // get output data views

  arrayView2d< real64 const, multifluid::USD_FLUID > totalDensity   = fluid.totalDensity();
  arrayView3d< real64 const, multifluid::USD_PHASE > phaseFraction  = fluid.phaseFraction();
  arrayView3d< real64 const, multifluid::USD_PHASE > phaseDensity   = fluid.phaseDensity();
  arrayView3d< real64 const, multifluid::USD_PHASE > phaseViscosity = fluid.phaseViscosity();

  // create kernel wrapper

  typename FLUID_TYPE::KernelWrapper kernelWrapper = fluid.createKernelWrapper();

  // set composition to user specified feed
  // it is more convenient to provide input in molar, so perform molar to mass conversion here

  GEOSX_ASSERT_EQ( NC, m_feed.size() );
  array2d< real64, compflow::LAYOUT_COMP > compositionValues( 1, NC );

  real64 sum = 0.0;
  for( localIndex i = 0; i < NC; ++i )
  {
    compositionValues[0][i] = m_feed[i] * fluid.componentMolarWeights()[i];
    sum += compositionValues[0][i];
  }
  for( localIndex i = 0; i < NC; ++i )
  {
    compositionValues[0][i] /= sum;
  }

  arrayView2d< real64 const, compflow::USD_COMP > const composition = compositionValues;

  // perform fluid update using table (P,T) and save resulting total density, etc.
  // note: column indexing should be kept consistent with output file header below.

  integer numSteps = m_numSteps;
  using ExecPolicy = typename FLUID_TYPE::exec_policy;
  forAll< ExecPolicy >( 1, [=]  GEOSX_HOST_DEVICE ( localIndex const ei )
  {
    for( integer n = 0; n <= numSteps; ++n )
    {
      kernelWrapper.update( ei, 0, table( n, PRES ), table( n, TEMP ), composition[0] );
      table( n, TEMP+1 ) = totalDensity( ei, 0 );

      for( integer p=0; p<NP; ++p )
      {
        table( n, TEMP+2+p ) = phaseFraction( ei, 0, p );
        table( n, TEMP+2+p+NP ) = phaseDensity( ei, 0, p );
        table( n, TEMP+2+p+2*NP ) = phaseViscosity( ei, 0, p );
      }
    }
  } );

}


void PVTDriver::outputResults()
{
  // TODO: improve file path output to grab command line -o directory
  //       for the moment, we just use the specified m_outputFile directly

  FILE * fp = fopen( m_outputFile.c_str(), "w" );

  fprintf( fp, "# column 1 = time\n" );
  fprintf( fp, "# column 2 = pressure\n" );
  fprintf( fp, "# column 3 = temperature\n" );
  fprintf( fp, "# column 4 = density\n" );
  fprintf( fp, "# columns %d-%d = phase fractions\n", 5, 4+m_numPhases );
  fprintf( fp, "# columns %d-%d = phase densities\n", 5+m_numPhases, 4+2*m_numPhases );
  fprintf( fp, "# columns %d-%d = phase viscosities\n", 5+2*m_numPhases, 4+3*m_numPhases );

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


void PVTDriver::compareWithBaseline()
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
                        PVTDriver,
                        string const &, dataRepository::Group * const )

} /* namespace geosx */
