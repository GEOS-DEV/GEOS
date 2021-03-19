/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TriaxialDriver.cpp
 */

#include "TriaxialDriver.hpp"

#include "managers/ProblemManager.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;


TriaxialDriver::TriaxialDriver( const string & name,
                                Group * const parent ):
  TaskBase( name, parent )
{
  enableLogLevelInput();

  registerWrapper( viewKeyStruct::solidMaterialNameString(), &m_solidMaterialName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Solid material to test" );

  registerWrapper( viewKeyStruct::modeString(), &m_mode ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Test mode [triaxial, volumetric, oedometer]" );

  registerWrapper( viewKeyStruct::strainFunctionString(), &m_strainFunctionName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Function controlling strain loading (role depends on test mode)" );

  registerWrapper( viewKeyStruct::stressFunctionString(), &m_stressFunctionName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Function controlling stress loading (role depends on test mode)" );

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


TriaxialDriver::~TriaxialDriver()
{}


void TriaxialDriver::postProcessInput()
{

  GEOSX_THROW_IF( m_mode != "triaxial" && m_mode != "volumetric" && m_mode != "oedometer",
                  "Test mode \'" << m_mode << "\' not recognized.",
                  InputError );

  // initialize table functions

  FunctionManager & functionManager = getGlobalState().getFunctionManager();

  TableFunction & strainFunction = functionManager.getGroup< TableFunction >( m_strainFunctionName );
  TableFunction & stressFunction = functionManager.getGroup< TableFunction >( m_stressFunctionName );

  strainFunction.initializeFunction();
  stressFunction.initializeFunction();

  // determine time increment

  ArrayOfArraysView< real64 > coordinates = strainFunction.getCoordinates();
  real64 const minTime = coordinates[0][0];
  real64 const maxTime = coordinates[0][coordinates.sizeOfArray( 0 )-1];
  real64 const dt = (maxTime-minTime) / m_numSteps;

  // resize data arrays

  localIndex const length = m_numSteps+1;
  m_table.resize( length, m_numColumns );

  // set time and axial strain columns

  for( localIndex n=0; n<length; ++n )
  {
    m_table( n, TIME ) = minTime + n*dt;
    m_table( n, EPS0 ) = strainFunction.evaluate( &m_table( n, TIME ) );
  }

  // initial stress is always isotropic, using stressFunction at t=tmin

  m_table( 0, SIG0 ) = stressFunction.evaluate( &m_table( 0, TIME ) ); // init stress
  m_table( 0, SIG1 ) = stressFunction.evaluate( &m_table( 0, TIME ) ); // init stress
  m_table( 0, SIG2 ) = stressFunction.evaluate( &m_table( 0, TIME ) ); // init stress

  // other columns depend on testing mode:
  //  * triaxial ...... specified axial strain and radial stress
  //  * oedometeter ... specified axial strain and zero radial strain (1D compression)
  //  * volumetric .... axial strain equals radial strain (isotropic compression)

  for( localIndex n=0; n<length; ++n )
  {
    if( m_mode == "triaxial" )
    {
      m_table( n, SIG1 ) = stressFunction.evaluate( &m_table( n, TIME ) );
      m_table( n, SIG2 ) = m_table( n, SIG1 );
    }
    else if( m_mode == "volumetric" )
    {
      m_table( n, EPS1 ) = m_table( n, EPS0 );
      m_table( n, EPS2 ) = m_table( n, EPS0 );
    }
  }
}


bool TriaxialDriver::execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                              real64 const GEOSX_UNUSED_PARAM( dt ),
                              integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                              integer const GEOSX_UNUSED_PARAM( eventCounter ),
                              real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                              DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  // this code only makes sense in serial

  GEOSX_THROW_IF( MpiWrapper::commRank() > 0, "Triaxial Driver should only be run in serial", std::runtime_error );

  // get the solid out of the constitutive manager.
  // for the moment it is of type SolidBase.

  ConstitutiveManager & constitutiveManager = getGlobalState()
                                                .getProblemManager()
                                                .getDomainPartition()
                                                .getConstitutiveManager();

  SolidBase & baseSolid = constitutiveManager.getGroup< SolidBase >( m_solidMaterialName );

  // depending on logLevel, print some useful info

  if( this->getLogLevel() > 0 )
  {
    GEOSX_LOG_RANK_0( "Launching Triaxial Driver" );
    GEOSX_LOG_RANK_0( "  Material .......... " << m_solidMaterialName );
    GEOSX_LOG_RANK_0( "  Type .............. " << baseSolid.getCatalogName() );
    GEOSX_LOG_RANK_0( "  Mode .............. " << m_mode );
    GEOSX_LOG_RANK_0( "  Strain Function ... " << m_strainFunctionName );
    GEOSX_LOG_RANK_0( "  Stress Function ... " << m_stressFunctionName );
    GEOSX_LOG_RANK_0( "  Steps ............. " << m_numSteps );
    GEOSX_LOG_RANK_0( "  Output ............ " << m_outputFile );
  }

  // create a dummy discretization with one quadrature point for
  // storing constitutive data

  conduit::Node node;
  dataRepository::Group rootGroup( "root", node );
  dataRepository::Group discretization( "discretization", &rootGroup );

  discretization.resize( 1 );   // one element
  baseSolid.allocateConstitutiveData( discretization, 1 );   // one quadrature point

  // set the initial stress state using the data table

  arrayView3d< real64, solid::STRESS_USD > stressArray = baseSolid.getStress();

  stressArray( 0, 0, 0 ) = m_table( 0, SIG0 );
  stressArray( 0, 0, 1 ) = m_table( 0, SIG1 );
  stressArray( 0, 0, 2 ) = m_table( 0, SIG2 );

  baseSolid.saveConvergedState();

  // pass the solid through the ConstitutivePassThru to downcast from the
  // base type to a known model type.  the lambda here then executes the
  // appropriate test driver. note that these calls will move data to device if available.

  ConstitutivePassThru< SolidBase >::execute( baseSolid, [&]( auto & selectedSolid )
  {
    using SOLID_TYPE = TYPEOFREF( selectedSolid );

    if( m_mode == "triaxial" )
    {
      runMixedControlTest< SOLID_TYPE >( selectedSolid, m_table );
    }
    else if( m_mode == "volumetric" || m_mode == "oedometer" )
    {
      runStrainControlTest< SOLID_TYPE >( selectedSolid, m_table );
    }
  } );

  // move table back to host for output
  m_table.move( LvArray::MemorySpace::CPU );

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


void TriaxialDriver::outputResults()
{
  FILE * fp = fopen( m_outputFile.c_str(), "w" );

  fprintf( fp, "# column 1 = time\n" );
  fprintf( fp, "# column 2 = axial_strain\n" );
  fprintf( fp, "# column 3 = radial_strain_1\n" );
  fprintf( fp, "# column 4 = radial_strain_2\n" );
  fprintf( fp, "# column 5 = axial_stress\n" );
  fprintf( fp, "# column 6 = radial_strain_1\n" );
  fprintf( fp, "# column 7 = radial_stress_2\n" );
  fprintf( fp, "# column 8 = newton_iter\n" );
  fprintf( fp, "# column 9 = residual_norm\n" );

  for( localIndex n=0; n<=m_numSteps; ++n )
  {
    for( localIndex col=0; col<m_numColumns; ++col )
    {
      fprintf( fp, "%.3e ", m_table( n, col ) );
    }
    fprintf( fp, "\n" );
  }
  fclose( fp );
}


void TriaxialDriver::compareWithBaseline()
{}

REGISTER_CATALOG_ENTRY( TaskBase,
                        TriaxialDriver,
                        string const &, dataRepository::Group * const )

} /* namespace geosx */
