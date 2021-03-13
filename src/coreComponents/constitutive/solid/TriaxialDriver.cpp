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

  registerWrapper( viewKeyStruct::outputString(), &m_outputFileName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Output file" );
}

TriaxialDriver::~TriaxialDriver()
{}

void TriaxialDriver::postProcessInput()
{
  // initialize table functions

  FunctionManager & functionManager = getGlobalState().getFunctionManager();

  TableFunction & strainFunction = functionManager.getGroup< TableFunction >( m_strainFunctionName );
  TableFunction & stressFunction = functionManager.getGroup< TableFunction >( m_stressFunctionName );

  strainFunction.initializeFunction(); //TODO: why not automatically initialized?
  stressFunction.initializeFunction(); //TODO: why not automatically initialized?

  // determine time increment

  ArrayOfArraysView< real64 > coordinates = strainFunction.getCoordinates();
  real64 const minTime = coordinates[0][0];
  real64 const maxTime = coordinates[0][coordinates.sizeOfArray( 0 )-1];
  real64 const dt = (maxTime-minTime) / m_numSteps;

  // resize data arrays

  localIndex const length = m_numSteps+1;

  m_time.resize( length );
  m_axialStrain.resize( length );
  m_axialStress.resize( length );
  m_radialStrain.resize( length );
  m_radialStress.resize( length );

  // set time array

  for( localIndex n=0; n<length; ++n )
  {
    m_time[n] = minTime + n*dt;
  }

  // set other arrays based on testing mode
  // initial stress is always isotropic, using stressFunction at t=tmin

  if( m_mode == "triaxial" ) // specified axial strain and radial stress
  {
    for( localIndex n=0; n<length; ++n )
    {
      m_axialStrain[n] = strainFunction.evaluate( &m_time[n] );
      m_radialStress[n] = stressFunction.evaluate( &m_time[n] );
    }
    m_axialStress[0] = stressFunction.evaluate( &m_time[0] ); // init stress
  }
  else if( m_mode == "volumetric" ) // specified axial strain = radial strain
  {
    for( localIndex n=0; n<length; ++n )
    {
      m_axialStrain[n] = strainFunction.evaluate( &m_time[n] );
      m_radialStrain[n] = m_axialStrain[n];
    }
    m_axialStress[0] = stressFunction.evaluate( &m_time[0] ); // init stress
    m_radialStress[0] = stressFunction.evaluate( &m_time[0] ); // init stress

  }
  else if( m_mode == "oedometer" ) // specified axial strain, zero radial strain
  {
    for( localIndex n=0; n<length; ++n )
    {
      m_axialStrain[n] = strainFunction.evaluate( &m_time[n] );
    }
    m_axialStress[0] = stressFunction.evaluate( &m_time[0] ); // init stress
    m_radialStress[0] = stressFunction.evaluate( &m_time[0] ); // init stress
  }
  else // unrecognized option
  {
    GEOSX_THROW( "Test mode \'" << m_mode << "\' not recognized.", InputError );
  }
}


void TriaxialDriver::outputResults()
{
  FILE * fp = fopen( m_outputFileName.c_str(), "w" );

  fprintf( fp, "# column 1 = time\n" );
  fprintf( fp, "# column 2 = axial strain\n" );
  fprintf( fp, "# column 3 = radial strain\n" );
  fprintf( fp, "# column 4 = axial stress\n" );
  fprintf( fp, "# column 5 = radial stress\n" );

  for( localIndex n=0; n<m_time.size(); ++n )
  {
    fprintf( fp,
             "%.4e %.4e %.4e %.4e %.4e\n",
             m_time[n],
             m_axialStrain[n],
             m_radialStrain[n],
             m_axialStress[n],
             m_radialStress[n] );
  }
  fclose( fp );
}


REGISTER_CATALOG_ENTRY( TaskBase,
                        TriaxialDriver,
                        string const &, dataRepository::Group * const )

} /* namespace geosx */
