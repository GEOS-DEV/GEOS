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
 * @file TriaxialDriver.cpp
 */

#include "TriaxialDriver.hpp"
#include "constitutiveDrivers/LogLevelsInfo.hpp"
#include "constitutive/solid/SolidBase.hpp"
#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/ConstitutiveManager.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

TriaxialDriver::TriaxialDriver( const string & name,
                                Group * const parent ):
  ConstitutiveDriver( name, parent )
{

  registerWrapper( viewKeyStruct::solidMaterialNameString(), &m_solidMaterialName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Solid material to test" );

  registerWrapper( viewKeyStruct::modeString(), &m_mode ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Test mode [stressControl, strainControl, mixedControl]" );

  registerWrapper( viewKeyStruct::axialFunctionString(), &m_axialFunctionName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Function controlling axial stress or strain (depending on test mode)" );

  registerWrapper( viewKeyStruct::radialFunctionString(), &m_radialFunctionName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Function controlling radial stress or strain (depending on test mode)" );

  registerWrapper( viewKeyStruct::initialStressString(), &m_initialStress ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Initial stress (scalar used to set an isotropic stress state)" );

  addLogLevel< logInfo::LogOutput >();
}

void TriaxialDriver::postInputInitialization()
{
  // initialize table functions

  FunctionManager & functionManager = FunctionManager::getInstance();

  TableFunction & axialFunction = functionManager.getGroup< TableFunction >( m_axialFunctionName );
  TableFunction & radialFunction = functionManager.getGroup< TableFunction >( m_radialFunctionName );

  axialFunction.initializeFunction();
  radialFunction.initializeFunction();

  // determine time increment
  string_array columnNames;
  getColumnNames( columnNames );
  integer const numCols = static_cast< integer >(columnNames.size());

  ArrayOfArraysView< real64 > coordinates = axialFunction.getCoordinates();
  real64 const minTime = coordinates[0][0];
  real64 const maxTime = coordinates[0][coordinates.sizeOfArray( 0 )-1];

  // Allocate the data
  allocateTable( numCols, minTime, maxTime );

  integer const numRows = m_table.size( 0 );

  // initial stress is always isotropic

  m_table( 0, SIG0 ) = m_initialStress;
  m_table( 0, SIG1 ) = m_initialStress;
  m_table( 0, SIG2 ) = m_initialStress;

  // preset certain columns depending on testing mode:
  //   mixedControl .... specified axial strain and radial stress
  //   strainControl ... specified axial and radial strain
  //   stressControl ... specified axial and radial stress

  for( integer n=0; n<numRows; ++n )
  {
    real64 axi = axialFunction.evaluate( &m_table( n, TIME ) );
    real64 rad = radialFunction.evaluate( &m_table( n, TIME ) );

    switch( m_mode )
    {
      case Mode::MixedControl:
        m_table( n, EPS0 ) = axi;
        m_table( n, SIG1 ) = rad;
        m_table( n, SIG2 ) = rad;
        break;

      case Mode::StrainControl:
        m_table( n, EPS0 ) = axi;
        m_table( n, EPS1 ) = rad;
        m_table( n, EPS2 ) = rad;
        break;

      case Mode::StressControl:
        m_table( n, SIG0 ) = axi;
        m_table( n, SIG1 ) = rad;
        m_table( n, SIG2 ) = rad;
        break;
    }
  }

  // double check the initial stress value is consistent with any function values that
  // may overwrite it.

  GEOS_THROW_IF( !isEqual( m_initialStress, m_table( 0, SIG0 ), 1e-6 ),
                 "Initial stress values indicated by initialStress and axialFunction(time=0) appear inconsistent",
                 InputError, getDataContext() );

  GEOS_THROW_IF( !isEqual( m_initialStress, m_table( 0, SIG1 ), 1e-6 ),
                 "Initial stress values indicated by initialStress and radialFunction(time=0) appear inconsistent",
                 InputError, getDataContext() );
}


template< typename SOLID_TYPE >
void TriaxialDriver::runStrainControlTest( SOLID_TYPE & solid, arrayView2d< real64 > const & table )
{
  typename SOLID_TYPE::KernelWrapper updates = solid.createKernelUpdates();
  integer const numRows = m_table.size( 0 );

  forAll< parallelDevicePolicy<> >( 1, [=]  GEOS_HOST_DEVICE ( integer const ei )
  {
    real64 stress[6] = {};
    real64 timeIncrement = 0;
    real64 strainIncrement[6] = {};
    real64 stiffness[6][6] = {{}};

    for( integer n = 1; n <= numRows; ++n )
    {
      strainIncrement[0] = table( n, EPS0 )-table( n-1, EPS0 );
      strainIncrement[1] = table( n, EPS1 )-table( n-1, EPS1 );
      strainIncrement[2] = table( n, EPS2 )-table( n-1, EPS2 );

      timeIncrement = table( n, TIME )-table( n-1, TIME );

      updates.smallStrainUpdate( ei, 0, timeIncrement, strainIncrement, stress, stiffness );
      updates.saveConvergedState ( ei, 0 );

      table( n, SIG0 ) = stress[0];
      table( n, SIG1 ) = stress[1];
      table( n, SIG2 ) = stress[2];

      table( n, ITER ) = 0;
    }
  } );
}


template< typename SOLID_TYPE >
void TriaxialDriver::runMixedControlTest( SOLID_TYPE & solid, arrayView2d< real64 > const & table )
{
  typename SOLID_TYPE::KernelWrapper updates = solid.createKernelUpdates();
  integer const numSteps = m_numSteps;
  integer const maxIter = m_maxIter;
  integer const maxCuts = m_maxCuts;
  real64 const newtonTol = m_newtonTol;

  forAll< parallelDevicePolicy<> >( 1, [=]  GEOS_HOST_DEVICE ( integer const ei )
  {
    real64 stress[6] = {};
    real64 timeIncrement = 0;
    real64 strainIncrement[6] = {};
    real64 deltaStrainIncrement = 0;
    real64 stiffness[6][6] = {{}};

    real64 scale = 0;
    for( integer n = 1; n <= numSteps; ++n )
    {
      scale += fabs( table( n, SIG0 )) + fabs( table( n, SIG1 )) + fabs( table( n, SIG2 ));
    }
    scale = 3 * numSteps / scale;

    for( integer n=1; n<=numSteps; ++n )
    {
      strainIncrement[0] = table( n, EPS0 )-table( n-1, EPS0 );
      strainIncrement[1] = 0;
      strainIncrement[2] = 0;

      timeIncrement = table( n, TIME )-table( n-1, TIME );

      real64 norm, normZero = 1e30;
      integer k = 0;
      integer cuts = 0;

      for(; k < maxIter; ++k )
      {
        updates.smallStrainUpdate( ei, 0, timeIncrement, strainIncrement, stress, stiffness );

        norm = scale*fabs( stress[1]-table( n, SIG1 ) );

        if( k == 0 )
        {
          normZero = norm;
        }

        if( norm < newtonTol ) // success
        {
          break;
        }
        else if( k > 0 && norm > normZero && cuts < maxCuts ) // backtrack by half delta
        {
          cuts++;
          deltaStrainIncrement *= 0.5;
          strainIncrement[1] += deltaStrainIncrement;
          strainIncrement[2]  = strainIncrement[1];
        }
        else // newton update
        {
          deltaStrainIncrement  = (stress[1]-table( n, SIG1 )) / (stiffness[1][1]+stiffness[1][2]);
          strainIncrement[1]   -= deltaStrainIncrement;
          strainIncrement[2]    = strainIncrement[1];
        }
      }

      updates.saveConvergedState ( ei, 0 );

      table( n, SIG0 ) = stress[0];
      table( n, EPS1 ) = table( n-1, EPS1 )+strainIncrement[1];
      table( n, EPS2 ) = table( n, EPS1 );

      table( n, ITER ) = k;
      table( n, NORM ) = norm;

      if( norm > newtonTol )
      {
        break;
      }
    }
  } );
}


template< typename SOLID_TYPE >
void TriaxialDriver::runStressControlTest( SOLID_TYPE & solid, arrayView2d< real64 > const & table )
{
  typename SOLID_TYPE::KernelWrapper updates = solid.createKernelUpdates();
  integer const numSteps = m_numSteps;
  integer const maxIter = m_maxIter;
  integer const maxCuts = m_maxCuts;
  real64 const newtonTol = m_newtonTol;

  forAll< parallelDevicePolicy<> >( 1, [=]  GEOS_HOST_DEVICE ( integer const ei )
  {
    real64 stress[6] = {};
    real64 timeIncrement = 0;
    real64 strainIncrement[6] = {};
    real64 deltaStrainIncrement[6] = {};
    real64 stiffness[6][6] = {{}};

    real64 resid[2] = {};
    real64 jacobian[2][2] = {{}};

    real64 scale = 0;
    for( integer n = 1; n <= numSteps; ++n )
    {
      scale += fabs( table( n, SIG0 )) + fabs( table( n, SIG1 )) + fabs( table( n, SIG2 ));
    }
    scale = 3 * numSteps / scale;

    for( integer n = 1; n <= numSteps; ++n )
    {
      //   std::cout<<"time step="<<n<<std::endl;
      strainIncrement[0] = 0;
      strainIncrement[1] = 0;
      strainIncrement[2] = 0;

      timeIncrement = table( n, TIME )-table( n-1, TIME );

      real64 norm, normZero = 1e30, det;
      integer k = 0;
      integer cuts = 0;

      for(; k < maxIter; ++k )
      {
        updates.smallStrainUpdate( ei, 0, timeIncrement, strainIncrement, stress, stiffness );

        resid[0] = scale * (stress[0]-table( n, SIG0 ));
        resid[1] = scale * (stress[1]-table( n, SIG1 ));

        norm = sqrt( resid[0]*resid[0] + resid[1]*resid[1] );
        //  std::cout<<"k= "<<k<<std::endl;
        // std::cout<<"norm ="<<norm<<std::endl;

        if( k == 0 )
        {
          normZero = norm;
        }

        if( norm < newtonTol ) // success
        {
          break;
        }
        else if( k > 0 && norm > normZero && cuts < maxCuts ) // backtrack by half delta
        {
          cuts++;
          deltaStrainIncrement[0] *= 0.5;
          deltaStrainIncrement[1] *= 0.5;
          strainIncrement[0] += deltaStrainIncrement[0];
          strainIncrement[1] += deltaStrainIncrement[1];
          strainIncrement[2]  = strainIncrement[1];
          //  std::cout<<"k="<<k<<" , cuts="<<cuts<<std::endl;
        }
        else // newton update
        {
          cuts = 0;
          jacobian[0][0] = scale*stiffness[0][0];
          jacobian[1][0] = scale*stiffness[1][0];
          jacobian[0][1] = scale*(stiffness[0][1] + stiffness[0][2]);
          jacobian[1][1] = scale*(stiffness[1][1] + stiffness[1][2]);

          det = jacobian[0][0]*jacobian[1][1]-jacobian[0][1]*jacobian[1][0];

          deltaStrainIncrement[0] = (jacobian[1][1]*resid[0]-jacobian[0][1]*resid[1] ) / det;
          deltaStrainIncrement[1] = (jacobian[0][0]*resid[1]-jacobian[1][0]*resid[0] ) / det;

          strainIncrement[0] -= deltaStrainIncrement[0];
          strainIncrement[1] -= deltaStrainIncrement[1];
          strainIncrement[2]  = strainIncrement[1];
        }
      }

      updates.saveConvergedState ( ei, 0 );

      table( n, EPS0 ) = table( n-1, EPS0 )+strainIncrement[0];
      table( n, EPS1 ) = table( n-1, EPS1 )+strainIncrement[1];
      table( n, EPS2 ) = table( n, EPS1 );

      table( n, ITER ) = k;
      table( n, NORM ) = norm;

      if( norm > newtonTol )
      {
        break;
      }
    }
  } );
}


bool TriaxialDriver::execute()
{
  // get the solid out of the constitutive manager.
  // for the moment it is of type SolidBase.

  SolidBase & baseSolid = getSolid();

  // depending on logLevel, print some useful info

  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "Launching Triaxial Driver" );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Material .......... " << m_solidMaterialName );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Type .............. " << baseSolid.getCatalogName() );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Mode .............. " << m_mode );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Axial Control ..... " << m_axialFunctionName );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Radial Control .... " << m_radialFunctionName );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Initial Stress .... " << m_initialStress );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Steps ............. " << m_numSteps );
  GEOS_LOG_LEVEL_RANK_0( logInfo::LogOutput, "  Output ............ " << (m_outputFile.empty() ? "<none>" : m_outputFile) );

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

  ConstitutivePassThruTriaxialDriver< SolidBase >::execute( baseSolid, [&]( auto & selectedSolid )
  {
    using SOLID_TYPE = TYPEOFREF( selectedSolid );

    switch( m_mode )
    {
      case Mode::MixedControl:
        runMixedControlTest< SOLID_TYPE >( selectedSolid, m_table );
        break;

      case Mode::StrainControl:
        runStrainControlTest< SOLID_TYPE >( selectedSolid, m_table );
        break;

      case Mode::StressControl:
        runStressControlTest< SOLID_TYPE >( selectedSolid, m_table );
        break;
    }
  } );

  return false;
}

void TriaxialDriver::getColumnNames( string_array & columnNames ) const
{
  columnNames.emplace_back( "time" );
  columnNames.emplace_back( "strain,axial" );
  columnNames.emplace_back( "strain,radial_1" );
  columnNames.emplace_back( "strain,radial_2" );
  columnNames.emplace_back( "stress,axial" );
  columnNames.emplace_back( "stress,radial_1" );
  columnNames.emplace_back( "stress,radial_2" );
  columnNames.emplace_back( "newton_iter" );
  columnNames.emplace_back( "residual_norm" );
}

bool TriaxialDriver::validateResults()
{
  for( integer n=0; n<m_numSteps; ++n )
  {
    if( m_table( n, NORM ) > m_newtonTol )
    {
      GEOS_LOG_RANK_0( "WARNING: Material driver failed to converge at loadstep " << n << "." );
      GEOS_LOG_RANK_0( "         This usually indicates the material has completely failed and/or the loading state is inadmissible." );
      GEOS_LOG_RANK_0( "         In rare cases, it may indicate a problem in the material model implementation." );

      for( integer col=EPS0; col<ITER; ++col )
      {
        m_table( n, col ) = 0;
      }
    }
  }
  return true;
}

#if 0
void TriaxialDriver::outputResults()
{
  // TODO: improve file path output to grab command line -o directory
  //       for the moment, we just use the specified m_outputFile directly

  FILE * fp = fopen( m_outputFile.c_str(), "w" );

  /*
     string const outputDir = OutputBase::getOutputDirectory();
     string const outputPath = joinPath( outputDir, m_outputFile );
     FILE * fp = fopen( outputPath.c_str(), "w" );
   */

  fprintf( fp, "# column 1 = time\n" );
  fprintf( fp, "# column 2 = axial_strain\n" );
  fprintf( fp, "# column 3 = radial_strain_1\n" );
  fprintf( fp, "# column 4 = radial_strain_2\n" );
  fprintf( fp, "# column 5 = axial_stress\n" );
  fprintf( fp, "# column 6 = radial_stress_1\n" );
  fprintf( fp, "# column 7 = radial_stress_2\n" );
  fprintf( fp, "# column 8 = newton_iter\n" );
  fprintf( fp, "# column 9 = residual_norm\n" );

  for( integer n=0; n<=m_numSteps; ++n )
  {
    for( integer col=0; col<m_numColumns; ++col )
    {
      fprintf( fp, "%.4e ", m_table( n, col ) );
    }
    fprintf( fp, "\n" );
  }
  fclose( fp );
}
#endif

SolidBase & TriaxialDriver::getSolid()
{
  return getConstitutiveManager().getGroup< SolidBase >( m_solidMaterialName );
}

SolidBase const & TriaxialDriver::getSolid() const
{
  return getConstitutiveManager().getGroup< SolidBase >( m_solidMaterialName );
}


REGISTER_CATALOG_ENTRY( TaskBase,
                        TriaxialDriver,
                        string const &, dataRepository::Group * const )

} /* namespace geos */
