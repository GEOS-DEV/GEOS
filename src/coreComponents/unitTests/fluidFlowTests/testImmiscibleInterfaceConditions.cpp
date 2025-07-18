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

#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "constitutive/fluid/twophaseimmisciblefluid/TwoPhaseImmiscibleFluid.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/ImmiscibleMultiphaseFlow.hpp"
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"
#include "constitutive/unitTests/FluidModelTest.hpp"
#include "constitutive/unitTests/FluidModelTest_impl.hpp"
#include "common/initializeEnvironment.hpp"
#include "unitTests/constitutiveTests/constitutiveTestHelpers.hpp"
#include "functions/FunctionManager.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureFields.hpp"

#include <conduit.hpp>

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::constitutive;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

TwoPhaseImmiscibleFluid * makeTwoPhaseImmiscibleFluid( TwoPhaseImmiscibleFluid & fluid )
{
  
  FunctionManager & functionManager = FunctionManager::getInstance();

  // 1D table with linear interpolation
  localIndex constexpr Naxis = 6;
  localIndex constexpr NaxisSingle = 1;

  real64_array densityCoordPhase0;
  // fill( densityCoordPhase0, Feed< Naxis >{ 0.22, 0.3, 0.5, 0.6, 0.8, 1.0 } );
  for (auto v : {0.0})
    densityCoordPhase0.emplace_back(v);
  real64_array densityValuesPhase0;
  // fill( densityValuesPhase0, Feed< Naxis >{ 0.00603, 0.04224, 0.04224, 0.22423, 0.31311, 0.40203 } );
  for (auto v : {1000.0})
  densityValuesPhase0.emplace_back(v);

  real64_array densityCoordPhase1;
  // fill( densityCoordPhase1, Feed< Naxis >{ 1.22, 1.3, 1.5, 1.6, 1.8, 2.0 } );
  for (auto v : {0.0})
    densityCoordPhase1.emplace_back(v);
  real64_array densityValuesPhase1;
  // fill( densityValuesPhase1, Feed< Naxis >{ 0.00603, 0.04224, 0.04224, 0.22423, 0.31311, 0.40203 } );
  for (auto v : {100.0})
  densityValuesPhase1.emplace_back(v);

  real64_array viscosityCoordPhase0;
  // fill( viscosityCoordPhase0, Feed< Naxis >{ 0.22, 0.3, 0.5, 0.6, 0.8, 1.0 } );
  for (auto v : {0.0})
  viscosityCoordPhase0.emplace_back(v);
  real64_array viscosityValuesPhase0;
  // fill( viscosityValuesPhase0, Feed< Naxis >{ 40203, 31311, 22423, 15011, 4224, 603 } );
  for (auto v : {0.001})
  viscosityValuesPhase0.emplace_back(v);

  real64_array viscosityCoordPhase1;
  // fill( viscosityCoordPhase1, Feed< NaxisSingle >{ 0.22 } );
  for (auto v : {0.0})
  viscosityCoordPhase1.emplace_back(v);

  real64_array viscosityValuesPhase1;
  // fill( viscosityValuesPhase1, Feed< NaxisSingle >{ 45 } );
  for (auto v : {0.001})
  viscosityValuesPhase1.emplace_back(v);

  TableFunction & table_density0 = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "densityTablePhase0" ) );
  array1d<real64_array> coords_density0;
  coords_density0.emplace_back(densityCoordPhase0);
  table_density0.setTableCoordinates( coords_density0, { units::Dimensionless } );
  table_density0.setTableValues( densityValuesPhase0, units::Dimensionless );
  table_density0.reInitializeFunction();

  table_density0.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  TableFunction & table_density1 = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "densityTablePhase1" ) );
  array1d<real64_array> coords_density1;
  coords_density1.emplace_back(densityCoordPhase1);
  table_density1.setTableCoordinates( coords_density1, { units::Dimensionless } );
  table_density1.setTableValues( densityValuesPhase1, units::Dimensionless );
  table_density1.reInitializeFunction();

  table_density1.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  TableFunction & table_viscosity0 = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "viscosityTablePhase0" ) );
  array1d<real64_array> coords_viscosity0;
  coords_viscosity0.emplace_back(viscosityCoordPhase0);
  table_viscosity0.setTableCoordinates( coords_viscosity0, { units::Dimensionless } );
  table_viscosity0.setTableValues( viscosityValuesPhase0, units::Dimensionless );
  table_viscosity0.reInitializeFunction();

  table_viscosity0.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  TableFunction & table_viscosity1 = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "viscosityTablePhase1" ) );
  array1d<real64_array> coords_viscosity1;
  coords_viscosity1.emplace_back(viscosityCoordPhase1);
  table_viscosity1.setTableCoordinates( coords_viscosity1, { units::Dimensionless } );
  table_viscosity1.setTableValues( viscosityValuesPhase1, units::Dimensionless );
  table_viscosity1.reInitializeFunction();

  table_viscosity1.setInterpolationMethod( TableFunction::InterpolationType::Linear );


  // createTable( "densityTablePhase0", densityCoordPhase0, densityValuesPhase0 );
  // createTable( "densityTablePhase1", densityCoordPhase1, densityValuesPhase1 );
  // createTable( "viscosityTablePhase0", viscosityCoordPhase0, viscosityValuesPhase0 );
  // createTable( "viscosityTablePhase1", viscosityCoordPhase1, viscosityValuesPhase1 );

  // 2) Set up the constitutive model

  string_array & phaseNames = fluid.getReference< string_array >( TwoPhaseImmiscibleFluid::viewKeyStruct::phaseNamesString() );
  phaseNames.emplace_back( "water" );
  phaseNames.emplace_back( "gas" );

  string_array & densityTableNames = fluid.getReference< string_array >( TwoPhaseImmiscibleFluid::viewKeyStruct::densityTableNamesString() );
  densityTableNames.emplace_back( "densityTablePhase0" );
  densityTableNames.emplace_back( "densityTablePhase1" );

  string_array & viscosityTableNames = fluid.getReference< string_array >( TwoPhaseImmiscibleFluid::viewKeyStruct::viscosityTableNamesString() );
  viscosityTableNames.emplace_back( "viscosityTablePhase0" );
  viscosityTableNames.emplace_back( "viscosityTablePhase1" );

  fluid.postInputInitializationRecursive();
  fluid.initialize(); // to test all the checks

  return &fluid;
}

CapillaryPressureBase & makeJFunctionCapPressureTwoPhase( string const & name, Group & parent )
{
  FunctionManager & functionManager = FunctionManager::getInstance();

  // 1) First, define the tables

  // 1D table, various interpolation methods
  localIndex const Naxis = 12;

  // Setup table
  array1d< real64_array > coordinates;
  coordinates.resize( 1 );
  coordinates[0].resize( Naxis );

  coordinates[0][0] = 0.0;
  coordinates[0][1] = 0.05;
  coordinates[0][2] = 0.15;
  coordinates[0][3] = 0.25;
  coordinates[0][4] = 0.35;
  coordinates[0][5] = 0.45;
  coordinates[0][6] = 0.55;
  coordinates[0][7] = 0.65;
  coordinates[0][8] = 0.75;
  coordinates[0][9] = 0.85;
  coordinates[0][10] = 0.95;
  coordinates[0][11] = 1.0;

  real64_array values( Naxis );
  values[0] = 4.331729359;
  values[1] = 3.523266264;
  values[2] = 2.677103439;
  values[3] = 2.356150157;
  values[4] = 2.166062360;
  values[5] = 2.034158727;
  values[6] = 1.934627222;
  values[7] = 1.855494313;
  values[8] = 1.790286970;
  values[9] = 1.735134860;
  values[10] = 1.687551617;
  values[11] = 1.666049754;


  // localIndex const Naxis = 2;

  // // Setup table
  // array1d< real64_array > coordinates;
  // coordinates.resize( 1 );
  // coordinates[0].resize( Naxis );

  // coordinates[0][0] = 0.0;
  // coordinates[0][1] = 1.0;

  // real64_array values( Naxis );
  // values[0] = 4.331729359;
  // values[1] = 1.666049754;


  TableFunction & table_w = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "water_jFunction" ) );
  table_w.setTableCoordinates( coordinates, { units::Dimensionless } );
  table_w.setTableValues( values, units::Dimensionless );
  table_w.reInitializeFunction();

  table_w.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  // 2) Then set up the constitutive model

  JFunctionCapillaryPressure & capPressure = parent.registerGroup< JFunctionCapillaryPressure >( name );

  string_array & phaseNames = capPressure.getReference< string_array >( CapillaryPressureBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "water"; phaseNames[1] = "gas";

  auto & waterTableName = capPressure.getReference< string >( JFunctionCapillaryPressure::viewKeyStruct::wettingNonWettingJFuncTableNameString() );
  waterTableName = "water_jFunction";

  auto & surfaceTension = capPressure.getReference< real64 >( JFunctionCapillaryPressure::viewKeyStruct::wettingNonWettingSurfaceTensionString() );
  //surfaceTension = 23.86955676433857e-3;
  surfaceTension = 0.02;

  auto & permeabilityDirection =
    capPressure.getReference< JFunctionCapillaryPressure::PermeabilityDirection >( JFunctionCapillaryPressure::viewKeyStruct::permeabilityDirectionString() );
  permeabilityDirection = JFunctionCapillaryPressure::PermeabilityDirection::XY;

  capPressure.postInputInitializationRecursive();
  capPressure.initialize(); // to test all the checks

  return capPressure;
}

CapillaryPressureBase & makeTableCapPressureTwoPhase( string const & name, Group & parent )
{
  FunctionManager & functionManager = FunctionManager::getInstance();

  // 1) First, define the tables

  // 1D table, various interpolation methods
  localIndex Naxis = 12;

  // Setup table
  array1d< real64_array > coordinates;
  coordinates.resize( 1 );
  coordinates[0].resize( Naxis );
  coordinates[0][0] = 0.0;
  coordinates[0][1] = 0.05;
  coordinates[0][2] = 0.15;
  coordinates[0][3] = 0.25;
  coordinates[0][4] = 0.35;
  coordinates[0][5] = 0.45;
  coordinates[0][6] = 0.55;
  coordinates[0][7] = 0.65;
  coordinates[0][8] = 0.75;
  coordinates[0][9] = 0.85;
  coordinates[0][10] = 0.95;
  coordinates[0][11] = 1.0;

  real64_array values( Naxis );
  values[0] = 129999.999994362;
  values[1] = 105737.126293106;
  values[2] = 80342.8419025813;
  values[3] = 70710.6781175791;
  values[4] = 65005.9326081244;
  values[5] = 61047.3583602598;
  values[6] = 58060.3075597394;
  values[7] = 55685.4412380058;
  values[8] = 53728.4965890932;
  values[9] = 52073.3206292211;
  values[10] = 50645.2947584729;
  values[11] = 50000.0000139914;


  TableFunction & table_w = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "water_pc" ) );
  table_w.setTableCoordinates( coordinates, { units::Dimensionless } );
  table_w.setTableValues( values, units::Pressure );
  table_w.reInitializeFunction();

  table_w.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  // 2) Then set up the constitutive model

  TableCapillaryPressure & capPressure = parent.registerGroup< TableCapillaryPressure >( name );

  string_array & phaseNames = capPressure.getReference< string_array >( CapillaryPressureBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "water"; phaseNames[1] = "gas";

  auto & waterTableName = capPressure.getReference< string >( TableCapillaryPressure::viewKeyStruct::wettingNonWettingCapPresTableNameString() );
  waterTableName = "water_pc";

  capPressure.postInputInitializationRecursive();
  capPressure.initialize(); // to test all the checks
  return capPressure;
}

CapillaryPressureBase & makeBrooksCoreyCapPressureTwoPhase1( string const & name, Group & parent )
{
  BrooksCoreyCapillaryPressure & capPressure = parent.registerGroup< BrooksCoreyCapillaryPressure >( name );

  string_array & phaseNames = capPressure.getReference< string_array >( CapillaryPressureBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "water"; phaseNames[1] = "gas";

  array1d< real64 > & phaseMinSat = capPressure.getReference< array1d< real64 > >( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseMinVolumeFractionString() );
  phaseMinSat.resize( 2 );
  phaseMinSat[0] = 0.0; phaseMinSat[1] = 0.0;

  array1d< real64 > & phaseCapPressureExpInv =
    capPressure.getReference< array1d< real64 > >( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseCapPressureExponentInvString() );
  phaseCapPressureExpInv.resize( 2 );
  phaseCapPressureExpInv[0] = 4; phaseCapPressureExpInv[1] = 1;

  array1d< real64 > & phaseEntryPressure = capPressure.getReference< array1d< real64 > >( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseEntryPressureString() );
  phaseEntryPressure.resize( 2 );
  phaseEntryPressure[0] = 1.5e5; phaseEntryPressure[1] = 0;

  real64 & capPressureEpsilon = capPressure.getReference< real64 >( BrooksCoreyCapillaryPressure::viewKeyStruct::capPressureEpsilonString() );
  capPressureEpsilon = 1.5e-8;

  capPressure.postInputInitializationRecursive();
  return capPressure;
}

CapillaryPressureBase & makeBrooksCoreyCapPressureTwoPhase2( string const & name, Group & parent )
{
  BrooksCoreyCapillaryPressure & capPressure = parent.registerGroup< BrooksCoreyCapillaryPressure >( name );

  string_array & phaseNames = capPressure.getReference< string_array >( CapillaryPressureBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "water"; phaseNames[1] = "gas";

  array1d< real64 > & phaseMinSat = capPressure.getReference< array1d< real64 > >( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseMinVolumeFractionString() );
  phaseMinSat.resize( 2 );
  phaseMinSat[0] = 0.0; phaseMinSat[1] = 0.0;

  array1d< real64 > & phaseCapPressureExpInv =
    capPressure.getReference< array1d< real64 > >( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseCapPressureExponentInvString() );
  phaseCapPressureExpInv.resize( 2 );
  phaseCapPressureExpInv[0] = 4; phaseCapPressureExpInv[1] = 1;

  array1d< real64 > & phaseEntryPressure = capPressure.getReference< array1d< real64 > >( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseEntryPressureString() );
  phaseEntryPressure.resize( 2 );
  phaseEntryPressure[0] = 1.0e5; phaseEntryPressure[1] = 0;

  real64 & capPressureEpsilon = capPressure.getReference< real64 >( BrooksCoreyCapillaryPressure::viewKeyStruct::capPressureEpsilonString() );
  capPressureEpsilon = 1e-8;

  capPressure.postInputInitializationRecursive();
  return capPressure;
}

RelativePermeabilityBase & makeBrooksCoreyRelPerm( string const & name, Group & parent )
{
  BrooksCoreyRelativePermeability & relPerm = parent.registerGroup< BrooksCoreyRelativePermeability >( name );

  string_array & phaseNames = relPerm.getReference< string_array >( RelativePermeabilityBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "water"; phaseNames[1] = "gas";

  array1d< real64 > & phaseMinSat = relPerm.getReference< array1d< real64 > >( BrooksCoreyRelativePermeability::viewKeyStruct::phaseMinVolumeFractionString() );
  phaseMinSat.resize( 2 );
  phaseMinSat[0] = 0.0; phaseMinSat[1] = 0.0;

  array1d< real64 > & phaseRelPermExp = relPerm.getReference< array1d< real64 > >( BrooksCoreyRelativePermeability::viewKeyStruct::phaseRelPermExponentString() );
  phaseRelPermExp.resize( 2 );
  phaseRelPermExp[0] = 1.0; phaseRelPermExp[1] = 1.0;

  array1d< real64 > & phaseRelPermMaxVal = relPerm.getReference< array1d< real64 > >( BrooksCoreyRelativePermeability::viewKeyStruct::phaseRelPermMaxValueString() );
  phaseRelPermMaxVal.resize( 2 );
  phaseRelPermMaxVal[0] = 1.0; phaseRelPermMaxVal[1] = 1.0;

  relPerm.postInputInitializationRecursive();
  return relPerm;
}

class ImmiscibleInterfaceConditionsTest : public FluidModelTest< TwoPhaseImmiscibleFluid, 2 >
{
public:
  ImmiscibleInterfaceConditionsTest(): state( std::make_unique< CommandLineOptions >( g_commandLineOptions )), 
  m_parent( "TestParentGroup", m_node )

  {}

protected:

  static real64 constexpr time = 0.0;
  static real64 constexpr dt = 1e4;
  static real64 constexpr eps = std::numeric_limits< real64 >::epsilon();

  GeosxState state;
  ImmiscibleMultiphaseFlow *solver;
  conduit::Node m_node;
  dataRepository::Group m_parent;
};

real64 constexpr ImmiscibleInterfaceConditionsTest::time;
real64 constexpr ImmiscibleInterfaceConditionsTest::dt;
real64 constexpr ImmiscibleInterfaceConditionsTest::eps;



TEST_F( ImmiscibleInterfaceConditionsTest, LocalNonlinearSolverConvergence )
{

  // using Base = FluidModelTest< TwoPhaseImmiscibleFluid, 2 >;
  createFluid( "fluid", [this]( TwoPhaseImmiscibleFluid & fluid ){
    makeTwoPhaseImmiscibleFluid( fluid );
 
  // getting constitutive models:
    RelativePermeabilityBase & relPerm = makeBrooksCoreyRelPerm( "relPerm" , this->m_parent);
    RelativePermeabilityBase * relPermPtr = &relPerm;
    // CapillaryPressureBase & capPressure0 =  makeJFunctionCapPressureTwoPhase( "capPressure0", this->m_parent );
    // CapillaryPressureBase * capPressurePtr0 = &capPressure0;

    // CapillaryPressureBase & capPressure1 =  makeTableCapPressureTwoPhase( "capPressure1", this->m_parent );
    // CapillaryPressureBase * capPressurePtr1 = &capPressure1;
    CapillaryPressureBase & capPressure0 =  makeBrooksCoreyCapPressureTwoPhase1( "capPressure0", this->m_parent );
    CapillaryPressureBase * capPressurePtr0 = &capPressure0;

    CapillaryPressureBase & capPressure1 =  makeBrooksCoreyCapPressureTwoPhase2( "capPressure1", this->m_parent );
    CapillaryPressureBase * capPressurePtr1 = &capPressure1;
    
  std::vector< RelativePermeabilityBase * > relPerms = {relPermPtr, relPermPtr};
  std::vector< CapillaryPressureBase * > capPressures = {capPressurePtr0, capPressurePtr1};
  std::vector< TwoPhaseImmiscibleFluid * > fluids = { &fluid, &fluid };
    // real64 uT = 3.2864545889999906e-05;
  
  real64 uT = 3.3e-5;
// real64 uT = 1e-17;
  stdVector< real64 > saturations = {0.2, 0.4};
  stdVector< real64 > pressures = {1e7, 1e7};
  stdVector< real64 > JFMultipliers = {45016.662822296035, 30011.108548197357};
  stdVector< real64 > transHats = {1.9738466000000002e-12, 4.4411548500000007e-12};
  stdVector< real64 > dTransHats_dP = {0.0, 0.0};
  stdVector< real64 > gravCoefHats = {490.5, 490.5};
  stdVector< real64 > gravCoefs = {465.97500000000002, 515.02499999999998};


  stdVector< real64 > phi = {0.0, 0.0};
  stdVector< real64 > grad_phi = {0.0, 0.0, 0.0, 0.0};

  std::ofstream outFile( "local_solver_results.csv" );


          // Write data to the file
          outFile << "Si";
          outFile << ",";
          outFile << "Sj";
          outFile << ",";
          outFile << "Fw_alpha";
          outFile << ",";
          outFile << "Fn_alpha";
          outFile << ",";
          outFile << "Pc_int";
          outFile << ",";
          outFile << "Residual";
          outFile << ",";
          outFile << "newton";
          outFile << std::endl;

          real64 const start_sat = 0.0;
          real64 const end_sat   = 1.0;
          real64 const dS        = 1e-2;
            // real64 Si = 0.98;
            // real64 Sj = 0.98;
          
          for( real64 Si = start_sat; Si <= end_sat + 1e-8; Si += dS )
          {
            for( real64 Sj = start_sat; Sj <= end_sat + 1e-8; Sj += dS )
            {
              saturations[0] = Si;
              saturations[1] = Sj;

    geos::immiscibleMultiphaseKernels::local_solver( uT, saturations, pressures, JFMultipliers, transHats, dTransHats_dP, gravCoefHats, gravCoefs,
      relPerms, capPressures, fluids, phi, grad_phi );


                  // Write data to the file
                  outFile << GEOS_FMT( "{:10.10e}", saturations[0] );
                  outFile << GEOS_FMT( ",{:10.10e}", saturations[1] );
                  outFile << GEOS_FMT( ",{:10.10e}", phi[0] );
                  outFile << GEOS_FMT( ",{:10.10e}", grad_phi[0] );
                  outFile << GEOS_FMT( ",{:10.10e}", grad_phi[1] );
                  outFile << GEOS_FMT( ",{:10.10e}", grad_phi[2] );
                  outFile << GEOS_FMT( ",{:10.10e}", grad_phi[3] );
                  outFile << std::endl;
                  phi[0] = 0;
                  grad_phi[0] = 0;
                  grad_phi[1] = 0;
                  grad_phi[2] = 0;
                  grad_phi[3] = 0;

  }
}

  outFile.close();

} );
}

int main( int argc, char * *argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}

// maybe needed later on
// TEST_F( CapillaryPressureTest, numericalDerivatives_jFunctionCapPressureTwoPhase )
// {
//   initialize( makeJFunctionCapPressureTwoPhase( "capPressure", m_parent ) );

//   // here, we have to apply a special treatment to this test
//   // to make sure that the J-function multiplier is initialized using initializeRockState
//   // this requires calling allocateConstitutiveData in advance (it will be called again later, in the "test" function)

//   // setup some values for porosity and permeability
//   array2d< real64 > porosity;
//   porosity.resize( 1, 1 );
//   porosity[0][0] = 0.13496794266569806;
//   array3d< real64 > permeability;
//   permeability.resize( 1, 1, 3 );
//   permeability[0][0][0] = 0.1722194e-15;
//   permeability[0][0][1] = 0.3423156e-15;
//   permeability[0][0][2] = 0.2324191e-15;

//   // initialize the J-function multiplier (done on GPU if GPU is available)
//   m_model->allocateConstitutiveData( m_parent, 1 );
//   m_model->initializeRockState( porosity.toViewConst(), permeability.toViewConst() );

//   // move the multiplier back to the CPU since the test is performed on the CPU
//   auto & jFuncMultiplier =
//     m_model->getReference< array2d< real64 > >( fields::cappres::jFuncMultiplier::key() );
//   jFuncMultiplier.move( hostMemorySpace, false );

//   // we are ready to proceed to the test

//   real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
//   real64 const tol = 1e-4;

//   real64 const start_sat = 0.3;
//   real64 const end_sat   = 0.9;
//   real64 const dS        = 1e-1;
//   array1d< real64 > sat( 2 );
//   sat[0] = start_sat; sat[1] = 1-sat[0];
//   while( sat[0] <= end_sat )
//   {
//     test( sat, eps, tol );
//     sat[0] += dS;
//     sat[1] = 1 - sat[0];
//   }

// }