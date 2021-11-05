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

// Source inclues
#include "constitutiveTestHelpers.hpp"
#include "functions/FunctionManager.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"

using namespace geosx;
using namespace geosx::testing;
using namespace geosx::constitutive;
using namespace geosx::constitutive::cappres;
using namespace geosx::dataRepository;

CapillaryPressureBase & makeBrooksCoreyCapPressureTwoPhase( string const & name, Group & parent )
{
  BrooksCoreyCapillaryPressure & capPressure = parent.registerGroup< BrooksCoreyCapillaryPressure >( name );

  string_array & phaseNames = capPressure.getReference< string_array >( CapillaryPressureBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "water"; phaseNames[1] = "oil";

  array1d< real64 > & phaseMinSat = capPressure.getReference< array1d< real64 > >( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseMinVolumeFractionString() );
  phaseMinSat.resize( 2 );
  phaseMinSat[0] = 0.1; phaseMinSat[1] = 0.05;

  array1d< real64 > & phaseCapPressureExpInv =
    capPressure.getReference< array1d< real64 > >( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseCapPressureExponentInvString() );
  phaseCapPressureExpInv.resize( 2 );
  phaseCapPressureExpInv[0] = 2; phaseCapPressureExpInv[1] = 4;

  array1d< real64 > & phaseEntryPressure = capPressure.getReference< array1d< real64 > >( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseEntryPressureString() );
  phaseEntryPressure.resize( 2 );
  phaseEntryPressure[0] = 1; phaseEntryPressure[1] = 1;

  real64 & capPressureEpsilon = capPressure.getReference< real64 >( BrooksCoreyCapillaryPressure::viewKeyStruct::capPressureEpsilonString() );
  capPressureEpsilon = 1e-4;

  capPressure.postProcessInputRecursive();
  return capPressure;
}


CapillaryPressureBase & makeBrooksCoreyCapPressureThreePhase( string const & name, Group & parent )
{
  BrooksCoreyCapillaryPressure & capPressure = parent.registerGroup< BrooksCoreyCapillaryPressure >( name );

  string_array & phaseNames = capPressure.getReference< string_array >( CapillaryPressureBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 3 );
  phaseNames[0] = "water"; phaseNames[1] = "oil"; phaseNames[2] = "gas";

  array1d< real64 > & phaseMinSat = capPressure.getReference< array1d< real64 > >( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseMinVolumeFractionString() );
  phaseMinSat.resize( 3 );
  phaseMinSat[0] = 0.04; phaseMinSat[1] = 0.02; phaseMinSat[2] = 0.1;

  array1d< real64 > & phaseCapPressureExpInv =
    capPressure.getReference< array1d< real64 > >( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseCapPressureExponentInvString() );
  phaseCapPressureExpInv.resize( 3 );
  phaseCapPressureExpInv[0] = 2; phaseCapPressureExpInv[1] = -3; phaseCapPressureExpInv[2] = 2.5;

  array1d< real64 > & phaseEntryPressure = capPressure.getReference< array1d< real64 > >( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseEntryPressureString() );
  phaseEntryPressure.resize( 3 );
  phaseEntryPressure[0] = 1; phaseEntryPressure[1] = -1; phaseEntryPressure[2] = 2;

  real64 & capPressureEpsilon = capPressure.getReference< real64 >( BrooksCoreyCapillaryPressure::viewKeyStruct::capPressureEpsilonString() );
  capPressureEpsilon = 1e-7;

  capPressure.postProcessInputRecursive();
  return capPressure;
}


CapillaryPressureBase & makeVanGenuchtenCapPressureTwoPhase( string const & name, Group & parent )
{
  VanGenuchtenCapillaryPressure & capPressure = parent.registerGroup< VanGenuchtenCapillaryPressure >( name );

  string_array & phaseNames = capPressure.getReference< string_array >( CapillaryPressureBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas";

  array1d< real64 > & phaseMinSat = capPressure.getReference< array1d< real64 > >( VanGenuchtenCapillaryPressure::viewKeyStruct::phaseMinVolumeFractionString() );
  phaseMinSat.resize( 2 );
  phaseMinSat[0] = 0.04; phaseMinSat[1] = 0.1;

  array1d< real64 > & phaseCapPressureExpInv = capPressure.getReference< array1d< real64 > >(
    VanGenuchtenCapillaryPressure::viewKeyStruct::phaseCapPressureExponentInvString() );
  phaseCapPressureExpInv.resize( 2 );
  phaseCapPressureExpInv[0] = 0.4; phaseCapPressureExpInv[1] = 0.5;

  array1d< real64 > & phaseCapPressureMultiplier = capPressure.getReference< array1d< real64 > >(
    VanGenuchtenCapillaryPressure::viewKeyStruct::phaseCapPressureMultiplierString() );
  phaseCapPressureMultiplier.resize( 2 );
  phaseCapPressureMultiplier[0] = 0.5; phaseCapPressureMultiplier[1] = 1;

  real64 & capPressureEpsilon = capPressure.getReference< real64 >( VanGenuchtenCapillaryPressure::viewKeyStruct::capPressureEpsilonString() );
  capPressureEpsilon = 1e-4;

  capPressure.postProcessInputRecursive();
  return capPressure;
}

CapillaryPressureBase & makeVanGenuchtenCapPressureThreePhase( string const & name, Group & parent )
{
  VanGenuchtenCapillaryPressure & capPressure = parent.registerGroup< VanGenuchtenCapillaryPressure >( name );

  string_array & phaseNames = capPressure.getReference< string_array >( CapillaryPressureBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 3 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas"; phaseNames[2] = "water";

  array1d< real64 > & phaseMinSat = capPressure.getReference< array1d< real64 > >( VanGenuchtenCapillaryPressure::viewKeyStruct::phaseMinVolumeFractionString() );
  phaseMinSat.resize( 3 );
  phaseMinSat[0] = 0.04; phaseMinSat[1] = 0.1; phaseMinSat[2] = 0.;

  array1d< real64 > & phaseCapPressureExpInv = capPressure.getReference< array1d< real64 > >(
    VanGenuchtenCapillaryPressure::viewKeyStruct::phaseCapPressureExponentInvString() );
  phaseCapPressureExpInv.resize( 3 );
  phaseCapPressureExpInv[0] = 0.33; phaseCapPressureExpInv[1] = 0.4; phaseCapPressureExpInv[ 2 ] = 0.5;

  array1d< real64 > & phaseCapPressureMultiplier = capPressure.getReference< array1d< real64 > >(
    VanGenuchtenCapillaryPressure::viewKeyStruct::phaseCapPressureMultiplierString() );
  phaseCapPressureMultiplier.resize( 3 );
  phaseCapPressureMultiplier[0] = 0.5; phaseCapPressureMultiplier[1] = 1; phaseCapPressureMultiplier[2] = 0.2;

  real64 & capPressureEpsilon = capPressure.getReference< real64 >( VanGenuchtenCapillaryPressure::viewKeyStruct::capPressureEpsilonString() );
  capPressureEpsilon = 1e-4;

  capPressure.postProcessInputRecursive();
  return capPressure;
}

CapillaryPressureBase & makeTableCapPressureTwoPhase( string const & name, Group & parent )
{
  FunctionManager & functionManager = FunctionManager::getInstance();

  // 1) First, define the tables

  // 1D table, various interpolation methods
  localIndex Naxis = 6;

  // Setup table
  array1d< real64_array > coordinates;
  coordinates.resize( 1 );
  coordinates[0].resize( Naxis );
  coordinates[0][0] = 0.22;
  coordinates[0][1] = 0.3;
  coordinates[0][2] = 0.5;
  coordinates[0][3] = 0.6;
  coordinates[0][4] = 0.8;
  coordinates[0][5] = 1.0;

  real64_array values( Naxis );
  values[0] = 48263.3;
  values[1] = 27579;
  values[2] = 20684.3;
  values[3] = 13789.5;
  values[4] = 6894.76;
  values[5] = 3204.28;

  TableFunction & table_w = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "water_pc" ) );
  table_w.setTableCoordinates( coordinates );
  table_w.setTableValues( values );
  table_w.reInitializeFunction();

  table_w.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  // 2) Then set up the constitutive model

  TableCapillaryPressure & capPressure = parent.registerGroup< TableCapillaryPressure >( name );

  string_array & phaseNames = capPressure.getReference< string_array >( CapillaryPressureBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "water"; phaseNames[1] = "gas";

  auto & waterTableName = capPressure.getReference< string >( TableCapillaryPressure::viewKeyStruct::wettingNonWettingCapPresTableNameString() );
  waterTableName = "water_pc";

  capPressure.postProcessInputRecursive();
  capPressure.initialize(); // to test all the checks
  return capPressure;
}

CapillaryPressureBase & makeTableCapPressureThreePhase( string const & name, Group & parent )
{
  FunctionManager & functionManager = FunctionManager::getInstance();

  // 1) First, define the tables

  // 1D table, various interpolation methods
  localIndex Naxis = 6;

  // 1.a) First pair of phases (ow)

  // Setup table
  array1d< real64_array > coordinates_w;
  coordinates_w.resize( 1 );
  coordinates_w[0].resize( Naxis );
  coordinates_w[0][0] = 0.22;
  coordinates_w[0][1] = 0.3;
  coordinates_w[0][2] = 0.5;
  coordinates_w[0][3] = 0.6;
  coordinates_w[0][4] = 0.8;
  coordinates_w[0][5] = 1.0;

  real64_array values_w( Naxis );
  values_w[0] = 48263.3;
  values_w[1] = 27579;
  values_w[2] = 20684.3;
  values_w[3] = 13789.5;
  values_w[4] = 6894.76;
  values_w[5] = 3294.76;

  TableFunction & table_ow_w = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "water_pc" ) );
  table_ow_w.setTableCoordinates( coordinates_w );
  table_ow_w.setTableValues( values_w );
  table_ow_w.reInitializeFunction();

  table_ow_w.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  // 1.a) Second pair of phases (og)

  array1d< real64_array > coordinates_g;
  coordinates_g.resize( 1 );
  coordinates_g[0].resize( Naxis );
  coordinates_g[0][0] = 0.0;
  coordinates_g[0][1] = 0.04;
  coordinates_g[0][2] = 0.1;
  coordinates_g[0][3] = 0.2;
  coordinates_g[0][4] = 0.7;
  coordinates_g[0][5] = 0.78;

  real64_array values_g( Naxis );
  values_g[0] = 0.0;
  values_g[1] = 1723.689;
  values_g[2] = 3447.38;
  values_g[3] = 6894.76;
  values_g[4] = 24131.7;
  values_g[5] = 26889.6;

  TableFunction & table_og_g = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "gas_pc" ) );
  table_og_g.setTableCoordinates( coordinates_g );
  table_og_g.setTableValues( values_g );
  table_og_g.reInitializeFunction();

  table_og_g.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  // 2) Then set up the constitutive model

  TableCapillaryPressure & capPressure = parent.registerGroup< TableCapillaryPressure >( name );

  string_array & phaseNames = capPressure.getReference< string_array >( CapillaryPressureBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 3 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas"; phaseNames[2] = "water";

  auto & waterTableName = capPressure.getReference< string >( TableCapillaryPressure::viewKeyStruct::wettingIntermediateCapPresTableNameString() );
  waterTableName = "water_pc";

  auto & gasTableName = capPressure.getReference< string >( TableCapillaryPressure::viewKeyStruct::nonWettingIntermediateCapPresTableNameString() );
  gasTableName = "gas_pc";

  capPressure.postProcessInputRecursive();
  capPressure.initialize(); // to test all the checks
  return capPressure;
}

class CapillaryPressureTest : public ConstitutiveTestBase< CapillaryPressureBase >
{
public:
  void test( arraySlice1d< real64 const > const sat, real64 const eps, real64 const tol )
  {
    arrayView3d< real64 const, USD_CAPPRES > phaseCapPressure;
    arrayView4d< real64 const, USD_CAPPRES_DS > dPhaseCapPressure_dPhaseVolFraction;
    testNumericalDerivatives( m_parent,
                              *m_model,
                              sat,
                              eps,
                              tol,
                              "phaseCapPressure",
                              [&phaseCapPressure] ( CapillaryPressureBase & relPerm )
    {
      phaseCapPressure = relPerm.phaseCapPressure();
      return phaseCapPressure[ 0 ][ 0 ];
    },
                              [&dPhaseCapPressure_dPhaseVolFraction] ( CapillaryPressureBase & relPerm )
    {
      dPhaseCapPressure_dPhaseVolFraction = relPerm.dPhaseCapPressure_dPhaseVolFraction();
      return dPhaseCapPressure_dPhaseVolFraction[ 0 ][ 0 ];
    }
                              );
  }
};

TEST_F( CapillaryPressureTest, numericalDerivatives_brooksCoreyCapPressureTwoPhase )
{
  initialize( makeBrooksCoreyCapPressureTwoPhase( "capPressure", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const start_sat = 0.4;
  real64 const end_sat   = 0.6;
  real64 const dS = 1e-1;
  array1d< real64 > sat( 2 );
  sat[0] = start_sat; sat[1] = 1.0-sat[0];
  while( sat[0] <= end_sat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = 1 - sat[0];
  }
}


TEST_F( CapillaryPressureTest, numericalDerivatives_brooksCoreyCapPressureThreePhase )
{
  initialize( makeBrooksCoreyCapPressureThreePhase( "capPressure", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const start_sat = 0.4;
  real64 const end_sat   = 0.6;
  real64 const dS = 1e-1;
  array1d< real64 > sat( 3 );
  sat[0] = start_sat;
  sat[1] = 0.5*(1-sat[0]);
  sat[2] = 1.0-sat[0]-sat[1];
  while( sat[0] <= end_sat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = 0.5 * ( 1-sat[0] );
    sat[2] = 1.0 - sat[0] - sat[1];
  }
}


TEST_F( CapillaryPressureTest, numericalDerivatives_vanGenuchtenCapPressureTwoPhase )
{
  initialize( makeVanGenuchtenCapPressureTwoPhase( "capPressure", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const start_sat = 0.4;
  real64 const end_sat   = 0.6;
  real64 const dS        = 1e-1;
  array1d< real64 > sat( 2 );
  sat[0] = start_sat; sat[1] = 1-sat[1];
  while( sat[0] <= end_sat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = 1 - sat[0];
  }

}


TEST_F( CapillaryPressureTest, numericalDerivatives_vanGenuchtenCapPressureThreePhase )
{
  initialize( makeVanGenuchtenCapPressureThreePhase( "capPressure", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const start_sat = 0.4;
  real64 const end_sat   = 0.6;
  real64 const dS        = 1e-1;
  array1d< real64 > sat( 3 );
  sat[0] = start_sat;
  sat[1] = 0.5*(1-sat[0]);
  sat[2] = 1.0-sat[0]-sat[1];
  while( sat[0] <= end_sat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = 0.5*(1-sat[0]);
    sat[2] = 1 - sat[0] - sat[1];
  }
}


TEST_F( CapillaryPressureTest, numericalDerivatives_tableCapPressureTwoPhase )
{
  initialize( makeTableCapPressureTwoPhase( "capPressure", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const start_sat = 0.25;
  real64 const end_sat   = 0.75;
  real64 const dS        = 1e-1;
  array1d< real64 > sat( 2 );
  sat[0] = start_sat; sat[1] = 1-sat[1];
  while( sat[0] <= end_sat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = 1 - sat[0];
  }

}


TEST_F( CapillaryPressureTest, numericalDerivatives_tableCapPressureThreePhase )
{
  initialize( makeTableCapPressureThreePhase( "capPressure", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const start_sat = 0.25;
  real64 const end_sat   = 0.75;
  real64 const dS        = 1e-1;
  array1d< real64 > sat( 3 );
  sat[0] = start_sat;
  sat[1] = 0.5*(1-sat[0]);
  sat[2] = 1.0-sat[0]-sat[1];
  while( sat[0] <= end_sat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = 0.5*(1-sat[0]);
    sat[2] = 1 - sat[0] - sat[1];
  }
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::GeosxState state( geosx::basicSetup( argc, argv ) );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
