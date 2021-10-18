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

// Source includes
#include "constitutiveTestHelpers.hpp"
#include "functions/FunctionManager.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"

using namespace geosx;
using namespace geosx::testing;
using namespace geosx::constitutive;
using namespace geosx::constitutive::relperm;
using namespace geosx::dataRepository;

RelativePermeabilityBase & makeBrooksCoreyRelPerm( string const & name, Group & parent )
{
  BrooksCoreyRelativePermeability & relPerm = parent.registerGroup< BrooksCoreyRelativePermeability >( name );

  string_array & phaseNames = relPerm.getReference< string_array >( RelativePermeabilityBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas";

  array1d< real64 > & phaseMinSat = relPerm.getReference< array1d< real64 > >( BrooksCoreyRelativePermeability::viewKeyStruct::phaseMinVolumeFractionString() );
  phaseMinSat.resize( 2 );
  phaseMinSat[0] = 0.1; phaseMinSat[1] = 0.15;

  array1d< real64 > & phaseRelPermExp = relPerm.getReference< array1d< real64 > >( BrooksCoreyRelativePermeability::viewKeyStruct::phaseRelPermExponentString() );
  phaseRelPermExp.resize( 2 );
  phaseRelPermExp[0] = 2.0; phaseRelPermExp[1] = 2.0;

  array1d< real64 > & phaseRelPermMaxVal = relPerm.getReference< array1d< real64 > >( BrooksCoreyRelativePermeability::viewKeyStruct::phaseRelPermMaxValueString() );
  phaseRelPermMaxVal.resize( 2 );
  phaseRelPermMaxVal[0] = 0.8; phaseRelPermMaxVal[1] = 0.9;

  relPerm.postProcessInputRecursive();
  return relPerm;
}

RelativePermeabilityBase & makeBrooksCoreyBakerRelPermTwoPhase( string const & name, Group & parent )
{
  BrooksCoreyBakerRelativePermeability & relPerm = parent.registerGroup< BrooksCoreyBakerRelativePermeability >( name );

  string_array & phaseNames = relPerm.getReference< string_array >( RelativePermeabilityBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "water"; phaseNames[1] = "oil";

  array1d< real64 > & phaseMinSat = relPerm.getReference< array1d< real64 > >( BrooksCoreyBakerRelativePermeability::viewKeyStruct::phaseMinVolumeFractionString() );
  phaseMinSat.resize( 2 );
  phaseMinSat[0] = 0.03; phaseMinSat[1] = 0.01;

  array1d< real64 > & waterOilRelPermExp = relPerm.getReference< array1d< real64 > >( BrooksCoreyBakerRelativePermeability::viewKeyStruct::waterOilRelPermExponentString() );
  waterOilRelPermExp.resize( 2 );
  waterOilRelPermExp[0] = 1.9; waterOilRelPermExp[1] = 3.95;

  array1d< real64 > & waterOilRelPermMaxVal =
    relPerm.getReference< array1d< real64 > >( BrooksCoreyBakerRelativePermeability::viewKeyStruct::waterOilRelPermMaxValueString() );
  waterOilRelPermMaxVal.resize( 2 );
  waterOilRelPermMaxVal[0] = 0.8; waterOilRelPermMaxVal[1] = 0.75;

  relPerm.postProcessInputRecursive();
  return relPerm;
}

RelativePermeabilityBase & makeBrooksCoreyBakerRelPermThreePhase( string const & name, Group & parent )
{
  BrooksCoreyBakerRelativePermeability & relPerm = parent.registerGroup< BrooksCoreyBakerRelativePermeability >( name );

  string_array & phaseNames = relPerm.getReference< string_array >( RelativePermeabilityBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 3 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas"; phaseNames[2] = "water";

  array1d< real64 > & phaseMinSat = relPerm.getReference< array1d< real64 > >( BrooksCoreyBakerRelativePermeability::viewKeyStruct::phaseMinVolumeFractionString() );
  phaseMinSat.resize( 3 );
  phaseMinSat[0] = 0.03; phaseMinSat[1] = 0.01; phaseMinSat[2] = 0.025;

  array1d< real64 > & waterOilRelPermExp = relPerm.getReference< array1d< real64 > >( BrooksCoreyBakerRelativePermeability::viewKeyStruct::waterOilRelPermExponentString() );
  waterOilRelPermExp.resize( 2 );
  waterOilRelPermExp[0] = 2.4; waterOilRelPermExp[1] = 1.5;

  array1d< real64 > & waterOilRelPermMaxVal =
    relPerm.getReference< array1d< real64 > >( BrooksCoreyBakerRelativePermeability::viewKeyStruct::waterOilRelPermMaxValueString() );
  waterOilRelPermMaxVal.resize( 2 );
  waterOilRelPermMaxVal[0] = 0.9; waterOilRelPermMaxVal[1] = 0.65;

  array1d< real64 > & gasOilRelPermExp = relPerm.getReference< array1d< real64 > >( BrooksCoreyBakerRelativePermeability::viewKeyStruct::gasOilRelPermExponentString() );
  gasOilRelPermExp.resize( 2 );
  gasOilRelPermExp[0] = 1.9; gasOilRelPermExp[1] = 3.95;

  array1d< real64 > & gasOilRelPermMaxVal = relPerm.getReference< array1d< real64 > >( BrooksCoreyBakerRelativePermeability::viewKeyStruct::gasOilRelPermMaxValueString() );
  gasOilRelPermMaxVal.resize( 2 );
  gasOilRelPermMaxVal[0] = 0.8; gasOilRelPermMaxVal[1] = 0.95;

  relPerm.postProcessInputRecursive();
  return relPerm;
}

RelativePermeabilityBase & makeVanGenuchtenBakerRelPermTwoPhase( string const & name, Group & parent )
{
  VanGenuchtenBakerRelativePermeability & relPerm = parent.registerGroup< VanGenuchtenBakerRelativePermeability >( name );

  string_array & phaseNames = relPerm.getReference< string_array >( RelativePermeabilityBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas";

  array1d< real64 > & phaseMinSat = relPerm.getReference< array1d< real64 > >( VanGenuchtenBakerRelativePermeability::viewKeyStruct::phaseMinVolumeFractionString() );
  phaseMinSat.resize( 2 );
  phaseMinSat[0] = 0.02; phaseMinSat[1] = 0.05;

  array1d< real64 > & gasOilRelPermExpInv =
    relPerm.getReference< array1d< real64 > >( VanGenuchtenBakerRelativePermeability::viewKeyStruct::gasOilRelPermExponentInvString() );
  gasOilRelPermExpInv.resize( 2 );
  gasOilRelPermExpInv[0] = 1.7; gasOilRelPermExpInv[1] = 2.15;

  array1d< real64 > & gasOilRelPermMaxVal = relPerm.getReference< array1d< real64 > >( VanGenuchtenBakerRelativePermeability::viewKeyStruct::gasOilRelPermMaxValueString() );
  gasOilRelPermMaxVal.resize( 2 );
  gasOilRelPermMaxVal[0] = 0.5; gasOilRelPermMaxVal[1] = 0.75;

  relPerm.postProcessInputRecursive();
  return relPerm;
}

RelativePermeabilityBase & makeVanGenuchtenBakerRelPermThreePhase( string const & name, Group & parent )
{
  VanGenuchtenBakerRelativePermeability & relPerm = parent.registerGroup< VanGenuchtenBakerRelativePermeability >( name );

  string_array & phaseNames = relPerm.getReference< string_array >( RelativePermeabilityBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 3 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas"; phaseNames[2] = "water";

  array1d< real64 > & phaseMinSat = relPerm.getReference< array1d< real64 > >( VanGenuchtenBakerRelativePermeability::viewKeyStruct::phaseMinVolumeFractionString() );
  phaseMinSat.resize( 3 );
  phaseMinSat[0] = 0.03; phaseMinSat[1] = 0.01; phaseMinSat[2] = 0.025;

  array1d< real64 > & waterOilRelPermExpInv = relPerm.getReference< array1d< real64 > >(
    VanGenuchtenBakerRelativePermeability::viewKeyStruct::waterOilRelPermExponentInvString() );
  waterOilRelPermExpInv.resize( 2 );
  waterOilRelPermExpInv[0] = 2.4; waterOilRelPermExpInv[1] = 2.5;

  array1d< real64 > & waterOilRelPermMaxVal =
    relPerm.getReference< array1d< real64 > >( VanGenuchtenBakerRelativePermeability::viewKeyStruct::waterOilRelPermMaxValueString() );
  waterOilRelPermMaxVal.resize( 2 );
  waterOilRelPermMaxVal[0] = 0.9; waterOilRelPermMaxVal[1] = 0.75;

  array1d< real64 > & gasOilRelPermExpInv =
    relPerm.getReference< array1d< real64 > >( VanGenuchtenBakerRelativePermeability::viewKeyStruct::gasOilRelPermExponentInvString() );
  gasOilRelPermExpInv.resize( 2 );
  gasOilRelPermExpInv[0] = 1.9; gasOilRelPermExpInv[1] = 3.95;

  array1d< real64 > & gasOilRelPermMaxVal = relPerm.getReference< array1d< real64 > >( BrooksCoreyBakerRelativePermeability::viewKeyStruct::gasOilRelPermMaxValueString() );
  gasOilRelPermMaxVal.resize( 2 );
  gasOilRelPermMaxVal[0] = 0.8; gasOilRelPermMaxVal[1] = 0.75;

  relPerm.postProcessInputRecursive();
  return relPerm;
}

RelativePermeabilityBase & makeTableRelPermTwoPhase( string const & name, Group & parent )
{
  FunctionManager & functionManager = FunctionManager::getInstance();

  // 1) First, define the tables (to values that matters for our use cases)

  // 1D table, various interpolation methods
  localIndex Naxis = 41;

  // Setup table
  array1d< real64_array > coordinates_w;
  coordinates_w.resize( 1 );
  coordinates_w[0].resize( Naxis );

  coordinates_w[0][0] = 3.000000e-1;
  coordinates_w[0][1] = 3.175000e-1;
  coordinates_w[0][2] = 3.350000e-1;
  coordinates_w[0][3] = 3.525000e-1;
  coordinates_w[0][4] = 3.700000e-1;
  coordinates_w[0][5] = 3.875000e-1;
  coordinates_w[0][6] = 4.050000e-1;
  coordinates_w[0][7] = 4.225000e-1;
  coordinates_w[0][8] = 4.400000e-1;
  coordinates_w[0][9] = 4.575000e-1;
  coordinates_w[0][10] = 4.750000e-1;
  coordinates_w[0][11] = 4.925000e-1;
  coordinates_w[0][12] = 5.100000e-1;
  coordinates_w[0][13] = 5.275000e-1;
  coordinates_w[0][14] = 5.450000e-1;
  coordinates_w[0][15] = 5.625000e-1;
  coordinates_w[0][16] = 5.800000e-1;
  coordinates_w[0][17] = 5.975000e-1;
  coordinates_w[0][18] = 6.150000e-1;
  coordinates_w[0][19] = 6.325000e-1;
  coordinates_w[0][20] = 6.500000e-1;
  coordinates_w[0][21] = 6.675000e-1;
  coordinates_w[0][22] = 6.850000e-1;
  coordinates_w[0][23] = 7.025000e-1;
  coordinates_w[0][24] = 7.200000e-1;
  coordinates_w[0][25] = 7.375000e-1;
  coordinates_w[0][26] = 7.550000e-1;
  coordinates_w[0][27] = 7.725000e-1;
  coordinates_w[0][28] = 7.900000e-1;
  coordinates_w[0][29] = 8.054617e-1;
  coordinates_w[0][30] = 8.209233e-1;
  coordinates_w[0][31] = 8.404617e-1;
  coordinates_w[0][32] = 8.600000e-1;
  coordinates_w[0][33] = 8.775000e-1;
  coordinates_w[0][34] = 8.950000e-1;
  coordinates_w[0][35] = 9.125000e-1;
  coordinates_w[0][36] = 9.300000e-1;
  coordinates_w[0][37] = 9.475000e-1;
  coordinates_w[0][38] = 9.650000e-1;
  coordinates_w[0][39] = 9.825000e-1;
  coordinates_w[0][40] = 1.000000;

  array1d< real64_array > coordinates_g;
  coordinates_g.resize( 1 );
  coordinates_g[0].resize( Naxis );

  coordinates_g[0][0] = 0.000000;
  coordinates_g[0][1] = 1.750000e-2;
  coordinates_g[0][2] = 3.500000e-2;
  coordinates_g[0][3] = 5.250000e-2;
  coordinates_g[0][4] = 7.000000e-2;
  coordinates_g[0][5] = 8.750000e-2;
  coordinates_g[0][6] = 1.050000e-1;
  coordinates_g[0][7] = 1.225000e-1;
  coordinates_g[0][8] = 1.400000e-1;
  coordinates_g[0][9] = 1.595383e-1;
  coordinates_g[0][10] = 1.790767e-1;
  coordinates_g[0][11] = 1.945383e-1;
  coordinates_g[0][12] = 2.100000e-1;
  coordinates_g[0][13] = 2.275000e-1;
  coordinates_g[0][14] = 2.450000e-1;
  coordinates_g[0][15] = 2.625000e-1;
  coordinates_g[0][16] = 2.800000e-1;
  coordinates_g[0][17] = 2.975000e-1;
  coordinates_g[0][18] = 3.150000e-1;
  coordinates_g[0][19] = 3.325000e-1;
  coordinates_g[0][20] = 3.500000e-1;
  coordinates_g[0][21] = 3.675000e-1;
  coordinates_g[0][22] = 3.850000e-1;
  coordinates_g[0][23] = 4.025000e-1;
  coordinates_g[0][24] = 4.200000e-1;
  coordinates_g[0][25] = 4.375000e-1;
  coordinates_g[0][26] = 4.550000e-1;
  coordinates_g[0][27] = 4.725000e-1;
  coordinates_g[0][28] = 4.900000e-1;
  coordinates_g[0][29] = 5.075000e-1;
  coordinates_g[0][30] = 5.250000e-1;
  coordinates_g[0][31] = 5.425000e-1;
  coordinates_g[0][32] = 5.600000e-1;
  coordinates_g[0][33] = 5.775000e-1;
  coordinates_g[0][34] = 5.950000e-1;
  coordinates_g[0][35] = 6.125000e-1;
  coordinates_g[0][36] = 6.300000e-1;
  coordinates_g[0][37] = 6.475000e-1;
  coordinates_g[0][38] = 6.650000e-1;
  coordinates_g[0][39] = 6.825000e-1;
  coordinates_g[0][40] = 7.000000e-1;

  real64_array values_w( Naxis );
  values_w[0] = 0.000000;
  values_w[1] = 1.069690e-7;
  values_w[2] = 1.523818e-6;
  values_w[3] = 7.304599e-6;
  values_w[4] = 2.242961e-5;
  values_w[5] = 5.398050e-5;
  values_w[6] = 1.113999e-4;
  values_w[7] = 2.068239e-4;
  values_w[8] = 3.554932e-4;
  values_w[9] = 5.762517e-4;
  values_w[10] = 8.921512e-4;
  values_w[11] = 1.331180e-3;
  values_w[12] = 1.927144e-3;
  values_w[13] = 2.720726e-3;
  values_w[14] = 3.760776e-3;
  values_w[15] = 5.105868e-3;
  values_w[16] = 6.826186e-3;
  values_w[17] = 9.005830e-3;
  values_w[18] = 1.174561e-2;
  values_w[19] = 1.516648e-2;
  values_w[20] = 1.941368e-2;
  values_w[21] = 2.466185e-2;
  values_w[22] = 3.112128e-2;
  values_w[23] = 3.904542e-2;
  values_w[24] = 4.874017e-2;
  values_w[25] = 6.057494e-2;
  values_w[26] = 7.499593e-2;
  values_w[27] = 9.254174e-2;
  values_w[28] = 1.138611e-1;
  values_w[29] = 1.364565e-1;
  values_w[30] = 1.632363e-1;
  values_w[31] = 2.042135e-1;
  values_w[32] = 2.547712e-1;
  values_w[33] = 3.097943e-1;
  values_w[34] = 3.755964e-1;
  values_w[35] = 4.536528e-1;
  values_w[36] = 5.451093e-1;
  values_w[37] = 6.502388e-1;
  values_w[38] = 7.674166e-1;
  values_w[39] = 8.909226e-1;
  values_w[40] = 1.000000;

  real64_array values_g( Naxis );
  values_g[0] = 0.000000;
  values_g[1] = 8.885248e-4;
  values_g[2] = 2.483741e-3;
  values_g[3] = 4.583224e-3;
  values_g[4] = 7.135315e-3;
  values_g[5] = 1.012132e-2;
  values_g[6] = 1.353719e-2;
  values_g[7] = 1.738728e-2;
  values_g[8] = 2.168159e-2;
  values_g[9] = 2.701850e-2;
  values_g[10] = 3.295183e-2;
  values_g[11] = 3.808925e-2;
  values_g[12] = 4.363513e-2;
  values_g[13] = 5.042783e-2;
  values_g[14] = 5.779578e-2;
  values_g[15] = 6.577020e-2;
  values_g[16] = 7.438478e-2;
  values_g[17] = 8.367565e-2;
  values_g[18] = 9.368138e-2;
  values_g[19] = 1.044429e-1;
  values_g[20] = 1.160032e-1;
  values_g[21] = 1.284076e-1;
  values_g[22] = 1.417029e-1;
  values_g[23] = 1.559376e-1;
  values_g[24] = 1.711607e-1;
  values_g[25] = 1.874214e-1;
  values_g[26] = 2.047679e-1;
  values_g[27] = 2.232459e-1;
  values_g[28] = 2.428968e-1;
  values_g[29] = 2.637550e-1;
  values_g[30] = 2.858446e-1;
  values_g[31] = 3.091747e-1;
  values_g[32] = 3.337331e-1;
  values_g[33] = 3.594782e-1;
  values_g[34] = 3.863263e-1;
  values_g[35] = 4.141347e-1;
  values_g[36] = 4.426735e-1;
  values_g[37] = 4.715782e-1;
  values_g[38] = 5.002513e-1;
  values_g[39] = 5.275887e-1;
  values_g[40] = 5.500000e-1;

  TableFunction & table_w = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "water_swg" ) );
  table_w.setTableCoordinates( coordinates_w );
  table_w.setTableValues( values_w );
  table_w.reInitializeFunction();

  table_w.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  TableFunction & table_g = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "gas_swg" ) );
  table_g.setTableCoordinates( coordinates_g );
  table_g.setTableValues( values_g );
  table_g.reInitializeFunction();

  table_g.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  // 2) Then set up the constitutive model

  auto & relPerm = parent.registerGroup< TableRelativePermeability >( name );

  auto & phaseNames = relPerm.getReference< string_array >( RelativePermeabilityBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "water"; phaseNames[1] = "gas";

  auto & waterOilTableNames = relPerm.getReference< array1d< string > >( TableRelativePermeability::viewKeyStruct::wettingNonWettingRelPermTableNamesString() );
  waterOilTableNames.resize( 2 );
  waterOilTableNames[0] = "water_swg"; waterOilTableNames[1] = "gas_swg";

  relPerm.postProcessInputRecursive();
  relPerm.initialize(); // to test all the checks
  return relPerm;
}

RelativePermeabilityBase & makeTableRelPermThreePhase( string const & name, Group & parent )
{
  FunctionManager & functionManager = FunctionManager::getInstance();

  // 1) First, define the tables

  // 1D table, various interpolation methods
  localIndex Naxis = 6;

  // 1.a) First pair of phases (ow)

  // Setup table
  array1d< real64_array > coordinates;
  coordinates.resize( 1 );
  coordinates[0].resize( Naxis );
  coordinates[0][0] = 0.0;
  coordinates[0][1] = 0.11;
  coordinates[0][2] = 0.23;
  coordinates[0][3] = 0.54;
  coordinates[0][4] = 0.85;
  coordinates[0][5] = 1.0;

  real64_array values( Naxis );
  for( localIndex i = 0; i < coordinates[0].size(); ++i )
  {
    values[i] = coordinates[0][i]*coordinates[0][i];
  }

  TableFunction & table_ow_w = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "water_swof" ) );
  table_ow_w.setTableCoordinates( coordinates );
  table_ow_w.setTableValues( values );
  table_ow_w.reInitializeFunction();

  table_ow_w.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  TableFunction & table_ow_o = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "oil_swof" ) );
  table_ow_o.setTableCoordinates( coordinates );
  table_ow_o.setTableValues( values );
  table_ow_o.reInitializeFunction();

  table_ow_o.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  // 1.a) Second pair of phases (og)

  coordinates[0].resize( Naxis );
  coordinates[0][0] = 0.0;
  coordinates[0][1] = 0.01;
  coordinates[0][2] = 0.23;
  coordinates[0][3] = 0.44;
  coordinates[0][4] = 0.83;
  coordinates[0][5] = 1.0;

  for( localIndex i = 0; i < coordinates[0].size(); ++i )
  {
    values[i] = coordinates[0][i]*coordinates[0][i]*coordinates[0][i];
  }

  TableFunction & table_og_g = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "gas_sgof" ) );
  table_og_g.setTableCoordinates( coordinates );
  table_og_g.setTableValues( values );
  table_og_g.reInitializeFunction();

  table_og_g.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  TableFunction & table_og_o = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "oil_sgof" ) );
  table_og_o.setTableCoordinates( coordinates );
  table_og_o.setTableValues( values );
  table_og_o.reInitializeFunction();

  table_og_o.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  // 2) Then set up the constitutive model

  auto & relPerm = parent.registerGroup< TableRelativePermeability >( name );

  auto & phaseNames = relPerm.getReference< string_array >( RelativePermeabilityBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 3 );
  phaseNames[0] = "oil"; phaseNames[1] = "water"; phaseNames[2] = "gas";

  auto & waterOilTableNames = relPerm.getReference< array1d< string > >( TableRelativePermeability::viewKeyStruct::wettingIntermediateRelPermTableNamesString() );
  waterOilTableNames.resize( 2 );
  waterOilTableNames[0] = "water_swof"; waterOilTableNames[1] = "oil_swof";

  auto & gasOilTableNames = relPerm.getReference< array1d< string > >( TableRelativePermeability::viewKeyStruct::nonWettingIntermediateRelPermTableNamesString() );
  gasOilTableNames.resize( 2 );
  gasOilTableNames[0] = "gas_sgof"; gasOilTableNames[1] = "oil_sgof";

  relPerm.postProcessInputRecursive();
  relPerm.initialize(); // to test all the checks
  return relPerm;
}

class RelPermTest : public ConstitutiveTestBase< RelativePermeabilityBase >
{
public:
  void test( arraySlice1d< real64 const > const sat, real64 const eps, real64 const tol )
  {
    arrayView3d< real64 const, USD_RELPERM > phaseRelPerm;
    arrayView4d< real64 const, USD_RELPERM_DS > dPhaseRelPerm_dPhaseVolFraction;
    testNumericalDerivatives( m_parent,
                              *m_model,
                              sat,
                              eps,
                              tol,
                              "phaseRelPerm",
                              [&phaseRelPerm] ( RelativePermeabilityBase & relPerm )
    {
      phaseRelPerm = relPerm.phaseRelPerm();
      return phaseRelPerm[ 0 ][ 0 ];
    },
                              [&dPhaseRelPerm_dPhaseVolFraction] ( RelativePermeabilityBase & relPerm )
    {
      dPhaseRelPerm_dPhaseVolFraction = relPerm.dPhaseRelPerm_dPhaseVolFraction();
      return dPhaseRelPerm_dPhaseVolFraction[ 0 ][ 0 ];
    }
                              );
  }
};

TEST_F( RelPermTest, numericalDerivatives_brooksCoreyRelPerm )
{
  initialize( makeBrooksCoreyRelPerm( "relPerm", m_parent ) );

  // TODO test over a range of values
  array1d< real64 > sat( 2 );
  sat[0] = 0.7; sat[1] = 0.3;

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  test( sat, eps, tol );
}

TEST_F( RelPermTest, numericalDerivatives_BrooksCoreyBakerRelPermTwoPhase )
{
  initialize( makeBrooksCoreyBakerRelPermTwoPhase( "relPerm", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  // TODO test over a range of values
  real64 const start_sat = 0.3;
  real64 const end_sat   = 0.7;
  real64 const dS = 1e-1;
  real64 const alpha = 0.4;
  array1d< real64 > sat( 2 );
  sat[0] = start_sat;
  sat[1] = alpha*(1.0-sat[0]);
  while( sat[0] <= end_sat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = 1-sat[0];
  }
}

TEST_F( RelPermTest, numericalDerivatives_BrooksCoreyBakerRelPermThreePhase )
{
  initialize( makeBrooksCoreyBakerRelPermThreePhase( "relPerm", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const start_sat = 0.3;
  real64 const end_sat   = 0.7;
  real64 const dS = 1e-1;
  real64 const alpha = 0.4;
  array1d< real64 > sat( 3 );
  sat[0] = start_sat;
  sat[1] = alpha*(1.0-sat[0]);
  sat[2] = (1-alpha)*(1.0-sat[0]);
  while( sat[0] <= end_sat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = alpha *(1-sat[0]);
    sat[2] = (1-alpha) *(1-sat[0]);
  }
}


TEST_F( RelPermTest, numericalDerivatives_VanGenuchtenBakerRelPermTwoPhase )
{
  initialize( makeVanGenuchtenBakerRelPermTwoPhase( "relPerm", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const start_sat = 0.3;
  real64 const end_sat   = 0.7;
  real64 const dS = 1e-1;
  real64 const alpha = 0.4;
  array1d< real64 > sat( 2 );
  sat[0] = start_sat;
  sat[1] = alpha*(1.0-sat[0]);
  while( sat[0] <= end_sat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = 1-sat[0];
  }
}

TEST_F( RelPermTest, numericalDerivatives_VanGenuchtenBakerRelPermThreePhase )
{
  initialize( makeVanGenuchtenBakerRelPermThreePhase( "relPerm", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const start_sat = 0.3;
  real64 const end_sat   = 0.7;
  real64 const dS = 1e-1;
  real64 const alpha = 0.4;
  array1d< real64 > sat( 3 );
  sat[0] = start_sat;
  sat[1] = alpha*(1.0-sat[0]);
  sat[2] = (1-alpha)*(1.0-sat[0]);
  while( sat[0] <= end_sat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = alpha *(1-sat[0]);
    sat[2] = (1-alpha) *(1-sat[0]);
  }
}

TEST_F( RelPermTest, numericalDerivatives_TableRelPermTwoPhase )
{
  initialize( makeTableRelPermTwoPhase( "relPerm", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const start_sat = 0.2;
  real64 const end_sat   = 0.6;
  real64 const dS = 1e-1;
  real64 const alpha = 0.4;
  array1d< real64 > sat( 2 );
  sat[0] = start_sat;
  sat[1] = alpha*(1.0-sat[0]);
  while( sat[0] <= end_sat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = 1-sat[0];
  }
}

TEST_F( RelPermTest, numericalDerivatives_TableRelPermThreePhase )
{
  initialize( makeTableRelPermThreePhase( "relPerm", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const start_sat = 0.3;
  real64 const end_sat   = 0.7;
  real64 const dS = 1e-1;
  real64 const alpha = 0.4;
  array1d< real64 > sat( 3 );
  sat[0] = start_sat;
  sat[1] = alpha*(1.0-sat[0]);
  sat[2] = (1-alpha)*(1.0-sat[0]);
  while( sat[0] <= end_sat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = alpha *(1-sat[0]);
    sat[2] = (1-alpha) *(1-sat[0]);
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
