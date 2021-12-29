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
#include "codingUtilities/UnitTestUtilities.hpp"
#include "constitutiveTestHelpers.hpp"
#include "functions/FunctionManager.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"

using namespace geosx;
using namespace geosx::testing;
using namespace geosx::constitutive;
using namespace geosx::constitutive::relperm;
using namespace geosx::dataRepository;


TableRelativePermeabilityHysteresis& makeTableRelPermHysteresisTwoPhase( string const & name, Group & parent )
{
  FunctionManager & functionManager = FunctionManager::getInstance();

  // 1) First, define the tables (to values that matters for our use cases)

  // 1D table, various interpolation methods
  localIndex Naxis1 = 17;
  localIndex Naxis2 = 20;

  // Setup table
  array1d< real64_array > coordinates_w;
  coordinates_w.resize( 1 );
  coordinates_w[0].resize( Naxis1 );

  coordinates_w[0][0] = 0.20000;
  coordinates_w[0][1] = 0.22000;
  coordinates_w[0][2] = 0.25000;
  coordinates_w[0][3] = 0.30000;
  coordinates_w[0][4] = 0.35000;
  coordinates_w[0][5] = 0.40000;
  coordinates_w[0][6] = 0.45000;
  coordinates_w[0][7] = 0.50000;
  coordinates_w[0][8] = 0.55000;
  coordinates_w[0][9] = 0.60000;
  coordinates_w[0][10] = 0.65000;
  coordinates_w[0][11] = 0.66000;
  coordinates_w[0][12] = 0.68000;
  coordinates_w[0][13] = 0.72000;
  coordinates_w[0][14] = 0.82000;
  coordinates_w[0][15] = 0.91000;
  coordinates_w[0][16] = 1.00000;

  array1d< real64_array > coordinates_g;
  coordinates_g.resize( 1 );
  coordinates_g[0].resize( Naxis2 );

  coordinates_g[0][0] = 0.000;
  coordinates_g[0][1] = 0.005;
  coordinates_g[0][2] = 0.010;
  coordinates_g[0][3] = 0.030;
  coordinates_g[0][4] = 0.050;
  coordinates_g[0][5] = 0.100;
  coordinates_g[0][6] = 0.150;
  coordinates_g[0][7] = 0.200;
  coordinates_g[0][8] = 0.250;
  coordinates_g[0][9] = 0.300;
  coordinates_g[0][10] = 0.350;
  coordinates_g[0][11] = 0.400;
  coordinates_g[0][12] = 0.450;
  coordinates_g[0][13] = 0.500;
  coordinates_g[0][14] = 0.550;
  coordinates_g[0][15] = 0.600;
  coordinates_g[0][16] = 0.650;
  coordinates_g[0][17] = 0.700;
  coordinates_g[0][18] = 0.750;
  coordinates_g[0][19] = 0.780;

  real64_array drainageValues_w( Naxis1 );

  drainageValues_w[0] =  0.00000;
  drainageValues_w[1] =  0.00000;
  drainageValues_w[2] =  0.00100;
  drainageValues_w[3] =  0.00300;
  drainageValues_w[4] =  0.01000;
  drainageValues_w[5] =  0.01800;
  drainageValues_w[6] =  0.03500;
  drainageValues_w[7] =  0.04000;
  drainageValues_w[8] =  0.05700;
  drainageValues_w[9] =  0.08800;
  drainageValues_w[10] = 0.14500;
  drainageValues_w[11] = 0.16000;
  drainageValues_w[12] = 0.19000;
  drainageValues_w[13] = 0.26300;
  drainageValues_w[14] = 0.45500;
  drainageValues_w[15] = 0.69200;
  drainageValues_w[16] = 1.;


  real64_array drainageValues_g( Naxis2 );

  drainageValues_g[0] = 0.00000;
  drainageValues_g[1] = 0.00000;
  drainageValues_g[2] = 0.00200;
  drainageValues_g[3] = 0.00700;
  drainageValues_g[4] = 0.01000;
  drainageValues_g[5] = 0.02000;
  drainageValues_g[6] = 0.04000;
  drainageValues_g[7] = 0.07500;
  drainageValues_g[8] = 0.12700;
  drainageValues_g[9] = 0.18000;
  drainageValues_g[10] = 0.24000;
  drainageValues_g[11] = 0.31000;
  drainageValues_g[12] = 0.37300;
  drainageValues_g[13] = 0.46000;
  drainageValues_g[14] = 0.55000;
  drainageValues_g[15] = 0.64000;
  drainageValues_g[16] = 0.73000;
  drainageValues_g[17] = 0.82500;
  drainageValues_g[18] = 0.92000;
  drainageValues_g[19] = 1.00000;

  real64_array imbibitionValues_w( Naxis1 );

  imbibitionValues_w[0] = 0;
  imbibitionValues_w[1] = 0;
  imbibitionValues_w[2] = 0.0156;
  imbibitionValues_w[3] = 0.0680;
  imbibitionValues_w[4] = 0.1409;
  imbibitionValues_w[5] = 0.2296;
  imbibitionValues_w[6] = 0.3317;
  imbibitionValues_w[7] = 0.4455;
  imbibitionValues_w[8] = 0.5700;
  imbibitionValues_w[9] = 0.7044;
  imbibitionValues_w[10] = 0.8479;
  imbibitionValues_w[11] = 0.8776;
  imbibitionValues_w[12] = 0.9382;
  imbibitionValues_w[13] = 0.9382;
  imbibitionValues_w[14] = 0.9382;
  imbibitionValues_w[15] = 0.9382;
  imbibitionValues_w[16] = 0.9382;

  real64_array imbibitionValues_g( Naxis2 );

  imbibitionValues_g[0] = 0.0000;
  imbibitionValues_g[1] = 0.0000;
  imbibitionValues_g[2] = 0.0000;
  imbibitionValues_g[3] = 0.0000;
  imbibitionValues_g[4] = 0.0000;
  imbibitionValues_g[5] = 0.0000;
  imbibitionValues_g[6] = 0.0000;
  imbibitionValues_g[7] = 0.0000;
  imbibitionValues_g[8] = 0.0000;
  imbibitionValues_g[9] = 0.0000;
  imbibitionValues_g[10] = 0.03361965;
  imbibitionValues_g[11] = 0.09509072;
  imbibitionValues_g[12] = 0.17469281;
  imbibitionValues_g[13] = 0.26895718;
  imbibitionValues_g[14] = 0.37587908;
  imbibitionValues_g[15] = 0.49410588;
  imbibitionValues_g[16] = 0.62264458;
  imbibitionValues_g[17] = 0.76072577;
  imbibitionValues_g[18] = 0.90773047;
  imbibitionValues_g[19] = 1.;


  TableFunction & drainageTable_w = dynamicCast< TableFunction & >(
    *functionManager.createChild( "TableFunction", "drainageWater_swg" ) );
  drainageTable_w.setTableCoordinates( coordinates_w );
  drainageTable_w.setTableValues( drainageValues_w );
  drainageTable_w.reInitializeFunction();

  drainageTable_w.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  TableFunction & drainageTable_g = dynamicCast< TableFunction & >(
    *functionManager.createChild( "TableFunction", "drainageGas_swg" ) );
  drainageTable_g.setTableCoordinates( coordinates_g );
  drainageTable_g.setTableValues( drainageValues_g );
  drainageTable_g.reInitializeFunction();

  drainageTable_g.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  TableFunction & imbibitionTable_w = dynamicCast< TableFunction & >(
    *functionManager.createChild( "TableFunction", "imbibitionWater_swg" ) );
  imbibitionTable_w.setTableCoordinates( coordinates_w );
  imbibitionTable_w.setTableValues( imbibitionValues_w );
  imbibitionTable_w.reInitializeFunction();

  imbibitionTable_w.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  TableFunction & imbibitionTable_g = dynamicCast< TableFunction & >(
    *functionManager.createChild( "TableFunction", "imbibitionGas_swg" ) );
  imbibitionTable_g.setTableCoordinates( coordinates_g );
  imbibitionTable_g.setTableValues( imbibitionValues_g );
  imbibitionTable_g.reInitializeFunction();

  imbibitionTable_g.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  // 2) Then set up the constitutive model

  auto & relPerm = parent.registerGroup< TableRelativePermeabilityHysteresis >( name );

  auto & phaseNames = relPerm.getReference< string_array >( RelativePermeabilityBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "water"; phaseNames[1] = "gas";

  using keys = TableRelativePermeabilityHysteresis::viewKeyStruct;

  auto & drainageWaterGasTableNames = relPerm.getReference< array1d< string > >( keys::drainageWettingNonWettingRelPermTableNamesString() );
  drainageWaterGasTableNames.resize( 2 );
  drainageWaterGasTableNames[0] = "drainageWater_swg"; drainageWaterGasTableNames[1] = "drainageGas_swg";

  auto & imbibitionWaterTableName = relPerm.getReference< string >( keys::imbibitionWettingRelPermTableNameString() );
  imbibitionWaterTableName = "imbibitionWater_swg";

  auto & imbibitionGasTableName = relPerm.getReference< string >( keys::imbibitionNonWettingRelPermTableNameString() );
  imbibitionGasTableName = "imbibitionGas_swg";

  relPerm.postProcessInputRecursive();
  relPerm.initialize(); // to test all the checks
  return relPerm;
}

//todo (jacques) revise values in accordance with what's above
/*RelativePermeabilityBase & makeTableRelPermThreePhase( string const & name, Group & parent )
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

  coordinates[0][0] = 0.20000;
  coordinates[0][1] = 0.22000;
  coordinates[0][2] = 0.25000;
  coordinates[0][3] = 0.30000;
  coordinates[0][4] = 0.35000;
  coordinates[0][5] = 0.40000;
  coordinates[0][6] = 0.45000;
  coordinates[0][7] = 0.50000;
  coordinates[0][8] = 0.55000;
  coordinates[0][9] = 0.60000;
  coordinates[0][10] = 0.65000;
  coordinates[0][11] = 0.66000;
  coordinates[0][12] = 0.68000;

  real64_array values( Naxis1 );

  values[0] = 0.00000;
  values[1] = 0.00000;
  values[2] = 0.00100;
  values[3] = 0.00300;
  values[4] = 0.01000;
  values[5] = 0.01800;
  values[6] = 0.03500;
  values[7] = 0.04000;
  values[8] = 0.57000;
  values[9] = 0.08800;
  values[10] = 0.14500;
  values[11] = 0.16000;
  values[12] = 0.19000;


  TableFunction & table_ow_w = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "water_swof" ) );
  table_ow_w.setTableCoordinates( coordinates );
  table_ow_w.setTableValues( values );
  table_ow_w.reInitializeFunction();

  table_ow_w.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  array1d< real64_array > coordinates_inv;
  coordinates_inv.resize( 1 );
  coordinates_inv[0].resize( Naxis );
  //reuse coordinates
  for( int i = 0; i <Naxis1 ; ++i )
  {
    coordinates_inv[0][i] = 1. - coordinates[0][Naxis1-1-i];
  }

  values[0] = 0.00000;
  values[1] = 0.00000;
  values[2] = 0.00700;
  values[3] = 0.04500;
  values[4] = 0.09500;
  values[5] = 0.16000;
  values[6] = 0.24000;
  values[7] = 0.33000;
  values[8] = 0.44000;
  values[9] = 0.57000;
  values[10] = 0.73000;
  values[11] = 0.85000;
  values[12] = 0.85000;

  TableFunction & table_ow_o = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "oil_swof" ) );
  table_ow_o.setTableCoordinates( coordinates_inv );
  table_ow_o.setTableValues( values );
  table_ow_o.reInitializeFunction();

  table_ow_o.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  // 1.a) Second pair of phases (og)
  localIndex Naxis2 = 19;
  coordinates[0].resize( Naxis2 );

  coordinates[0][0] = 0.0;
  coordinates[0][1] = 0.005;
  coordinates[0][2] = 0.010;
  coordinates[0][3] = 0.030;
  coordinates[0][4] = 0.050;
  coordinates[0][5] = 0.100;
  coordinates[0][6] = 0.150;
  coordinates[0][7] = 0.200;
  coordinates[0][8] = 0.250;
  coordinates[0][9] = 0.300;
  coordinates[0][10] = 0.350;
  coordinates[0][11] = 0.400;
  coordinates[0][12] = 0.450;
  coordinates[0][13] = 0.500;
  coordinates[0][14] = 0.550;
  coordinates[0][15] = 0.600;
  coordinates[0][16] = 0.650;
  coordinates[0][17] = 0.700;
  coordinates[0][18] = 0.750;
  coordinates[0][19] = 0.780;

  //v1
  values.resize( Naxis2 );
  values[0] = 0.00000;
  values[1] = 0.00000;
  values[2] = 0.00200;
  values[3] = 0.00700;
  values[4] = 0.01000;
  values[5] = 0.02000;
  values[6] = 0.04000;
  values[7] = 0.07500;
  values[8] = 0.12700;
  values[9] = 0.18000;
  values[10] = 0.24000;
  values[11] = 0.31000;
  values[12] = 0.37300;
  values[13] = 0.46000;
  values[14] = 0.55000;
  values[15] = 0.64000;
  values[16] = 0.73000;
  values[17] = 0.82500;
  values[18] = 0.92000;
  values[19] = 1.00000;

  TableFunction & table_og_g = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "gas_sgof" ) );
  table_og_g.setTableCoordinates( coordinates );
  table_og_g.setTableValues( values );
  table_og_g.reInitializeFunction();

  table_og_g.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  coordinates_inv[0].resize( Naxis2 );
  //reuse coordinates
  for( int i = 0; i <Naxis2 ; ++i )
  {
    coordinates_inv[0][i] = 1. - coordinates[0][Naxis2-1-i];
  }


  values[0] = 0.00000;
  values[1] = 0.00330;
  values[2] = 0.00870;
  values[3] = 0.01420;
  values[4] = 0.01960;
  values[5] = 0.02500;
  values[6] = 0.03500;
  values[7] = 0.05000;
  values[8] = 0.07500;
  values[9] = 0.14000;
  values[10] = 0.22000;
  values[11] = 0.33000;
  values[12] = 0.45000;
  values[13] = 0.58750;
  values[14] = 0.72500;
  values[15] = 0.86250;
  values[16] = 0.91750;
  values[17] = 0.97250;
  values[18] = 1.00000;
  values[19] = 1.00000;

  TableFunction & table_og_o = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "oil_sgof" ) );
  table_og_o.setTableCoordinates( coordinates_inv );
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
}*/


class KilloughHysteresisTest : public ConstitutiveTestBase< TableRelativePermeabilityHysteresis >
{
public:

};

template< typename TBL_WRAPPER >
void testValuesAgainstPreviousImplementation( TBL_WRAPPER const & relpermTblWrapper,
                                              real64 const sat_nw,
                                              real64 const & oldImplValue_w,
                                              real64 const & oldImplValue_nw,
                                              real64 const relTol )
{
  static localIndex const MAX_PHASE = 2;

  static localIndex const ipWetting = 0;
  static localIndex const ipNonWetting = 1;

  stackArray1d< real64, MAX_PHASE > phaseVolFraction( MAX_PHASE );
  stackArray1d< real64, MAX_PHASE > phaseMaxHistoricalVolFraction( MAX_PHASE );
  stackArray1d< real64, MAX_PHASE > phaseMinHistoricalVolFraction( MAX_PHASE );
  stackArray1d< real64, MAX_PHASE > phaseRelPerm( MAX_PHASE );
  stackArray2d< real64, MAX_PHASE > dPhaseRelPerm_dPhaseVolFrac( MAX_PHASE, MAX_PHASE );
/*
  relpermTblWrapper.computeImbibitionWettingRelPerm( drainageKernel,
                                                     imbibitionKernel,
                                                     jerauld_a,
                                                     jerauld_b,
                                                     landParam,
                                                     drainageKrEndPoint, //krwd(smxd)
                                                     imbibitionKrEndPoint, //krwi(smxi)
                                                     phaseVolFraction,//s
                                                     phaseMinHistoricalVolFraction,//shy
                                                     imbibitionMinPhaseWettingVolFraction,//swc
                                                     drainagePhaseMaxVolFraction,//smxd
                                                     imbibitionPhaseMaxVolFraction,//smxi
                                                     phaseRelPerm,// kr_w(s) out intend
                                                     dPhaseRelPerm_dPhaseVolFrac );// dkr_w_dS out intend
  checkRelativeError( phaseRelPerm, oldImplValue, relTol );

  relpermTblWrapper.computeImbibitionNonWettingRelPerm( drainageRelPermKernelWrapper,
                                                        imbibitionRelPermKernelWrapper,
                                                        jerauld_a,
                                                        jerauld_b,
                                                        landParam,
                                                        imbibitionPhaseMinVolFraction, //sncri
                                                        drainagePhaseMinVolFraction,//sncrd
                                                        phaseMaxVolFraction,//smax
                                                        phaseVolFraction, //sg
                                                        phaseMaxHistoricalVolFraction, //shy
                                                        phaseRelPerm,// kr_nw(s) out intend
                                                        dPhaseRelPerm_dPhaseVolFrac ); //dkr_nw_dS ou intend
  checkRelativeError( phaseRelPerm, oldImplValue, relTol );*/

  relpermTblWrapper.computeTwoPhase( ipWetting,
                   ipNonWetting,
                   phaseVolFraction,
                   phaseMaxHistoricalVolFraction,
                   phaseMinHistoricalVolFraction,
                   phaseRelPerm,
                   dPhaseRelPerm_dPhaseVolFrac );

    checkRelativeError( phaseRelPerm[ipWetting], oldImplValue_w, relTol );
    checkRelativeError( phaseRelPerm[ipNonWetting], oldImplValue_nw, relTol );
//    checkRelativeError( dPhaseRelPerm_dPhaseVolFrac[ip], dOldImplValue_dPhaseVolFrac[ip], relTol );
}

//template< typename TBL_WRAPPER >
//void testValuesAgainstPreviousImplementation( TBL_WRAPPER const & relpermTblWrapper,
//                                              real64 const sat_w,
//                                              real64 const sat_nw,
//                                              real64 const & oldImplValue_w,
//                                              real64 const & oldImplValue_nw,
//                                              real64 const relTol )
//{
//   relpermTblWrapper.computeTwoPhase( ipWetting,
//                   ipNonWetting,
//                   phaseVolFraction,
//                   phaseMaxHistoricalVolFraction,
//                   phaseMinHistoricalVolFraction,
//                   phaseRelPerm,
//                   dPhaseRelPerm_dPhaseVolFrac );
//
//
//   relpermTblWrapper.computeThreePhase( ipWetting,
//                     ipInter,
//                     ipNonWetting,
//                     phaseVolFraction,
//                     phaseMaxHistoricalVolFraction,
//                     phaseMinHistoricalVolFraction,
//                     phaseRelPerm,
//                     dPhaseRelPerm_dPhaseVolFrac );
//
//
//  checkRelativeError( phaseRelPerm[ip], oldImplValue[ip], relTol );
//  checkRelativeError( dPhaseRelPerm_dPhaseVolFrac[ip], dOldImplValue_dPhaseVolFrac[ip], relTol );
//
//}


TEST_F( KilloughHysteresisTest, KilloughTwoPhaseHysteresisTest )
{
//  real64 const Sd0_[3] = { 0.1, 0.1, 0.1 }; // todo (jacques) check if relevant
  real64 const shy[3] = { 0.5, 0.75, 0.78 }; // todo (jacques) check if relevant

//  // snrt for shy values all in gas sat
//  real64 const Sncrt_[3] = { 0.243, 0.295, .3 }; // shoudl be in [sncrd,sncri]
//  real64 const Swcrt_[3] = { 0.75652999, 0.71058803, 0.7063253 }; //should be in [smxi,smxd]
//  // o(1e-2) is on the consistency check Sncrt_ = 1-Swcrt_
//  real64 const krwe_[3] = {};
//
//  //assert on carac from data
//  real64 const Swc = .22;
//  real64 const Swmxd = 1.;
//  real64 const Swmxi = .7; // independent but data req. consistent here
//  real64 const Sncri = .3;
//  real64 const Sncrd = .0;
//  real64 const Snmx = 1.;
//
// //compute Land coeff
//  real64 const Cw = (Swmxd - Swc)/(Swmxd - Swmxi) - 1.;
//  real64 const Cnw = (Snmx - Sncrd)/(Sncri - Sncrd) - 1;
//
//  //have default value for Jerauld coeffs
//  real64 const jerauld_a = .1;
//  real64 const jerauld_b = .0;
//  real64 const A_ = 1 + jerauld_a * ( Snmx - Shy );
//  //and alpha2
//  real64 const alpha_2 = 1.;

/// old impl data -- saved data

  localIndex const ninc = 15;
  real64 const sg_inc[] = { 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6,
                            0.65, 0.7, 0.75, 0.78 };

  localIndex const ndsc = 12;
  localIndex const off[3] = { 0, 1, 7 };
  real64 const sg_dsc[] = { 0.78, 0.73, 0.68, 0.63, 0.58, 0.53, 0.48, 0.43, 0.38, 0.33, 0.28, 0.23 };

  real64 const drainage_w_values[] = {
    0.66591, 0.53419, 0.41681, 0.32093, 0.23426, 0.16754, 0.11563, 0.07573, 0.04626, 0.02564,
    0.01229, 0.00463, 0.00108, 0.00006, 0.00000, 0.00000, 0.00000, 0.00000 };

  real64 const scanning_w_i[][ndsc] = {
    { 0.000000, 0.033542, 0.089103, 0.161722, 0.247921, 0.345789, 0.454057, 0.571801, 0.698309, 0.832241, 0.974551,
      1.127816 },
    { 0.000057, 0.012572, 0.061368, 0.131894, 0.218619, 0.318701, 0.430554, 0.553054, 0.685347, 0.826014, 0.975993,
      1.137386 },
    { 0.046258, 0.046258, 0.046258, 0.046258, 0.046258, 0.079032, 0.194995, 0.360339, 0.562881, 0.797637, 1.061088,
      1.346964 } };


//krgd

  real64 const drainage_g_values[] = { 0.02000, 0.04000, 0.07500, 0.12700, 0.18000, 0.24000, 0.31000,
                                       0.37300, 0.46000, 0.55000, 0.64000, 0.73000, 0.82500, 0.92000, 1.00000 };


  real64 const scanning_g_i[][ndsc] = { //imbibtion curve for Shy = .78
    { 1.000000, 0.848929, 0.705493, 0.571229, 0.446815, 0.333110, 0.231251, 0.142852, 0.070502, 0.020172, 0.000000,
      0.000000 },
    { 0.920000, 0.845672, 0.702498, 0.568606, 0.444687, 0.331435, 0.229982, 0.141933, 0.069856, 0.019582, 0.000000,
      0.000000 },
    { //imbibtion curve for Shy = .5
      0.460000, 0.460000, 0.460000, 0.460000, 0.460000, 0.390993, 0.270810, 0.166261, 0.080132, 0.018934, 0.000000,
      0.00000 } };


  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const relTol = 5e-5;

//saved cycle
  initialize( makeTableRelPermHysteresisTwoPhase("relPerm", m_parent) );
//  TableRelativePermeability::KernelWrapper // ?
  auto relpermTblWrapper = m_model->createKernelWrapper();
  // all sat are nonwetting sat
  localIndex counter = 0;
  for( localIndex count = 0; count < 3; ++counter )
  {
    //drainage
    for( localIndex isat = 0; isat < ninc; ++isat )
    {
      if( shy[isat] < sg_inc[isat] ) break;

      testValuesAgainstPreviousImplementation( relpermTblWrapper,
                                               sg_inc[isat], drainage_w_values[isat], drainage_g_values[isat],
                                               relTol );
//
//      testNumericalDerivatives( relpermTblWrapper, ds[isat], xxxx, comp,
//      true, eps, relTol );

    }
    //imbibition
    for( localIndex isat = off[counter]; isat < ndsc; ++isat )
    {
      testValuesAgainstPreviousImplementation( relpermTblWrapper,
                                               sg_dsc[isat], scanning_w_i[isat][counter], scanning_g_i[isat][counter],
                                               relTol );

//      testNumericalDerivatives( relpermTblWrapper, ds[isat], xxxx, comp,
//      true, eps, relTol );
    }
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
