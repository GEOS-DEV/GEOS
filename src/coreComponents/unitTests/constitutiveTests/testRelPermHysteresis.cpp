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
#include "common/DataTypes.hpp"
#include "constitutiveTestHelpers.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"

using namespace geos;
using namespace geos::testing;
using namespace geos::constitutive;
using namespace geos::constitutive::relperm;
using namespace geos::dataRepository;

// The point of this unit test is to carefully trace drainage-imbibition cycles and
// verify that hysteresis works well for relative permeability and ultimately, capillary pressure.

// For now, we verify that relative permeability hysteresis implementation matches the expected relative permeability hysteresis values
// by doing an imbibition cycles for the data provided "Reservoir simulation with history-dependent saturation functions",
// Killough, J. E., Society of Petroleum Engineers Journal, 16(01), 37-48 (1976).

TableRelativePermeabilityHysteresis & makeTableRelPermHysteresisTwoPhase( string const & name, Group & parent )
{
  // 1) First, define the tables (to values that matters for our use cases)

  // First, define the saturation ranges

  // Water phase saturation, first column of Table 2
  array1d< real64_array > coordinates_dw;
// Swc = 0.22
// consistent with Swmaxd = 1-Sgcrd = 1-0 = 1
  geos::testing::fillArray( coordinates_dw,
                            {.22, .25, .3, .35, .4, .45, .5, .55, .6, .65, .66, .68, .72, .82, .91, 1.} );

  array1d< real64_array > coordinates_iw;
// Swc = 0.22
// consistent with Swmaxi = 1-Sgcri = 1-0.3 = 0.7
  geos::testing::fillArray( coordinates_iw, {.22, .25, .3, .35, .4, .45, .5, .55, .6, .65, .66, .7} );

  // Gas phase saturation, fifth column of Table 2
  array1d< real64_array > coordinates_dg;
  // Sgcrd = 0.0
  // consistent with Swc = 0.22
  geos::testing::fillArray( coordinates_dg,
                            {0.000, 0.010, 0.030, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450,
                             0.500,
                             0.550, 0.600, 0.650, 0.700, 0.750, 0.780} );

  array1d< real64_array > coordinates_ig;
  // Sgcri = 0.30;
  // consistent with Swc = 0.22
  geos::testing::fillArray( coordinates_ig,
                            {0.300, 0.350, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.78} );


  // then define the bounding drainage and imbibibition relative permeability

  // Water phase drainage relperm
  real64_array drainageValues_w;
  geos::testing::fillArray( drainageValues_w, {0.00000, 0.00100, 0.00300, 0.01000, 0.01800, 0.03500, 0.04000, 0.05700,
                                               0.08800, 0.14500, 0.16000, 0.19000, 0.26300, 0.45500, 0.69200, 1.} );
  // Gas phase drainage relperm, seventh column of Table 2
  real64_array drainageValues_g;
  geos::testing::fillArray( drainageValues_g, {0.00000, 0.00200, 0.00700, 0.01000, 0.02000, 0.04000, 0.07500,
                                               0.12700, 0.18000, 0.24000, 0.31000, 0.37300, 0.46000, 0.55000,
                                               0.64000, 0.73000, 0.82500, 0.92000, 1.00000} );

  real64_array imbibitionValues_w;
  geos::testing::fillArray( imbibitionValues_w, {0, 0.0156, 0.0680, 0.1409, 0.2296, 0.3317, 0.4455, 0.5700,
                                                 0.7044, 0.8479, 0.8776, 0.9382} );

  real64_array imbibitionValues_g;
  geos::testing::fillArray( imbibitionValues_g, {0.0000, 0.03361965, 0.09509072, 0.17469281, 0.26895718,
                                                 0.37587908, 0.49410588, 0.62264458, 0.76072577, 0.90773047, 1.} );

  initializeTable( "drainageWater_swg",
                   coordinates_dw,
                   drainageValues_w );
  initializeTable( "drainageGas_swg",
                   coordinates_dg,
                   drainageValues_g );
  initializeTable( "imbibitionWater_swg",
                   coordinates_iw,
                   imbibitionValues_w );
  initializeTable( "imbibitionGas_swg",
                   coordinates_ig,
                   imbibitionValues_g );

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

class KilloughHysteresisTest : public ConstitutiveTestBase< TableRelativePermeabilityHysteresis >
{
public:

};

template< typename TBL_WRAPPER >
void testValuesAgainstReference( TBL_WRAPPER const & relpermTblWrapper,
                                 real64 const & sat_nw,
                                 real64 const & shy,
                                 real64 const & refRelPerm_w,
                                 real64 const & refRelPerm_nw,
                                 real64 const & relTol )
{
  integer const numPhases = 2;

  integer const ipWetting = 0;
  integer const ipNonWetting = 1;

  StackArray< real64, 2, constitutive::RelativePermeabilityBase::MAX_NUM_PHASES,
              compflow::LAYOUT_PHASE > phaseVolFraction( 1, numPhases );
  phaseVolFraction[0][0] = 1. - sat_nw;
  phaseVolFraction[0][1] = sat_nw;

  StackArray< real64, 2, constitutive::RelativePermeabilityBase::MAX_NUM_PHASES,
              compflow::LAYOUT_PHASE > phaseMaxHistoricalVolFraction( 1, numPhases );
  phaseMaxHistoricalVolFraction[0][0] = 1.;
  phaseMaxHistoricalVolFraction[0][1] = shy;

  StackArray< real64, 2, constitutive::RelativePermeabilityBase::MAX_NUM_PHASES,
              compflow::LAYOUT_PHASE > phaseMinHistoricalVolFraction( 1, numPhases );
  phaseMinHistoricalVolFraction[0][0] = 1. - shy;
  phaseMinHistoricalVolFraction[0][1] = 0.;

  StackArray< real64, 3, constitutive::RelativePermeabilityBase::MAX_NUM_PHASES,
              relperm::LAYOUT_RELPERM > phaseTrappedVolFrac( 1, 1, numPhases );

  StackArray< real64, 3, constitutive::RelativePermeabilityBase::MAX_NUM_PHASES,
              relperm::LAYOUT_RELPERM > phaseRelPerm( 1, 1, numPhases );

  StackArray< real64, 4, constitutive::RelativePermeabilityBase::MAX_NUM_PHASES *constitutive::RelativePermeabilityBase::MAX_NUM_PHASES,
              relperm::LAYOUT_RELPERM_DS > dPhaseRelPerm_dPhaseVolFrac( 1, 1, numPhases, numPhases );

  relpermTblWrapper.computeTwoPhase( ipWetting,
                                     ipNonWetting,
                                     phaseVolFraction[0],
                                     phaseMaxHistoricalVolFraction[0],
                                     phaseMinHistoricalVolFraction[0],
                                     phaseTrappedVolFrac[0][0],
                                     phaseRelPerm[0][0],
                                     dPhaseRelPerm_dPhaseVolFrac[0][0] );

  checkRelativeError( phaseRelPerm[0][0][ipWetting], refRelPerm_w, relTol );
  checkRelativeError( phaseRelPerm[0][0][ipNonWetting], refRelPerm_nw, relTol );
  //TODO check phaseTrappedVolFraction
}

TEST_F( KilloughHysteresisTest, KilloughTwoPhaseHysteresisTest )
{
  real64 const shy[3] = { 0.5, 0.75, 0.78 };

  localIndex const nIncreasingGasSat = 15;
  real64 const increasingGasSat[] = { 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.78 };

  localIndex const nDecreasingGasSat = 12;
  localIndex const offset[3] = { 6, 1, 0 };
  real64 const decreasingGasSat[] = { 0.78, 0.73, 0.68, 0.63, 0.58, 0.53, 0.48, 0.43, 0.38, 0.33, 0.28, 0.23 };

  /// Reference data
  real64 const drainageRelPerm_w_values[] =
  { 0.66567, 0.53400, 0.41660, 0.32060, 0.22650, 0.14500, 0.08800, 0.05700, 0.04, 0.035,
    0.018, 0.01, 0.003, 0.001, 0.00000, 0.00000, 0.00000, 0.00000 };

  real64 const drainageRelPerm_g_values[] =
  { 0.02000, 0.04000, 0.07500, 0.12700, 0.18000, 0.24000, 0.31000,
    0.37300, 0.46000, 0.55000, 0.64000, 0.73000, 0.82500, 0.92000, 1.00000 };

  real64 const scanningRelPerm_w_i[][nDecreasingGasSat] =
  { { 0.046258, 0.046258, 0.046258, 0.046258, 0.046258, 0.046258, 0.0616975, 0.172378, 0.337489, 0.541488, 0.778040, 0.910729 },
    { 0.000057, 0.011931, 0.062327, 0.136442, 0.228608, 0.336255, 0.4569098, 0.589277, 0.732427, 0.879957, 0.935521, 0.935521 },
    { 0.0, 0.036560, 0.097156, 0.176380, 0.270440, 0.377219, 0.4953000, 0.623760, 0.761800, 0.892749, 0.938200, 0.938200 } };

  real64 const scanningRelPerm_g_i[][nDecreasingGasSat] =
  { { 0.460000, 0.460000, 0.460000, 0.460000, 0.460000, 0.460000, 0.407517, 0.285259, 0.178848, 0.090686, 0.025844, 0.000000 },
    { 0.920000, 0.860285, 0.716467, 0.581413, 0.456038, 0.341595, 0.238939, 0.149497, 0.075697, 0.022777, 0.000000, 0.000000 },
    { 1.000000, 0.848929, 0.705493, 0.571229, 0.446815, 0.333110, 0.231251, 0.142852, 0.070502, 0.020172, 0.000000, 0.000000 } };

  real64 const relTol = 5e-5;

  // saved cycle
  initialize( makeTableRelPermHysteresisTwoPhase( "relPerm", m_parent ) );
  auto relpermTblWrapper = m_model->createKernelWrapper();

  // all saturations are nonwetting saturations
  for( integer count = 0; count < 3; ++count )
  {
    // drainage
    for( integer iSat = 0; iSat < nIncreasingGasSat; ++iSat )
    {
      if( increasingGasSat[iSat] > shy[count] )
      {
        break; // exit as the drainage bounding curve is not scanned anymore
      }
      testValuesAgainstReference( relpermTblWrapper,
                                  increasingGasSat[iSat], increasingGasSat[iSat],
                                  drainageRelPerm_w_values[iSat], drainageRelPerm_g_values[iSat],
                                  relTol );
    }

    // imbibition
    for( integer iSat = offset[count]; iSat < nDecreasingGasSat; ++iSat )
    {
      testValuesAgainstReference( relpermTblWrapper,
                                  decreasingGasSat[iSat], shy[count],
                                  scanningRelPerm_w_i[count][iSat], scanningRelPerm_g_i[count][iSat],
                                  relTol );
    }
  }
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::GeosxState state( geos::basicSetup( argc, argv ) );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
