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

using namespace geosx;
using namespace geosx::testing;
using namespace geosx::constitutive;
using namespace geosx::constitutive::cappres;
using namespace geosx::dataRepository;

CommandLineOptions g_commandLineOptions;
char const *xmlInput =
        "<Problem>\n"
        "  <Mesh>\n"
        "    <InternalMesh name=\"mesh\"\n"
        "                  elementTypes=\"{C3D8}\" \n"
        "                  xCoords=\"{0, 1}\"\n"
        "                  yCoords=\"{0, 1}\"\n"
        "                  zCoords=\"{0, 1}\"\n"
        "                  nx=\"{1}\"\n"
        "                  ny=\"{1}\"\n"
        "                  nz=\"{1}\"\n"
        "                  cellBlockNames=\"{cb1}\"/>\n"
        "  </Mesh>\n"
        "  <ElementRegions>\n"
        "    <CellElementRegion name=\"region\" cellBlocks=\"{cb1}\" materialList=\"{ cappres }\" />\n"
        "  </ElementRegions>\n"
        "  <Constitutive>\n"
        "    <TableCapillaryPressureHysteresis\n"
        "        name=\"cappres\"\n"
        "        phaseNames=\"{ gas, water }\"\n"
        "        KilloughModelName=\"KilloughHyst\"\n"
        "        drainageWettingNonWettingCapPressureTableName=\"dummyA\"\n"
        "        imbibitionWettingNonWettingCapPressureTableName=\"dummyB\"\n"
        "        />\n"
        "    <KilloughHysteresis\n"
        "        name=\"KilloughHyst\" />\n"
        " </Constitutive>\n"
        " <Functions>\n"
        "    <TableFunction name=\"dummyA\"\n"
        "                   coordinates=\"{0.0, 1.0}\"\n"
        "                   values=\"{1.0, 0.}\"/>\n"
        "    <TableFunction name=\"dummyB\"\n"
        "                   coordinates=\"{0.0, 1.0}\"\n"
        "                   values=\"{1.0, 0.}\"/>\n"
        "</Functions>\n"
        "</Problem>";


TableCapillaryPressureHysteresis& makeTableCapPresHysteresisTwoPhase(string const &name, Group &constbase)
{
    // 1) First, define the tables
    array1d< real64_array > coordinates_dw;
    real64_array drainageValues_w;
    geosx::testing::fill_array(coordinates_dw,{0.2,0.278,0.356,0.433,0.511,0.589,0.667,0.744,0.822,0.9});
    geosx::testing::fill_array(drainageValues_w,{2236.068, 1897.367, 1677.051, 1519.109, 1398.757, 1303.117,
       1224.745, 1159.001, 1102.822, 1054.093});
    array1d< real64_array > coordinates_iw;
    real64_array imbibitionValues_w;
    geosx::testing::fill_array(coordinates_iw,{0.2,0.233,0.267,0.3,0.333,0.367,0.4,0.433,0.467,0.5 });
    geosx::testing::fill_array(imbibitionValues_w,{2236.068, 1897.367, 1677.051, 1519.109, 1398.757, 1303.117,
                                   1224.745, 1159.001, 1102.822, 1054.093});

    initializeTable( "drainageWater_swg",
                     coordinates_dw,
                     drainageValues_w );

    initializeTable( "imbibitionWater_swg",
                     coordinates_iw,
                     imbibitionValues_w );

// 2) Then set up the constitutive model


    auto& capPres = constbase.getGroupByPath("/Problem/domain/Constitutive").getGroup<TableCapillaryPressureHysteresis>(name);

    using keys = TableCapillaryPressureHysteresis::viewKeyStruct;

    auto & drainageWaterGasTableNames = capPres.getReference< string >(keys::drainageWettingNonWettingCapPresTableNameString() );
    drainageWaterGasTableNames = "drainageWater_swg";

    auto & imbibitionWaterGasTableName = capPres.getReference< string >(keys::imbibitionWettingNonWettingCapPresTableNameString());
    imbibitionWaterGasTableName = "imbibitionWater_swg";

    capPres.postProcessInputRecursive();
    capPres.initialize(); // to test all the checks
    return capPres;
}

//TableCapillaryPressureHysteresis& makeTableCapPresHysteresisThreePhase( string const & name, Group & parent )
//{
//
//}


/// class

class KilloughHysteresisTest : public ::testing::Test
{
public:

    KilloughHysteresisTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
    {}

protected:
     void SetUp() override {
         setupProblemFromXML(state.getProblemManager(), xmlInput);
         DomainPartition &domain = state.getProblemManager().getDomainPartition();

    }

     GeosxState state;

};


/// templated comparator

template< typename TBL_WRAPPER >
void testValuesAgainstReference( TBL_WRAPPER const & cappresTblWrapper,
                                 real64 const & sat_w,
                                 real64 const & shy,
                                 real64 const & refCapPres_w,
                                 real64 const & relTol )
{
    integer const numPhases = 2;

    integer const ipWetting = 0;
    integer const ipNonWetting = 1;

    StackArray< real64, 2, constitutive::CapillaryPressureBase::MAX_NUM_PHASES,
            compflow::LAYOUT_PHASE > phaseVolFraction( 1, numPhases );
    phaseVolFraction[0][0] = sat_w;
    phaseVolFraction[0][1] = 1. - sat_w;

    StackArray< real64, 2, constitutive::CapillaryPressureBase::MAX_NUM_PHASES,
            compflow::LAYOUT_PHASE > phaseMaxHistoricalVolFraction( 1, numPhases );
    phaseMaxHistoricalVolFraction[0][0] = 1.;
    phaseMaxHistoricalVolFraction[0][1] = 1. - shy;

    StackArray< real64, 2, constitutive::CapillaryPressureBase::MAX_NUM_PHASES,
            compflow::LAYOUT_PHASE > phaseMinHistoricalVolFraction( 1, numPhases );
    phaseMinHistoricalVolFraction[0][0] = shy;
    phaseMinHistoricalVolFraction[0][1] = 0.;

    StackArray< real64, 3, constitutive::CapillaryPressureBase::MAX_NUM_PHASES,
            relperm::LAYOUT_RELPERM > phaseTrappedVolFrac( 1, 1, numPhases );

    StackArray< real64, 3, constitutive::CapillaryPressureBase::MAX_NUM_PHASES,
            relperm::LAYOUT_RELPERM > phaseCapPres(1, 1, numPhases );

    StackArray< real64, 4, constitutive::CapillaryPressureBase::MAX_NUM_PHASES * constitutive::CapillaryPressureBase::MAX_NUM_PHASES,
            relperm::LAYOUT_RELPERM_DS > dPhaseCapPres_dPhaseVolFrac(1, 1, numPhases, numPhases );

    cappresTblWrapper.computeTwoPhase(ipWetting,
                                      ipNonWetting,
                                      phaseVolFraction[0],
                                      phaseMaxHistoricalVolFraction[0],
                                      phaseMinHistoricalVolFraction[0],
                                      phaseTrappedVolFrac[0][0],
                                      phaseCapPres[0][0],
                                      dPhaseCapPres_dPhaseVolFrac[0][0] );

    checkRelativeError(phaseCapPres[0][0][ipWetting], refCapPres_w, relTol );
    //TODO check phaseTrappedVolFraction
}


/// gtest magicwand
TEST_F( KilloughHysteresisTest, KilloughTwoPhaseHysteresisTest )
{
    //define eps to avoid falling in the out of range case
    real64 const eps = 0.00001;
    real64 const shy[3] = { 0.2 + eps, 0.3, 0.4 };
    // which correspond resp to 0.4, 0.36, 0.32 trapped sat

    localIndex const nDecreasingWaterSat = 11;
    real64 const decreasingWaterSat[] = { 0.9 - eps, 0.8, 0.7, 0.65, 0.6, 0.55, 0.5, 0.4, 0.3, 0.25, 0.2 + eps };

    localIndex const nIncreasingWaterSat = 11;
    localIndex const offset[3] = { 0, 2, 3 };
    real64 const increasingWaterSat[] = { 0.2 + eps , 0.25, 0.3 , 0.4 , 0.5 , 0.55, 0.6 , 0.65, 0.7 , 0.8 , 0.9 - eps  };

    /// Reference data
    real64 const drainageCapPres_w_values[] =
            { 1054.093, 1118.667, 1196.569, 1241.826, 1292.065, 1350.937, 1415.730 , 1586.798, 1835.227 , 2018.952, 2236.068  };

    real64 const scanningCapPres_w_i[][nIncreasingWaterSat] =
            { {2236.068, 1929.0456, 1651.251, 1305.843, 1100.016, 1082.189, 1070.461, 1063.032, 1058.364, 1054.239, 1054.093},
              {1835.226, 1835.226, 1835.226, 1364.184, 1119.252, 1090.195, 1072.584, 1062.267, 1056.403, 1054.093, 1054.093},
              { 1586.7984285714285, 1586.7984285714285, 1586.7984285714285, 1586.7984285714285, 1168.670, 1107.583, 1075.298, 1058.607, 1054.093, 1054.093, 1054.093} };

    real64 const relTol = 5e-5;

    // saved cycle
    TableCapillaryPressureHysteresis&  table = makeTableCapPresHysteresisTwoPhase( "cappres",
                                                                                   state.getProblemManager().getDomainPartition().getConstitutiveManager() );
//    initialize( table );
    auto cappresTblWrapper = table.createKernelWrapper();

    // all saturations are wetting saturations
    auto ncycles = 3;
    for( integer count = 0; count < ncycles; ++count )
    {
        // drainage
        for(integer iSat = 0; iSat < nDecreasingWaterSat; ++iSat )
        {
            if(decreasingWaterSat[iSat] < shy[count] )
            {
                break; // exit as the drainage bounding curve is not scanned anymore
            }
            testValuesAgainstReference(cappresTblWrapper,
                                       decreasingWaterSat[iSat], decreasingWaterSat[iSat],
                                       drainageCapPres_w_values[iSat],
                                       relTol );
        }

        // imbibition
        for( integer iSat = offset[count]; iSat < nIncreasingWaterSat; ++iSat )
        {
            testValuesAgainstReference(cappresTblWrapper,
                                       increasingWaterSat[iSat], shy[count],
                                       scanningCapPres_w_i[count][iSat],
                                       relTol );
        }
    }
}

/// main
int main( int argc, char * * argv )
{
    ::testing::InitGoogleTest( &argc, argv );
    g_commandLineOptions = *geosx::basicSetup( argc, argv );
    int const result = RUN_ALL_TESTS();

    geosx::basicCleanup();

    return result;
}

