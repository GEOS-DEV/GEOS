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

#include <functional>
#include "codingUtilities/UnitTestUtilities.hpp"
#include "common/DataTypes.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/GeosxState.hpp"
#include "functions/FunctionManager.hpp"

#include "physicsSolvers/fluidFlow/wells/PipeFlowTableFunction.hpp"

using namespace geos;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

char const * xmlInput =
  R"xml(
<Problem>
<Functions>
<PipeFlowTableFunction
name="1"
wellHeadPressure = "{ 2000000.000 , 6000000.000 , 10000000.000 }"
gasFractionType = "gor"
gfr = "{ 550.000 }"
waterFractionType = "wct"
wfr = "{ 0.000 , 0.400 , 0.800 }"
rateType = "liq"
rate = "{ 0.006 , 0.012 , 0.023 , 0.035 , 0.046 , 0.058 }"
bottomHolePressure = "{ 12886800.000 , 10877900.000 , 9446100.000 , 9395970.000 , 9689360.000 , 10053300.000 , 19004500.000 , 18928100.000 , 19006100.000 , 19124500.000
                       , 19304000.000 , 19595800.000 , 23433800.000 , 23473300.000 , 23620700.000 , 23853100.000 , 24169600.000 , 24568500.000 , 15858800.000 , 15440700.000
                       , 14512800.000 , 14066000.000 , 14120300.000 , 14381000.000 , 21296900.000 , 21267300.000 , 21357200.000 , 21571800.000 , 21882700.000 , 22284400.000
                       , 25605100.000 , 25642400.000 , 25781500.000 , 26002200.000 , 26302800.000 , 26682400.000 , 19488700.000 , 19341100.000 , 19341600.000 , 19527300.000
                       , 19840900.000 , 20197400.000 , 24098100.000 , 24107700.000 , 24232500.000 , 24434100.000 , 24719800.000 , 25081700.000 , 28250600.000 , 28285800.000
                       , 28416200.000 , 28623000.000 , 28904900.000 , 29261400.000 }"
/>
</Functions>
</Problem>
  )xml";

auto makeLine = []( double slope, double intercept )
{
  return [slope, intercept]( double x ) { return slope * x + intercept; };
};
auto makeInverseLine = []( double slope, double intercept )
{
  return [slope, intercept]( double y ) { return (y - intercept) / slope; };
};
void testfindIPR_VFPIntersection(   )
{
  // Define the rates as functions of bottom hole pressure
  auto oilRate =  makeLine( -2.29129742E-08, 4.66732788E-01 );
  auto gasRate =  makeLine( -4.46909867E-09, 9.10293585E-02 );
  auto waterRate =  makeLine( -5.379855053E-12, 1.096008746E-04 );
  // well BHP as a function of liquid rate
  //auto wellBHPFromLiq = makeInverseLine( -4.474478524E-09, 9.113895942E-02 );
  auto wellBHPFromLiq =  makeInverseLine( -2.29129742E-08, 4.66732788E-01 );

  const double whpMin = 7.80e6;
  const double bhpMin = 1.90e7;

  // Get the flow table function
  std::string tableName="1";
  FunctionManager & functionManager = FunctionManager::getInstance();
  PipeFlowTableFunction & m_flowTable =  functionManager.getGroup< PipeFlowTableFunction >( tableName );
  m_flowTable.postInputInitialization();
  integer flowTableSolveState;
  // Liquid constraint is used to find intersection of IPR and VLP
  const array1d< real64 > & liquidRates = m_flowTable.getRates();
  array1d< real64 > phaseRates( 3 );
  double tableBHP0;
  for( integer i=0; i <  liquidRates.size(); ++i )
  {
    // get well ipr
    double liqRate0 = liquidRates[i];
    double wellbhp = wellBHPFromLiq( liqRate0 );
    phaseRates[0] = oilRate( wellbhp );
    phaseRates[1] = gasRate( wellbhp );
    phaseRates[2] = waterRate( wellbhp );

    std::cout << phaseRates << std::endl;
    ;

    m_flowTable.calculateBHP( phaseRates, whpMin, tableBHP0, flowTableSolveState );
    std::cout << wellbhp << "," << tableBHP0 << "," << whpMin << "," << liqRate0 << "," <<  phaseRates[0] << "," <<  phaseRates[1] << " ," <<  phaseRates[2] << std::endl;
  }
  const integer nRates = liquidRates.size();

  array1d< real64 > phaseRates0( 3 );
  //double tableBHP0;
  double liqRate0 = liquidRates[nRates-2];
  double wellbhp0 = wellBHPFromLiq( liqRate0 );
  phaseRates0[0] = oilRate( wellbhp0 );
  phaseRates0[1] = gasRate( wellbhp0 );
  phaseRates0[2] = waterRate( wellbhp0 );
  m_flowTable.calculateBHP( phaseRates0, whpMin, tableBHP0, flowTableSolveState );
  double q0 = wellbhp0 - tableBHP0;

  array1d< real64 > phaseRates1( 3 );
  double tableBHP1;
  double liqRate1 = liquidRates[nRates-1];
  double wellbhp1 = wellBHPFromLiq( liqRate1 );
  phaseRates1[0] = oilRate( wellbhp1 );
  phaseRates1[1] = gasRate( wellbhp1 );
  phaseRates1[2] = waterRate( wellbhp1 );
  m_flowTable.calculateBHP( phaseRates1, whpMin, tableBHP1, flowTableSolveState );
  double q1 = wellbhp1 - tableBHP1;
  ;

  integer const maxIters=100;
  real64 const tol = 1e-6;
  integer iter = 0;
  std::ofstream of;
  of.open( "fl.csv" );

  of << " error,wellbhp,tablebhp,whp,orate,grate,wrate "  << std::endl;
  real64 wellbhpN;

  while( iter < maxIters ) //&& std::abs( wellbhp1 - tableBHP )/wellbhp1 > tol )
  {

    if( true )
    {
      wellbhp1 = 0.5 *(tableBHP1+wellbhp1);
      phaseRates1[0] = oilRate( wellbhp1 );
      phaseRates1[1] = gasRate( wellbhp1 );
      phaseRates1[2] = waterRate( wellbhp1 );
      m_flowTable.calculateBHP( phaseRates1, whpMin, tableBHP1, flowTableSolveState );
      if( std::abs( wellbhp1-tableBHP1 )/wellbhp1 < tol )
      {
        break;
      }
      of << std::abs( wellbhp1-tableBHP1 ) << ","<< wellbhp1 << "," << tableBHP1<< ","<< whpMin <<  "," << phaseRates1[0] << "," << phaseRates1[1] << " ," << phaseRates1[2] <<
        std::endl;
      std::cout << std::abs( wellbhp1-tableBHP1 ) << ","<<wellbhp1 << "," << tableBHP1<< ","<< whpMin <<  "," << phaseRates1[0] << "," << phaseRates1[1] << " ," <<
        phaseRates1[2] << std::endl;
    }
    else
    {
      real64 liqrate = liqRate1 - 0.5*(q1* (liqRate1-liqRate0)/(q1-q0));
      if( std::abs( liqrate-liqRate1 ) < tol )
      {
        break;
      }
      liqRate0=liqRate1;
      q0=q1;
      liqRate1=liqrate;
      wellbhp1 = wellBHPFromLiq( liqRate1 );
      phaseRates1[0] = oilRate( wellbhp1 );
      phaseRates1[1] = gasRate( wellbhp1 );
      phaseRates1[2] = liqRate1-phaseRates1[0];

      wellbhp0=wellbhp1;
      tableBHP0 = tableBHP1;
      m_flowTable.calculateBHP( phaseRates1, whpMin, tableBHP1, flowTableSolveState );
      q1 = wellbhp1 - tableBHP1;
      of << std::abs( q1 ) << "," << std::abs( liqrate-liqRate1 ) << ","<< wellbhp1 << "," << tableBHP1<< ","<< whpMin <<  "," << phaseRates1[0] << "," << phaseRates1[1] << " ," << phaseRates1[2] <<
        std::endl;
      std::cout << std::abs( q1 ) << "," <<  std::abs( liqrate-liqRate1 ) << ","<<wellbhp1 << "," << tableBHP1<< ","<< whpMin <<  "," << phaseRates1[0] << "," << phaseRates1[1] << " ," <<
        phaseRates1[2] << std::endl;

      ++iter;
    }



  }
  of.close();

}

class VFPTabelTest : public ::testing::Test
{
public:

  VFPTabelTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:

  void SetUp() override
  {
    ProblemManager & problemManager = state.getProblemManager();
    xmlWrapper::xmlDocument xmlDocument;
    xmlWrapper::xmlResult xmlResult = xmlDocument.loadString( xmlInput );
    if( !xmlResult )
    {
      GEOS_LOG_RANK_0( "XML parsed with errors!" );
      GEOS_LOG_RANK_0( "Error description: " << xmlResult.description());
      GEOS_LOG_RANK_0( "Error offset: " << xmlResult.offset );
    }

    int mpiSize = MpiWrapper::commSize( MPI_COMM_GEOS );

    dataRepository::Group & commandLine =
      problemManager.getGroup< dataRepository::Group >( problemManager.groupKeys.commandLine );

    commandLine.registerWrapper< integer >( problemManager.viewKeys.xPartitionsOverride.key() ).
      setApplyDefaultValue( mpiSize );

    xmlWrapper::xmlNode xmlProblemNode = xmlDocument.getChild( dataRepository::keys::ProblemManager );
    problemManager.processInputFileRecursive( xmlDocument, xmlProblemNode );

    DomainPartition & domain = problemManager.getDomainPartition();
  }


  GeosxState state;

};

TEST_F( VFPTabelTest, ipr_vfp_intersect )
{
  testfindIPR_VFPIntersection();
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
