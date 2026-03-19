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

// Source includes
#include "common/format/table/TableData.hpp"
#include "common/format/table/TableMpiComponents.hpp"
#include "common/initializeEnvironment.hpp"
#include "common/MpiWrapper.hpp"
// TPL includes
#include <gtest/gtest.h>
#include <gtest/gtest-spi.h>

using namespace geos;

class MpiTestScope
{
public:

  MpiTestScope( int argc, char * argv[] )
  {
    ::testing::InitGoogleTest( &argc, argv );
    geos::setupEnvironment( argc, argv );
  }

  ~MpiTestScope()
  {
    geos::cleanupEnvironment();
  }

};

TEST( testMpiTables, testDifferentRankData )
{
  struct TestCase
  {
    stdVector< stdVector< std::pair< integer, real64 > > > m_ranksValues;
    string m_expectedResult;
  };

  stdVector< TestCase > const testCases =
  {
    {
      { // m_ranksValues: in this test, rank 2 has no value
        { {1, 0.502} },
        { {2, 0.624}, {3, 0.791} },
        {},
        { {4, 0.243}, {5, 0.804}, {6, 0.302} },
      },
      "\n" // m_expectedResult
      "-------------------------------------------\n"
      "|  Summary of negative pressure elements  |\n"
      "|-----------------------------------------|\n"
      "|    Global Id     |    pressure [Pa]     |\n"
      "|------------------|----------------------|\n"
      "|------------Rank 0, 1 values-------------|\n"
      "|               1  |               0.502  |\n"
      "|------------Rank 1, 2 values-------------|\n"
      "|               2  |               0.624  |\n"
      "|               3  |               0.791  |\n"
      "|------------Rank 3, 3 values-------------|\n"
      "|               4  |               0.243  |\n"
      "|               5  |               0.804  |\n"
      "|               6  |               0.302  |\n"
      "-------------------------------------------\n"
    },
    { // m_ranksValues: in this test, rank 0 has no value
      {
        {},
        { {4, 0.243}, {5, 0.804}, {6, 0.302} },
        { {1, 0.502} },
        { {2, 0.624}, {3, 0.791} },
      },
      "\n" // m_expectedResult
      "-------------------------------------------\n"
      "|  Summary of negative pressure elements  |\n"
      "|-----------------------------------------|\n"
      "|    Global Id     |    pressure [Pa]     |\n"
      "|------------------|----------------------|\n"
      "|------------Rank 1, 3 values-------------|\n"
      "|               4  |               0.243  |\n"
      "|               5  |               0.804  |\n"
      "|               6  |               0.302  |\n"
      "|------------Rank 2, 1 values-------------|\n"
      "|               1  |               0.502  |\n"
      "|------------Rank 3, 2 values-------------|\n"
      "|               2  |               0.624  |\n"
      "|               3  |               0.791  |\n"
      "-------------------------------------------\n"
    },
  };
  for( TestCase const & testCase: testCases )
  {
    int const rankId = MpiWrapper::commRank();
    int const nbRanks = MpiWrapper::commSize();
    if( nbRanks > 1 )
    {
      ASSERT_EQ( nbRanks, 4 );

      TableLayout const layout = TableLayout().
                                   setTitle( "Summary of negative pressure elements" ).
                                   addColumns( { "Global Id", "pressure [Pa]" } ).
                                   setDefaultHeaderAlignment( TableLayout::Alignment::left );
      TableData data;
      auto const & rankTestData = testCase.m_ranksValues[rankId];

      TableMpiLayout mpiLayout;
      mpiLayout.m_separatorBetweenRanks = true;

      if( !rankTestData.empty() )
      {
        mpiLayout.m_rankTitle = GEOS_FMT( "Rank {}, {} values", rankId, rankTestData.size() );
        for( auto const & [id, value] : rankTestData )
        {
          data.addRow( id, value );
        }
      }

      TableTextMpiFormatter const formatter = TableTextMpiFormatter( layout, mpiLayout );
      std::ostringstream oss;
      formatter.toStream( oss, data );
      if( rankId == 0 )
      {
        EXPECT_STREQ( testCase.m_expectedResult.data(),
                      oss.str().data() );
      }
    }
  }
}

TEST( testMpiTables, testSortingMethod )
{
  struct TestCase
  {
    stdVector< stdVector< std::pair< integer, real64 > > > m_ranksValues;
    string m_expectedResult;
  };

  TestCase const testCase =
  {
    {   // m_ranksValues: in this test, rank 2 has no value
      { {1, 0.502} },
      { {2, 0.624}, {3, 0.791} },
      {},
      { {4, 0.243}, {5, 0.804}, {6, 0.302} },
    },
    "\n"   // m_expectedResult
    "-------------------------------------------\n"
    "|  Summary of negative pressure elements  |\n"
    "|-----------------------------------------|\n"
    "|    Global Id     |    pressure [Pa]     |\n"
    "|------------------|----------------------|\n"
    "|               1  |               0.502  |\n"
    "|               2  |               0.624  |\n"
    "|               3  |               0.791  |\n"
    "|               4  |               0.243  |\n"
    "|               5  |               0.804  |\n"
    "|               6  |               0.302  |\n"
    "-------------------------------------------\n"
  };

  int const rankId = MpiWrapper::commRank();
  int const nbRanks = MpiWrapper::commSize();
  if( nbRanks > 1 )
  {
    ASSERT_EQ( nbRanks, 4 );

    TableLayout const layout = TableLayout().
                                 setTitle( "Summary of negative pressure elements" ).
                                 addColumns( { "Global Id", "pressure [Pa]" } ).
                                 setDefaultHeaderAlignment( TableLayout::Alignment::left );
    TableData data;
    auto const & rankTestData = testCase.m_ranksValues[rankId];

    TableMpiLayout mpiLayout;
    mpiLayout.m_separatorBetweenRanks = true;

    if( !rankTestData.empty() )
    {
      mpiLayout.m_rankTitle = GEOS_FMT( "Rank {}, {} values", rankId, rankTestData.size() );
      for( auto const & [id, value] : rankTestData )
      {
        data.addRow( id, value );
      }
    }

    TableTextMpiFormatter formatter = TableTextMpiFormatter( layout, mpiLayout );
    formatter.setSortingFunc( []( std::vector< TableData::CellData > const & row1,
                                  std::vector< TableData::CellData > const & row2 ) {
      return tabledatasorting::positiveNumberStringComp( row1[0].value, row2[0].value );
    } );

    std::ostringstream oss;
    formatter.toStream( oss, data );
    if( rankId == 0 )
    {
      std::cout << "ma boula " ""<< oss.str()<<std::endl;
      EXPECT_STREQ( testCase.m_expectedResult.data(),
                    oss.str().data() );
    }
  }

}

TEST( testMpiTables, testCompPositiveValueTable )
{
  EXPECT_FALSE ( tabledatasorting::positiveNumberStringComp( "123", "45" ));
  EXPECT_TRUE( tabledatasorting::positiveNumberStringComp( "45", "123" ));
  EXPECT_FALSE( tabledatasorting::positiveNumberStringComp( "42", "42" ));
  EXPECT_FALSE( tabledatasorting::positiveNumberStringComp( "0", "0" ));
  EXPECT_FALSE ( tabledatasorting::positiveNumberStringComp( "9", "1" ));
  EXPECT_FALSE ( tabledatasorting::positiveNumberStringComp( "10000000", "9999999" ));
  EXPECT_TRUE( tabledatasorting::positiveNumberStringComp( "9999999", "10000000" ));

  EXPECT_FALSE ( tabledatasorting::positiveNumberStringComp( "10.5", "9.99" ));
  EXPECT_TRUE( tabledatasorting::positiveNumberStringComp( "9.99", "10.5" ));

  EXPECT_FALSE ( tabledatasorting::positiveNumberStringComp( "10", "9.999" ));
  EXPECT_TRUE( tabledatasorting::positiveNumberStringComp( "9.999", "10" ));

  EXPECT_TRUE( tabledatasorting::positiveNumberStringComp( "1.2", "1.9" ));
  EXPECT_FALSE ( tabledatasorting::positiveNumberStringComp( "1.9", "1.2" ));

  EXPECT_FALSE( tabledatasorting::positiveNumberStringComp( "1.5", "1.50" ));
  EXPECT_FALSE( tabledatasorting::positiveNumberStringComp( "1.50", "1.5" ));
  EXPECT_FALSE( tabledatasorting::positiveNumberStringComp( "1.500", "1.5" ));
  EXPECT_FALSE( tabledatasorting::positiveNumberStringComp( "1.5", "1.500" ));

  EXPECT_FALSE( tabledatasorting::positiveNumberStringComp( "1.51", "1.510" ));
  EXPECT_FALSE ( tabledatasorting::positiveNumberStringComp( "1.51", "1.509" ));

  EXPECT_FALSE( tabledatasorting::positiveNumberStringComp( "3.14", "3.14" ));
  EXPECT_FALSE( tabledatasorting::positiveNumberStringComp( "0.001", "0.001" ));
  EXPECT_FALSE( tabledatasorting::positiveNumberStringComp( "100.0", "100.0" ));

  EXPECT_FALSE( tabledatasorting::positiveNumberStringComp( "5", "5.0" ));
  EXPECT_FALSE( tabledatasorting::positiveNumberStringComp( "5.0", "5" ));

  EXPECT_FALSE ( tabledatasorting::positiveNumberStringComp( "5595", "5155" ));
  EXPECT_TRUE( tabledatasorting::positiveNumberStringComp( "5155", "5595" ));


}

int main( int argc, char * * argv )
{
  int r;
  {
    MpiTestScope testScope{ argc, argv };
    r = RUN_ALL_TESTS();
  }
  return r;
}
