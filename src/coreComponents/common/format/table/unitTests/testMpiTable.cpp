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
  stdVector< stdVector< std::pair< integer, real64 > > > const testValues = {
    { {1, 0.502} },
    { {2, 0.624}, {3, 0.791} },
    {},
    { {4, 0.243}, {5, 0.804}, {6, 0.302} },
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
    auto const & rankTestData = testValues[rankId];

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

    TableTextMpiOutput const formatter = TableTextMpiOutput( layout, mpiLayout );
    std::ostringstream oss;
    formatter.toStream( oss, data );
    if( rankId == 0 )
    {
      EXPECT_STREQ( "\n"
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
                    "-------------------------------------------\n",
                    oss.str().data() );
    }
  }
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
