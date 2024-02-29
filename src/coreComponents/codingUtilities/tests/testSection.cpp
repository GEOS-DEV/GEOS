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

#include "../../dataRepository/Group.hpp"
// TPL includes
#include "codingUtilities/Section.hpp"
#include <gtest/gtest.h>

using namespace geos;

TEST( sectionTable, sectionClass )
{
  std::vector< string > testSectionOutput;

  std::ostringstream oss;

  Section section1;
  section1.setName( "section name" );
  section1.begin( oss );
  testSectionOutput.push_back( oss.str() );
  oss.clear();
  oss.str( "" );
  EXPECT_EQ( testSectionOutput[0], "" );

  section1.end( oss );
  testSectionOutput.push_back( oss.str() );
  oss.clear();
  oss.str( "" );
  EXPECT_EQ( testSectionOutput[1], "" );

  Section section2;
  section2.setName( "section name" );
  section2.addDescription( "description name" );
  section2.begin( oss );
  testSectionOutput.push_back( oss.str() );
  oss.clear();
  oss.str( "" );
  EXPECT_EQ( testSectionOutput[2], "" );

  Section section3;
  section3.setName( "section name" );
  section3.addDescription( "description name 1" );
  section3.addDescription( "description name 2" );
  section3.begin( oss );
  testSectionOutput.push_back( oss.str() );
  oss.clear();
  oss.str( "" );
  EXPECT_EQ( testSectionOutput[3], "" );

  Section section4;
  section4.setName( "section name" );
  section4.addDescription( "description name 1" );
  section4.addDescription( "description name 2" );
  section4.setMinWidth( 100 );
  section4.begin( oss );
  testSectionOutput.push_back( oss.str() );
  oss.clear();
  oss.str( "" );
  EXPECT_EQ( testSectionOutput[4], "" );
  section4.end( oss );
  testSectionOutput.push_back( oss.str() );
  oss.clear();
  oss.str( "" );
  EXPECT_EQ( testSectionOutput[5], "" );


}

int main( int argc, char * * argv )
{
  testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();;
}
