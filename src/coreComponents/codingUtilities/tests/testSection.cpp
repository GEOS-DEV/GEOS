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

TEST( testSection, sectionClass )
{

  std::ostringstream oss;

  //testing section format with only title
  {
    Section section;
    section.setName( "section name" );
    section.begin( oss );
    EXPECT_EQ( oss.str(),
               "\n##############################\n"
               "##  Section : section name  ##\n"
               "##############################\n\n"
               );
    oss.clear();
    oss.str( "" );
    section.end( oss );
    EXPECT_EQ( oss.str(),
               "\n##    End : section name    ##\n"
               "##############################\n\n"
               );
    oss.clear();
    oss.str( "" );
  }

  //testing section format with  title and one description
  {
    Section section;
    section.setName( "section name" );
    section.addDescription( "description name" );
    section.begin( oss );
    EXPECT_EQ( oss.str(),
               "\n##############################\n"
               "##  Section : section name  ##\n"
               "##############################\n"
               "##  description name        ##\n\n"
               );
    oss.clear();
    oss.str( "" );
  }

  //testing section format with title and multiple description with min width
  {
    Section section;
    section.setName( "section name" );
    section.addDescription( "description name 1" );
    section.addDescription( "description name 2" );
    section.setMinWidth( 100 );
    section.begin( oss );
    EXPECT_EQ( oss.str(),
               "\n####################################################################################################\n"
               "##                                     Section : section name                                     ##\n"
               "####################################################################################################\n"
               "##  description name 1                                                                            ##\n"
               "##  description name 2                                                                            ##\n\n"
               );
    oss.clear();
    oss.str( "" );
    section.end( oss );
    EXPECT_EQ( oss.str(),
               "\n##                                       End : section name                                       ##\n"
               "####################################################################################################\n\n"
               );
    oss.clear();
    oss.str( "" );
  }

}

int main( int argc, char * * argv )
{
  testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();;
}
