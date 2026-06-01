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

#include "common/DataTypes.hpp"
#include "common/format/LogPart.hpp"
#include "common/initializeEnvironment.hpp"
#include <gtest/gtest.h>

using namespace geos;

TEST( testLogPart, sectionWithTitle )
{
  std::ostringstream oss;
  LogPart logPart( "section name", true );
  logPart.begin( oss );
  EXPECT_EQ( oss.str(),
             "\n"
             "####################################################################################################\n"
             "##                                          section name                                          ##\n"
             "####################################################################################################\n\n"
             );
  oss.clear();
  oss.str( "" );

  logPart.end( oss );
  EXPECT_EQ( oss.str(),
             "\n"
             "##                                      End of section name                                       ##\n"
             "####################################################################################################\n\n"
             );
  oss.clear();
}

TEST( testLogPart, sectionWithTitleAndOneDescription )
{
  std::ostringstream oss;
  LogPart logPart( "section name", true );
  logPart.addDescription( "description name" );
  logPart.begin( oss );
  EXPECT_EQ( oss.str(),
             "\n"
             "####################################################################################################\n"
             "##                                          section name                                          ##\n"
             "####################################################################################################\n"
             "##  description name                                                                              ##\n\n"
             );
  oss.clear();
}

TEST( testLogPart, sectionWithSetWidth )
{
  std::ostringstream oss;
  LogPart logPart( "section name", true );
  logPart.addDescription( "description name 1" );
  logPart.addDescription( "description name 2" );
  logPart.setMinWidth( 100 );
  logPart.begin( oss );
  EXPECT_EQ( oss.str(),
             "\n"
             "####################################################################################################\n"
             "##                                          section name                                          ##\n"
             "####################################################################################################\n"
             "##  description name 1                                                                            ##\n"
             "##  description name 2                                                                            ##\n\n"
             );
  oss.clear();
  oss.str( "" );

  logPart.end( oss );
  std::cout <<oss.str() << std::endl;
  EXPECT_EQ( oss.str(),
             "\n"
             "##                                      End of section name                                       ##\n"
             "####################################################################################################\n\n"
             );
  oss.clear();
}

TEST( testLogPart, sectionMultipleDescriptions )
{
  std::ostringstream oss;
  LogPart logPart( "TIMESTEP START", true );
  logPart.addDescription( "- Time", "00h08m20s out of 2d, 21h26m40s (0% completed)", "500 s / 250000 s" );
  logPart.addDescription( "- Delta Time", "00h16m40s (1000 s)" );
  logPart.addDescription( "Description test" );
  logPart.setMinWidth( 70 );
  logPart.begin( oss );
  EXPECT_EQ ( oss.str(),
              "\n"
              "####################################################################################################\n"
              "##                                         TIMESTEP START                                         ##\n"
              "####################################################################################################\n"
              "##  - Time       : 00h08m20s out of 2d, 21h26m40s (0% completed)                                  ##\n"
              "##                 500 s / 250000 s                                                               ##\n"
              "##  - Delta Time : 00h16m40s (1000 s)                                                             ##\n"
              "##  Description test                                                                              ##\n\n"
              );
  oss.clear();
  oss.str( "" );

  logPart.end( oss );
  std::cout <<oss.str() << std::endl;

  EXPECT_EQ( oss.str(),
             "\n"
             "##                                     End of TIMESTEP START                                      ##\n"
             "####################################################################################################\n\n"
             );
  oss.clear();
}

TEST( testLogPart, sectionEndDescription )
{
  std::ostringstream oss;
  LogPart logPart( "TIMESTEP START", true );
  logPart.addEndDescription( "test end description" );
  logPart.setMinWidth( 70 );
  logPart.begin( oss );
  oss.clear();
  oss.str( "" );

  logPart.end( oss );

  EXPECT_EQ( oss.str(),
             "\n"
             "##  test end description                                                                          ##\n"
             "####################################################################################################\n"
             "##                                     End of TIMESTEP START                                      ##\n"
             "####################################################################################################\n\n"
             );
  oss.clear();
}

TEST( testLogPart, valuesMultiLines )
{

  std::ostringstream oss;
  LogPart logPart( "TIMESTEP START", true );
  logPart.addDescription( "dummy name\ndummy name long", "long dummy values, long dummy values1, long dummy values2, long dummy values3" );
  logPart.addDescription( "dummy long very long for test crash", "long dummy values", "long dummy values", "long dummy values", "long dummy values" );
  logPart.addDescription( "long very long\nwith many \nlien return \nlong dummy name", "small dummy value" );
  logPart.addDescription( "another test\nwith second part of name very long", "small dummy value" );
  logPart.addDescription( "dummy name, long dummy values, long dummy values, long dummy values, long dummy values" );

  logPart.addEndDescription( "dummy name", "long dummy end values, long dummy end values, long dummy end values, long dummy end values" );
  logPart.addEndDescription( "dummy name", "long dummy end values", "long dummy end values", "long dummy end values", "long dummy end values" );
  logPart.addEndDescription( "dummy name", "small dummy end value" );
  logPart.addEndDescription( "Ceci est un timestep extremement long 10h00h00545ms ( 1h 30 s en heure)" );
  logPart.setMaxWidth( 60 );
  logPart.begin( oss );
  EXPECT_EQ( oss.str(),
             "\n"
             "############################################################\n"
             "##                     TIMESTEP START                     ##\n"
             "############################################################\n"
             "##  dummy name                          : long dummy      ##\n"
             "##  dummy name long                     : values, long    ##\n"
             "##                                        dummy values1,  ##\n"
             "##                                        long dummy      ##\n"
             "##                                        values2, long   ##\n"
             "##                                        dummy values3   ##\n"
             "##  dummy long very long for test crash : long dummy      ##\n"
             "##                                        values          ##\n"
             "##                                        long dummy      ##\n"
             "##                                        values          ##\n"
             "##                                        long dummy      ##\n"
             "##                                        values          ##\n"
             "##                                        long dummy      ##\n"
             "##                                        values          ##\n"
             "##  long very long                      : small dummy     ##\n"
             "##  with many                           : value           ##\n"
             "##  lien return                                           ##\n"
             "##  long dummy name                                       ##\n"
             "##  another test                        : small dummy     ##\n"
             "##  with second part of name very long  : value           ##\n"
             "##  dummy name, long dummy values, long dummy values,     ##\n"
             "##  long dummy values, long dummy values                  ##\n\n"
             );

  oss.clear();
  oss.str( "" );

  logPart.end( oss );
  EXPECT_EQ( oss.str(),
             "\n"
             "##  dummy name : long dummy end values, long dummy end    ##\n"
             "##               values, long dummy end values, long      ##\n"
             "##               dummy end values                         ##\n"
             "##  dummy name : long dummy end values                    ##\n"
             "##               long dummy end values                    ##\n"
             "##               long dummy end values                    ##\n"
             "##               long dummy end values                    ##\n"
             "##  dummy name : small dummy end value                    ##\n"
             "##  Ceci est un timestep extremement long 10h00h00545ms   ##\n"
             "##  ( 1h 30 s en heure)                                   ##\n"
             "############################################################\n"
             "##                 End of TIMESTEP START                  ##\n"
             "############################################################\n\n"
             );
  oss.clear();
}

TEST( testLogPart, multiLineWithExtraSpace )
{
  std::ostringstream oss;
  LogPart logPart( "TIMESTEP", true );
  logPart.addDescription( "- Time", "00h00m00s out of 2y, 269d, 12h21m36s (0% completed), 0 s / 86400000 s" );
  logPart.addDescription( "- Delta Time", "00h00m00s (0.001 s)" );
  logPart.addDescription( "- Cycle", "0" );
  logPart.setMaxWidth( 60 );
  logPart.begin( oss );

  EXPECT_EQ( oss.str(),
             "\n"
             "############################################################\n"
             "##                        TIMESTEP                        ##\n"
             "############################################################\n"
             "##  - Time       : 00h00m00s out of 2y, 269d, 12h21m36s   ##\n"
             "##                 (0% completed), 0 s / 86400000 s       ##\n"
             "##  - Delta Time : 00h00m00s (0.001 s)                    ##\n"
             "##  - Cycle      : 0                                      ##\n\n" );
  oss.clear();
  oss.str( "" );
}

int main( int argc, char * * argv )
{
  testing::InitGoogleTest( &argc, argv );

  int const result = RUN_ALL_TESTS();

  return result;
}
