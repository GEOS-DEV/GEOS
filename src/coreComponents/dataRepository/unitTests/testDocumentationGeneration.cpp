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

//include LogLevel
#include "dataRepository/Group.hpp"
#include "dataRepository/LogLevelsRegistry.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <gtest/gtest-spi.h>


using namespace geos;
using namespace dataRepository;


TEST( testDocumentationGeneration, testLogLevelDocumentation )
{
  // Keep track of log levels & descriptions
  std::unique_ptr< LogLevelsRegistry > logLevelsRegistry = std::make_unique< LogLevelsRegistry >();

  {
    logLevelsRegistry->addEntry( 0, "description1" );
    logLevelsRegistry->addEntry( 1, "description1" );
    logLevelsRegistry->addEntry( 2, "description1" );
    logLevelsRegistry->addEntry( 3, "description1" );

    string const descriptionBuild = logLevelsRegistry->buildLogLevelDescription();
    EXPECT_EQ( descriptionBuild,
               "Sets the level of information to write in the standard output (the console typically).\n"
               "Information output from lower logLevels is added with the desired log level\n"
               "0\n"
               " - description1\n"
               "1\n"
               " - description1\n"
               "2\n"
               " - description1\n"
               "3\n"
               " - description1"
               );
  }

  { // duplication
    logLevelsRegistry->addEntry( 0, "description1" );
    logLevelsRegistry->addEntry( 1, "description1" );
    logLevelsRegistry->addEntry( 2, "description1" );
    logLevelsRegistry->addEntry( 3, "description1" );
    string const duplicationBuild = logLevelsRegistry->buildLogLevelDescription();

    EXPECT_EQ( duplicationBuild,
               "Sets the level of information to write in the standard output (the console typically).\n"
               "Information output from lower logLevels is added with the desired log level\n"
               "0\n"
               " - description1\n"
               "1\n"
               " - description1\n"
               "2\n"
               " - description1\n"
               "3\n"
               " - description1"
               );
  }

  {
    logLevelsRegistry->addEntry( 0, "description2" );
    logLevelsRegistry->addEntry( 0, "description1" );
    logLevelsRegistry->addEntry( 0, "description3" );
    logLevelsRegistry->addEntry( 0, "description4" );
    string const duplicationBuild = logLevelsRegistry->buildLogLevelDescription();

    EXPECT_EQ( duplicationBuild,
               "Sets the level of information to write in the standard output (the console typically).\n"
               "Information output from lower logLevels is added with the desired log level\n"
               "0\n"
               " - description1\n"
               " - description2\n"
               " - description3\n"
               " - description4\n"
               "1\n"
               " - description1\n"
               "2\n"
               " - description1\n"
               "3\n"
               " - description1"
               );
  }

}

int main( int argc, char * * argv )
{
  testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
