/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "gtest/gtest.h"

#include "tests/wellFileNames.hpp"

#include "fileIO/las/LASFile.hpp"

#include "SetSignalHandling.hpp"

#include "managers/initialization.hpp"

#include "stackTrace.hpp"

using namespace geosx;

/*!
 * @brief This test load a LAS file, save it, load it again and compare the contents
 * with the original one
 */
TEST( LASImport, testXML )
{
  // Load and save the LAS file
  LASFile lasFile;
  lasFile.Load(lasFilePath);
  lasFile.Save("save.las");

  // Reload it
  LASFile lasFile_saved;
  lasFile_saved.Load( "save.las" );

  // Check if the information sections are the same
  integer informationSectionIndex = 0;
  lasFile_saved.forInformationSections([&]( auto & informationSection )
  {
    informationSection->forLines([&]( auto & line )
    {
      EXPECT_EQ( line.GetUnit(),
          lasFile.GetInformationSection( informationSectionIndex).GetLine( line.GetKeyword() ).GetUnit() );
      EXPECT_EQ( line.GetDescription(),
          lasFile.GetInformationSection( informationSectionIndex).GetLine( line.GetKeyword() ).GetDescription() );
      EXPECT_EQ( line.GetData(),
          lasFile.GetInformationSection( informationSectionIndex).GetLine( line.GetKeyword() ).GetData() );
    });
    informationSectionIndex++;
  });

  // Check if the log sections are the same
  integer logSectionIndex = 0;
  real64 tolerance = 1e-5;
  lasFile_saved.forLogSections([&]( auto & logSection )
  {
    for( integer logIndex = 0; logIndex < logSection.NbLogs() ; logIndex++ )
    {
      auto log1 = logSection.GetLog( logIndex );
      auto log2 = lasFile.GetLogSection(logSectionIndex).GetLog( logIndex );
      for( integer entryIndex = 0; entryIndex < logSection.LogSize(); entryIndex++ )
      {
        EXPECT_LE( ( log1[entryIndex] - log2[entryIndex] ) / log2[entryIndex], tolerance);
      }
    }
  });

}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  geosx::basicSetup( argc, argv);

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
