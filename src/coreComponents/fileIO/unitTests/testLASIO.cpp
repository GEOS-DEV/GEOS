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

#include "stackTrace.hpp"

using namespace geosx;
using namespace geosx::dataRepository;

/*!
 * @brief This test load a LAS file, save it, load it again and compare the contents
 * with the original one
 */
TEST( PAMELAImport, testXML )
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
      GEOS_ERROR_IF( line.GetUnit() != lasFile.GetInformationSection( informationSectionIndex).GetLine( line.GetKeyword() ).GetUnit(),
                     "Mismatch between the unit of section " << informationSection->GetName() << ", keyword "<< line.GetKeyword() );
      GEOS_ERROR_IF( line.GetDescription() != lasFile.GetInformationSection( informationSectionIndex).GetLine( line.GetKeyword() ).GetDescription(),
                     "Mismatch between the description of section " << informationSection->GetName() << ", keyword "<< line.GetKeyword() );
      GEOS_ERROR_IF( line.GetData() != lasFile.GetInformationSection( informationSectionIndex).GetLine( line.GetKeyword() ).GetData(),
                     "Mismatch between the data of section " << informationSection->GetName() << ", keyword "<< line.GetKeyword() );
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
        GEOS_ERROR_IF( ( log1[entryIndex] - log2[entryIndex] ) / log2[entryIndex] > tolerance,
                       "Mismatch between logs");
      }
    }
  });

}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

#ifdef GEOSX_USE_MPI

  MPI_Init(&argc,&argv);

  MPI_Comm_dup( MPI_COMM_WORLD, &MPI_COMM_GEOSX );

  logger::InitializeLogger(MPI_COMM_GEOSX);
#else
  logger::InitializeLogger():
#endif

  cxx_utilities::setSignalHandling(cxx_utilities::handler1);

  int const result = RUN_ALL_TESTS();

  logger::FinalizeLogger();

#ifdef GEOSX_USE_MPI
  MPI_Comm_free( &MPI_COMM_GEOSX );
  MPI_Finalize();
#endif

  return result;
}
