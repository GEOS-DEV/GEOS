/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
#include "gtest/gtest.h"

#include "fileIO/xmlWrapper.hpp"

#include "tests/meshFileNames.hpp"

#include "meshUtilities/MeshManager.hpp"

#include "SetSignalHandling.hpp"

#include "stackTrace.hpp"

using namespace geosx;
using namespace geosx::dataRepository;

TEST( PAMELAImport, testXML )
{
  MeshManager meshManager("mesh", nullptr);

  std::stringstream inputStream;
  inputStream <<
  "<?xml version=\"1.0\" ?>" <<
  "  <Mesh xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"geos_v0.0.xsd\">" <<
  "  <PAMELAMeshGenerator name=\"ToyModel\" " <<
  "  file=\"" <<gmshFilePath.c_str()<< "\"/>"<<
  "</Mesh>";
  const string inputString = inputStream.str();

  xmlWrapper::xmlDocument xmlDocument;
  xmlWrapper::xmlResult xmlResult = xmlDocument.load_buffer( inputString.c_str(), inputString.size() );
  if (!xmlResult)
  {
    GEOS_LOG_RANK_0("XML parsed with errors!");
    GEOS_LOG_RANK_0("Error description: " << xmlResult.description());
    GEOS_LOG_RANK_0("Error offset: " << xmlResult.offset);
  }

  xmlWrapper::xmlNode xmlMeshNode = xmlDocument.child("Mesh");
  meshManager.ProcessInputFileRecursive( xmlMeshNode );
  meshManager.PostProcessInputRecursive();

  auto domain = std::unique_ptr< DomainPartition >( new DomainPartition( "domain", nullptr ) );
  meshManager.GenerateMeshes( domain.get() );

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
