/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wused-but-marked-unused"
#endif

#include "gtest/gtest.h"

#ifdef __clang__
#define __null nullptr
#endif

#include "SetSignalHandling.hpp"
#include "stackTrace.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "managers/ProblemManager.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"

#include "DofManager.hpp"

using namespace geosx;

namespace
{
int global_argc;
char** global_argv;
}

class DofManagerTest : public ::testing::Test
{
protected:

  static void SetUpTestCase()
  {
    problemManager.InitializePythonInterpreter();
    problemManager.ParseCommandLineInput( global_argc, global_argv );
    problemManager.ParseInputFile();
    problemManager.ProblemSetup();
  }

  static void TearDownTestCase()
  {
  }

  static ProblemManager problemManager;
};

ProblemManager DofManagerTest::problemManager( "ProblemManager", nullptr );

TEST_F(DofManagerTest, TestOne)
{
  DomainPartition * const domain = problemManager.getDomainPartition();

  DofManager dofManager;
  dofManager.setMesh( domain, 0, 0 );

  string_array displacementRegion;
  displacementRegion.push_back( "region1" );
  displacementRegion.push_back( "region3" );
  displacementRegion.push_back( "region4" );

  string_array pressureRegion;
  pressureRegion.push_back( "region1" );
  pressureRegion.push_back( "region2" );
  pressureRegion.push_back( "region3" );

  string_array testRegion1;
  testRegion1.push_back( "region1" );

  string_array testRegion2;
  testRegion2.push_back( "region2" );
  testRegion2.push_back( "region4" );

  string_array testRegion3;
  testRegion3.push_back( "region3" );

  dofManager.addField( "displacement", DofManager::Location::Node, DofManager::Connectivity::Elem, 1,
                       displacementRegion );
  //dofManager.addField( "displacement", DofManager::Location::Node, DofManager::Connectivity::Elem, 1);
  dofManager.addField( "pressure", DofManager::Location::Elem, DofManager::Connectivity::Face, pressureRegion );
  dofManager.addField( "massmatrix", DofManager::Location::Elem, DofManager::Connectivity::None, 2, testRegion3 );

  /*
   dofManager.addField( "facenode", DofManager::Location::Node, DofManager::Connectivity::Face, testRegion2 );
   dofManager.addField( "elemface", DofManager::Location::Face, DofManager::Connectivity::Elem );
   dofManager.addField( "nodeelem", DofManager::Location::Elem, DofManager::Connectivity::Node, pressureRegion );
   dofManager.addField( "nodeface", DofManager::Location::Face, DofManager::Connectivity::Node, testRegion2 );
   */

  dofManager.addCoupling( "displacement", "pressure", DofManager::Connectivity::Elem, testRegion3, false );
  dofManager.addCoupling( "massmatrix", "pressure", DofManager::Connectivity::Elem );
  /*
   dofManager.addCoupling( "facenode", "pressure", DofManager::Connectivity::Node, testRegion1 );
   dofManager.addCoupling( "elemface", "nodeelem", DofManager::Connectivity::Node, testRegion3, false );
   dofManager.addCoupling( "nodeface", "nodeelem", DofManager::Connectivity::Face, testRegion1 );
   dofManager.addCoupling( "displacement", "facenode", DofManager::Connectivity::Face );
   dofManager.addCoupling( "pressure", "nodeface", DofManager::Connectivity::Face, testRegion1 );
   dofManager.addCoupling( "nodeface", "displacement", DofManager::Connectivity::Elem );
   */

  ParallelMatrix pattern;

  dofManager.getSparsityPattern( pattern, "displacement", "displacement" );
  dofManager.printParallelMatrix( pattern, "displacement" );

  dofManager.getSparsityPattern( pattern, "pressure", "pressure" );
  dofManager.printParallelMatrix( pattern, "pressure" );

  dofManager.getSparsityPattern( pattern, "displacement", "pressure" );
  dofManager.printParallelMatrix( pattern, "coupling.mtx" );

  dofManager.getSparsityPattern( pattern, "pressure", "displacement" );
  dofManager.printParallelMatrix( pattern, "coupling_empty.mtx" );

  dofManager.getSparsityPattern( pattern, "massmatrix", "massmatrix" );
  dofManager.printParallelMatrix( pattern, "massmatrix.mtx" );

  dofManager.getSparsityPattern( pattern );
  dofManager.printParallelMatrix( pattern, "global" );

  int mpiRank;
  MPI_Comm_rank( MPI_COMM_GEOSX, &mpiRank );

  globalIndex_array indices;
  dofManager.getIndices( indices, DofManager::Connectivity::Elem, 1, 0, 10, "displacement" );
  if( indices.size() > 0 )
    std::cout << mpiRank << " - " << "displacement" << " " << indices << std::endl;

  dofManager.getIndices( indices, DofManager::Connectivity::Face, 30, "pressure" );
  if( indices.size() > 0 )
    std::cout << mpiRank << " - " << "pressure" << " " << indices << std::endl;

  std::cout << mpiRank << " - numGlobalDofs - " << dofManager.numGlobalDofs() << std::endl;
  std::cout << mpiRank << " - numGlobalDofs(" "displacement" ") - " << dofManager.numGlobalDofs( "displacement" )
            << std::endl;
  std::cout << mpiRank << " - numGlobalDofs(" "pressure" ") - " << dofManager.numGlobalDofs( "pressure" ) << std::endl;
  std::cout << mpiRank << " - numLocalDofs - " << dofManager.numLocalDofs() << std::endl;
  std::cout << mpiRank << " - numLocalDofs(" "displacement" ") - " << dofManager.numLocalDofs( "displacement" )
            << std::endl;
  std::cout << mpiRank << " - numLocalDofs(" "pressure" ") - " << dofManager.numLocalDofs( "pressure" ) << std::endl;

  dofManager.printConnectivityMatrix();

  dofManager.cleanUp();
}

int main( int argc, char** argv )
{
  ::testing::InitGoogleTest( &argc, argv );
#ifdef GEOSX_USE_MPI
  MPI_Init( &argc, &argv );
  MPI_Comm_dup( MPI_COMM_WORLD, &MPI_COMM_GEOSX );
  logger::InitializeLogger( MPI_COMM_GEOSX );
#else
  logger::InitializeLogger():
#endif
  cxx_utilities::setSignalHandling( cxx_utilities::handler1 );

  global_argc = argc;
  global_argv = new char*[static_cast<unsigned int>( global_argc )];
  for( int i = 0 ; i < argc ; ++i )
  {
    global_argv[i] = argv[i];
  }

  int const result = RUN_ALL_TESTS();

  delete[] global_argv;
  logger::FinalizeLogger();
#ifdef GEOSX_USE_MPI
  MPI_Comm_free( &MPI_COMM_GEOSX );
  MPI_Finalize();
#endif

  return result;
}

#ifdef __clang__
#pragma clang diagnostic pop
#endif
