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

  string_array oneRegion;
  oneRegion.push_back( "region2" );

  string_array twoRegion;
  twoRegion.push_back( "region1" );
  twoRegion.push_back( "region2" );

  //- dofManager.addField("displacement", DofManager::Location::Node, DofManager::Connectivity::Elem, 1, twoRegion);
  dofManager.addField( "pressure", DofManager::Location::Elem, DofManager::Connectivity::Face, 1, twoRegion );
  //dofManager.addField("face_node", DofManager::Location::Node, DofManager::Connectivity::Face, 1, twoRegion);
  //dofManager.addField("node_face", DofManager::Location::Face, DofManager::Connectivity::Node, 1, twoRegion);
  //dofManager.addField("massmatrix", DofManager::Location::Node, DofManager::Connectivity::None, 4, twoRegion);
  //- dofManager.addCoupling("displacement", "pressure", DofManager::Connectivity::Elem, true);
  //SparsityPattern const & pattern = dofManager.getSparsityPattern("displacement", "pressure");
  //SparsityPattern const & pattern = dofManager.getSparsityPattern("displacement", "displacement");
  //SparsityPattern const & pattern = dofManager.getSparsityPattern("displacement");
  //SparsityPattern const & pattern = dofManager.getSparsityPattern("", "displacement");
  //SparsityPattern const & pattern = dofManager.getSparsityPattern("pressure", "pressure");
  //SparsityPattern const & pattern = dofManager.getSparsityPattern("pressure");
  //SparsityPattern const & pattern = dofManager.getSparsityPattern("", "pressure");
  SparsityPattern const & pattern = dofManager.getSparsityPattern();

  int mpiRank;
  MPI_Comm_rank( MPI_COMM_GEOSX, &mpiRank );
  dofManager.printSparsityPattern( pattern, "global_pattern_" + std::to_string( mpiRank ) + ".csr" );

  globalIndex_array indices;
  //dofManager.getIndices(indices, DofManager::Connectivity::Elem, 1, 0, 10, "displacement");
  dofManager.getIndices( indices, DofManager::Connectivity::Face, 10, "pressure" );
  if( indices.size() > 0 )
    std::cout << indices << std::endl;

  std::cout << mpiRank << " - numGlobalDofs - " << dofManager.numGlobalDofs() << std::endl;
  //std::cout << mpiRank << " - numGlobalDofs(""pressure"") - " << dofManager.numGlobalDofs("pressure") << std::endl;
  std::cout << mpiRank << " - numLocalDofs - " << dofManager.numLocalDofs() << std::endl;

  /*
   dofManager.addField("acceleration", DofManager::Location::Node, DofManager::Connectivity::Elem, 3);
   dofManager.addField("flux", DofManager::Location::Face, DofManager::Connectivity::Elem, 1);
   dofManager.addField("pressure", DofManager::Location::Elem, DofManager::Connectivity::Face, 1, oneRegion);
   dofManager.addField("composition", DofManager::Location::Elem, DofManager::Connectivity::Face, 2, oneRegion);
   dofManager.addCoupling("displacement", "pressure", DofManager::Connectivity::Elem, true);
   */

  dofManager.printCoupling();
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
