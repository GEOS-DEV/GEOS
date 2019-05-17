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

#include <numeric>

#include "SetSignalHandling.hpp"
#include "stackTrace.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "managers/ProblemManager.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "createConnLocPattern.hpp"

#include "DofManager.hpp"

using namespace geosx;

namespace
{
int global_argc;
char** global_argv;
}

class DofManagerTest : public ::testing::Test
{
public:
  void setTimer( double &time )
  {
    time = MPI_Wtime();
  }
  void getElapsedTime( double &time )
  {
    time = MPI_Wtime() - time;
  }

protected:

  static void SetUpTestCase()
  {
    problemManager = new ProblemManager("Problem", nullptr);

    problemManager->InitializePythonInterpreter();
    problemManager->ParseCommandLineInput( global_argc, global_argv );
    problemManager->ParseInputFile();
    problemManager->ProblemSetup();
  }

  static void TearDownTestCase()
  {
    delete problemManager;
    problemManager = nullptr;
  }

  static ProblemManager * problemManager;
};

ProblemManager * DofManagerTest::problemManager = nullptr;

TEST_F(DofManagerTest, TestOne)
{
  DomainPartition * const domain = problemManager->getDomainPartition();

  DofManager dofManager;
  dofManager.setMesh( domain, 0, 0 );

  const bool printPattern = false;

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

  // ParallelMatrix userPattern;
  // createLocLocPattern( domain, 0, 0, userPattern );

  ParallelMatrix connLocInput;
  createConnLocPattern( domain, 0, 0, 1, connLocInput );

  double timeAddField, timeAddCoupling, timeGetSingleSparsityPattern, timeGetGlobalSparsityPattern, timeGetIndices;

  setTimer( timeAddField );

  dofManager.addField( "displacement", DofManager::Location::Node, DofManager::Connectivity::Elem, 1,
                       displacementRegion );
  //dofManager.addField( "displacement", DofManager::Location::Node, DofManager::Connectivity::Elem, 1);
  dofManager.addField( "pressure", DofManager::Location::Elem, DofManager::Connectivity::Face, pressureRegion );
  dofManager.addField( "massmatrix", DofManager::Location::Elem, DofManager::Connectivity::None, 2, testRegion3 );
  dofManager.addField( "user-defined", connLocInput, DofManager::Connectivity::Face );

  /*
   dofManager.addField( "facenode", DofManager::Location::Node, DofManager::Connectivity::Face, testRegion2 );
   dofManager.addField( "elemface", DofManager::Location::Face, DofManager::Connectivity::Elem );
   dofManager.addField( "nodeelem", DofManager::Location::Elem, DofManager::Connectivity::Node, pressureRegion );
   dofManager.addField( "nodeface", DofManager::Location::Face, DofManager::Connectivity::Node, testRegion2 );
   */

  getElapsedTime( timeAddField );
  setTimer( timeAddCoupling );

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

  getElapsedTime( timeAddCoupling );
  setTimer( timeGetSingleSparsityPattern );

  ParallelMatrix pattern;

  dofManager.getSparsityPattern( pattern, "displacement", "displacement" );
  if( printPattern )
    dofManager.printParallelMatrix( pattern, "displacement" );

  dofManager.getSparsityPattern( pattern, "pressure", "pressure" );
  if( printPattern )
    dofManager.printParallelMatrix( pattern, "pressure" );

  dofManager.getSparsityPattern( pattern, "massmatrix", "massmatrix" );
  if( printPattern )
    dofManager.printParallelMatrix( pattern, "massmatrix.mtx" );

  dofManager.getSparsityPattern( pattern, "displacement", "pressure" );
  if( printPattern )
    dofManager.printParallelMatrix( pattern, "coupling1.mtx" );

  dofManager.getSparsityPattern( pattern, "pressure", "displacement" );
  if( printPattern )
    dofManager.printParallelMatrix( pattern, "coupling1_empty.mtx" );

  dofManager.getSparsityPattern( pattern, "pressure", "massmatrix" );
  if( printPattern )
    dofManager.printParallelMatrix( pattern, "coupling2.mtx" );

  dofManager.getSparsityPattern( pattern, "massmatrix", "pressure" );
  if( printPattern )
    dofManager.printParallelMatrix( pattern, "coupling2_transp.mtx" );

  dofManager.getSparsityPattern( pattern, "user-defined", "user-defined" );
  if( printPattern )
    dofManager.printParallelMatrix( pattern, "user-defined.mtx" );

  getElapsedTime( timeGetSingleSparsityPattern );
  setTimer( timeGetGlobalSparsityPattern );

  dofManager.getSparsityPattern( pattern );
  if( printPattern )
    dofManager.printParallelMatrix( pattern, "global" );

  // Create the permutation collecting all unknowns belonging to each process
  ParallelMatrix permutation;
  dofManager.createPermutation( permutation );

  ParallelMatrix permutedPattern;
  // Apply the permutation
  dofManager.permuteSparsityPattern( pattern, permutation, permutedPattern );
  if( printPattern )
    dofManager.printParallelMatrix( permutedPattern, "permutatedGlobal" );

  getElapsedTime( timeGetGlobalSparsityPattern );
  setTimer( timeGetIndices );

  int mpiRank;
  MPI_Comm_rank( MPI_COMM_GEOSX, &mpiRank );

  globalIndex_array indices;
  dofManager.getIndices( indices, DofManager::Connectivity::Elem, 1, 0, 10, "displacement" );
  if( indices.size() > 0 )
    std::cout << mpiRank << " - " << "displacement" << " " << indices << std::endl;

  dofManager.getIndices( indices, DofManager::Connectivity::Face, 30, "pressure" );
  if( indices.size() > 0 )
    std::cout << mpiRank << " - " << "pressure" << " " << indices << std::endl;

  getElapsedTime( timeGetIndices );

  std::cout << mpiRank << " - numGlobalDofs - " << dofManager.numGlobalDofs() << std::endl;
  std::cout << mpiRank << " - numGlobalDofs(" "displacement" ") - " << dofManager.numGlobalDofs( "displacement" )
            << std::endl;
  std::cout << mpiRank << " - numGlobalDofs(" "pressure" ") - " << dofManager.numGlobalDofs( "pressure" ) << std::endl;
  std::cout << mpiRank << " - numLocalDofs - " << dofManager.numLocalDofs() << std::endl;
  std::cout << mpiRank << " - numLocalDofs(" "displacement" ") - " << dofManager.numLocalDofs( "displacement" )
            << std::endl;
  std::cout << mpiRank << " - numLocalDofs(" "pressure" ") - " << dofManager.numLocalDofs( "pressure" ) << std::endl;

  dofManager.printConnectivityMatrix();

  // Sum up all timings
  int mpiSize;
  MPI_Comm_size( MPI_COMM_GEOSX, &mpiSize );
  array1d<double> timesLocal( 5 ), timesSum( 5 );
  timesLocal[0] = timeAddField;
  timesLocal[1] = timeAddCoupling;
  timesLocal[2] = timeGetSingleSparsityPattern;
  timesLocal[3] = timeGetGlobalSparsityPattern;
  timesLocal[4] = timeGetIndices;

  MPI_Allreduce( timesLocal.data(), timesSum.data(), 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_GEOSX );

  timeAddField = timesSum[0];
  timeAddCoupling = timesSum[1];
  timeGetSingleSparsityPattern = timesSum[2];
  timeGetGlobalSparsityPattern = timesSum[3];
  timeGetIndices = timesSum[4];

  if( mpiRank == 0 )
  {
    double totalTime = std::accumulate( timesSum.begin(), timesSum.end(), 0.0 );
    std::cout << "TIMING" << std::endl;
    std::cout << "addField: " << std::fixed << std::setprecision( 0 ) << std::trunc( timeAddField * 1e3 ) << " [ms] -- "
              << std::fixed
              << std::setprecision( 2 ) << timeAddField / totalTime * 100.0 << "%" << std::endl;
    std::cout << "addCoupling: " << std::fixed << std::setprecision( 0 ) << std::trunc( timeAddCoupling * 1e3 )
              << " [ms] -- "
              << std::fixed << std::setprecision( 2 ) << timeAddCoupling / totalTime * 100.0 << "%"
              << std::endl;
    std::cout << "getSparsityPattern(local): " << std::fixed << std::setprecision( 0 )
              << std::trunc( timeGetSingleSparsityPattern * 1e3 )
              << " [ms] -- " << std::fixed << std::setprecision( 2 )
              << timeGetSingleSparsityPattern / totalTime * 100.0
              << "%" << std::endl;
    std::cout << "getSparsityPattern(global): " << std::fixed << std::setprecision( 0 )
              << std::trunc( timeGetGlobalSparsityPattern * 1e3 )
              << " [ms] -- " << std::fixed << std::setprecision( 2 )
              << timeGetGlobalSparsityPattern / totalTime * 100.0
              << "%" << std::endl;
    std::cout << "getIndices: " << std::fixed << std::setprecision( 0 ) << std::trunc( timeGetIndices * 1e3 )
              << " [ms] -- "
              << std::fixed << std::setprecision( 2 ) << timeGetIndices / totalTime * 100.0 << "%"
              << std::endl;
    std::cout << "TOTAL: " << std::fixed << std::setprecision( 0 ) << std::trunc( totalTime * 1e3 ) << " [ms] -- "
              << std::fixed
              << std::setprecision( 2 ) << 100.0 << "%" << std::endl;
  }
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
