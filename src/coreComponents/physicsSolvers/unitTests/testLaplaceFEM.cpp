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

/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
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
#include "managers/EventManager.hpp"
#include "managers/DomainPartition.hpp"
#include "meshUtilities/MeshManager.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/SimpleSolvers/LaplaceFEM.hpp"

using namespace geosx;
using namespace geosx::systemSolverInterface;

namespace
{
int global_argc;
char** global_argv;
int mpiRank = 0;
int mpiSize = 1;
}

class LaplaceFEMTest : public ::testing::Test
{
protected:

  static void SetUpTestCase()
  {
    string const inputStream =
    "<?xml version=\"1.0\" ?>"
    "<Problem xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"geos_v0.0.xsd\">"
    "  <Solvers>"
    "    <LaplaceFEM name=\"laplace\""
    "                discretization=\"FE1\""
    "                timeIntegrationOption=\"SteadyState\""
    "                fieldName=\"Temperature\""
    "                verboseLevel=\"0\""
    "                targetRegions=\"Region1\">"
    "    </LaplaceFEM>"
    "  </Solvers>"
    "  <Mesh>"
    "    <InternalMesh name=\"mesh1\""
    "                  elementTypes=\"C3D8\""
    "                  xCoords=\"0, 1\""
    "                  yCoords=\"0, 1\""
    "                  zCoords=\"0, 1\""
    "                  nx=\"10\""
    "                  ny=\"10\""
    "                  nz=\"10\""
    "                  cellBlockNames=\"cb1\" />"
    "  </Mesh>"
    "  <Events maxTime=\"1.0\">"
    "    <!-- This event is applied every cycle, and overrides the solver time-step request -->"
    "    <PeriodicEvent name=\"solverApplications\""
    "                   forceDt=\"1.0\""
    "                   target=\"/Solvers/laplace\" />"
    "    <!-- This event is applied every 1.0s.  The targetExactTimestep flag allows this event"
    "    to request a dt modification to match an integer multiple of the timeFrequency. -->"
    "    <PeriodicEvent name=\"outputs\""
    "                   timeFrequency=\"1.0\""
    "                   targetExactTimestep=\"1\""
    "                   target=\"/Outputs/siloOutput\" />"
    "  </Events>"
    "  <NumericalMethods>"
    "    <BasisFunctions>"
    "      <LagrangeBasis3 name=\"linearBasis\" degree=\"1\" />"
    "    </BasisFunctions>"
    "    <QuadratureRules>"
    "      <GaussQuadrature3 name=\"gaussian\" degree=\"2\" />"
    "    </QuadratureRules>"
    "    <FiniteElements>"
    "      <FiniteElementSpace name=\"FE1\" basis=\"linearBasis\" quadrature=\"gaussian\" />"
    "    </FiniteElements>"
    "  </NumericalMethods>"
    "  <ElementRegions>"
    "    <ElementRegion name=\"Region1\" cellBlocks=\"cb1\" materialList=\"shale\" />"
    "  </ElementRegions>"
    "  <Constitutive>"
    "    <LinearElasticIsotropic name=\"granite\""
    "                            defaultDensity=\"2700\""
    "                            defaultBulkModulus=\"5.5556e9\""
    "                            defaultShearModulus=\"4.16667e9\" />"
    "    <LinearElasticIsotropic name=\"shale\""
    "                            defaultDensity=\"2700\""
    "                            defaultBulkModulus=\"5.5556e9\""
    "                            defaultShearModulus=\"4.16667e9\" />"
    "  </Constitutive>"
    "  <FieldSpecifications>"
    "    <FieldSpecification name=\"sourceTerm\""
    "                        fieldName=\"Temperature\""
    "                        objectPath=\"nodeManager\""
    "                        scale=\"1000.0\""
    "                        setNames=\"source\" />"
    "    <FieldSpecification name=\"sinkTerm\""
    "                        fieldName=\"Temperature\""
    "                        objectPath=\"nodeManager\""
    "                        scale=\"0.0\""
    "                        setNames=\"sink\" />"
    "  </FieldSpecifications>"
    "  <Functions>"
    "    <TableFunction name=\"timeFunction\""
    "                   inputVarNames=\"time\""
    "                   coordinates=\"0.0 1.0e-6 2.0e-6 1.0e9\""
    "                   values=\"0.0 1.0 1.0 1.0\" />"
    "    <SymbolicFunction name=\"spaceFunction\""
    "                      inputVarNames=\"ReferencePosition\""
    "                      variableNames=\"x y z\""
    "                      expression=\"sqrt(pow(x,2)+pow(y,2)+pow(z,2))\" />"
    "  </Functions>"
    "  <Outputs>"
    "    <Silo name=\"siloOutput\" parallelThreads=\"32\" plotFileRoot=\"plot\" />"
    "  </Outputs>"
    "  <Geometry>"
    "    <Box name=\"source\" xMin=\"-0.01, -0.01, -0.01\" xMax=\"+0.01, +1.01, +1.01\" />"
    "    <Box name=\"sink\"   xMin=\"+0.99, -0.01, -0.01\" xMax=\"+1.01, +1.01, +1.01\" />"
    "  </Geometry>"
    "</Problem>";

    xmlWrapper::xmlDocument xmlDocument;
    xmlWrapper::xmlResult xmlResult = xmlDocument.load_buffer( inputStream.c_str(), inputStream.size() );
    if (!xmlResult)
    {
      GEOS_LOG_RANK_0("XML parsed with errors!");
      GEOS_LOG_RANK_0("Error description: " << xmlResult.description());
      GEOS_LOG_RANK_0("Error offset: " << xmlResult.offset);
    }

    dataRepository::ManagedGroup * commandLine =
      problemManager.GetGroup<dataRepository::ManagedGroup>( problemManager.groupKeys.commandLine );
    commandLine->RegisterViewWrapper<integer>( problemManager.viewKeys.zPartitionsOverride.Key() )->
      setApplyDefaultValue(mpiSize);

    xmlWrapper::xmlNode xmlProblemNode = xmlDocument.child( "Problem" );
    problemManager.InitializePythonInterpreter();
    problemManager.ProcessInputFileRecursive( xmlProblemNode );

    // The objects in domain are handled separately for now
    DomainPartition * domain  = problemManager.getDomainPartition();
    constitutive::ConstitutiveManager *
    constitutiveManager = domain->GetGroup<constitutive::ConstitutiveManager>( problemManager.groupKeys.constitutiveManager );
    xmlWrapper::xmlNode topLevelNode = xmlProblemNode.child(constitutiveManager->getName().c_str());
    constitutiveManager->ProcessInputFileRecursive( topLevelNode );
    constitutiveManager->PostProcessInputRecursive();

    // Open mesh levels
    MeshManager * meshManager = problemManager.GetGroup<MeshManager>( problemManager.groupKeys.meshManager );
    meshManager->GenerateMeshLevels(domain);

    ElementRegionManager * elementManager = domain->getMeshBody(0)->getMeshLevel(0)->getElemManager();
    topLevelNode = xmlProblemNode.child( elementManager->getName().c_str() );
    elementManager->ProcessInputFileRecursive( topLevelNode );
    elementManager->PostProcessInputRecursive();

    problemManager.ProblemSetup();

    /////////////////////////////
    /////////////////////////////
    /*
    problemManager.InitializePythonInterpreter();
    problemManager.ParseCommandLineInput( global_argc, global_argv );
    problemManager.ParseInputFile();
    problemManager.ProblemSetup();
    */
    /////////////////////////////
    /////////////////////////////

    solver = problemManager.GetPhysicsSolverManager().GetGroup<LaplaceFEM>( "laplace" );
  }

  static void TearDownTestCase()
  {

  }

  static ProblemManager problemManager;
  static LaplaceFEM * solver;
};

ProblemManager LaplaceFEMTest::problemManager( "Problem", nullptr );
LaplaceFEM * LaplaceFEMTest::solver = nullptr;

TEST_F(LaplaceFEMTest, laplaceSolverCheckSolution)
{
  real64 const eps = sqrt(std::numeric_limits<real64>::epsilon());

  string const fieldName = "Temperature";
  real64 const time = 1.0;
  real64 const dt = 1.0;
  real64 const scalingFactor = 1.0;
  int const cycleNumber = 0;

  DomainPartition * domain = problemManager.getDomainPartition();

  // Create and solve the problem
  LaplaceFEM laplaceFEM( fieldName, domain );
  EpetraBlockSystem * system = solver->getLinearSystemRepository();
  solver->ImplicitStepSetup( time, dt, domain, system );
  solver->AssembleSystem( domain, system, time, dt );
  solver->ApplyBoundaryConditions( domain, system, time, dt );
  solver->SolverStep( time, dt, cycleNumber, domain );
  solver->ApplySystemSolution( system, scalingFactor, domain );

  // Get matrix and matrix size
  Epetra_FECrsMatrix const * const matrix = system->GetMatrix( BlockIDs::dummyScalarBlock,
                                                               BlockIDs::dummyScalarBlock );
  real64 const matrixSize3 = std::pow( static_cast<real64>( matrix->NumGlobalRows64() ), 1.0/3.0 );
  real64 const tol = 4.0 * std::pow( matrixSize3, 2 ) * eps;

  // Get solution
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();
  real64_array & fieldVar = nodeManager->getReference<real64_array>( fieldName );

  // Compute relative error
  real64 xMin = 0.0;
  real64 xMax = 0.0;
  // TODO read them on input
  real64 const vMin = 1000.0;
  real64 const vMax = 0.0;

  // Compute domain bounds (x direction)
  r1_array const & referencePosition = nodeManager->getReference<r1_array>(dataRepository::keys::referencePositionString);
  localIndex const numNodes = nodeManager->size();
  for( localIndex a = 0 ; a < numNodes ; ++a )
  {
    R1Tensor nodePosition;
    nodePosition = referencePosition[a];
    if( xMin > nodePosition[0] )
    {
      xMin = nodePosition[0];
    }
    if( xMax < nodePosition[0] )
    {
      xMax = nodePosition[0];
    }
  }

  // Compute xMax and xMin across ranks
  real64_array gather;
  CommunicationTools::allGather( xMax, gather );
  xMax = *std::max_element( gather.begin(), gather.end() );
  CommunicationTools::allGather( xMin, gather );
  xMin = *std::min_element( gather.begin(), gather.end() );

  // Compute true solution and error
  real64 const slope = ( vMax - vMin ) / ( xMax - xMin );
  real64 error = 0.0;
  real64 normSol = 0.0;
  for( localIndex a = 0 ; a < numNodes ; ++a )
  {
    R1Tensor nodePosition;
    nodePosition = referencePosition[a];
    real64 refVal = slope * ( nodePosition[0] - xMin ) + vMin;
    error += std::pow( fieldVar[a] - refVal, 2 );
    normSol += std::pow( refVal, 2 );
  }

  // Gather errors across ranks
  CommunicationTools::allGather( error, gather );
  error = 0.0;
  for( localIndex p = 0 ; p < mpiSize ; ++p )
  {
    error += gather[p];
  }
  error = std::sqrt( error );

  // Gather solution norms across ranks
  CommunicationTools::allGather( normSol, gather );
  normSol = 0.0;
  for( localIndex p = 0 ; p < mpiSize ; ++p )
  {
    normSol += gather[p];
  }
  normSol = std::sqrt( normSol );

  // Compute and check relative error
  error /= normSol;
  if( mpiRank == 0)
  {
    std::cout << "Relative error: " << error << std::endl;
  }
  EXPECT_NEAR( error, 0.0, tol );
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

#ifdef GEOSX_USE_MPI
  MPI_Init(&argc,&argv);

  MPI_Comm_dup( MPI_COMM_WORLD, &MPI_COMM_GEOSX );

  mpiRank = CommunicationTools::MPI_Rank( MPI_COMM_GEOSX );
  mpiSize = CommunicationTools::MPI_Size( MPI_COMM_GEOSX );

  logger::InitializeLogger(MPI_COMM_GEOSX);
#else
  logger::InitializeLogger():
#endif

  cxx_utilities::setSignalHandling(cxx_utilities::handler1);

  global_argc = argc;
  global_argv = new char*[static_cast<unsigned int>(global_argc)];
  for( int i=0 ; i<argc ; ++i )
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
