/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2020-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// using some utility classes from the following unit test
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"

#include "common/DataTypes.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mainInterface/GeosxState.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/wavePropagation/WaveSolverBase.hpp"
#include "physicsSolvers/wavePropagation/AcousticWaveEquationSEM.hpp"

#include <gtest/gtest.h>

using namespace geosx;
using namespace geosx::dataRepository;
using namespace geosx::testing;

CommandLineOptions g_commandLineOptions;

// This unit test checks the interpolation done to extract seismic traces from a wavefield 

char const * xmlInput = 
  "<?xml version=\"1.0\" ?>\n"
  "\n"
  "<Problem>\n"
  "  <Solvers>\n"
  "    <!-- define the solver -->\n"
  "    <!-- define the source coordinates -->\n"
  "    <!-- define the time source frequency -->\n"
  "    <!-- define the receiver coordinates -->\n"
  "    <AcousticSEM\n"
  "      name=\"acousticSolver\"\n"
  "      cflFactor=\"0.25\"\n"
  "      discretization=\"FE1\"\n"
  "      targetRegions=\"{ Region }\"\n"
  "      sourceCoordinates=\"{ { 55, 55, 55 },\n"
  "                           { 10, 10, 14 } }\"\n"
  "      timeSourceFrequency=\"5.0\"\n"
  "      receiverCoordinates=\"{ { 5, 5, 11 },\n"
  "                             { 5, 50, 11 },\n"
  "                             { 5, 95, 11 } }\"/>\n"
  "  </Solvers>\n"
  "\n"
  "  <!-- hexahedral mesh generated internally by GEOSX -->\n"
  "  <Mesh>\n"
  "    <InternalMesh\n"
  "      name=\"mesh\"\n"
  "      elementTypes=\"{ C3D8 }\"\n"
  "      xCoords=\"{ 0, 101 }\"\n"
  "      yCoords=\"{ 0, 101 }\"\n"
  "      zCoords=\"{ 0, 101 }\"\n"
  "      nx=\"{ 10 }\"\n"
  "      ny=\"{ 10 }\"\n"
  "      nz=\"{ 10 }\"\n"
  "      cellBlockNames=\"{ cb }\"/>\n"
  "  </Mesh>\n"
  "\n"
  "  <Events\n"
  "    maxTime=\"0.2\">\n"
  "    <!-- control the timestepping here with forceDt -->\n"
  "    <PeriodicEvent\n"
  "      name=\"solverApplications\"\n"
  "      forceDt=\"0.005\"\n"
  "      target=\"/Solvers/acousticSolver\"/>\n"
  "\n"
  "    <!-- generate an output that can be read from VTK -->\n"
  "    <PeriodicEvent\n"
  "      name=\"vtk\"\n"
  "      timeFrequency=\"0.1\"\n"
  "      targetExactTimestep=\"0\"\n"
  "      target=\"/Outputs/vtkOutput\"/>\n"
  "\n"
  "    <!-- two events to output pressure in an hdf5 file -->\n"
  "    <PeriodicEvent\n"
  "      name=\"timeHistoryCollection\"\n"
  "      timeFrequency=\"0.005\"\n"
  "      target=\"/Tasks/pressureCollection\"/>\n"
  "\n"
  "    <PeriodicEvent\n"
  "      name=\"timeHistoryOutput\"\n"
  "      timeFrequency=\"0.05\"\n"
  "      targetExactTimestep=\"0\"\n"
  "      target=\"/Outputs/timeHistoryOutput\"/>\n"
  "\n"
  "    <!-- restart event -->\n"
  "    <PeriodicEvent\n"
  "      name=\"restarts\"\n"
  "      timeFrequency=\"0.1\"\n"
  "      target=\"/Outputs/restartOutput\"/>\n"
  "  </Events>\n"
  "\n"
  "  <NumericalMethods>\n"
  "    <FiniteElements>\n"
  "      <FiniteElementSpace\n"
  "        name=\"FE1\"\n"
  "        order=\"1\"/>\n"
  "    </FiniteElements>\n"
  "  </NumericalMethods>\n"
  "\n"
  "  <ElementRegions>\n"
  "    <CellElementRegion\n"
  "      name=\"Region\"\n"
  "      cellBlocks=\"{ cb }\"\n"
  "      materialList=\"{ nullModel }\"/>\n"
  "  </ElementRegions>\n"
  "\n"
  "  <Constitutive>\n"
  "    <NullModel\n"
  "      name=\"nullModel\"/>\n"
  "  </Constitutive>\n"
  "\n"
  "  <FieldSpecifications>\n"
  "    <!-- 1) The initial pressure field -->\n"
  "    <FieldSpecification\n"
  "      name=\"initialPressure\"\n"
  "      initialCondition=\"1\"\n"
  "      setNames=\"{ all }\"\n"
  "      objectPath=\"nodeManager\"\n"
  "      fieldName=\"pressure_n\"\n"
  "      scale=\"0.0\"/>\n"
  "\n"
  "    <FieldSpecification\n"
  "      name=\"initialPressure\"\n"
  "      initialCondition=\"1\"\n"
  "      setNames=\"{ all }\"\n"
  "      objectPath=\"nodeManager\"\n"
  "      fieldName=\"pressure_nm1\"\n"
  "      scale=\"0.0\"/>\n"
  "\n"
  "    <!-- 2) The velocity in the domain -->\n"
  "    <FieldSpecification\n"
  "      name=\"cellVelocity\"\n"
  "      initialCondition=\"1\"\n"
  "      objectPath=\"ElementRegions/Region/elementSubRegions/cb\"\n"
  "      fieldName=\"mediumVelocity\"\n"
  "      scale=\"1500\"\n"
  "      setNames=\"{ all }\"/>\n"
  "  </FieldSpecifications>\n"
  "\n"
  "  <!-- collect the pressure values at the nodes -->\n"
  "  <Tasks>\n"
  "    <PackCollection\n"
  "      name=\"pressureCollection\"\n"
  "      objectPath=\"nodeManager\"\n"
  "      fieldName=\"pressure_np1\"/>\n"
  "  </Tasks>\n"
  "\n"
  "  <Outputs>\n"
  "    <!-- output all the mesh values registered with a plot level LEVEL_0, LEVEL_1, LEVEL_2, LEVEL_3   -->\n"
  "    <VTK\n"
  "      name=\"vtkOutput\"\n"
  "      plotLevel=\"3\"/>\n"
  "\n"
  "    <!-- output the pressure values to a file named pressure_history.hdf5  -->\n"
  "    <TimeHistory\n"
  "      name=\"timeHistoryOutput\"\n"
  "      sources=\"{ /Tasks/pressureCollection }\"\n"
  "      filename=\"pressure_history\"/>\n"
  "\n"
  "    <Restart\n"
  "      name=\"restartOutput\"/>\n"
  "  </Outputs>\n"
  "</Problem>\n";

class AcousticWaveEquationSEMTest : public ::testing::Test
{
public:

  AcousticWaveEquationSEMTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:

  void SetUp() override
  {
    setupProblemFromXML( state.getProblemManager(), xmlInput );
    propagator = &state.getProblemManager().getPhysicsSolverManager().getGroup< AcousticWaveEquationSEM >( "acousticSolver" );

    DomainPartition & domain = state.getProblemManager().getDomainPartition();

    propagator->setupSystem( domain,
                         propagator->getDofManager(),
                         propagator->getLocalMatrix(),
                         propagator->getSystemRhs(),
                         propagator->getSystemSolution() );

    propagator->implicitStepSetup( time, dt, domain );
  }

  static real64 constexpr time = 0.0;
  static real64 constexpr dt = 1e-3;
  static real64 constexpr eps = std::numeric_limits< real64 >::epsilon();

  GeosxState state;
  AcousticWaveEquationSEM  * propagator;
};

real64 constexpr AcousticWaveEquationSEMTest::time;
real64 constexpr AcousticWaveEquationSEMTest::dt;
real64 constexpr AcousticWaveEquationSEMTest::eps;

TEST_F( AcousticWaveEquationSEMTest, SeismoTrace )
{

  DomainPartition & domain = state.getProblemManager().getDomainPartition();
 
  state.applyInitialConditions(); 
  state.run();

}



int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geosx::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
