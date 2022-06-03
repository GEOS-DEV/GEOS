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

#include "constitutive/fluid/SingleFluidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseFVM.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "unitTests/fluidFlowTests/testSingleFlowUtils.hpp"

using namespace geosx;
using namespace geosx::dataRepository;
using namespace geosx::constitutive;
using namespace geosx::testing;

CommandLineOptions g_commandLineOptions;

// Sphinx start after input XML

char const * xmlInput =
  "<Problem>\n"
    "<Solvers>\n"
      "<SinglePhaseFVM\n"
        "name=\"singleflow\"\n"
        "logLevel=\"1\"\n"
        "discretization=\"fluidTPFA\"\n"
        "temperature=\"368.15\"\n"
        "isThermal=\"1\"\n"
        "targetRegions=\"{ region }\">\n"
        "<NonlinearSolverParameters\n"
          "newtonTol=\"1.0e-6\"\n"
          "newtonMaxIter=\"100\"/>\n"
        "<LinearSolverParameters\n"
          "solverType=\"gmres\"\n"
          "krylovTol=\"1.0e-10\"/>\n"
      "</SinglePhaseFVM>\n"
    "</Solvers>\n"
    "<Mesh>\n"
      "<InternalMesh\n"
        "name=\"mesh\"\n"
        "elementTypes=\"{ C3D8 }\"\n"
        "xCoords=\"{ 0, 20 }\"\n"
        "yCoords=\"{ 0, 1 }\"\n"
        "zCoords=\"{ 0, 1 }\"\n"
        "nx=\"{ 5 }\"\n"
        "ny=\"{ 1 }\"\n"
        "nz=\"{ 1 }\"\n"
        "cellBlockNames=\"{ cb }\"/>\n"
    "</Mesh>\n"
    "<Geometry>\n"
      "<Box\n"
        "name=\"sink\"\n"
        "xMin=\"{ -0.01, -0.01, -0.01 }\"\n"
        "xMax=\"{ 4.01, 1.01, 1.01 }\"/>\n"
      "<Box\n"
        "name=\"source\"\n"
        "xMin=\"{ -0.01, -0.01, -0.01 }\"\n"
        "xMax=\"{ 4.01, 1.01, 1.01 }\"/>\n"
    "</Geometry>\n"
    "<Events\n"
      "maxTime=\"1000\">\n"
      "<PeriodicEvent\n"
        "name=\"solverApplications\"\n"
        "maxEventDt=\"1000\"\n"
        "target=\"/Solvers/singleflow\"/>\n"
    "</Events>\n"
    "<NumericalMethods>\n"
      "<FiniteVolume>\n"
        "<TwoPointFluxApproximation\n"
          "name=\"fluidTPFA\"/>\n"
      "</FiniteVolume>\n"
    "</NumericalMethods>\n"
    "<ElementRegions>\n"
      "<CellElementRegion\n"
        "name=\"region\"\n"
        "cellBlocks=\"{ cb }\"\n"
        "materialList=\"{ water, rock, thermalCond }\"/>\n"
    "</ElementRegions>\n"
    "<Constitutive>\n"
      "<CompressibleSolidConstantPermeability\n"
        "name=\"rock\"\n"
        "solidModelName=\"nullSolid\"\n"
        "porosityModelName=\"rockPorosity\"\n"
        "permeabilityModelName=\"rockPerm\"\n"
        "solidInternalEnergyModelName=\"rockInternalEnergy\"/>\n"
      "<NullModel\n"
        "name=\"nullSolid\"/>\n"
      "<PressurePorosity\n"
        "name=\"rockPorosity\"\n"
        "defaultReferencePorosity=\"0.05\"\n"
        "referencePressure=\"0.0\"\n"
        "compressibility=\"1.0e-9\"/>\n"
      "<SolidInternalEnergy\n"
        "name=\"rockInternalEnergy\"\n"
        "volumetricHeatCapacity=\"1.95e6\"\n"
        "referenceTemperature=\"368.15\"\n"
        "referenceInternalEnergy=\"0\"/>\n"
      "<ConstantPermeability\n"
        "name=\"rockPerm\"\n"
        "permeabilityComponents=\"{ 1.0e-13, 1.0e-13, 1.0e-13 }\"/>\n"
      "<CompressibleSinglePhaseFluid\n"
        "name=\"water\"\n"
        "isThermal=\"1\"\n"
        "defaultDensity=\"1000\"\n"
        "defaultViscosity=\"0.001\"\n"
        "referencePressure=\"0.0\"\n"
        "referenceTemperature=\"0.0\"\n"
        "compressibility=\"5e-10\"\n"
        "thermalExpansionCoeff=\"7e-4\"\n"
        "viscosibility=\"0.0\"\n"
        "volumetricHeatCapacity=\"4.5e3\"/>\n"
      "<SinglePhaseConstantThermalConductivity\n"
        "name=\"thermalCond\"\n"
        "thermalConductivityComponents=\"{ 0.6, 0.6, 0.6 }\"/>\n"
    "</Constitutive>\n"
    "<FieldSpecifications>\n"
      "<FieldSpecification\n"
        "name=\"initialPressure\"\n"
        "initialCondition=\"1\"\n"
        "setNames=\"{ all }\"\n"
        "objectPath=\"ElementRegions/region/cb\"\n"
        "fieldName=\"pressure\"\n"
        "scale=\"9e6\"/>\n"
      "<FieldSpecification\n"
        "name=\"initialTemperature\"\n"
        "initialCondition=\"1\"\n"
        "setNames=\"{ all }\"\n"
        "objectPath=\"ElementRegions/region/cb\"\n"
        "fieldName=\"temperature\"\n"
        "scale=\"368.15\"/>\n"
      "<FieldSpecification\n"
        "name=\"sinkPressure\"\n"
        "setNames=\"{ sink }\"\n"
        "objectPath=\"ElementRegions/region/cb\"\n"
        "fieldName=\"pressure\"\n"
        "scale=\"7e6\"/>\n"
      "<FieldSpecification\n"
        "name=\"sinkTemperature\"\n"
        "setNames=\"{ sink }\"\n"
        "objectPath=\"ElementRegions/region/cb\"\n"
        "fieldName=\"temperature\"\n"
        "scale=\"368.15\"/>\n"
      "<FieldSpecification\n"
        "name=\"sourcePressure\"\n"
        "setNames=\"{ source }\"\n"
        "objectPath=\"ElementRegions/region/cb\"\n"
        "fieldName=\"pressure\"\n"
        "scale=\"1.45e7\"/>\n"
      "<FieldSpecification\n"
        "name=\"sourceTemperature\"\n"
        "setNames=\"{ source }\"\n"
        "objectPath=\"ElementRegions/region/cb\"\n"
        "fieldName=\"temperature\"\n"
        "scale=\"300.15\"/>\n"
    "</FieldSpecifications>\n"
  "</Problem>\n";

// Sphinx end before input XML

template< typename LAMBDA >
void testNumericalJacobian( SinglePhaseFVM< SinglePhaseBase > & solver,
                            DomainPartition & domain,
                            real64 const perturbParameter,
                            real64 const relTol,
                            LAMBDA assembleFunction )
{
  CRSMatrix< real64, globalIndex > const & jacobian = solver.getLocalMatrix();
  array1d< real64 > residual( jacobian.numRows() );

  // assemble the analytical residual
  solver.resetStateToBeginningOfStep( domain );

  residual.zero();
  jacobian.zero();

  assembleFunction( jacobian.toViewConstSizes(), residual.toView() );
  residual.move( LvArray::MemorySpace::host, false );

  // copy the analytical residual
  array1d< real64 > residualOrig( residual );

  // create the numerical jacobian
  jacobian.move( LvArray::MemorySpace::host );
  CRSMatrix< real64, globalIndex > jacobianFD( jacobian );
  jacobianFD.zero();

  // fill jacobian FD
  fillCellCenteredNumericalJacobian( solver,
                                     domain,
                                     true,
                                     perturbParameter,
                                     residual.toView(),
                                     residualOrig.toView(),
                                     jacobian.toView(),
                                     jacobianFD.toView(),
                                     assembleFunction );

  // assemble the analytical jacobian
  solver.resetStateToBeginningOfStep( domain );

  residual.zero();
  jacobian.zero();
  assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

  compareLocalMatrices( jacobian.toViewConst(), jacobianFD.toViewConst(), relTol, 1e-6 );
}

class ThermalSinglePhaseFlowTest : public ::testing::Test
{
public:

  ThermalSinglePhaseFlowTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:

  void SetUp() override
  {

    setupProblemFromXML( state.getProblemManager(), xmlInput );
    solver = &state.getProblemManager().getPhysicsSolverManager().getGroup< SinglePhaseFVM< SinglePhaseBase > >( "singleflow" );

    DomainPartition & domain = state.getProblemManager().getDomainPartition();

    solver->setupSystem( domain,
                         solver->getDofManager(),
                         solver->getLocalMatrix(),
                         solver->getSystemRhs(),
                         solver->getSystemSolution() );

    solver->implicitStepSetup( time, dt, domain );
  }

  static real64 constexpr time = 0.0;
  static real64 constexpr dt = 1e4;
  static real64 constexpr eps = std::numeric_limits< real64 >::epsilon();

  GeosxState state;
  SinglePhaseFVM< SinglePhaseBase > * solver;
};

real64 constexpr ThermalSinglePhaseFlowTest::time;
real64 constexpr ThermalSinglePhaseFlowTest::dt;
real64 constexpr ThermalSinglePhaseFlowTest::eps;

TEST_F( ThermalSinglePhaseFlowTest, derivativeNumericalCheck_mobility )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 5e-2; // 5% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testMobilityNumericalDerivatives( *solver, domain, true, perturb, tol );
}

TEST_F( ThermalSinglePhaseFlowTest, jacobianNumericalCheck_flux )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 2e-2; // 2% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    // The first input parameter denotes t_n, which is unused. Just input something here. 
    solver->assembleFluxTerms( 0.0, dt, domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}

#if 1
TEST_F( ThermalSinglePhaseFlowTest, jacobianNumericalCheck_accumulationBalance )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 2e-2; // 2% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->assembleAccumulationTerms( domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}
#endif

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geosx::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
