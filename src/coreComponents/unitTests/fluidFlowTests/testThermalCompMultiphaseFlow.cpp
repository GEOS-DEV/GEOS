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

#include "constitutive/fluid/MultiFluidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"

using namespace geosx;
using namespace geosx::dataRepository;
using namespace geosx::constitutive;
using namespace geosx::constitutive::multifluid;
using namespace geosx::testing;

CommandLineOptions g_commandLineOptions;

// Sphinx start after input XML

char const * pvtLiquid = "DensityFun PhillipsBrineDensity 1e6 7.5e7 5e5 295.15 370.15 25 0\n"
                         "ViscosityFun PhillipsBrineViscosity 0\n"
                         "EnthalpyFun BrineEnthalpy 1e6 7.5e7 5e5 295.15 370.15 25 0\n";

char const * pvtGas = "DensityFun SpanWagnerCO2Density 1e6 7.5e7 5e5 295.15 370.15 25\n"
                      "ViscosityFun FenghourCO2Viscosity 1e6 7.5e7 5e5 295.15 370.15 25\n"
                      "EnthalpyFun CO2Enthalpy 1e6 7.5e7 5e5 295.15 370.15 25\n";

char const * co2flash = "FlashModel CO2Solubility  1e6 7.5e7 5e5 295.15 370.15 25 0";

char const * xmlInput =
  "<Problem>\n"
  "<Solvers>\n"
  "<CompositionalMultiphaseFVM\n"
  "name=\"compflow\"\n"
  "logLevel=\"1\"\n"
  "discretization=\"fluidTPFA\"\n"
  "temperature=\"368.15\"\n"
  "useMass=\"1\"\n"
  "isThermal=\"1\"\n"
  "initialDt=\"1000\"\n"
  "maxCompFractionChange=\"0.5\"\n"
  "targetRegions=\"{ region }\">\n"
  "<NonlinearSolverParameters\n"
  "newtonTol=\"1.0e-6\"\n"
  "newtonMaxIter=\"100\"\n"
  "lineSearchAction=\"None\"\n"
  "maxTimeStepCuts=\"5\"/>\n"
  "<LinearSolverParameters\n"
  "directParallel=\"0\"/>\n"
  "</CompositionalMultiphaseFVM>\n"
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
  "target=\"/Solvers/compflow\"/>\n"
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
  "materialList=\"{ fluid, rock, relperm, thermalCond }\"/>\n"
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
  "defaultReferencePorosity=\"0.2\"\n"
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
  "<CO2BrinePhillipsThermalFluid\n"
  "name=\"fluid\"\n"
  "phaseNames=\"{ gas, water }\"\n"
  "componentNames=\"{ co2, water }\"\n"
  "componentMolarWeight=\"{ 44e-3, 18e-3 }\"\n"
  "phasePVTParaFiles=\"{ pvtgas.txt, pvtliquid.txt }\"\n"
  "flashModelParaFile=\"co2flash.txt\"/>\n"
  "<BrooksCoreyRelativePermeability\n"
  "name=\"relperm\"\n"
  "phaseNames=\"{ gas, water }\"\n"
  "phaseMinVolumeFraction=\"{ 0.0, 0.0 }\"\n"
  "phaseRelPermExponent=\"{ 1.5, 1.5 }\"\n"
  "phaseRelPermMaxValue=\"{ 0.9, 0.9 }\"/>\n"
  "<MultiPhaseConstantThermalConductivity\n"
  "name=\"thermalCond\"\n"
  "phaseNames=\"{ gas, water }\"\n"
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
  "name=\"initialComposition_co2\"\n"
  "initialCondition=\"1\"\n"
  "setNames=\"{ all }\"\n"
  "objectPath=\"ElementRegions/region/cb\"\n"
  "fieldName=\"globalCompFraction\"\n"
  "component=\"0\"\n"
  "scale=\"0.1\"/>\n"
  "<FieldSpecification\n"
  "name=\"initialComposition_water\"\n"
  "initialCondition=\"1\"\n"
  "setNames=\"{ all }\"\n"
  "objectPath=\"ElementRegions/region/cb\"\n"
  "fieldName=\"globalCompFraction\"\n"
  "component=\"1\"\n"
  "scale=\"0.9\"/>\n"
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
  "name=\"sinkTermComposition_co2\"\n"
  "setNames=\"{ sink }\"\n"
  "objectPath=\"ElementRegions/region/cb\"\n"
  "fieldName=\"globalCompFraction\"\n"
  "component=\"0\"\n"
  "scale=\"0.1\"/>\n"
  "<FieldSpecification\n"
  "name=\"sinkTermComposition_water\"\n"
  "setNames=\"{ sink }\"\n"
  "objectPath=\"ElementRegions/region/cb\"\n"
  "fieldName=\"globalCompFraction\"\n"
  "component=\"1\"\n"
  "scale=\"0.9\"/>\n"
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
  "<FieldSpecification\n"
  "name=\"sourceTermComposition_co2\"\n"
  "setNames=\"{ source }\"\n"
  "objectPath=\"ElementRegions/region/cb\"\n"
  "fieldName=\"globalCompFraction\"\n"
  "component=\"0\"\n"
  "scale=\"0.9\"/>\n"
  "<FieldSpecification\n"
  "name=\"sourceTermComposition_water\"\n"
  "setNames=\"{ source }\"\n"
  "objectPath=\"ElementRegions/region/cb\"\n"
  "fieldName=\"globalCompFraction\"\n"
  "component=\"1\"\n"
  "scale=\"0.1\"/>\n"
  "</FieldSpecifications>\n"
  "</Problem>\n";

// Sphinx end before input XML

template< typename LAMBDA >
void testNumericalJacobian( CompositionalMultiphaseFVM & solver,
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

void writeTableToFile( string const & filename, char const * str )
{
  std::ofstream os( filename );
  ASSERT_TRUE( os.is_open() );
  os << str;
  os.close();
}

void removeFile( string const & filename )
{
  int const ret = std::remove( filename.c_str() );
  ASSERT_TRUE( ret == 0 );
}

class ThermalCompositionalMultiphaseFlowTest : public ::testing::Test
{
public:

  ThermalCompositionalMultiphaseFlowTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {
    writeTableToFile( pvtLiquidFilename, pvtLiquid );
    writeTableToFile( pvtGasFilename, pvtGas );
    writeTableToFile( co2flashFilename, co2flash );
  }

  ~ThermalCompositionalMultiphaseFlowTest() override
  {
    removeFile( pvtLiquidFilename );
    removeFile( pvtGasFilename );
    removeFile( co2flashFilename );
  }

protected:

  void SetUp() override
  {

    setupProblemFromXML( state.getProblemManager(), xmlInput );
    solver = &state.getProblemManager().getPhysicsSolverManager().getGroup< CompositionalMultiphaseFVM >( "compflow" );

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
  CompositionalMultiphaseFVM * solver;

  string const pvtLiquidFilename = "pvtliquid.txt";
  string const pvtGasFilename = "pvtgas.txt";
  string const co2flashFilename = "co2flash.txt";
};

real64 constexpr ThermalCompositionalMultiphaseFlowTest::time;
real64 constexpr ThermalCompositionalMultiphaseFlowTest::dt;
real64 constexpr ThermalCompositionalMultiphaseFlowTest::eps;

TEST_F( ThermalCompositionalMultiphaseFlowTest, derivativeNumericalCheck_composition )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-4;

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testCompositionNumericalDerivatives( *solver, domain, perturb, tol );
}

TEST_F( ThermalCompositionalMultiphaseFlowTest, derivativeNumericalCheck_phaseVolumeFraction )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 5e-2; // 5% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();
  testPhaseVolumeFractionNumericalDerivatives( *solver, domain, true, perturb, tol );
}

TEST_F( ThermalCompositionalMultiphaseFlowTest, derivativeNumericalCheck_phaseMobility )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 5e-2; // 5% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testPhaseMobilityNumericalDerivatives( *solver, domain, true, perturb, tol );
}

TEST_F( ThermalCompositionalMultiphaseFlowTest, jacobianNumericalCheck_flux )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->assembleFluxTerms( dt, domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}

#if 0
TEST_F( ThermalCompositionalMultiphaseFlowTest, jacobianNumericalCheck_accumulationVolumeBalance )
{
  real64 const perturb = sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->assembleAccumulationAndVolumeBalanceTerms( domain, solver->getDofManager(), localMatrix, localRhs );
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
