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
#include "mainInterface/initialization.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/GeosxState.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseHybridFVM.hpp"
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"

using namespace geosx;
using namespace geosx::dataRepository;
using namespace geosx::constitutive;
using namespace geosx::testing;

CommandLineOptions g_commandLineOptions;

char const * xmlInput =
  "<Problem>\n"
  "  <Solvers gravityVector=\"{ 0.0, 0.0, -9.81 }\">\n"
  "    <CompositionalMultiphaseHybridFVM name=\"compflow\"\n"
  "                                 logLevel=\"0\"\n"
  "                                 discretization=\"fluidHM\"\n"
  "                                 targetRegions=\"{Region}\"\n"
  "                                 temperature=\"297.15\"\n"
  "                                 useMass=\"1\">\n"
  "                                 \n"
  "      <NonlinearSolverParameters newtonTol=\"1.0e-6\"\n"
  "                                 newtonMaxIter=\"2\"/>\n"
  "      <LinearSolverParameters solverType=\"gmres\"\n"
  "                              krylovTol=\"1.0e-10\"/>\n"
  "    </CompositionalMultiphaseHybridFVM>\n"
  "  </Solvers>\n"
  "  <Mesh>\n"
  "    <InternalMesh name=\"mesh1\"\n"
  "                  elementTypes=\"{C3D8}\" \n"
  "                  xCoords=\"{0, 1}\"\n"
  "                  yCoords=\"{0, 1}\"\n"
  "                  zCoords=\"{0, 10}\"\n"
  "                  nx=\"{1}\"\n"
  "                  ny=\"{1}\"\n"
  "                  nz=\"{4}\"\n"
  "                  cellBlockNames=\"{cb1}\"/>\n"
  "  </Mesh>\n"
  "  <NumericalMethods>\n"
  "    <FiniteVolume>\n"
  "      <HybridMimeticDiscretization name=\"fluidHM\"\n"
  "                                   innerProductType=\"beiraoDaVeigaLipnikovManzini\"/>\n"
  "    </FiniteVolume>\n"
  "  </NumericalMethods>\n"
  "  <ElementRegions>\n"
  "    <CellElementRegion name=\"Region\" cellBlocks=\"{cb1}\" materialList=\"{fluid1, rock, relperm}\" />\n"
  "  </ElementRegions>\n"
  "  <Constitutive>\n"
  "    <CompositionalMultiphaseFluid name=\"fluid1\"\n"
  "                                  phaseNames=\"{oil, gas}\"\n"
  "                                  equationsOfState=\"{PR, PR}\"\n"
  "                                  componentNames=\"{N2, C10, C20, H2O}\"\n"
  "                                  componentCriticalPressure=\"{34e5, 25.3e5, 14.6e5, 220.5e5}\"\n"
  "                                  componentCriticalTemperature=\"{126.2, 622.0, 782.0, 647.0}\"\n"
  "                                  componentAcentricFactor=\"{0.04, 0.443, 0.816, 0.344}\"\n"
  "                                  componentMolarWeight=\"{28e-3, 134e-3, 275e-3, 18e-3}\"\n"
  "                                  componentVolumeShift=\"{0, 0, 0, 0}\"\n"
  "                                  componentBinaryCoeff=\"{ {0, 0, 0, 0},\n"
  "                                                          {0, 0, 0, 0},\n"
  "                                                          {0, 0, 0, 0},\n"
  "                                                          {0, 0, 0, 0} }\"/>\n"
  "    <CompressibleSolidConstantPermeability name=\"rock\"\n"
  "        solidModelName=\"nullSolid\"\n"
  "        porosityModelName=\"rockPorosity\"\n"
  "        permeabilityModelName=\"rockPerm\"/>\n"
  "   <NullModel name=\"nullSolid\"/> \n"
  "   <PressurePorosity name=\"rockPorosity\"\n"
  "                     defaultReferencePorosity=\"0.05\"\n"
  "                     referencePressure = \"0.0\"\n"
  "                     compressibility=\"1.0e-9\"/>\n"
  "  <ConstantPermeability name=\"rockPerm\"\n"
  "                        permeabilityComponents=\"{2.0e-16, 2.0e-16, 2.0e-16}\"/> \n"
  "    <BrooksCoreyRelativePermeability name=\"relperm\"\n"
  "                                     phaseNames=\"{oil, gas}\"\n"
  "                                     phaseMinVolumeFraction=\"{0.1, 0.15}\"\n"
  "                                     phaseRelPermExponent=\"{2.0, 2.0}\"\n"
  "                                     phaseRelPermMaxValue=\"{0.8, 0.9}\"/>\n"
  "  </Constitutive>\n"
  "  <FieldSpecifications>\n"
  "    <FieldSpecification name=\"initialPressure\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region/cb1\"\n"
  "               fieldName=\"pressure\"\n"
  "               functionName=\"initialPressureFunc\"\n"
  "               scale=\"5e6\"/>\n"
  "    <FieldSpecification name=\"initialFacePressure\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"faceManager\"\n"
  "               fieldName=\"facePressure\"\n"
  "               functionName=\"initialFacePressureFunc\"\n"
  "               scale=\"5e6\"/>\n"
  "    <FieldSpecification name=\"initialComposition_N2\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region/cb1\"\n"
  "               fieldName=\"globalCompFraction\"\n"
  "               component=\"0\"\n"
  "               scale=\"0.099\"/>\n"
  "    <FieldSpecification name=\"initialComposition_C10\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region/cb1\"\n"
  "               fieldName=\"globalCompFraction\"\n"
  "               component=\"1\"\n"
  "               scale=\"0.3\"/>\n"
  "    <FieldSpecification name=\"initialComposition_C20\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region/cb1\"\n"
  "               fieldName=\"globalCompFraction\"\n"
  "               component=\"2\"\n"
  "               scale=\"0.6\"/>\n"
  "    <FieldSpecification name=\"initialComposition_H20\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region/cb1\"\n"
  "               fieldName=\"globalCompFraction\"\n"
  "               component=\"3\"\n"
  "               scale=\"0.001\"/>\n"
  "  </FieldSpecifications>\n"
  "  <Functions>\n"
  "    <TableFunction name=\"initialPressureFunc\"\n"
  "                   inputVarNames=\"{elementCenter}\"\n"
  "                   coordinates=\"{0.0, 2.0, 4.0, 6.0, 8.0, 10.0 }\"\n"
  "                   values=\"{ 1.0, 0.5, 0.2, 3.0, 2.1, 1.0 }\"/>\n"
  "    <TableFunction name=\"initialFacePressureFunc\"\n"
  "                   inputVarNames=\"{faceCenter}\"\n"
  "                   coordinates=\"{0.0, 2.0, 4.0, 6.0, 8.0, 10.0 }\"\n"
  "                   values=\"{ 2.0, 0.1, 0.2, 2.1, 2.1, 0.1 }\"/>\n"
  "  </Functions>"
  "</Problem>";

template< typename LAMBDA >
void testNumericalJacobian( CompositionalMultiphaseHybridFVM & solver,
                            DomainPartition & domain,
                            real64 const perturbParameter,
                            real64 const relTol,
                            LAMBDA assembleFunction )
{
  CRSMatrix< real64, globalIndex > const & jacobian = solver.getLocalMatrix();
  array1d< real64 > residual( jacobian.numRows() );
  DofManager const & dofManager = solver.getDofManager();

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

  // fill jacobian FD for cell centered variables
  fillCellCenteredNumericalJacobian( solver,
                                     domain,
                                     false,
                                     perturbParameter,
                                     residual.toView(),
                                     residualOrig.toView(),
                                     jacobian.toView(),
                                     jacobianFD.toView(),
                                     assembleFunction );

  solver.forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                       MeshLevel & mesh,
                                                                       arrayView1d< string const > const & )
  {

    FaceManager & faceManager = mesh.getFaceManager();

    // get the face-based pressure
    arrayView1d< real64 > const & facePres =
      faceManager.getExtrinsicData< extrinsicMeshData::flow::facePressure >();
    facePres.move( LvArray::MemorySpace::host, false );

    string const faceDofKey = dofManager.getKey( CompositionalMultiphaseHybridFVM::viewKeyStruct::faceDofFieldString() );

    arrayView1d< globalIndex const > const & faceDofNumber =
      faceManager.getReference< array1d< globalIndex > >( faceDofKey );
    faceDofNumber.move( LvArray::MemorySpace::host );

    arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();
    faceGhostRank.move( LvArray::MemorySpace::host );

    for( localIndex iface = 0; iface < faceManager.size(); ++iface )
    {
      if( faceGhostRank[iface] >= 0 )
      {
        continue;
      }

      solver.resetStateToBeginningOfStep( domain );

      facePres.move( LvArray::MemorySpace::host, true ); // to get the correct facePres after reset
      real64 const dFP = perturbParameter * ( facePres[iface] + perturbParameter );
      facePres[iface] += dFP;
#if defined(GEOSX_USE_CUDA)
      facePres.move( LvArray::MemorySpace::cuda, false );
#endif


      residual.zero();
      jacobian.zero();
      assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

      fillNumericalJacobian( residual.toViewConst(),
                             residualOrig.toViewConst(),
                             faceDofNumber[iface],
                             dFP,
                             jacobianFD.toViewConstSizes() );
    }

    // assemble the analytical jacobian
    solver.resetStateToBeginningOfStep( domain );

    residual.zero();
    jacobian.zero();
    assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

    compareLocalMatrices( jacobian.toViewConst(), jacobianFD.toViewConst(), relTol );

  } );
}

class CompositionalMultiphaseHybridFlowTest : public ::testing::Test
{
public:

  CompositionalMultiphaseHybridFlowTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:

  void SetUp() override
  {
    setupProblemFromXML( state.getProblemManager(), xmlInput );
    solver = &state.getProblemManager().getPhysicsSolverManager().getGroup< CompositionalMultiphaseHybridFVM >( "compflow" );

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
  CompositionalMultiphaseHybridFVM * solver;
};

real64 constexpr CompositionalMultiphaseHybridFlowTest::time;
real64 constexpr CompositionalMultiphaseHybridFlowTest::dt;
real64 constexpr CompositionalMultiphaseHybridFlowTest::eps;


TEST_F( CompositionalMultiphaseHybridFlowTest, jacobianNumericalCheck_flux )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 5e-3; // 10% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->assembleFluxTerms( dt, domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geosx::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
