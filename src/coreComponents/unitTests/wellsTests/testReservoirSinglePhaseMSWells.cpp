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


#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"

#include "common/DataTypes.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/multiphysics/SinglePhaseReservoirAndWells.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseFVM.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWell.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWellExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseExtrinsicData.hpp"

using namespace geosx;
using namespace geosx::dataRepository;
using namespace geosx::constitutive;
using namespace geosx::testing;

CommandLineOptions g_commandLineOptions;

char const * xmlInput =
  "<Problem>\n"
  "  <Solvers gravityVector=\"{ 0.0, 0.0, -9.81 }\">\n"
  "    <SinglePhaseReservoir name=\"reservoirSystem\"\n"
  "               flowSolverName=\"singlePhaseFlow\"\n"
  "               wellSolverName=\"singlePhaseWell\"\n"
  "               logLevel=\"1\"\n"
  "               targetRegions=\"{Region1,wellRegion1,wellRegion2}\">\n"
  "      <NonlinearSolverParameters newtonMaxIter=\"40\"/>\n"
  "      <LinearSolverParameters solverType=\"direct\"\n"
  "                              logLevel=\"2\"/>\n"
  "    </SinglePhaseReservoir>\n"
  "    <SinglePhaseFVM name=\"singlePhaseFlow\"\n"
  "                             logLevel=\"1\"\n"
  "                             discretization=\"singlePhaseTPFA\"\n"
  "                             targetRegions=\"{Region1}\">\n"
  "    </SinglePhaseFVM>\n"
  "    <SinglePhaseWell name=\"singlePhaseWell\"\n"
  "                     logLevel=\"1\"\n"
  "                     targetRegions=\"{wellRegion1,wellRegion2}\">\n"
  "        <WellControls name=\"wellControls1\"\n"
  "                      type=\"producer\"\n"
  "                      referenceElevation=\"2\"\n"
  "                      control=\"BHP\"\n"
  "                      targetBHP=\"5e5\"\n"
  "                      targetTotalRate=\"1e-3\"/>\n"
  "        <WellControls name=\"wellControls2\"\n"
  "                      type=\"injector\"\n"
  "                      referenceElevation=\"2\"\n"
  "                      control=\"totalVolRate\" \n"
  "                      targetBHP=\"2e7\"\n"
  "                      targetTotalRate=\"1e-4\"/>\n"
  "    </SinglePhaseWell>\n"
  "  </Solvers>\n"
  "  <Mesh>\n"
  "    <InternalMesh name=\"mesh1\"\n"
  "                  elementTypes=\"{C3D8}\" \n"
  "                  xCoords=\"{0, 5}\"\n"
  "                  yCoords=\"{0, 1}\"\n"
  "                  zCoords=\"{0, 1}\"\n"
  "                  nx=\"{3}\"\n"
  "                  ny=\"{1}\"\n"
  "                  nz=\"{1}\"\n"
  "                  cellBlockNames=\"{cb1}\"/>\n"
  "    <InternalWell name=\"well_producer1\"\n"
  "                  wellRegionName=\"wellRegion1\"\n"
  "                  wellControlsName=\"wellControls1\"\n"
  "                  meshName=\"mesh1\"\n"
  "                  polylineNodeCoords=\"{ {4.5, 0,  2  },\n"
  "                                         {4.5, 0,  0.5} }\"\n"
  "                  polylineSegmentConn=\"{ {0, 1} }\"\n"
  "                  radius=\"0.1\"\n"
  "                  numElementsPerSegment=\"1\">\n"
  "        <Perforation name=\"producer1_perf1\"\n"
  "                     distanceFromHead=\"1.45\"/> \n"
  "    </InternalWell> \n"
  "    <InternalWell name=\"well_injector1\"\n"
  "                  wellRegionName=\"wellRegion2\"\n"
  "                  wellControlsName=\"wellControls2\"\n"
  "                  meshName=\"mesh1\"\n"
  "                  polylineNodeCoords=\"{ {0.5, 0, 2  },\n"
  "                                         {0.5, 0, 0.5} }\"\n"
  "                  polylineSegmentConn=\"{ {0, 1} }\"\n"
  "                  radius=\"0.1\"\n"
  "                  numElementsPerSegment=\"1\">\n"
  "        <Perforation name=\"injector1_perf1\"\n"
  "                     distanceFromHead=\"1.45\"/> \n"
  "    </InternalWell>\n"
  "  </Mesh>\n"
  "  <NumericalMethods>\n"
  "    <FiniteVolume>\n"
  "      <TwoPointFluxApproximation name=\"singlePhaseTPFA\"/>\n"
  "    </FiniteVolume>\n"
  "  </NumericalMethods>\n"
  "  <ElementRegions>\n"
  "    <CellElementRegion name=\"Region1\"\n"
  "                       cellBlocks=\"{cb1}\"\n"
  "                       materialList=\"{water, rock}\"/>\n"
  "    <WellElementRegion name=\"wellRegion1\"\n"
  "                       materialList=\"{water}\"/> \n"
  "    <WellElementRegion name=\"wellRegion2\"\n"
  "                       materialList=\"{water}\"/> \n"
  "  </ElementRegions>\n"
  "  <Constitutive>\n"
  "    <CompressibleSinglePhaseFluid name=\"water\"\n"
  "                                  defaultDensity=\"1000\"\n"
  "                                  defaultViscosity=\"0.001\"\n"
  "                                  referencePressure=\"0.0\"\n"
  "                                  referenceDensity=\"1000\"\n"
  "                                  compressibility=\"5e-10\"\n"
  "                                  referenceViscosity=\"0.001\"\n"
  "                                  viscosibility=\"0.0\"/>\n"
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
  "  </Constitutive>\n"
  "  <FieldSpecifications>\n"
  "    <FieldSpecification name=\"initialPressure\"\n"
  "                        initialCondition=\"1\"\n"
  "                        setNames=\"{all}\"\n"
  "                        objectPath=\"ElementRegions/Region1/cb1\"\n"
  "                        fieldName=\"pressure\"\n"
  "                        scale=\"5e6\"/>\n"
  "  </FieldSpecifications>\n"
  "</Problem>";

template< typename LAMBDA >
void testNumericalJacobian( SinglePhaseReservoirAndWells< SinglePhaseBase > & solver,
                            DomainPartition & domain,
                            real64 const perturbParameter,
                            real64 const relTol,
                            LAMBDA && assembleFunction )
{
  SinglePhaseWell & wellSolver = *solver.wellSolver();
  SinglePhaseFVM< SinglePhaseBase > & flowSolver = dynamicCast< SinglePhaseFVM< SinglePhaseBase > & >( *solver.reservoirSolver() );

  CRSMatrix< real64, globalIndex > const & jacobian = solver.getLocalMatrix();
  array1d< real64 > residual( jacobian.numRows() );
  DofManager const & dofManager = solver.getDofManager();

  // assemble the analytical residual
  solver.resetStateToBeginningOfStep( domain );

  residual.zero();
  jacobian.zero();

  assembleFunction( jacobian.toViewConstSizes(), residual.toView() );
  residual.move( hostMemorySpace, false );

  // copy the analytical residual
  array1d< real64 > residualOrig( residual );

  // create the numerical jacobian
  jacobian.move( hostMemorySpace );
  CRSMatrix< real64, globalIndex > jacobianFD( jacobian );
  jacobianFD.zero();

  string const & resDofKey  = dofManager.getKey( wellSolver.resElementDofName() );
  string const & wellDofKey = dofManager.getKey( wellSolver.wellElementDofName() );

  // at this point we start assembling the finite-difference block by block

  ////////////////////////////////////////////////
  // Step 1) Compute the terms in J_RR and J_WR //
  ////////////////////////////////////////////////

  domain.forMeshBodies( [&] ( MeshBody & meshBody )
  {
    meshBody.forMeshLevels( [&] ( MeshLevel & mesh )
    {
      ElementRegionManager & elemManager = mesh.getElemManager();

      for( localIndex er = 0; er < elemManager.numRegions(); ++er )
      {
        ElementRegionBase & elemRegion = elemManager.getRegion( er );
        elemRegion.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
        {
          // get the dof numbers and ghosting information
          arrayView1d< globalIndex const > const & dofNumber =
            subRegion.getReference< array1d< globalIndex > >( resDofKey );

          // get the primary variables on reservoir elements
          arrayView1d< real64 > const & pres =
            subRegion.getExtrinsicData< extrinsicMeshData::well::pressure >();
          pres.move( hostMemorySpace, false );

          // a) compute all the derivatives wrt to the pressure in RESERVOIR elem ei
          for( localIndex ei = 0; ei < subRegion.size(); ++ei )
          {
            {
              solver.resetStateToBeginningOfStep( domain );

              // here is the perturbation in the pressure of the element
              real64 const dP = perturbParameter * (pres[ei] + perturbParameter);
              pres.move( hostMemorySpace, true );
              pres[ei] += dP;

              // after perturbing, update the pressure-dependent quantities in the reservoir
              flowSolver.forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                       MeshLevel & mesh2,
                                                                       arrayView1d< string const > const & regionNames2 )
              {
                mesh2.getElemManager().forElementSubRegions( regionNames2,
                                                             [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion2 )
                {
                  flowSolver.updateFluidState( subRegion2 );
                } );
              } );

              wellSolver.updateState( domain );

              residual.zero();
              jacobian.zero();
              assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

              fillNumericalJacobian( residual.toViewConst(),
                                     residualOrig.toViewConst(),
                                     dofNumber[ei],
                                     dP,
                                     jacobianFD.toViewConstSizes() );
            }
          }
        } );
      }
    } );
  } );

  /////////////////////////////////////////////////
  // Step 2) Compute the terms in J_RW and J_WW //
  /////////////////////////////////////////////////

  // loop over the wells
  wellSolver.forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                           MeshLevel & mesh,
                                                           arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                        [&]( localIndex const,
                                                                             WellElementSubRegion & subRegion )
    {

      // get the degrees of freedom and ghosting information
      arrayView1d< globalIndex const > const & wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );

      // get the primary variables on well elements
      arrayView1d< real64 > const & wellElemPressure =
        subRegion.getExtrinsicData< extrinsicMeshData::well::pressure >();
      wellElemPressure.move( hostMemorySpace, false );

      arrayView1d< real64 > const & connRate =
        subRegion.getExtrinsicData< extrinsicMeshData::well::connectionRate >();
      connRate.move( hostMemorySpace, false );

      // a) compute all the derivatives wrt to the pressure in WELL elem iwelem
      for( localIndex iwelem = 0; iwelem < subRegion.size(); ++iwelem )
      {
        {
          solver.resetStateToBeginningOfStep( domain );

          // here is the perturbation in the pressure of the well element
          real64 const dP = perturbParameter * ( wellElemPressure[iwelem] + perturbParameter );
          wellElemPressure.move( hostMemorySpace, true );
          wellElemPressure[iwelem] += dP;

          // after perturbing, update the pressure-dependent quantities in the well
          wellSolver.updateState( domain );

          residual.zero();
          jacobian.zero();
          assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

          // consider mass balance eq lid in RESERVOIR elems and WELL elems
          //      this is computing J_RW and J_WW
          fillNumericalJacobian( residual.toViewConst(),
                                 residualOrig.toViewConst(),
                                 wellElemDofNumber[iwelem] + singlePhaseWellKernels::ColOffset::DPRES,
                                 dP,
                                 jacobianFD.toViewConstSizes() );
        }
      }

      // b) compute all the derivatives wrt to the connection in WELL elem iwelem
      for( localIndex iwelem = 0; iwelem < subRegion.size(); ++iwelem )
      {
        {
          solver.resetStateToBeginningOfStep( domain );

          // here is the perturbation in the pressure of the well element
          real64 const dRate = perturbParameter * ( connRate[iwelem] + perturbParameter );
          connRate.move( hostMemorySpace, true );
          connRate[iwelem] += dRate;

          // after perturbing, update the rate-dependent quantities in the well (well controls)
          wellSolver.updateState( domain );

          residual.zero();
          jacobian.zero();
          assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

          // consider mass balance eq lid in RESERVOIR elems and WELL elems
          //      this is computing J_RW and J_WW
          fillNumericalJacobian( residual.toViewConst(),
                                 residualOrig.toViewConst(),
                                 wellElemDofNumber[iwelem] + singlePhaseWellKernels::ColOffset::DRATE,
                                 dRate,
                                 jacobianFD.toViewConstSizes() );
        }
      }
    } );
  } );

  // assemble the analytical jacobian
  solver.resetStateToBeginningOfStep( domain );

  residual.zero();
  jacobian.zero();
  assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

  compareLocalMatrices( jacobian.toViewConst(), jacobianFD.toViewConst(), relTol );
}

class SinglePhaseReservoirSolverTest : public ::testing::Test
{
public:

  SinglePhaseReservoirSolverTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:

  void SetUp() override
  {
    setupProblemFromXML( state.getProblemManager(), xmlInput );
    solver = &state.getProblemManager().getPhysicsSolverManager().getGroup< SinglePhaseReservoirAndWells< SinglePhaseBase > >( "reservoirSystem" );

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
  SinglePhaseReservoirAndWells< SinglePhaseBase > * solver;
};

real64 constexpr SinglePhaseReservoirSolverTest::time;
real64 constexpr SinglePhaseReservoirSolverTest::dt;
real64 constexpr SinglePhaseReservoirSolverTest::eps;

TEST_F( SinglePhaseReservoirSolverTest, jacobianNumericalCheck_Perforation )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->assembleCouplingTerms( time, dt, domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}

TEST_F( SinglePhaseReservoirSolverTest, jacobianNumericalCheck_Flux )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->wellSolver()->assembleFluxTerms( time, dt, domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}

TEST_F( SinglePhaseReservoirSolverTest, jacobianNumericalCheck_PressureRel )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->wellSolver()->assemblePressureRelations( domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}

TEST_F( SinglePhaseReservoirSolverTest, jacobianNumericalCheck_Accum )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->wellSolver()->assembleAccumulationTerms( domain, solver->getDofManager(), localMatrix, localRhs );
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
