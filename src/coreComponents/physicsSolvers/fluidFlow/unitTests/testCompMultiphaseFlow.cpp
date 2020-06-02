/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "constitutive/fluid/MultiFluidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/initialization.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "managers/ProblemManager.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/unitTests/testCompFlowUtils.hpp"

using namespace geosx;
using namespace geosx::dataRepository;
using namespace geosx::constitutive;
using namespace geosx::testing;

char const * xmlInput =
  "<Problem>\n"
  "  <Solvers gravityVector=\"0.0, 0.0, -9.81\">\n"
  "    <CompositionalMultiphaseFlow name=\"compflow\"\n"
  "                                 logLevel=\"0\"\n"
  "                                 discretization=\"fluidTPFA\"\n"
  "                                 targetRegions=\"{Region2}\"\n"
  "                                 fluidNames=\"{fluid1}\"\n"
  "                                 solidNames=\"{rock}\"\n"
  "                                 relPermNames=\"{relperm}\"\n"
  "                                 capPressureNames=\"{cappressure}\"\n"
  "                                 temperature=\"297.15\"\n"
  "                                 useMass=\"1\">\n"
  "                                 \n"
  "      <NonlinearSolverParameters newtonTol=\"1.0e-6\"\n"
  "                                 newtonMaxIter=\"2\"/>\n"
  "      <LinearSolverParameters solverType=\"gmres\"\n"
  "                              krylovTol=\"1.0e-10\"/>\n"
  "    </CompositionalMultiphaseFlow>\n"
  "  </Solvers>\n"
  "  <Mesh>\n"
  "    <InternalMesh name=\"mesh1\"\n"
  "                  elementTypes=\"{C3D8}\" \n"
  "                  xCoords=\"{0, 3}\"\n"
  "                  yCoords=\"{0, 1}\"\n"
  "                  zCoords=\"{0, 1}\"\n"
  "                  nx=\"{3}\"\n"
  "                  ny=\"{1}\"\n"
  "                  nz=\"{1}\"\n"
  "                  cellBlockNames=\"{cb1}\"/>\n"
  "  </Mesh>\n"
  "  <Geometry>\n"
  "    <Box name=\"source\" xMin=\"-0.01, -0.01, -0.01\" xMax=\"1.01, 1.01, 1.01\"/>\n"
  "    <Box name=\"sink\"   xMin=\"1.99, -0.01, -0.01\" xMax=\"3.01, 1.01, 1.01\"/>\n"
  "  </Geometry>\n"
  "  <NumericalMethods>\n"
  "    <FiniteVolume>\n"
  "      <TwoPointFluxApproximation name=\"fluidTPFA\"\n"
  "                                 fieldName=\"pressure\"\n"
  "                                 coefficientName=\"permeability\"/>\n"
  "    </FiniteVolume>\n"
  "  </NumericalMethods>\n"
  "  <ElementRegions>\n"
  "    <CellElementRegion name=\"Region2\" cellBlocks=\"{cb1}\" materialList=\"{fluid1, rock, relperm, cappressure}\" />\n"
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
  "    <PoreVolumeCompressibleSolid name=\"rock\"\n"
  "                                 referencePressure=\"0.0\"\n"
  "                                 compressibility=\"1e-9\"/>\n"
  "    <BrooksCoreyRelativePermeability name=\"relperm\"\n"
  "                                     phaseNames=\"{oil, gas}\"\n"
  "                                     phaseMinVolumeFraction=\"{0.1, 0.15}\"\n"
  "                                     phaseRelPermExponent=\"{2.0, 2.0}\"\n"
  "                                     phaseRelPermMaxValue=\"{0.8, 0.9}\"/>\n"
  "    <BrooksCoreyCapillaryPressure name=\"cappressure\"\n"
  "                                  phaseNames=\"{oil, gas}\"\n"
  "                                  phaseMinVolumeFraction=\"{0.2, 0.05}\"\n"
  "                                  phaseCapPressureExponentInv=\"{4.25, 3.5}\"\n"
  "                                  phaseEntryPressure=\"{0., 1e8}\"\n"
  "                                  capPressureEpsilon=\"0.0\"/> \n"
  "  </Constitutive>\n"
  "  <FieldSpecifications>\n"
  "    <FieldSpecification name=\"permx\"\n"
  "               component=\"0\"\n"
  "               initialCondition=\"1\"  \n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region2/cb1\"\n"
  "               fieldName=\"permeability\"\n"
  "               scale=\"2.0e-16\"/>\n"
  "    <FieldSpecification name=\"permy\"\n"
  "               component=\"1\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region2/cb1\"\n"
  "               fieldName=\"permeability\"\n"
  "               scale=\"2.0e-16\"/>\n"
  "    <FieldSpecification name=\"permz\"\n"
  "               component=\"2\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region2/cb1\"\n"
  "               fieldName=\"permeability\"\n"
  "               scale=\"2.0e-16\"/>\n"
  "    <FieldSpecification name=\"referencePorosity\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region2/cb1\"\n"
  "               fieldName=\"referencePorosity\"\n"
  "               scale=\"0.05\"/>\n"
  "    <FieldSpecification name=\"initialPressure\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region2/cb1\"\n"
  "               fieldName=\"pressure\"\n"
  "               functionName=\"initialPressureFunc\"\n"
  "               scale=\"5e6\"/>\n"
  "    <FieldSpecification name=\"initialComposition_N2\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region2/cb1\"\n"
  "               fieldName=\"globalCompFraction\"\n"
  "               component=\"0\"\n"
  "               scale=\"0.099\"/>\n"
  "    <FieldSpecification name=\"initialComposition_C10\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region2/cb1\"\n"
  "               fieldName=\"globalCompFraction\"\n"
  "               component=\"1\"\n"
  "               scale=\"0.3\"/>\n"
  "    <FieldSpecification name=\"initialComposition_C20\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region2/cb1\"\n"
  "               fieldName=\"globalCompFraction\"\n"
  "               component=\"2\"\n"
  "               scale=\"0.6\"/>\n"
  "    <FieldSpecification name=\"initialComposition_H20\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region2/cb1\"\n"
  "               fieldName=\"globalCompFraction\"\n"
  "               component=\"3\"\n"
  "               scale=\"0.001\"/>\n"
  "  </FieldSpecifications>\n"
  "  <Functions>\n"
  "    <TableFunction name=\"initialPressureFunc\"\n"
  "                   inputVarNames=\"{elementCenter}\"\n"
  "                   coordinates=\"{0.0, 3.0}\"\n"
  "                   values=\"{1.0, 0.5}\"/>\n"
  "  </Functions>"
  "</Problem>";

void testCompositionNumericalDerivatives( CompositionalMultiphaseFlow & solver,
                                          DomainPartition & domain,
                                          real64 const perturbParameter,
                                          real64 const relTol )
{
  localIndex const NC = solver.numFluidComponents();

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  solver.forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                         ElementSubRegionBase & subRegion )
  {
    SCOPED_TRACE( subRegion.getParent()->getParent()->getName() + "/" + subRegion.getName() );

    string const & fluidName = solver.fluidModelNames()[targetIndex];
    Group const * const constitutiveGroup = subRegion.GetConstitutiveModels();
    MultiFluidBase const & fluid = *constitutiveGroup->GetGroup< MultiFluidBase >( fluidName );
    arrayView1d< string const > const & components = fluid.componentNames();

    arrayView2d< real64 > & compDens =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::globalCompDensityString );

    arrayView2d< real64 > & dCompDens =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::deltaGlobalCompDensityString );

    arrayView2d< real64 > & compFrac =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::globalCompFractionString );

    arrayView3d< real64 > & dCompFrac_dCompDens =
      subRegion.getReference< array3d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

    // reset the solver state to zero out variable updates
    solver.ResetStateToBeginningOfStep( &domain );

    // make a copy of unperturbed values of component fractions
    array2d< real64 > compFracOrig( subRegion.size(), NC );
    for( localIndex ei = 0; ei < subRegion.size(); ++ei )
    {
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        compFracOrig[ei][ic] = compFrac[ei][ic];
      }
    }

    // update component density and check derivatives
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      // reset the solver state to zero out variable updates (resetting the whole domain is overkill...)
      solver.ResetStateToBeginningOfStep( &domain );

      // perturb a single component density in each cell
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        real64 const dRho = perturbParameter * ( compDens[ei][jc] + perturbParameter );
        dCompDens[ei][jc] = dRho;
      }

      // recompute component fractions
      solver.UpdateComponentFraction( subRegion );

      // check values in each cell
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        SCOPED_TRACE( "Element " + std::to_string( ei ) );

        auto dZ_dRho = invertLayout( dCompFrac_dCompDens[ei], NC, NC );
        string var = "compDens[" + components[jc] + "]";

        checkDerivative( compFrac[ei].toSliceConst(),
                         compFracOrig[ei].toSliceConst(),
                         dZ_dRho[jc].toSliceConst(),
                         dCompDens[ei][jc],
                         relTol,
                         "compFrac",
                         var,
                         components );
      }
    }
  } );
}


void testPhaseVolumeFractionNumericalDerivatives( CompositionalMultiphaseFlow & solver,
                                                  DomainPartition & domain,
                                                  real64 const perturbParameter,
                                                  real64 const relTol )
{
  localIndex const NC = solver.numFluidComponents();
  localIndex const NP = solver.numFluidPhases();

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  solver.forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                         ElementSubRegionBase & subRegion )
  {
    SCOPED_TRACE( subRegion.getParent()->getParent()->getName() + "/" + subRegion.getName() );

    string const & fluidName = solver.fluidModelNames()[targetIndex];
    Group const * const constitutiveGroup = subRegion.GetConstitutiveModels();
    MultiFluidBase const & fluid = *constitutiveGroup->GetGroup< MultiFluidBase >( fluidName );
    arrayView1d< string const > const & components = fluid.componentNames();
    arrayView1d< string const > const & phases = fluid.phaseNames();

    arrayView1d< real64 > & pres =
      subRegion.getReference< array1d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::pressureString );

    arrayView1d< real64 > & dPres =
      subRegion.getReference< array1d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::deltaPressureString );

    arrayView2d< real64 > & compDens =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::globalCompDensityString );

    arrayView2d< real64 > & dCompDens =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::deltaGlobalCompDensityString );

    arrayView2d< real64 > & phaseVolFrac =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::phaseVolumeFractionString );

    arrayView2d< real64 > & dPhaseVolFrac_dPres =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::dPhaseVolumeFraction_dPressureString );

    arrayView3d< real64 > & dPhaseVolFrac_dCompDens =
      subRegion.getReference< array3d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

    // reset the solver state to zero out variable updates
    solver.ResetStateToBeginningOfStep( &domain );

    // make a copy of unperturbed values of component fractions
    array2d< real64 > phaseVolFracOrig( subRegion.size(), NP );
    for( localIndex ei = 0; ei < subRegion.size(); ++ei )
    {
      for( localIndex ip = 0; ip < NP; ++ip )
      {
        phaseVolFracOrig[ei][ip] = phaseVolFrac[ei][ip];
      }
    }

    // update pressure and check derivatives
    {
      // perturb pressure in each cell
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        real64 const dP = perturbParameter * ( pres[ei] + perturbParameter );
        dPres[ei] = dP;
      }

      // recompute component fractions
      solver.UpdateState( subRegion, targetIndex );

      // check values in each cell
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        SCOPED_TRACE( "Element " + std::to_string( ei ) );

        checkDerivative( phaseVolFrac[ei].toSliceConst(),
                         phaseVolFracOrig[ei].toSliceConst(),
                         dPhaseVolFrac_dPres[ei].toSliceConst(),
                         dPres[ei],
                         relTol,
                         "phaseVolFrac",
                         "Pres",
                         phases );
      }
    }

    // update component density and check derivatives
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      // reset the solver state to zero out variable updates (resetting the whole domain is overkill...)
      solver.ResetStateToBeginningOfStep( &domain );

      // perturb a single component density in each cell
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        real64 const dRho = perturbParameter * ( compDens[ei][jc] + perturbParameter );
        dCompDens[ei][jc] = dRho;
      }

      // recompute component fractions
      solver.UpdateState( subRegion, targetIndex );

      // check values in each cell
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        SCOPED_TRACE( "Element " + std::to_string( ei ) );

        auto dS_dRho = invertLayout( dPhaseVolFrac_dCompDens[ei], NP, NC );
        string var = "compDens[" + components[jc] + "]";

        checkDerivative( phaseVolFrac[ei].toSliceConst(),
                         phaseVolFracOrig[ei].toSliceConst(),
                         dS_dRho[jc].toSliceConst(),
                         dCompDens[ei][jc],
                         relTol,
                         "phaseVolFrac",
                         var,
                         phases );
      }
    }
  } );
}

void testPhaseMobilityNumericalDerivatives( CompositionalMultiphaseFlow & solver,
                                            DomainPartition & domain,
                                            real64 const perturbParameter,
                                            real64 const relTol )
{
  localIndex const NC = solver.numFluidComponents();
  localIndex const NP = solver.numFluidPhases();

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  solver.forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                         ElementSubRegionBase & subRegion )
  {
    SCOPED_TRACE( subRegion.getParent()->getName() + "/" + subRegion.getName() );

    string const & fluidName = solver.fluidModelNames()[targetIndex];
    Group const * const constitutiveGroup = subRegion.GetConstitutiveModels();
    MultiFluidBase const & fluid = *constitutiveGroup->GetGroup< MultiFluidBase >( fluidName );
    arrayView1d< string const > const & components = fluid.componentNames();
    arrayView1d< string const > const & phases = fluid.phaseNames();

    arrayView1d< real64 > & pres =
      subRegion.getReference< array1d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::pressureString );

    arrayView1d< real64 > & dPres =
      subRegion.getReference< array1d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::deltaPressureString );

    arrayView2d< real64 > & compDens =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::globalCompDensityString );

    arrayView2d< real64 > & dCompDens =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::deltaGlobalCompDensityString );

    arrayView2d< real64 > & phaseMob =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::phaseMobilityString );

    arrayView2d< real64 > & dPhaseMob_dPres =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::dPhaseMobility_dPressureString );

    arrayView3d< real64 > & dPhaseMob_dCompDens =
      subRegion.getReference< array3d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::dPhaseMobility_dGlobalCompDensityString );

    // reset the solver state to zero out variable updates
    solver.ResetStateToBeginningOfStep( &domain );

    // make a copy of unperturbed values of component fractions
    array2d< real64 > phaseVolFracOrig( subRegion.size(), NP );
    for( localIndex ei = 0; ei < subRegion.size(); ++ei )
    {
      for( localIndex ip = 0; ip < NP; ++ip )
      {
        phaseVolFracOrig[ei][ip] = phaseMob[ei][ip];
      }
    }

    // update pressure and check derivatives
    {
      // perturb pressure in each cell
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        real64 const dP = perturbParameter * ( pres[ei] + perturbParameter );
        dPres[ei] = dP;
      }

      // recompute component fractions
      solver.UpdateState( subRegion, targetIndex );

      // check values in each cell
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        SCOPED_TRACE( "Element " + std::to_string( ei ) );

        checkDerivative( phaseMob[ei].toSliceConst(),
                         phaseVolFracOrig[ei].toSliceConst(),
                         dPhaseMob_dPres[ei].toSliceConst(),
                         dPres[ei],
                         relTol,
                         "phaseVolFrac",
                         "Pres",
                         phases );
      }
    }

    // update component density and check derivatives
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      // reset the solver state to zero out variable updates (resetting the whole domain is overkill...)
      solver.ResetStateToBeginningOfStep( &domain );

      // perturb a single component density in each cell
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        real64 const dRho = perturbParameter * ( compDens[ei][jc] + perturbParameter );
        dCompDens[ei][jc] = dRho;
      }

      // recompute component fractions
      solver.UpdateState( subRegion, targetIndex );

      // check values in each cell
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        SCOPED_TRACE( "Element " + std::to_string( ei ) );

        auto dS_dRho = invertLayout( dPhaseMob_dCompDens[ei], NP, NC );
        string var = "compDens[" + components[jc] + "]";

        checkDerivative( phaseMob[ei].toSliceConst(),
                         phaseVolFracOrig[ei].toSliceConst(),
                         dS_dRho[jc].toSliceConst(),
                         dCompDens[ei][jc],
                         relTol,
                         "phaseMob",
                         var,
                         phases );
      }
    }
  } );
}

class CompositionalMultiphaseFlowTest : public ::testing::Test
{
public:

  CompositionalMultiphaseFlowTest()
    : problemManager( std::make_unique< ProblemManager >( "Problem", nullptr ) )
  {}

protected:

  void SetUp() override
  {
    setupProblemFromXML( *problemManager, xmlInput );
    solver = problemManager->GetPhysicsSolverManager().GetGroup< CompositionalMultiphaseFlow >( "compflow" );
  }

  std::unique_ptr< ProblemManager > problemManager;
  CompositionalMultiphaseFlow * solver;
};

TEST_F( CompositionalMultiphaseFlowTest, derivativeNumericalCheck_composition )
{
  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const tol = 1e-4;

  DomainPartition * domain = problemManager->getDomainPartition();

  testCompositionNumericalDerivatives( *solver, *domain, eps, tol );
}

TEST_F( CompositionalMultiphaseFlowTest, derivativeNumericalCheck_phaseVolumeFraction )
{
  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const tol = 5e-2; // 5% error margin

  DomainPartition * domain = problemManager->getDomainPartition();

  testPhaseVolumeFractionNumericalDerivatives( *solver, *domain, eps, tol );
}

TEST_F( CompositionalMultiphaseFlowTest, derivativeNumericalCheck_phaseMobility )
{
  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const tol = 5e-2; // 5% error margin

  DomainPartition * domain = problemManager->getDomainPartition();

  testPhaseMobilityNumericalDerivatives( *solver, *domain, eps, tol );
}

/*
 * Accumulation numerical test not passing due to some numerical catastrophic cancellation
 * happenning in the kernel for the particular set of initial conditions we're running.
 * The test should be re-enabled and fixed at some point.
 */
#if 0
TEST_F( CompositionalMultiphaseFlowTest, jacobianNumericalCheck_accumulation )
{
  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-1; // 10% error margin

  real64 const time = 0.0;
  real64 const dt = 1e4;

  DomainPartition * domain = problemManager->getDomainPartition();

  solver->ImplicitStepSetup( time,
                             dt,
                             domain,
                             solver->getDofManager(),
                             solver->getSystemMatrix(),
                             solver->getSystemRhs(),
                             solver->getSystemSolution() );

  testNumericalJacobian( *solver, *domain, eps, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    MeshLevel const & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
    solver->AssembleAccumulationTerms( time, dt, mesh, solver->getDofManager(), localMatrix, localRhs );
  } );
}
#endif

TEST_F( CompositionalMultiphaseFlowTest, jacobianNumericalCheck_flux )
{
  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const tol = 1e-1; // 10% error margin

  real64 const time = 0.0;
  real64 const dt = 1e4;

  DomainPartition * domain = problemManager->getDomainPartition();

  solver->ImplicitStepSetup( time,
                             dt,
                             domain,
                             solver->getDofManager(),
                             solver->getSystemMatrix(),
                             solver->getSystemRhs(),
                             solver->getSystemSolution() );

  testNumericalJacobian( *solver, *domain, eps, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    MeshLevel const & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
    NumericalMethodsManager const & numericalMethodManager = domain->getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = *fvManager.getFluxApproximation( solver->getDiscretization() );

    solver->AssembleFluxTerms( time, dt, mesh, fluxApprox, solver->getDofManager(), localMatrix, localRhs );
  } );
}


TEST_F( CompositionalMultiphaseFlowTest, jacobianNumericalCheck_volumeBalance )
{
  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const tol = 1e-1; // 10% error margin

  real64 const time = 0.0;
  real64 const dt = 1e4;

  DomainPartition * domain = problemManager->getDomainPartition();

  solver->ImplicitStepSetup( time,
                             dt,
                             domain,
                             solver->getDofManager(),
                             solver->getSystemMatrix(),
                             solver->getSystemRhs(),
                             solver->getSystemSolution() );


  testNumericalJacobian( *solver, *domain, eps, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    MeshLevel const & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
    solver->AssembleVolumeBalanceTerms( time, dt, mesh, solver->getDofManager(), localMatrix, localRhs );
  } );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  geosx::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
