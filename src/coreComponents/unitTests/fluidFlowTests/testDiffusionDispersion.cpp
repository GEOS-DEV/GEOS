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

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"

#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseFVMKernels.hpp"

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::constitutive;
using namespace geos::constitutive::multifluid;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

// Sphinx start after input XML

char const * xmlInput =
  R"xml(
  <Problem>
    <Solvers gravityVector="{ 0.0, 0.0, -9.81 }">
      <CompositionalMultiphaseFVM name="compflow"
                                   logLevel="0"
                                   discretization="fluidTPFA"
                                   targetRegions="{region}"
                                   temperature="297.15"
                                   useMass="1">

        <NonlinearSolverParameters newtonTol="1.0e-6"
                                   newtonMaxIter="2"/>
        <LinearSolverParameters solverType="gmres"
                                krylovTol="1.0e-10"/>
      </CompositionalMultiphaseFVM>
    </Solvers>
    <Mesh>
      <InternalMesh name="mesh"
                    elementTypes="{C3D8}"
                    xCoords="{0, 3}"
                    yCoords="{0, 1}"
                    zCoords="{0, 1}"
                    nx="{3}"
                    ny="{1}"
                    nz="{1}"
                    cellBlockNames="{cb1}"/>
    </Mesh>
    <NumericalMethods>
      <FiniteVolume>
        <TwoPointFluxApproximation name="fluidTPFA"/>
      </FiniteVolume>
    </NumericalMethods>
    <ElementRegions>
      <CellElementRegion name="region" cellBlocks="{cb1}" materialList="{fluid, rock, relperm}" />
    </ElementRegions>
    <Constitutive>
      <CompositionalMultiphaseFluid name="fluid"
                                    phaseNames="{oil, gas}"
                                    equationsOfState="{PR, PR}"
                                    componentNames="{N2, C10, C20, H2O}"
                                    componentCriticalPressure="{34e5, 25.3e5, 14.6e5, 220.5e5}"
                                    componentCriticalTemperature="{126.2, 622.0, 782.0, 647.0}"
                                    componentAcentricFactor="{0.04, 0.443, 0.816, 0.344}"
                                    componentMolarWeight="{28e-3, 134e-3, 275e-3, 18e-3}"
                                    componentVolumeShift="{0, 0, 0, 0}"
                                    componentBinaryCoeff="{ {0, 0, 0, 0},
                                                            {0, 0, 0, 0},
                                                            {0, 0, 0, 0},
                                                            {0, 0, 0, 0} }"/>
      <CompressibleSolidConstantPermeability name="rock"
          solidModelName="nullSolid"
          porosityModelName="rockPorosity"
          permeabilityModelName="rockPerm"/>
     <NullModel name="nullSolid"/>
     <PressurePorosity name="rockPorosity"
                       defaultReferencePorosity="0.05"
                       referencePressure = "0.0"
                       compressibility="1.0e-9"/>
      <BrooksCoreyRelativePermeability name="relperm"
                                       phaseNames="{oil, gas}"
                                       phaseMinVolumeFraction="{0.1, 0.15}"
                                       phaseRelPermExponent="{2.0, 2.0}"
                                       phaseRelPermMaxValue="{0.8, 0.9}"/>
    <ConstantPermeability name="rockPerm"
                          permeabilityComponents="{2.0e-16, 2.0e-16, 2.0e-16}"/>
    <ConstantDiffusion
      name="diffusion"
      phaseNames="{ oil, gas }"
      defaultPhaseDiffusivityMultipliers="{ 1, 1 }"
      diffusivityComponents="{ 1e-5, 1e-5, 1e-5 }"/>

    <LinearIsotropicDispersion
      name="dispersion"
      phaseNames="{ oil, gas }"
      longitudinalDispersivity="1e-5" />
    </Constitutive>
    <FieldSpecifications>
      <FieldSpecification name="initialPressure"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="ElementRegions/region/cb1"
                 fieldName="pressure"
                 functionName="initialPressureFunc"
                 scale="5e6"/>
      <FieldSpecification name="initialComposition_N2"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="ElementRegions/region/cb1"
                 fieldName="globalCompFraction"
                 component="0"
                 scale="0.099"/>
      <FieldSpecification name="initialComposition_C10"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="ElementRegions/region/cb1"
                 fieldName="globalCompFraction"
                 component="1"
                 scale="0.3"/>
      <FieldSpecification name="initialComposition_C20"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="ElementRegions/region/cb1"
                 fieldName="globalCompFraction"
                 component="2"
                 scale="0.6"/>
      <FieldSpecification name="initialComposition_H20"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="ElementRegions/region/cb1"
                 fieldName="globalCompFraction"
                 component="3"
                 scale="0.001"/>
    </FieldSpecifications>
    <Functions>
      <TableFunction name="initialPressureFunc"
                     inputVarNames="{elementCenter}"
                     coordinates="{0.0, 3.0}"
                     values="{1.0, 0.5}"/>
    </Functions>
  </Problem>
  )xml";

class CompositionalMultiphaseFlowTest : public ::testing::Test
{
public:

  CompositionalMultiphaseFlowTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:

  void SetUp() override
  {
    setupProblemFromXML( state.getProblemManager(), xmlInput );
    solver = &state.getProblemManager().getPhysicsSolverManager().getGroup< CompositionalMultiphaseFVM >( "compflow" );


    DomainPartition & domain = state.getProblemManager().getDomainPartition();

    //a bit ugly but good enough for now
    auto & meshBodies = domain.getMeshBodies();
    MeshBody & meshBody = meshBodies.getGroup< MeshBody >( "mesh" );
    const char * lvl = "Level0";
    MeshLevel & meshLevel = meshBody.getMeshLevel( lvl );
    ElementRegionManager & elemRegionManager = meshLevel.getElemManager();
    ElementRegionBase & elemRegion = elemRegionManager.getRegion( "region" );
    ElementSubRegionBase & elemSubRegion = elemRegion.getSubRegion( "cb1" );
    domain.getConstitutiveManager().hangConstitutiveRelation( "diffusion", &elemSubRegion, 1 );
    domain.getConstitutiveManager().hangConstitutiveRelation( "dispersion", &elemSubRegion, 1 );

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
};

real64 constexpr CompositionalMultiphaseFlowTest::time;
real64 constexpr CompositionalMultiphaseFlowTest::dt;
real64 constexpr CompositionalMultiphaseFlowTest::eps;



TEST_F( CompositionalMultiphaseFlowTest, sameDispersionTest )
{
  DomainPartition & domain = state.getProblemManager().getDomainPartition();
  FluxApproximationBase const & fluxApprox = domain.getNumericalMethodManager().getFiniteVolumeManager().getFluxApproximation( "fluidTPFA" );
  string const & elemDofKey = solver->getDofManager().getKey( CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString() );

  solver->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                       MeshLevel const & mesh,
                                                                       arrayView1d< string const > const & ) {
    fluxApprox.forAllStencils( mesh, [&] ( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();


      using kernelType = isothermalCompositionalMultiphaseFVMKernels::DiffusionDispersionFaceBasedAssemblyKernel< 4, 5, typename TYPEOFREF( stencil ) ::KernelWrapper >;
      typename kernelType::CompFlowAccessors compFlowAccessors( mesh.getElemManager(), solver->getName() );
      typename kernelType::MultiFluidAccessors multiFluidAccessors( mesh.getElemManager(), solver->getName() );
      typename kernelType::DiffusionAccessors diffusionAccessors( mesh.getElemManager(), solver->getName() );
      typename kernelType::DispersionAccessors dispersionAccessors( mesh.getElemManager(), solver->getName() );
      typename kernelType::PorosityAccessors porosityAccessors( mesh.getElemManager(), solver->getName() );
      ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
        mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
      dofNumberAccessor.setName( solver->getName() + "/accessors/" + elemDofKey );

      ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > globalCellDimAccessor =
        mesh.getElemManager().constructArrayViewAccessor< real64, 2 >( CellElementSubRegion::viewKeyStruct::globalCellDimString() );


      BitFlags< isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernelFlags > kernelFlags;

      CRSMatrix< real64, globalIndex > const & jacobian = solver->getLocalMatrix();
      array1d< real64 > residual( jacobian.numRows() );

      kernelType kernel(
        2,
        solver->getDofManager().rankOffset(),
        stencilWrapper,
        dofNumberAccessor,
        globalCellDimAccessor,
        compFlowAccessors,
        multiFluidAccessors,
        diffusionAccessors,
        dispersionAccessors,
        porosityAccessors,
        dt,
        jacobian.toViewConstSizes(),
        residual.toView(),
        kernelFlags
        );

      forAll< parallelDevicePolicy<> >( stencilWrapper.size(), [=] GEOS_HOST_DEVICE ( localIndex const iconn ) {
        typename kernelType::StackVariables stack( kernel.stencilSize( iconn ),
                                                   kernel.numPointsInFlux( iconn ));
        kernel.setup( iconn, stack );
        kernel.computeDiffusionFlux( iconn, stack );
        real64 diffusiveFlux[4], dispersiveFlux[4];
        real64 xp_diffusiveFlux[2][2][4] = { {{ 0.0032541334546537454, -0.00094419582575865383, -0.0018912528276302915, -1.4646839777785416e-06 },
          { -0.0032541334546537454, 0.00094419582575865383, 0.0018912528276302915, 1.4646839777785416e-06 }},
          {{ 0.0028595396719695945, -0.00079184397572525838, -0.0015861104962527795, -2.0324312559781912e-07},
            { -0.0028595396719695945, 0.00079184397572525838, 0.0015861104962527795, 2.0324312559781912e-07} }};

        real64 xp_dispersiveFlux[2][2][4] = {{{ 0.16156107701739061, -0.053505727994398372, -0.10720112367369174, 3.2864329379105197e-05 },
          { -0.15505281010808311, 0.051617336342881071, 0.10341861801843115, -3.5793697334662278e-05 }},
          {{ 0.16176374826400264, -0.053509594669204075, -0.10740925205848552, 0.00021669857431804405},
            {-0.15604466892006347, 0.051925906717753559, 0.10423703106597997, -0.00021710506056923972}} };                                           //assuming
                                                                                                                                                     // initial
                                                                                                                                                     // velocity
                                                                                                                                                     // of
                                                                                                                                                     // 1

        for( int k = 0; k < stack.numConnectedElems; ++k )
        {

          for( int ic = 0; ic < 4; ++ic )
          {
            integer const eqIndex0 = k * 4 + ic;
            diffusiveFlux[ic]  = stack.localFlux[eqIndex0];
            EXPECT_EQ( diffusiveFlux[ic], xp_diffusiveFlux[iconn][k][ic] );                                  //expect equality under these
                                                                                                             // conditions
          }
        }
        //do something with velocity or access saturation

        kernel.computeDispersionFlux( iconn, stack );
        for( int k = 0; k < stack.numConnectedElems; ++k )
        {
          for( int ic = 0; ic < 4; ++ic )
          {
            integer const eqIndex0 = k * 4 + ic;
            dispersiveFlux[ic] = stack.localFlux[eqIndex0] - diffusiveFlux[ic];
            EXPECT_EQ( dispersiveFlux[ic], xp_dispersiveFlux[iconn][k][ic] );                                  //expect equality under these
                                                                                                               // conditions
          }
        }


        kernel.complete( iconn, stack );

      } );
    } );
  } );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
