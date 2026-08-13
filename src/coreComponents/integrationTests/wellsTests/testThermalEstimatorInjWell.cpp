/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "integrationTests/fluidFlowTests/testCompFlowUtils.hpp"

#include "common/DataTypes.hpp"
#include "mainInterface/initialization.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/multiphysics/CompositionalMultiphaseReservoirAndWells.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWell.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/wells/kernels/CompositionalMultiphaseWellKernels.hpp"

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::constitutive;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

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
char const * co2flash = "FlashModel CO2Solubility   1e5 7.5e7 5e5 283.15 414.15 10 0\n";
char const * pvtLiquid = "DensityFun PhillipsBrineDensity 1e5 7.5e7 5e5 283.15 414.15 10 0\n"
                         "ViscosityFun PhillipsBrineViscosity 0\n"
                         "EnthalpyFun BrineEnthalpy 1e5 7.5e7 5e5 283.15 414.15 10 0\n";

char const * pvtGas = "DensityFun SpanWagnerCO2Density 1e5 7.5e7 5e5 283.15 414.15 10\n"
                      "ViscosityFun FenghourCO2Viscosity 1e5 7.5e7 5e5 283.15 414.15 10\n"
                      "EnthalpyFun CO2Enthalpy 1e5 7.5e7 5e5 283.15 414.15 10\n";
char const * xmlInput =
  R"xml(
<?xml version="1.0" ?>
<Problem>
  <Solvers>
    <CompositionalMultiphaseReservoir
      name="reservoirSystem"
      flowSolverName="compflow"
      wellSolverName="compositionalMultiphaseWell"
      logLevel="4"
      initialDt="1e4"
      targetRegions="{ region, injwell }">
      <NonlinearSolverParameters
        newtonTol="1.0e-8"
        lineSearchAction="None"
        newtonMaxIter="40">
      </NonlinearSolverParameters>
      <LinearSolverParameters
        directParallel="1">
      </LinearSolverParameters>
    </CompositionalMultiphaseReservoir>
    <CompositionalMultiphaseFVM
      name="compflow"
      logLevel="4"
      discretization="fluidTPFA"
      temperature="368.15"
      useMass="0"
      isThermal="1"
      initialDt="1e4"
      targetRelativePressureChangeInTimeStep="1"
      targetRelativeTemperatureChangeInTimeStep="1"
      targetPhaseVolFractionChangeInTimeStep="1"
      maxCompFractionChange="0.5"
      targetRegions="{ region }">
    </CompositionalMultiphaseFVM>
    <WellManager
      name="compositionalMultiphaseWell"
      targetRegions="{ injwell }"
      isThermal="1"
      logLevel="1"
      initialDt="1e4"
      useMass="0">
      <LinearSolverParameters
        directParallel="0">
      </LinearSolverParameters>
      <NonlinearSolverParameters
        newtonTol="1.0e-8"
        lineSearchAction="None"
        newtonMaxIter="20">
      </NonlinearSolverParameters>
      <CompositionalMultiphaseWell
        name="WC_CO2_INJ"
              type="injector"
              writeCSV="1"
        useSurfaceConditions="1"
        surfacePressure="101325"
        surfaceTemperature="288.706"
        control="totalVolRate"
        enableCrossflow="0"
        estimateWellSolution="0">
        <MaximumBHPConstraint
          name="maxbhp"
          targetBHP="1e8"
          referenceElevation="-2738"/>
        <InjectionVolumeRateConstraint
          name="maxvolrateinj"
          volumeRate="292"
          injectionStream="{ 0.999, 0.001 }"
          injectionTemperature="306"/>
              </CompositionalMultiphaseWell>
      </WellManager>

  </Solvers>
  <Mesh>
    <InternalMesh
      name="mesh1"
      elementTypes="{ C3D8 }"
      xCoords="{ 0, 100 }"
      yCoords="{ 0, 100 }"
      zCoords="{ 0, 1 }"
      nx="{ 2 }"
      ny="{ 2 }"
      nz="{ 1 }"
      cellBlockNames="{ cb }">
      <InternalWell
        name="inj11"
        wellRegionName="injwell"
        wellControlsName="WC_CO2_INJ"
        logLevel="1"
        polylineNodeCoords="{{5.0,5.0,1.01},{5.0,5.0,0.01}} "
        polylineSegmentConn="{{0,1}} "
        radius="0.1"
        numElementsPerSegment="2">
        <Perforation
          name="inj1_perf1"
          distanceFromHead="0.75">
        </Perforation>
      </InternalWell>
    </InternalMesh>
  </Mesh>
  <Events
    maxTime="1.5e5">
    <PeriodicEvent
      name="outputs"
      timeFrequency="2.5e4"
      target="/Outputs/vtkOutput">
    </PeriodicEvent>
    <PeriodicEvent
      name="solverApplications"
      maxEventDt="2.5e4"
      target="/Solvers/reservoirSystem">
    </PeriodicEvent>
    <PeriodicEvent
      name="restarts"
      timeFrequency="7.5e5"
      target="/Outputs/sidreRestart">
    </PeriodicEvent>
  </Events>
   <NumericalMethods>
    <FiniteVolume>
      <TwoPointFluxApproximation
        name="fluidTPFA">
      </TwoPointFluxApproximation>
    </FiniteVolume>
  </NumericalMethods>
  <ElementRegions>
    <CellElementRegion
      name="region"
      cellBlocks="{ cb }"
      materialList="{ fluid, rock, relperm, thermalCond, diffusion }">
    </CellElementRegion>
    <WellElementRegion
      name="injwell"
      materialList="{ fluid, relperm }">
    </WellElementRegion>
  </ElementRegions> 
  <Constitutive>
    <CompressibleSolidConstantPermeability
      name="rock"
      solidModelName="nullSolid"
      porosityModelName="rockPorosity"
      permeabilityModelName="rockPerm"
      solidInternalEnergyModelName="rockInternalEnergy">
    </CompressibleSolidConstantPermeability>
    <NullModel
      name="nullSolid">
    </NullModel>
    <PressurePorosity
      name="rockPorosity"
      defaultReferencePorosity="0.2"
      referencePressure="0.0"
      compressibility="1.0e-9">
    </PressurePorosity>
    <SolidInternalEnergy
      name="rockInternalEnergy"
      referenceVolumetricHeatCapacity="1.95e6"
      referenceTemperature="368.15"
      referenceInternalEnergy="0">
    </SolidInternalEnergy>
    <ConstantPermeability
      name="rockPerm"
      permeabilityComponents="{ 1.0e-13, 1.0e-13, 1.0e-13 }">
    </ConstantPermeability>
    <CO2BrinePhillipsThermalFluid
      name="fluid"
      logLevel="1"
      phaseNames="{ gas, water }"
      componentNames="{ co2, water }"
      componentMolarWeight="{ 44e-3, 18e-3 }"
      phasePVTParaFiles="{  pvtgas.txt,  pvtliquid.txt }"
      flashModelParaFile="co2flash.txt">
    </CO2BrinePhillipsThermalFluid>
    <BrooksCoreyRelativePermeability
      name="relperm"
      phaseNames="{ gas, water }"
      phaseMinVolumeFraction="{ 0.0, 0.0 }"
      phaseRelPermExponent="{ 1.5, 1.5 }"
      phaseRelPermMaxValue="{ 0.9, 0.9 }">
    </BrooksCoreyRelativePermeability>
    <MultiPhaseConstantThermalConductivity
      name="thermalCond"
      phaseNames="{ gas, water }"
      thermalConductivityComponents="{ 0.6, 0.6, 0.6 }">
    </MultiPhaseConstantThermalConductivity>
    <ConstantDiffusion
      name="diffusion"
      phaseNames="{ gas, water }"
      defaultPhaseDiffusivityMultipliers="{ 20, 1 }"
      diffusivityComponents="{ 1e-9, 1e-9, 1e-9 }">
    </ConstantDiffusion>
  </Constitutive>
  <FieldSpecifications>
    <FieldSpecification
      name="initialPressure"
      initialCondition="1"
      setNames="{ all }"
      objectPath="ElementRegions/region/cb"
      fieldName="pressure"
      scale="9e6">
    </FieldSpecification>
    <FieldSpecification
      name="initialTemperature"
      initialCondition="1"
      setNames="{ all }"
      objectPath="ElementRegions/region/cb"
      fieldName="temperature"
      scale="368.15">
    </FieldSpecification>
    <FieldSpecification
      name="initialComposition_co2"
      initialCondition="1"
      setNames="{ all }"
      objectPath="ElementRegions/region/cb"
      fieldName="globalCompFraction"
      component="0"
      scale="0.005">
    </FieldSpecification>
    <FieldSpecification
      name="initialComposition_water"
      initialCondition="1"
      setNames="{ all }"
      objectPath="ElementRegions/region/cb"
      fieldName="globalCompFraction"
      component="1"
      scale="0.995">
    </FieldSpecification>
  </FieldSpecifications>
  <Outputs>
    <VTK
      name="vtkOutput">
    </VTK>
    <Restart
      name="sidreRestart">
    </Restart>
  </Outputs>
</Problem>
)xml";

template< typename T, typename COL_INDEX >
void printCompareLocalMatrices( CRSMatrixView< T const, COL_INDEX const > const & matrix1,
                                CRSMatrixView< T const, COL_INDEX const > const & matrix2, std::string const & testName )
{
  matrix1.move( hostMemorySpace, false );
  matrix2.move( hostMemorySpace, false );

  std::ofstream omat1( testName+".csv" );


  std::vector< std::vector< double > > fmat1( matrix1.numRows(), std::vector< double >( matrix1.numRows(), 0.0 ));
  std::vector< std::vector< double > > fmat2( matrix2.numRows(), std::vector< double >( matrix2.numRows(), 0.0 ));

  for( localIndex i = 0; i < matrix1.numRows(); ++i )
  {
    arraySlice1d< globalIndex const > indices1 = matrix1.getColumns( i );
    arraySlice1d< globalIndex const > indices2 = matrix2.getColumns( i );
    arraySlice1d< double const > values1  = matrix1.getEntries( i );
    arraySlice1d< double const > values2  = matrix2.getEntries( i );
    for( integer j=0; j<indices1.size(); j++ )
    {
      fmat1[i][indices1[j]] = values1[j];
    }
    for( integer j=0; j<indices2.size(); j++ )
    {
      fmat2[i][indices2[j]] = values2[j];
    }
  }
  for( integer i=0; i<matrix1.numRows(); i++ )
  {
    for( integer j=0; j<matrix1.numColumns(); j++ )
    {
      omat1 << "," << fmat1[i][j];
    }
    omat1 << "\n";
  }
  omat1 << "\n";
  for( integer i=0; i<matrix2.numRows(); i++ )
  {
    for( integer j=0; j<matrix2.numColumns(); j++ )
    {
      omat1 << "," << fmat2[i][j];

    }
    omat1 << "\n";

  }
  omat1.close();
}


void printResiduals( array1d< real64 > const & rsd1,
                     array1d< real64 > const & rsd2, std::string const & testName )
{
  std::ofstream omat1( testName+".csv" );

  for( integer i=0; i<rsd1.size(); i++ )
  {
    omat1 << i <<"," << rsd1[i] << "," << rsd2[i] << "\n";
  }

  omat1.close();
}
template< typename LAMBDA >
void testWellEstimatorNumericalJacobian( CompositionalMultiphaseReservoirAndWells< CompositionalMultiphaseBase > & solver,
                                         DomainPartition & domain,
                                         real64 const perturbParameter,
                                         real64 const time_n,
                                         real64 const relTol, bool diag_check, std::string const & testName,
                                         LAMBDA && assembleFunction )
{
  GEOS_UNUSED_VAR( testName );
  WellManager & wellSolver = *solver.wellSolver();

  wellSolver.forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                                           MeshLevel & mesh,
                                                                           string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                        [&]( localIndex const,
                                                                             WellElementSubRegion & subRegion )
    {
      WellControls & wellControls = wellSolver.getWellControls( subRegion );
      CompositionalMultiphaseWell * compWell = dynamic_cast< CompositionalMultiphaseWell * >(&wellControls);
      compWell->setWellState( 1 );
      compWell->initializeWell( domain, domain.getMeshBodies(), meshBodyName, mesh, subRegion, time_n );
    } );
  } );

  CompositionalMultiphaseFVM & flowSolver = dynamicCast< CompositionalMultiphaseFVM & >( *solver.reservoirSolver() );

  localIndex const NC = flowSolver.numFluidComponents();

  CRSMatrix< real64, globalIndex > const & jacobian = solver.getLocalMatrix();
  array1d< real64 > residual( jacobian.numRows() );
  DofManager const & dofManager = solver.getDofManager();

  // assemble the analytical residual
  solver.resetStateToBeginningOfStep( domain );
  string const resDofKey  = dofManager.getKey( wellSolver.resElementDofName() );
  string const wellDofKey = dofManager.getKey( wellSolver.wellElementDofName() );
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
  ////////////////////////////////////////////////
  // Step 1) Compute the terms in J_RR and J_WR //
  ////////////////////////////////////////////////
#if 1
  domain.forMeshBodies( [&] ( MeshBody & meshBody )
  {
    bool processMesh = true;
    meshBody.forMeshLevels( [&] ( MeshLevel & mesh )
    {
      if( !processMesh )
        return;
      processMesh = false;
      ElementRegionManager & elemManager = mesh.getElemManager();
      for( localIndex er = 0; er < elemManager.numRegions(); ++er )
      {
        ElementRegionBase & elemRegion = elemManager.getRegion( er );
        elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const, CellElementSubRegion & subRegion )
        {
          // get the degrees of freedom and ghosting information
          arrayView1d< globalIndex const > const & dofNumber =
            subRegion.getReference< array1d< globalIndex > >( resDofKey );

          // get the primary variables on the reservoir elements
          arrayView1d< real64 > const & pres =
            subRegion.getField< fields::flow::pressure >();
          pres.move( hostMemorySpace, false );

          arrayView1d< real64 > const & temperature =
            subRegion.getField< fields::flow::temperature >();
          temperature.move( hostMemorySpace, false );

          arrayView2d< real64, compflow::USD_COMP > const & compDens =
            subRegion.getField< fields::flow::globalCompDensity >();
          compDens.move( hostMemorySpace, false );

          // a) compute all the derivatives wrt to the pressure in RESERVOIR elem ei
          for( localIndex ei = 0; ei < subRegion.size(); ++ei )
          {
#if 1
            {
              solver.resetStateToBeginningOfStep( domain );

              // here is the perturbation in the pressure of the element
              real64 const dP = perturbParameter * (pres[ei] + perturbParameter);
              pres.move( hostMemorySpace, true );
              pres[ei] += dP;

              // after perturbing, update the pressure-dependent quantities in the reservoir
              flowSolver.forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                                       MeshLevel & mesh2,
                                                                                       string_array const & regionNames2 )
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
#endif
#if 1
            {
              solver.resetStateToBeginningOfStep( domain );

              // here is the perturbation in the temperature of the element
              real64 const dT = 1;//perturbParameter * (temperature[ei] + perturbParameter);
              temperature.move( hostMemorySpace, true );
              temperature[ei] += dT;

              // after perturbing, update the temperature-dependent quantities in the reservoir
              flowSolver.forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                                       MeshLevel & mesh2,
                                                                                       string_array const & regionNames2 )
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
                                     dofNumber[ei]+NC+1,
                                     dT,
                                     jacobianFD.toViewConstSizes() );
            }
#endif
#if 1
            real64 totalDensity = 0.0;
            for( localIndex ic = 0; ic < NC; ++ic )
            {
              totalDensity += compDens[ei][ic];
            }

            for( localIndex jc = 0; jc < NC; ++jc )
            {
              solver.resetStateToBeginningOfStep( domain );

              real64 const dRho = perturbParameter * totalDensity;
              compDens.move( hostMemorySpace, true );
              compDens[ei][jc] += dRho;

              flowSolver.forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                                       MeshLevel & mesh2,
                                                                                       string_array const & regionNames2 )
              {
                mesh2.getElemManager().forElementSubRegions( regionNames2,
                                                             [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion2 )
                {
                  flowSolver.updateFluidState( subRegion2 );
                } );
              } );

              residual.zero();
              jacobian.zero();
              assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

              fillNumericalJacobian( residual.toViewConst(),
                                     residualOrig.toViewConst(),
                                     dofNumber[ei] + jc + 1,
                                     dRho,
                                     jacobianFD.toViewConstSizes() );
            }
#endif
          }

          return;
        } );
      }
      return;
    } );
  } );

#endif
  // at this point we start assembling the finite-difference block by block


  /////////////////////////////////////////////////
  // Step 2) Compute the terms in J_RW and J_WW //
  /////////////////////////////////////////////////

  // loop over the wells
  if( 1 )
    wellSolver.forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                             MeshLevel & mesh,
                                                                             string_array const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                          [&]( localIndex const,
                                                                               WellElementSubRegion & subRegion )
      {
        // get the degrees of freedom, ghosting info and next well elem index
        arrayView1d< globalIndex const > const & wellElemDofNumber =
          subRegion.getReference< array1d< globalIndex > >( wellDofKey );

        // get the primary variables on the well elements
        arrayView1d< real64 > const & wellElemPressure =
          subRegion.getField< fields::well::pressure >();
        wellElemPressure.move( hostMemorySpace, false );

        arrayView1d< real64 > const & wellElemTemperature =
          subRegion.getField< fields::well::temperature >();
        wellElemTemperature.move( hostMemorySpace, false );

        arrayView2d< real64, compflow::USD_COMP > const & wellElemCompDens =
          subRegion.getField< fields::well::globalCompDensity >();
        wellElemCompDens.move( hostMemorySpace, false );

        arrayView1d< real64 > const & connRate =
          subRegion.getField< fields::well::connectionRate >();
        connRate.move( hostMemorySpace, false );

        // a) compute all the derivatives wrt to the pressure in WELL elem iwelem
        for( localIndex iwelem = 0; iwelem < subRegion.size(); ++iwelem )
        {

#if 1
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

            fillNumericalJacobian( residual.toViewConst(),
                                   residualOrig.toViewConst(),
                                   wellElemDofNumber[iwelem] + compositionalMultiphaseWellKernels::ColOffset::DPRES,
                                   dP,
                                   jacobianFD.toViewConstSizes() );
          }
#endif
#if 1
          real64 wellElemTotalDensity = 0.0;
          for( localIndex ic = 0; ic < NC; ++ic )
          {
            wellElemTotalDensity += wellElemCompDens[iwelem][ic];
          }
          for( localIndex jc = 0; jc < NC; ++jc )
          {
            solver.resetStateToBeginningOfStep( domain );

            real64 const dRho = perturbParameter * wellElemTotalDensity;
            wellElemCompDens.move( hostMemorySpace, true );
            wellElemCompDens[iwelem][jc] += dRho;

            wellSolver.updateState( domain );

            residual.zero();
            jacobian.zero();
            assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

            fillNumericalJacobian( residual.toViewConst(),
                                   residualOrig.toViewConst(),
                                   wellElemDofNumber[iwelem] + compositionalMultiphaseWellKernels::ColOffset::DCOMP + jc,
                                   dRho,
                                   jacobianFD.toViewConstSizes() );
          }
          {
            solver.resetStateToBeginningOfStep( domain );
            residual.zero();
            jacobian.zero();
            if( diag_check || iwelem > 0 )
            {
              // here is the perturbation in the temperature of the well element
              real64 const dT = perturbParameter * ( wellElemTemperature[iwelem] + perturbParameter );
              wellElemTemperature.move( hostMemorySpace, true );
              wellElemTemperature[iwelem] += dT;

              // after perturbing, update the pressure-dependent quantities in the well
              wellSolver.updateState( domain );

              assembleFunction( jacobian.toViewConstSizes(), residual.toView() );
              fillNumericalJacobian( residual.toViewConst(),
                                     residualOrig.toViewConst(),
                                     wellElemDofNumber[iwelem] + compositionalMultiphaseWellKernels::ColOffset::DCOMP + NC+1,
                                     dT,
                                     jacobianFD.toViewConstSizes() );
              if( iwelem == 1 )
              {
                real64 dRdX = 0.0;
                localIndex rowIndex = wellElemDofNumber[0] + compositionalMultiphaseWellKernels::ColOffset::DCOMP + NC+1;;
                for( integer ider=0; ider< 3; ider++ )
                {
                  globalIndex colIndex = wellElemDofNumber[0]+ ider;
                  setNumericalJacobianValue( rowIndex, colIndex, dRdX, jacobianFD.toViewConstSizes() );
                }
                globalIndex colIndex = wellElemDofNumber[1]+3;
                setNumericalJacobianValue( rowIndex, colIndex, dRdX, jacobianFD.toViewConstSizes() );

              }
            }
            else
            {
              localIndex rowIndex = wellElemDofNumber[iwelem] + compositionalMultiphaseWellKernels::ColOffset::DCOMP + NC+1;;
              globalIndex colIndex = wellElemDofNumber[iwelem] + compositionalMultiphaseWellKernels::ColOffset::DCOMP + NC+1;;
              setNumericalJacobianValue( rowIndex, colIndex, 1.0, jacobianFD.toViewConstSizes() );
            }
          }
#endif
        }

#if 1
        // b) compute all the derivatives wrt to the connection in WELL elem
        // iwelem
        for( localIndex iwelem = 0; iwelem < subRegion.size(); ++iwelem )
        {
          {
            solver.resetStateToBeginningOfStep( domain );

            // here is the perturbation in the rate of the well element
            real64 const dRate = perturbParameter * ( connRate[iwelem] + perturbParameter );
            connRate.move( hostMemorySpace, true );
            connRate[iwelem] += dRate;

            wellSolver.updateState( domain );

            residual.zero();
            jacobian.zero();
            assembleFunction( jacobian.toViewConstSizes(), residual.toView() );
            fillNumericalJacobian( residual.toViewConst(),
                                   residualOrig.toViewConst(),
                                   wellElemDofNumber[iwelem] + compositionalMultiphaseWellKernels::ColOffset::DCOMP + NC,
                                   dRate,
                                   jacobianFD.toViewConstSizes() );
          }
        }
#endif
      } );
    } );

  // assemble the analytical jacobian
  solver.resetStateToBeginningOfStep( domain );

  residual.zero();
  jacobian.zero();
  assembleFunction( jacobian.toViewConstSizes(), residual.toView() );
  //printCompareLocalMatrices( jacobian.toViewConst(), jacobianFD.toViewConst(), testName );
  compareLocalMatrices( jacobian.toViewConst(), jacobianFD.toViewConst(), relTol );
}

class CompositionalMultiphaseReservoirSolverTest : public ::testing::Test
{
public:

  CompositionalMultiphaseReservoirSolverTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:

  void SetUp() override
  {
    setupProblemFromXML( state.getProblemManager(), xmlInput );
    solver = &state.getProblemManager().getPhysicsSolverManager().getGroup< CompositionalMultiphaseReservoirAndWells< CompositionalMultiphaseBase > >( "reservoirSystem" );

    DomainPartition & domain = state.getProblemManager().getDomainPartition();

    solver->setupSystem( domain,
                         solver->getDofManager(),
                         solver->getLocalMatrix(),
                         solver->getSystemRhs(),
                         solver->getSystemSolution() );

    solver->implicitStepSetup( time, dt, domain );
    WellManager & wellSolver = *solver->wellSolver();
    wellSolver.forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const & meshBodyName,
                                                                            MeshLevel & meshLevel,
                                                                            string_array const & regionNames )
    {
      ElementRegionManager & elementRegionManager = meshLevel.getElemManager();
      elementRegionManager.forElementRegions< WellElementRegion >( regionNames,
                                                                   [&]( localIndex const,
                                                                        WellElementRegion & region )
      {
        WellElementSubRegion & subRegion = region.getGroup( ElementRegionBase::viewKeyStruct::elementSubRegions() )
                                             .getGroup< WellElementSubRegion >( region.getSubRegionName() );
        WellControls & wellControls = wellSolver.getWellControls( subRegion );
        wellControls.initializeWell( domain, domain.getMeshBodies(), meshBodyName, meshLevel, subRegion, time );
      } );
    } );
  }

  static real64 constexpr time = 0.0;
  static real64 constexpr dt = 1e4;
  static real64 constexpr eps = std::numeric_limits< real64 >::epsilon();

  GeosxState state;
  CompositionalMultiphaseReservoirAndWells< CompositionalMultiphaseBase > * solver;
};

real64 constexpr CompositionalMultiphaseReservoirSolverTest::time;
real64 constexpr CompositionalMultiphaseReservoirSolverTest::dt;
real64 constexpr CompositionalMultiphaseReservoirSolverTest::eps;



#if 0
// There are a few terms that cause this test to fail, these are dt terms of cell connected to cell with boundary condtion.
// The could be zeroed out in the jacobian assembly and then this test should pass.
// Otherwise the test is good, uncomment out printCompareLocalMatrices and look at FD and computed derivatives
TEST_F( CompositionalMultiphaseReservoirSolverTest, jacobianNumericalCheck_System )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain =  state.getProblemManager().getDomainPartition();

  testWellEstimatorNumericalJacobian( *solver, domain, perturb, time, tol, false, "Check_System",
                                      [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                            arrayView1d< real64 > const & localRhs )
  {

    DofManager const & dofManager = solver->getDofManager();
    WellManager * wellSolver = solver->wellSolver();
    wellSolver->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                              MeshLevel & mesh,
                                                                              string_array const & regionNames )
    {
      ElementRegionManager & elemManager = mesh.getElemManager();
      elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                [&]( localIndex const,
                                                                     WellElementSubRegion & subRegion )
      {
        // call assemble to fill the matrix and the rhs
        WellControls & wellControls = wellSolver->getWellControls( subRegion );
        CompositionalMultiphaseWell * compWell = dynamic_cast< CompositionalMultiphaseWell * >(&wellControls);
        compWell->assembleSystem( time,
                                  dt,
                                  0,
                                  elemManager,
                                  subRegion,
                                  dofManager,
                                  localMatrix,
                                  localRhs );

      } );
    } );

    solver->assembleCouplingTerms( time, dt, domain, solver->getDofManager(), localMatrix, localRhs );

  } );
}
#endif
TEST_F( CompositionalMultiphaseReservoirSolverTest, jacobianNumericalCheck_Accum )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain =  state.getProblemManager().getDomainPartition();

  testWellEstimatorNumericalJacobian( *solver, domain, perturb, time, tol, false, "Check_Accum",
                                      [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                            arrayView1d< real64 > const & localRhs )
  {

    //DofManager const & dofManager = solver->wellSolver()->getDofManager();
    WellManager * wellSolver = solver->wellSolver();
    wellSolver->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                              MeshLevel & mesh,
                                                                              string_array const & regionNames )
    {
      ElementRegionManager & elemManager = mesh.getElemManager();
      elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                [&]( localIndex const,
                                                                     WellElementSubRegion & subRegion )
      {
        // call assemble to fill the matrix and the rhs
        WellControls & wellControls = wellSolver->getWellControls( subRegion );
        wellControls.assembleWellAccumulationTerms( time, dt, subRegion, solver->getDofManager(), localMatrix, localRhs );
      } );
    } );
  } );
}
#if 1
TEST_F( CompositionalMultiphaseReservoirSolverTest, jacobianNumericalCheck_PressureRelation )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain =  state.getProblemManager().getDomainPartition();

  testWellEstimatorNumericalJacobian( *solver, domain, perturb, time, tol, true, "Check_PressureRelation",
                                      [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                            arrayView1d< real64 > const & localRhs )
  {

    //DofManager const & dofManager = solver->wellSolver()->getDofManager();
    DofManager const & dofManager = solver->getDofManager();
    WellManager * wellSolver = solver->wellSolver();
    wellSolver->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                              MeshLevel & mesh,
                                                                              string_array const & regionNames )
    {
      ElementRegionManager & elemManager = mesh.getElemManager();
      elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                [&]( localIndex const,
                                                                     WellElementSubRegion & subRegion )
      {
        // call assemble to fill the matrix and the rhs
        WellControls & wellControls = wellSolver->getWellControls( subRegion );
        //wellControls.assembleWellConstraintTerms( time, dt, subRegion, dofManager, localMatrix.toViewConstSizes(), localRhs );
        wellControls.assembleWellPressureRelations( time, dt, subRegion, dofManager, localMatrix.toViewConstSizes(), localRhs );

      } );
    } );
  } );
}
#endif
#if 1
TEST_F( CompositionalMultiphaseReservoirSolverTest, jacobianNumericalCheck_Flux )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain =  state.getProblemManager().getDomainPartition();

  testWellEstimatorNumericalJacobian( *solver, domain, perturb, time, tol, true, "Check_Flux",
                                      [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                            arrayView1d< real64 > const & localRhs )
  {

    //DofManager const & dofManager = solver->wellSolver()->getDofManager();
    WellManager * wellSolver = solver->wellSolver();
    wellSolver->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                              MeshLevel & mesh,
                                                                              string_array const & regionNames )
    {
      ElementRegionManager & elemManager = mesh.getElemManager();
      elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                [&]( localIndex const,
                                                                     WellElementSubRegion & subRegion )
      {
        // call assemble to fill the matrix and the rhs
        WellControls & wellControls = wellSolver->getWellControls( subRegion );
        wellControls.assembleWellFluxTerms( time, dt, subRegion, solver->getDofManager(), localMatrix, localRhs );
      } );
    } );
  } );
}
#endif
int main( int argc, char * * argv )
{
  writeTableToFile( "co2flash.txt", co2flash );
  writeTableToFile( "pvtliquid.txt", pvtLiquid );
  writeTableToFile( "pvtgas.txt", pvtGas );
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  removeFile( "co2flash.txt" );
  removeFile( "pvtliquid.txt" );
  removeFile( "pvtgas.txt" );

  return result;
}
