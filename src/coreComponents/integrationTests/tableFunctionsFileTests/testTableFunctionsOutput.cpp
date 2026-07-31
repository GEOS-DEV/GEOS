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

#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "codingUtilities/Parsing.hpp"
#include "physicsSolvers/PhysicsSolverBase.hpp"
#include "physicsSolvers/SolverStatistics.hpp"

#include <gtest/gtest.h>
#include <filesystem>


using namespace geos;

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
char const * co2flash = "FlashModel CO2Solubility  1e6 7.5e7 5e5 299.15 369.15 10 0";
char const * pvtLiquid = "DensityFun PhillipsBrineDensity 1e6 7.5e7 5e5 299.15 369.15 10 0\n"
                         "ViscosityFun PhillipsBrineViscosity 0\n"
                         "EnthalpyFun BrineEnthalpy 1e6 7.5e7 5e5 299.15 369.15 10 0\n";

char const * pvtGas = "DensityFun SpanWagnerCO2Density 1e6 7.5e7 5e5 299.15 369.15 10\n"
                      "ViscosityFun FenghourCO2Viscosity 1e6 7.5e7 5e5 299.15 369.15 10\n"
                      "EnthalpyFun CO2Enthalpy 1e6 7.5e7 5e5 299.15 369.15 10\n";
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
      initialDt="1e3"
      targetRegions="{ region, injwell }">
      <NonlinearSolverParameters
         logLevel="1"
        newtonTol="1.0e-5"
        lineSearchAction="None"
        newtonMaxIter="40"/>
      <LinearSolverParameters
        logLevel="4"
        directParallel="0"/>
    </CompositionalMultiphaseReservoir>

    <CompositionalMultiphaseFVM
      name="compflow"
      logLevel="4"
      discretization="fluidTPFA"
      temperature="368.15"
      useMass="1"
      initialDt="1e3"
      targetRelativePressureChangeInTimeStep="1"
      targetRelativeTemperatureChangeInTimeStep="1"
      targetPhaseVolFractionChangeInTimeStep="1"      
      maxCompFractionChange="0.5"
      targetRegions="{ region }">
    </CompositionalMultiphaseFVM>

    <WellManager
      name="compositionalMultiphaseWell"
      targetRegions="{ injwell }"
      logLevel="1"
      useMass="1">
    <CompositionalMultiphaseWell
      name="WC_CO2_INJ"
      logLevel="2"
      type="injector"
      enableCrossflow="0"
      useSurfaceConditions="1"
      control="totalVolRate"
      surfacePressure="1.45e7"
      surfaceTemperature="300.15">
      <MaximumBHPConstraint
        name="maxbhp"
        targetBHP="1.45e7"
        referenceElevation="-0.01"/>
      <InjectionVolumeRateConstraint
        name="maxvolrateinj"
        volumeRate="0.001"
        injectionStream="{ 1.0, 0.0 }"
        injectionTemperature="300.15"/>
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
        name="inj1"
        wellRegionName="injwell"
        wellControlsName="WC_CO2_INJ"
        logLevel="1"
        polylineNodeCoords="{ { 5.0, 5.0, 1.01 },
                              { 5.0, 5.0, -0.01 } }"
        polylineSegmentConn="{ { 0, 1 } }"
        radius="0.1"
        numElementsPerSegment="1">
        <Perforation
          name="injector1_perf1"
          distanceFromHead="0.75"/>
      </InternalWell>
    </InternalMesh>

  </Mesh>

  <Geometry>
    <Box
      name="sink"
      xMin="{ 49.99, 49.99, -0.01 }"
      xMax="{ 101.01, 101.01, 1.01 }"/>
  </Geometry>

  <Events
    maxTime="1.5e5">
    <PeriodicEvent
      name="outputs"
      timeFrequency="2.5e4"
      target="/Outputs/vtkOutput"/>

    <PeriodicEvent
      name="solverApplications"
      maxEventDt="2.5e4"
      target="/Solvers/coupledFlowAndWells"/>

    <PeriodicEvent
      name="restarts"
      timeFrequency="7.5e5"
      target="/Outputs/sidreRestart"/>
  </Events>

  <NumericalMethods>
    <FiniteVolume>
      <TwoPointFluxApproximation
        name="fluidTPFA"/>
    </FiniteVolume>
  </NumericalMethods>

  <ElementRegions>
    <CellElementRegion
      name="region"
      cellBlocks="{ cb }"
      materialList="{ fluid, rock, relperm }"/>
    <WellElementRegion
      name="injwell"
      materialList="{ fluid, relperm }"/>
  </ElementRegions>

  <Constitutive>

    <CompressibleSolidConstantPermeability
      name="rock"
      solidModelName="nullSolid"
      porosityModelName="rockPorosity"
      permeabilityModelName="rockPerm"/>
    <NullModel
      name="nullSolid"/>
    <PressurePorosity
      name="rockPorosity"
      defaultReferencePorosity="0.2"
      referencePressure="0.0"
      compressibility="1.0e-9"/>
    <ConstantPermeability
      name="rockPerm"
      permeabilityComponents="{ 1.0e-13, 1.0e-13, 1.0e-13 }"/>

    <CO2BrinePhillipsFluid
      name="fluid"
      logLevel="1"
      phaseNames="{ gas, water }"
      componentNames="{ co2, water }"
      componentMolarWeight="{ 44e-3, 18e-3 }"
      phasePVTParaFiles="{ pvtgas.txt, pvtliquid.txt }"
      flashModelParaFile="co2flash.txt"/>

    <BrooksCoreyRelativePermeability
      name="relperm"
      phaseNames="{ gas, water }"
      phaseMinVolumeFraction="{ 0.0, 0.0 }"
      phaseRelPermExponent="{ 1.5, 1.5 }"
      phaseRelPermMaxValue="{ 0.9, 0.9 }"/>
    
  </Constitutive>

  <FieldSpecifications>

    <FieldSpecification
      name="initialPressure"
      initialCondition="1"
      setNames="{ all }"
      objectPath="ElementRegions/region/cb"
      fieldName="pressure"
      scale="9e6"/>
    <FieldSpecification
      name="initialTemperature"
      initialCondition="1"
      setNames="{ all }"
      objectPath="ElementRegions/region/cb"
      fieldName="temperature"
      scale="368.15"/>
    <FieldSpecification
      name="initialComposition_co2"
      initialCondition="1"
      setNames="{ all }"
      objectPath="ElementRegions/region/cb"
      fieldName="globalCompFraction"
      component="0"
      scale="0.005"/>
    <FieldSpecification
      name="initialComposition_water"
      initialCondition="1"
      setNames="{ all }"
      objectPath="ElementRegions/region/cb"
      fieldName="globalCompFraction"
      component="1"
      scale="0.995"/>

    <FieldSpecification
      name="sinkPressure"
      setNames="{ sink }"       
      objectPath="ElementRegions/region/cb"
      fieldName="pressure"
      scale="7e6"/>
     <FieldSpecification
      name="sinkTermComposition_co2"
      setNames="{ sink }"
      objectPath="ElementRegions/region/cb"
      fieldName="globalCompFraction"
      component="0"
      scale="0.005"/>
    <FieldSpecification
      name="sinkTermComposition_water"
      setNames="{ sink }"
      objectPath="ElementRegions/region/cb"
      fieldName="globalCompFraction"
      component="1"
      scale="0.995"/>

  </FieldSpecifications>

  <Outputs>
    <VTK
      name="vtkOutput"/>

    <Restart
      name="sidreRestart"/>
  </Outputs>
</Problem>
  )xml";

CommandLineOptions g_commandLineOptions;

TEST( testTableFunctionsOutput, testOutputFiles )
{

  std::filesystem::path f1{"fluid_phaseModel1_PhillipsBrineDensity_table.csv"};
  std::filesystem::path f2{"fluid_phaseModel2_SpanWagnerCO2Density_table.csv"};
  std::filesystem::path f3{"fluid_phaseModel2_FenghourCO2Viscosity_table.csv"};
  std::filesystem::path f4{ "fluid_CO2Solubility_co2Dissolution_table.csv" };
  std::filesystem::path f5{ "fluid_phaseModel1_PhillipsBrineViscosity_table.csv" };
  // setup
  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
  ProblemManager & problem = state.getProblemManager();
  problem.parseInputString( xmlInput );
  problem.problemSetup();
  problem.applyInitialConditions();

  EXPECT_TRUE( std::filesystem::exists( f1 ));
  EXPECT_TRUE( std::filesystem::exists( f2 ));
  EXPECT_TRUE( std::filesystem::exists( f3 ));
  EXPECT_TRUE( std::filesystem::exists( f4 ));
  EXPECT_FALSE( std::filesystem::exists( f5 ));

  std::filesystem::remove( "fluid_phaseModel1_PhillipsBrineDensity_table.csv" );
  std::filesystem::remove( "fluid_phaseModel2_SpanWagnerCO2Density_table.csv" );
  std::filesystem::remove( "fluid_phaseModel2_FenghourCO2Viscosity_table.csv" );
  std::filesystem::remove( "fluid_CO2Solubility_co2Dissolution_table.csv" );
  std::filesystem::remove( "fluid_phaseModel1_PhillipsBrineViscosity_table.csv" );
}

int main( int argc, char * * argv )
{
  writeTableToFile( "co2flash.txt", co2flash );
  writeTableToFile( "pvtliquid.txt", pvtLiquid );
  writeTableToFile( "pvtgas.txt", pvtGas );
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  basicCleanup();
  return result;
}
