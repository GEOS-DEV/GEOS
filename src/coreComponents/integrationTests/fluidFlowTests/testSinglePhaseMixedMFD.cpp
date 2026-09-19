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

/**
 * @file testSinglePhaseMixedMFD.cpp
 *
 * Jacobian of the mixed MFD solver against finite differences of its residual, on a hybrid mesh of hexahedra,
 * pyramids and tetrahedra, isothermal and thermal, with the three classifications (TPFA everywhere, consistent
 * everywhere, adaptive) and the three kinds of boundary faces (Dirichlet, Neumann with a value, Robin).
 * The thermal system is then checked after the row operation that removes the enthalpy reference.
 * Exact solutions: linear steady states with Neumann and Robin conditions, and the hydrostatic equilibrium on
 * distorted cells.
 */

#include "integrationTests/fluidFlowTests/testSingleFlowUtils.hpp"

#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "mixedMimetic/MixedMimeticBoundaryConditions.hpp"
#include "mixedMimetic/MixedMimeticFields.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseMixedMFD.hpp"

#include <gtest/gtest.h>

#include <fstream>
#include <iomanip>

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

namespace
{

real64 constexpr time_n = 0.0;
real64 constexpr dt = 1.0e4;
real64 constexpr specificHeat = 4000.0;

struct Setup
{
  bool isThermal;
  string consistencyTolerance;     ///< 1.0e+20 TPFA everywhere, 1.0e-20 consistent everywhere, 0.1 adaptive
  real64 referenceTemperature;
  bool temperatureDependentFluid;  ///< thermal expansion and temperature-dependent viscosity
};

string generateXml( string const & meshFile, Setup const & setup )
{
  std::ostringstream os;
  os <<
    R"xml(
  <Problem>
    <Solvers gravityVector="{ 0.0, 0.0, -9.81 }">
      <SinglePhaseMixedMFD name="flow" logLevel="0" discretization="mixedMFD" targetRegions="{ domain }")xml"
     << ( setup.isThermal ? R"xml( isThermal="1" temperature="350.0")xml" : "" ) <<
    R"xml(>
        <NonlinearSolverParameters newtonTol="1.0e-8" newtonMaxIter="10"/>
        <LinearSolverParameters directParallel="0"/>
      </SinglePhaseMixedMFD>
    </Solvers>
    <NumericalMethods>
      <MixedMimetic>
        <MixedMimeticDiscretization name="mixedMFD" innerProductType="RT" adaptiveConsistency="1"
                                    consistencyTolerance=")xml" << setup.consistencyTolerance <<
    R"xml(" nominalGradient="{ 1.0, 1.0, 1.0 }"/>
      </MixedMimetic>
    </NumericalMethods>
    <Mesh>
      <VTKMesh name="mesh" partitionRefinement="0" useGlobalIds="0" file=")xml" << meshFile <<
    R"xml("/>
    </Mesh>
    <Geometry>
      <Box name="west" xMin="{ -0.01, -0.01, -0.01 }" xMax="{ 0.01, 1.01, 1.01 }"/>
      <Box name="east" xMin="{ 0.99, -0.01, -0.01 }" xMax="{ 1.01, 1.01, 1.01 }"/>
      <Box name="top" xMin="{ -0.01, -0.01, 0.99 }" xMax="{ 1.01, 1.01, 1.01 }"/>
    </Geometry>
    <ElementRegions>
      <CellElementRegion name="domain" cellBlocks="{ * }" materialList=")xml"
     << ( setup.isThermal ? "{ rock, fluid, thermalCond }" : "{ rock, fluid }" ) <<
    R"xml("/>
    </ElementRegions>
    <Constitutive>
      <NullModel name="nullSolid"/>
      <PressurePorosity name="rockPorosity" defaultReferencePorosity="0.2" referencePressure="1.0e7" compressibility="1.0e-9"/>
      <ConstantPermeability name="rockPerm" permeabilityComponents="{ 1.0e-13, 2.0e-13, 0.5e-13 }"/>)xml";
  if( setup.isThermal )
  {
    os <<
      R"xml(
      <CompressibleSolidConstantPermeability name="rock" solidModelName="nullSolid" porosityModelName="rockPorosity"
                                             permeabilityModelName="rockPerm" solidInternalEnergyModelName="rockInternalEnergy"/>
      <SolidInternalEnergy name="rockInternalEnergy" referenceVolumetricHeatCapacity="2.0e6" referenceTemperature="300.0"
                           referenceInternalEnergy="0.0"/>
      <SinglePhaseThermalConductivity name="thermalCond" defaultThermalConductivityComponents="{ 2.0, 3.0, 1.0 }"
                                      thermalConductivityGradientComponents="{ 0, 0, 0 }" referenceTemperature="350.0"/>
      <ThermalCompressibleSinglePhaseFluid name="fluid" defaultDensity="1000" defaultViscosity="0.001" referencePressure="1.0e7"
                                           compressibility="1.0e-9" viscosibility="1.0e-9" specificHeatCapacity="4000"
                                           referenceInternalEnergy="1.0e6" referenceTemperature=")xml"
       << setup.referenceTemperature << "\" thermalExpansionCoeff=\""
       << ( setup.temperatureDependentFluid ? "3.0e-4" : "0.0" ) << "\" temperatureViscosityCoefficient=\""
       << ( setup.temperatureDependentFluid ? "1.0e-3" : "0.0" ) << "\"/>";
  }
  else
  {
    os <<
      R"xml(
      <CompressibleSolidConstantPermeability name="rock" solidModelName="nullSolid" porosityModelName="rockPorosity"
                                             permeabilityModelName="rockPerm"/>
      <CompressibleSinglePhaseFluid name="fluid" defaultDensity="1000" defaultViscosity="0.001" referencePressure="1.0e7"
                                    referenceDensity="1000" referenceViscosity="0.001" compressibility="1.0e-9"
                                    viscosibility="1.0e-9"/>)xml";
  }
  os <<
    R"xml(
    </Constitutive>
    <FieldSpecifications>
      <FieldSpecification name="initialPressure" initialCondition="1" setNames="{ all }" objectPath="ElementRegions"
                          fieldName="pressure" scale="1.0e7"/>
      <FieldSpecification name="westPressure" setNames="{ west }" objectPath="faceManager" fieldName="bcPressure" scale="1.1e7"/>
      <FieldSpecification name="eastPressure" setNames="{ east }" objectPath="faceManager" fieldName="bcPressure" scale="0.9e7"/>
      <FieldSpecification name="topMassFlux" setNames="{ top }" objectPath="faceManager" fieldName="bcMassFlux" scale="1.0e-4"/>)xml";
  if( setup.isThermal )
  {
    os <<
      R"xml(
      <FieldSpecification name="initialTemperature" initialCondition="1" setNames="{ all }" objectPath="ElementRegions"
                          fieldName="temperature" scale="350.0"/>
      <FieldSpecification name="westTemperature" setNames="{ west }" objectPath="faceManager" fieldName="bcTemperature" scale="380.0"/>
      <FieldSpecification name="eastCoefficient" setNames="{ east }" objectPath="faceManager" fieldName="bcHeatTransferCoefficient"
                          scale="5.0"/>
      <FieldSpecification name="eastAmbient" setNames="{ east }" objectPath="faceManager" fieldName="bcAmbientTemperature"
                          scale="300.0"/>
      <FieldSpecification name="topHeatFlux" setNames="{ top }" objectPath="faceManager" fieldName="bcHeatFlux" scale="20.0"/>
      <FieldSpecification name="westEnthalpy" setNames="{ west }" objectPath="faceManager" fieldName="bcEnthalpy" scale=")xml"
       << specificHeat * ( 380.0 - setup.referenceTemperature ) << "\"/>";
  }
  os <<
    R"xml(
    </FieldSpecifications>
  </Problem>
  )xml";
  return os.str();
}

/**
 * @brief Problem set up from the XML, with a non-equilibrium state stored as the beginning of the step.
 */
class MixedMFDProblem
{
public:

  MixedMFDProblem( Setup const & setup ):
    m_setup( setup ),
    m_state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {
    string const meshFile = string( TEST_BINARY_DIR ) + "/hex_pyr_tet_nested_mixed.vtu";
    setupProblemFromXML( m_state.getProblemManager(), generateXml( meshFile, setup ).c_str() );
    m_solver = &m_state.getProblemManager().getPhysicsSolverManager().getGroup< SinglePhaseMixedMFD >( "flow" );

    m_solver->setupSystem( domain(), m_solver->getDofManager(), m_solver->getLocalMatrix(),
                           m_solver->getSystemRhs(), m_solver->getSystemSolution() );
    m_solver->implicitStepSetup( time_n, dt, domain() );

    // fluxes, pressures and temperatures vary from one entity to the next: no term of the residual vanishes
    // and the upwind directions are invariant under the perturbations
    setState();
    m_solver->updateState( domain() );
    m_solver->implicitStepSetup( time_n, dt, domain() );
  }

  DomainPartition & domain() { return m_state.getProblemManager().getDomainPartition(); }

  SinglePhaseMixedMFD & solver() { return *m_solver; }

  /// residual and jacobian of the step, without the boundary conditions of the cells
  void assemble( CRSMatrixView< real64, globalIndex const > const & matrix, arrayView1d< real64 > const & rhs )
  {
    m_solver->assembleSystem( time_n, dt, domain(), m_solver->getDofManager(), matrix, rhs );
  }

  template< typename LAMBDA >
  void forSubRegions( LAMBDA && lambda )
  {
    m_solver->forDiscretizationOnMeshTargets( domain().getMeshBodies(), [&]( string const &, MeshLevel & mesh,
                                                                             string_array const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const, ElementSubRegionBase & subRegion )
      {
        lambda( subRegion );
      } );
    } );
  }

  FaceManager & faceManager()
  {
    FaceManager * result = nullptr;
    m_solver->forDiscretizationOnMeshTargets( domain().getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & )
    {
      result = &mesh.getFaceManager();
    } );
    return *result;
  }

private:

  void setState()
  {
    real64 const tRef = m_setup.referenceTemperature;
    bool const isThermal = m_setup.isThermal;
    forSubRegions( [&]( ElementSubRegionBase & subRegion )
    {
      arrayView2d< real64 const > const center = subRegion.getElementCenter();
      arrayView1d< real64 > const pres = subRegion.getField< fields::flow::pressure >();
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        real64 const s = std::sin( 7.0 * center[ei][0] + 3.0 * center[ei][1] + 5.0 * center[ei][2] );
        pres[ei] = 1.0e7 + 2.0e5 * s - 1.5e5 * center[ei][0];
        if( isThermal )
        {
          real64 const temp = 350.0 + 15.0 * std::cos( 4.0 * center[ei][0] - 6.0 * center[ei][2] );
          subRegion.getField< fields::flow::temperature >()[ei] = temp;
          // the enthalpy does not satisfy the equation of state: the closure equations have a nonzero residual
          subRegion.getField< fields::mixedMimetic::enthalpy >()[ei] = specificHeat * ( temp - tRef ) + 500.0 * s;
        }
      }
    } );

    FaceManager & faces = faceManager();
    arrayView2d< real64 const > const faceCenter = faces.faceCenter();
    arrayView1d< real64 > const massFlux = faces.getField< fields::mixedMimetic::faceMassFlux >();
    for( localIndex kf = 0; kf < faces.size(); ++kf )
    {
      real64 const s = std::sin( 11.0 * faceCenter[kf][0] + 13.0 * faceCenter[kf][1] + 17.0 * faceCenter[kf][2] );
      massFlux[kf] = ( s >= 0.0 ? 1.0 : -1.0 ) * 1.0e-3 * ( 1.0 + LvArray::math::abs( s ) );
      if( isThermal )
      {
        faces.getField< fields::mixedMimetic::faceHeatFlux >()[kf] = 3.0 * std::cos( 5.0 * faceCenter[kf][0] + 2.0 * faceCenter[kf][2] );
      }
    }
  }

  Setup const m_setup;
  GeosxState m_state;
  SinglePhaseMixedMFD * m_solver;
};

/**
 * @brief Column of the finite-difference jacobian for one perturbed value.
 */
void fillColumn( MixedMFDProblem & problem,
                 real64 & value,
                 real64 const perturbation,
                 globalIndex const dofIndex,
                 arrayView1d< real64 > const & residual,
                 arrayView1d< real64 const > const & residualOrig,
                 CRSMatrix< real64, globalIndex > & jacobian,
                 CRSMatrix< real64, globalIndex > & jacobianFD )
{
  problem.solver().resetStateToBeginningOfStep( problem.domain() );
  value += perturbation;
  problem.solver().updateState( problem.domain() );

  residual.zero();
  jacobian.zero();
  problem.assemble( jacobian.toViewConstSizes(), residual );
  fillNumericalJacobian( residual.toViewConst(), residualOrig, dofIndex, perturbation, jacobianFD.toViewConstSizes() );
}

void testNumericalJacobian( Setup const & setup, real64 const relTol )
{
  MixedMFDProblem problem( setup );
  SinglePhaseMixedMFD & solver = problem.solver();
  DofManager const & dofManager = solver.getDofManager();
  real64 const perturb = std::sqrt( std::numeric_limits< real64 >::epsilon() );

  CRSMatrix< real64, globalIndex > & jacobian = solver.getLocalMatrix();
  array1d< real64 > residual( jacobian.numRows() );

  solver.resetStateToBeginningOfStep( problem.domain() );
  solver.updateState( problem.domain() );
  residual.zero();
  jacobian.zero();
  problem.assemble( jacobian.toViewConstSizes(), residual.toView() );
  array1d< real64 > const residualOrig( residual );

  CRSMatrix< real64, globalIndex > jacobianFD( jacobian );
  jacobianFD.zero();

  // face unknowns
  FaceManager & faces = problem.faceManager();
  // the residual is linear in the face fluxes: a perturbation of 10 % keeps the upwind directions and resolves the
  // small couplings of the inner product, which a perturbation of sqrt( eps ) loses against residuals of 1.0e7 Pa
  real64 const fluxPerturb = 0.1;
  stdVector< std::pair< string, real64 > > faceFields = { { fields::mixedMimetic::faceMassFlux::key(), 0.0 } };
  if( setup.isThermal )
  {
    faceFields.emplace_back( fields::mixedMimetic::faceHeatFlux::key(), 1.0 );
  }
  for( auto const & faceField : faceFields )
  {
    arrayView1d< real64 > const flux = faces.getReference< array1d< real64 > >( faceField.first );
    arrayView1d< globalIndex const > const dofNumber = faces.getReference< array1d< globalIndex > >( dofManager.getKey( faceField.first ) );
    for( localIndex kf = 0; kf < faces.size(); ++kf )
    {
      if( dofNumber[kf] < 0 )
      {
        continue;
      }
      real64 const scale = LvArray::math::abs( flux[kf] );
      fillColumn( problem, flux[kf], fluxPerturb * ( scale + faceField.second ), dofNumber[kf],
                  residual.toView(), residualOrig.toViewConst(), jacobian, jacobianFD );
    }
  }

  // cell unknowns
  string const elemDofKey = dofManager.getKey( SinglePhaseMixedMFD::viewKeyStruct::elemDofFieldString() );
  problem.forSubRegions( [&]( ElementSubRegionBase & subRegion )
  {
    arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( elemDofKey );
    stdVector< arrayView1d< real64 > > unknowns = { subRegion.getField< fields::flow::pressure >() };
    if( setup.isThermal )
    {
      unknowns.push_back( subRegion.getField< fields::flow::temperature >() );
      unknowns.push_back( subRegion.getField< fields::mixedMimetic::enthalpy >() );
    }
    for( localIndex ei = 0; ei < subRegion.size(); ++ei )
    {
      for( std::size_t c = 0; c < unknowns.size(); ++c )
      {
        real64 const scale = LvArray::math::abs( unknowns[c][ei] );
        fillColumn( problem, unknowns[c][ei], perturb * ( scale + 1.0 ), dofNumber[ei] + c,
                    residual.toView(), residualOrig.toViewConst(), jacobian, jacobianFD );
      }
    }
  } );

  // analytical jacobian at the unperturbed state
  solver.resetStateToBeginningOfStep( problem.domain() );
  solver.updateState( problem.domain() );
  residual.zero();
  jacobian.zero();
  problem.assemble( jacobian.toViewConstSizes(), residual.toView() );

  compareLocalMatrices( jacobian.toViewConst(), jacobianFD.toViewConst(), relTol, 1.0e-6 );
}

/**
 * @brief System of the thermal step before and after the boundary conditions, which carry the row operation.
 */
struct ReducedSystem
{
  CRSMatrix< real64, globalIndex > jacobian;
  array1d< real64 > residual;
  CRSMatrix< real64, globalIndex > reducedJacobian;
  array1d< real64 > reducedResidual;
  array1d< globalIndex > cellDof;
  array1d< real64 > cellEnthalpy;
};

ReducedSystem assembleReducedSystem( Setup const & setup )
{
  MixedMFDProblem problem( setup );
  SinglePhaseMixedMFD & solver = problem.solver();
  CRSMatrix< real64, globalIndex > & jacobian = solver.getLocalMatrix();
  array1d< real64 > residual( jacobian.numRows() );

  ReducedSystem result;
  residual.zero();
  jacobian.zero();
  problem.assemble( jacobian.toViewConstSizes(), residual.toView() );
  result.jacobian = jacobian;
  result.residual = residual;

  solver.applyBoundaryConditions( time_n, dt, problem.domain(), solver.getDofManager(), jacobian.toViewConstSizes(), residual.toView() );
  result.reducedJacobian = jacobian;
  result.reducedResidual = residual;

  string const elemDofKey = solver.getDofManager().getKey( SinglePhaseMixedMFD::viewKeyStruct::elemDofFieldString() );
  problem.forSubRegions( [&]( ElementSubRegionBase & subRegion )
  {
    arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( elemDofKey );
    arrayView1d< real64 const > const enthalpy = subRegion.getField< fields::mixedMimetic::enthalpy >();
    for( localIndex ei = 0; ei < subRegion.size(); ++ei )
    {
      result.cellDof.emplace_back( dofNumber[ei] );
      result.cellEnthalpy.emplace_back( enthalpy[ei] );
    }
  } );
  return result;
}

} // namespace

TEST( SinglePhaseMixedMFDJacobian, isothermal_tpfa )
{
  testNumericalJacobian( { false, "1.0e+20", 350.0, true }, 1.0e-3 );
}

TEST( SinglePhaseMixedMFDJacobian, isothermal_consistent )
{
  testNumericalJacobian( { false, "1.0e-20", 350.0, true }, 1.0e-3 );
}

TEST( SinglePhaseMixedMFDJacobian, isothermal_adaptive )
{
  testNumericalJacobian( { false, "0.1", 350.0, true }, 1.0e-3 );
}

TEST( SinglePhaseMixedMFDJacobian, thermal_tpfa )
{
  testNumericalJacobian( { true, "1.0e+20", 350.0, true }, 1.0e-3 );
}

TEST( SinglePhaseMixedMFDJacobian, thermal_consistent )
{
  testNumericalJacobian( { true, "1.0e-20", 350.0, true }, 1.0e-3 );
}

TEST( SinglePhaseMixedMFDJacobian, thermal_adaptive )
{
  testNumericalJacobian( { true, "0.1", 350.0, true }, 1.0e-3 );
}

// the field specifications of the deck end up in the right descriptor of the right operator on every face:
// flow Dirichlet west and east, Neumann 1.0e-4 on top; heat Dirichlet west, Robin east, Neumann 20 on top;
// Neumann 0 (no-flow, adiabatic) on the other boundary faces and no condition inside
TEST( SinglePhaseMixedMFDBoundary, descriptors )
{
  using namespace geos::mixedMimeticBoundary;

  MixedMFDProblem problem( { true, "0.1", 350.0, true } );
  FaceManager & faces = problem.faceManager();
  arrayView2d< real64 const > const center = faces.faceCenter();
  arrayView2d< localIndex const > const elemList = faces.elementList();

  arrayView1d< integer const > const flowType = faces.getField< fields::mixedMimetic::flowBoundaryType >();
  arrayView1d< real64 const > const flowValue = faces.getField< fields::mixedMimetic::flowBoundaryValue >();
  arrayView1d< integer const > const heatType = faces.getField< fields::mixedMimetic::heatBoundaryType >();
  arrayView1d< real64 const > const heatValue = faces.getField< fields::mixedMimetic::heatBoundaryValue >();
  arrayView1d< real64 const > const heatCoefficient = faces.getField< fields::mixedMimetic::heatBoundaryCoefficient >();
  arrayView1d< integer const > const isEnthalpyBcFace = faces.getField< fields::mixedMimetic::isEnthalpyBcFace >();

  localIndex numWest = 0, numEast = 0, numTop = 0, numHomogeneousNeumann = 0, numInterior = 0;
  for( localIndex kf = 0; kf < faces.size(); ++kf )
  {
    SCOPED_TRACE( GEOS_FMT( "face {} at ( {}, {}, {} )", kf, center[kf][0], center[kf][1], center[kf][2] ) );
    bool const onBoundary = elemList[kf][0] < 0 || elemList[kf][1] < 0;
    if( !onBoundary )
    {
      ++numInterior;
      EXPECT_EQ( flowType[kf], BoundaryType::interior );
      EXPECT_EQ( heatType[kf], BoundaryType::interior );
    }
    else if( center[kf][0] < 0.01 )
    {
      ++numWest;
      EXPECT_EQ( flowType[kf], BoundaryType::dirichlet );
      EXPECT_DOUBLE_EQ( flowValue[kf], 1.1e7 );
      EXPECT_EQ( heatType[kf], BoundaryType::dirichlet );
      EXPECT_DOUBLE_EQ( heatValue[kf], 380.0 );
      EXPECT_EQ( isEnthalpyBcFace[kf], 1 );
    }
    else if( center[kf][0] > 0.99 )
    {
      ++numEast;
      EXPECT_EQ( flowType[kf], BoundaryType::dirichlet );
      EXPECT_DOUBLE_EQ( flowValue[kf], 0.9e7 );
      EXPECT_EQ( heatType[kf], BoundaryType::robin );
      EXPECT_DOUBLE_EQ( heatValue[kf], 300.0 );
      EXPECT_DOUBLE_EQ( heatCoefficient[kf], 5.0 );
      EXPECT_EQ( isEnthalpyBcFace[kf], 0 );
    }
    else if( center[kf][2] > 0.99 )
    {
      ++numTop;
      EXPECT_EQ( flowType[kf], BoundaryType::neumann );
      EXPECT_DOUBLE_EQ( flowValue[kf], 1.0e-4 );
      EXPECT_EQ( heatType[kf], BoundaryType::neumann );
      EXPECT_DOUBLE_EQ( heatValue[kf], 20.0 );
    }
    else
    {
      ++numHomogeneousNeumann;
      EXPECT_EQ( flowType[kf], BoundaryType::neumann );
      EXPECT_DOUBLE_EQ( flowValue[kf], 0.0 );
      EXPECT_EQ( heatType[kf], BoundaryType::neumann );
      EXPECT_DOUBLE_EQ( heatValue[kf], 0.0 );
    }
  }

  // every kind is present: the jacobian tests above do exercise them
  EXPECT_GT( numWest, 0 );
  EXPECT_GT( numEast, 0 );
  EXPECT_GT( numTop, 0 );
  EXPECT_GT( numHomogeneousNeumann, 0 );
  EXPECT_GT( numInterior, 0 );
}

// energy row <- energy row - h_K * mass row, every other row untouched
TEST( SinglePhaseMixedMFDEnthalpyReference, rowOperation )
{
  ReducedSystem const system = assembleReducedSystem( { true, "0.1", 350.0, true } );
  CRSMatrixView< real64 const, globalIndex const > const J = system.jacobian.toViewConst();
  CRSMatrixView< real64 const, globalIndex const > const Jr = system.reducedJacobian.toViewConst();

  array1d< integer > isEnergyRow( J.numRows() );
  isEnergyRow.zero();
  for( localIndex k = 0; k < system.cellDof.size(); ++k )
  {
    localIndex const massRow = system.cellDof[k];
    localIndex const energyRow = massRow + 1;
    real64 const h = system.cellEnthalpy[k];
    isEnergyRow[energyRow] = 1;

    // the pattern of the mass row is contained in the energy row
    ASSERT_EQ( J.numNonZeros( massRow ), J.numNonZeros( energyRow ) );
    for( localIndex j = 0; j < J.numNonZeros( energyRow ); ++j )
    {
      ASSERT_EQ( J.getColumns( massRow )[j], J.getColumns( energyRow )[j] );
      real64 const expected = J.getEntries( energyRow )[j] - h * J.getEntries( massRow )[j];
      EXPECT_NEAR( Jr.getEntries( energyRow )[j], expected, 1.0e-12 * ( LvArray::math::abs( expected ) + 1.0 ) );
    }
    real64 const expectedResidual = system.residual[energyRow] - h * system.residual[massRow];
    EXPECT_NEAR( system.reducedResidual[energyRow], expectedResidual, 1.0e-12 * ( LvArray::math::abs( expectedResidual ) + 1.0 ) );
  }

  for( localIndex i = 0; i < J.numRows(); ++i )
  {
    if( isEnergyRow[i] == 1 )
    {
      continue;
    }
    EXPECT_DOUBLE_EQ( system.reducedResidual[i], system.residual[i] );
    for( localIndex j = 0; j < J.numNonZeros( i ); ++j )
    {
      EXPECT_DOUBLE_EQ( Jr.getEntries( i )[j], J.getEntries( i )[j] );
    }
  }
}

// a shift of the enthalpy reference changes h by a constant and adds that constant times the mass row to the
// energy row: the reduced system does not see it. The density and the viscosity do not depend on the temperature
// here, so that the reference temperature of the fluid only shifts its enthalpy and its internal energy
TEST( SinglePhaseMixedMFDEnthalpyReference, referenceInvariance )
{
  ReducedSystem const systemA = assembleReducedSystem( { true, "0.1", 350.0, false } );
  ReducedSystem const systemB = assembleReducedSystem( { true, "0.1", 0.0, false } );

  // the unreduced energy rows do differ
  real64 maxDifference = 0.0;
  for( localIndex k = 0; k < systemA.cellDof.size(); ++k )
  {
    localIndex const energyRow = systemA.cellDof[k] + 1;
    maxDifference = LvArray::math::max( maxDifference, LvArray::math::abs( systemA.residual[energyRow] - systemB.residual[energyRow] ) );
  }
  EXPECT_GT( maxDifference, 1.0 );

  compareLocalMatrices( systemA.reducedJacobian.toViewConst(), systemB.reducedJacobian.toViewConst(), 1.0e-7, 1.0e-6 );
  for( localIndex i = 0; i < systemA.reducedResidual.size(); ++i )
  {
    real64 const scale = LvArray::math::abs( systemA.reducedResidual[i] ) + 1.0;
    EXPECT_NEAR( systemA.reducedResidual[i], systemB.reducedResidual[i], 1.0e-7 * scale ) << "row " << i;
  }
}

//======================== exact steady states in one dimension ======================
// Domain [0,1] x [0,0.1] x [0,0.1] with 5 cells along x, Dirichlet value at x = 0 and the tested condition at x = 1.
// The steady solutions are linear in x: both inner products reproduce them to round-off in a single time step
// much larger than the diffusion time:
//   Neumann  sigma f = g |f|               x( s ) = x_0 - g s / c
//   Robin    sigma f = alpha |f| ( x_f - g )   x( 1 ) = ( c x_0 + alpha g ) / ( c + alpha ),  linear in between
// with c the conductance per unit area: rho k / mu for the flow operator, the conductivity for the heat operator.

namespace
{

enum class BarCase : integer { flowNeumann, heatNeumann, heatRobin, thermalViscosity };

real64 constexpr barArea = 0.01;
real64 constexpr barPressure = 1.0e7;
real64 constexpr barMassFlux = 1.0e-3;
real64 constexpr barMobility = 1000.0 * 1.0e-13 / 1.0e-3;
real64 constexpr barTemperature = 400.0;
real64 constexpr barConductivity = 2.0;
real64 constexpr barHeatFlux = 30.0;
real64 constexpr barAlpha = 4.0;
real64 constexpr barAmbient = 300.0;
real64 constexpr barPressureDrop = 1.0e4;
real64 constexpr barTemperatureDrop = 100.0;
real64 constexpr barViscosityCoefficient = 5.0e-3;

string generateBarXml( BarCase const barCase, string const & consistencyTolerance )
{
  bool const isThermal = barCase != BarCase::flowNeumann;
  std::ostringstream os;
  os <<
    R"xml(
  <Problem>
    <Solvers gravityVector="{ 0.0, 0.0, 0.0 }">
      <SinglePhaseMixedMFD name="flow" logLevel="0" discretization="mixedMFD" targetRegions="{ domain }")xml"
     << ( isThermal ? R"xml( isThermal="1" temperature="350.0")xml" : "" ) <<
    R"xml(>
        <NonlinearSolverParameters newtonTol="1.0e-10" newtonMaxIter="10"/>
        <LinearSolverParameters directParallel="0"/>
      </SinglePhaseMixedMFD>
    </Solvers>
    <NumericalMethods>
      <MixedMimetic>
        <MixedMimeticDiscretization name="mixedMFD" innerProductType="RT" adaptiveConsistency="1"
                                    consistencyTolerance=")xml" << consistencyTolerance <<
    R"xml(" nominalGradient="{ 1.0, 1.0, 1.0 }"/>
      </MixedMimetic>
    </NumericalMethods>
    <Mesh>
      <InternalMesh name="mesh" elementTypes="{ C3D8 }" xCoords="{ 0, 1 }" yCoords="{ 0, 0.1 }" zCoords="{ 0, 0.1 }"
                    nx="{ 5 }" ny="{ 1 }" nz="{ 1 }" cellBlockNames="{ block }"/>
    </Mesh>
    <Geometry>
      <Box name="west" xMin="{ -0.001, -0.001, -0.001 }" xMax="{ 0.001, 0.101, 0.101 }"/>
      <Box name="east" xMin="{ 0.999, -0.001, -0.001 }" xMax="{ 1.001, 0.101, 0.101 }"/>
    </Geometry>
    <ElementRegions>
      <CellElementRegion name="domain" cellBlocks="{ * }" materialList=")xml"
     << ( isThermal ? "{ rock, fluid, thermalCond }" : "{ rock, fluid }" ) <<
    R"xml("/>
    </ElementRegions>
    <Constitutive>
      <NullModel name="nullSolid"/>
      <PressurePorosity name="rockPorosity" defaultReferencePorosity="0.2" referencePressure="1.0e7" compressibility="0.0"/>
      <ConstantPermeability name="rockPerm" permeabilityComponents="{ 1.0e-13, 1.0e-13, 1.0e-13 }"/>)xml";
  if( isThermal )
  {
    os <<
      R"xml(
      <CompressibleSolidConstantPermeability name="rock" solidModelName="nullSolid" porosityModelName="rockPorosity"
                                             permeabilityModelName="rockPerm" solidInternalEnergyModelName="rockInternalEnergy"/>
      <SolidInternalEnergy name="rockInternalEnergy" referenceVolumetricHeatCapacity="2.0e6" referenceTemperature="350.0"
                           referenceInternalEnergy="0.0"/>
      <SinglePhaseThermalConductivity name="thermalCond" defaultThermalConductivityComponents="{ 2.0, 2.0, 2.0 }"
                                      thermalConductivityGradientComponents="{ 0, 0, 0 }" referenceTemperature="350.0"/>
      <ThermalCompressibleSinglePhaseFluid name="fluid" defaultDensity="1000" defaultViscosity="0.001" referencePressure="1.0e7"
                                           referenceTemperature="350.0" compressibility="0.0" thermalExpansionCoeff="0.0"
                                           viscosibility="0.0" referenceInternalEnergy="1.0e6" temperatureViscosityCoefficient=")xml"
       << ( barCase == BarCase::thermalViscosity ? barViscosityCoefficient : 0.0 ) << "\" specificHeatCapacity=\""
       << ( barCase == BarCase::thermalViscosity ? 0.0 : 4000.0 ) << "\"/>";
  }
  else
  {
    os <<
      R"xml(
      <CompressibleSolidConstantPermeability name="rock" solidModelName="nullSolid" porosityModelName="rockPorosity"
                                             permeabilityModelName="rockPerm"/>
      <CompressibleSinglePhaseFluid name="fluid" defaultDensity="1000" defaultViscosity="0.001" referencePressure="1.0e7"
                                    referenceDensity="1000" referenceViscosity="0.001" compressibility="0.0" viscosibility="0.0"/>)xml";
  }
  os <<
    R"xml(
    </Constitutive>
    <FieldSpecifications>
      <FieldSpecification name="initialPressure" initialCondition="1" setNames="{ all }" objectPath="ElementRegions"
                          fieldName="pressure" scale="1.0e7"/>
      <FieldSpecification name="westPressure" setNames="{ west }" objectPath="faceManager" fieldName="bcPressure" scale="1.0e7"/>)xml";
  if( isThermal )
  {
    os <<
      R"xml(
      <FieldSpecification name="initialTemperature" initialCondition="1" setNames="{ all }" objectPath="ElementRegions"
                          fieldName="temperature" scale="350.0"/>
      <FieldSpecification name="westTemperature" setNames="{ west }" objectPath="faceManager" fieldName="bcTemperature" scale="400.0"/>)xml";
  }
  switch( barCase )
  {
    case BarCase::flowNeumann:
      os << R"xml(
      <FieldSpecification name="eastMassFlux" setNames="{ east }" objectPath="faceManager" fieldName="bcMassFlux" scale="1.0e-3"/>)xml";
      break;
    case BarCase::heatNeumann:
      os << R"xml(
      <FieldSpecification name="eastHeatFlux" setNames="{ east }" objectPath="faceManager" fieldName="bcHeatFlux" scale="30.0"/>)xml";
      break;
    case BarCase::thermalViscosity:
      os <<
        R"xml(
      <FieldSpecification name="eastPressure" setNames="{ east }" objectPath="faceManager" fieldName="bcPressure" scale="0.999e7"/>
      <FieldSpecification name="eastTemperature" setNames="{ east }" objectPath="faceManager" fieldName="bcTemperature" scale="300.0"/>)xml";
      break;
    case BarCase::heatRobin:
      os <<
        R"xml(
      <FieldSpecification name="eastCoefficient" setNames="{ east }" objectPath="faceManager" fieldName="bcHeatTransferCoefficient"
                          scale="4.0"/>
      <FieldSpecification name="eastAmbient" setNames="{ east }" objectPath="faceManager" fieldName="bcAmbientTemperature"
                          scale="300.0"/>)xml";
      break;
  }
  os <<
    R"xml(
    </FieldSpecifications>
  </Problem>
  )xml";
  return os.str();
}

/**
 * @brief Solve one time step, then compare the cell values and the face fluxes with the exact steady solution.
 * @param value0 the Dirichlet value at x = 0
 * @param outwardFlux the exact outward flux per unit area at x = 1
 * @param conductance the conductance per unit area
 */
void testBar( BarCase const barCase,
              string const & consistencyTolerance,
              real64 const value0,
              real64 const outwardFlux,
              real64 const conductance )
{
  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
  setupProblemFromXML( state.getProblemManager(), generateBarXml( barCase, consistencyTolerance ).c_str() );
  DomainPartition & domain = state.getProblemManager().getDomainPartition();
  SinglePhaseMixedMFD & solver = state.getProblemManager().getPhysicsSolverManager().getGroup< SinglePhaseMixedMFD >( "flow" );

  // time step much larger than the diffusion time: the accumulation term is below round-off
  solver.execute( 0.0, 1.0e15, 0, 0, 0, domain );

  bool const isFlow = barCase == BarCase::flowNeumann;
  MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();
  CellElementSubRegion const & subRegion = mesh.getElemManager().getRegion( 0 ).getSubRegion< CellElementSubRegion >( 0 );
  arrayView2d< real64 const > const center = subRegion.getElementCenter();
  arrayView1d< real64 const > const value = isFlow ? subRegion.getField< fields::flow::pressure >().toViewConst()
                                                   : subRegion.getField< fields::flow::temperature >().toViewConst();
  real64 const drop = LvArray::math::abs( outwardFlux / conductance );
  for( localIndex ei = 0; ei < subRegion.size(); ++ei )
  {
    real64 const exact = value0 - outwardFlux * center[ei][0] / conductance;
    EXPECT_NEAR( value[ei], exact, 1.0e-7 * drop ) << "cell at x = " << center[ei][0];
  }

  // the flux unknown of every face normal to x equals the exact flux, the others vanish
  FaceManager const & faces = mesh.getFaceManager();
  arrayView2d< real64 const > const faceNormal = faces.faceNormal();
  arrayView1d< real64 const > const flux = isFlow ? faces.getField< fields::mixedMimetic::faceMassFlux >().toViewConst()
                                                  : faces.getField< fields::mixedMimetic::faceHeatFlux >().toViewConst();
  for( localIndex kf = 0; kf < faces.size(); ++kf )
  {
    bool const normalToBar = LvArray::math::abs( faceNormal[kf][0] ) > 0.5;
    real64 const exact = normalToBar ? LvArray::math::abs( outwardFlux ) * barArea : 0.0;
    EXPECT_NEAR( LvArray::math::abs( flux[kf] ), exact, 1.0e-7 * LvArray::math::abs( outwardFlux ) * barArea ) << "face " << kf;
  }
}

/**
 * @brief Flow driven by a pressure drop with mu = mu_ref ( 1 - beta ( T - T_ref ) ) and a fluid of zero heat capacity.
 *
 * Without advection T = T_0 - dT x exactly, so mu( x ) = mu_0 + mu' x is linear. The mass flux per unit area m is
 * uniform and p( x ) = p_0 - ( m / ( rho k ) ) int_0^x mu. The two-point relations integrate mu by the midpoint rule,
 * exact for a linear function over whole cells: m = rho k dp / int_0^1 mu is exact, and the half cell that ends at a
 * cell centre gives p_K = p( x_K ) - m mu' h^2 / ( 8 rho k ).
 */
void testTemperatureDependentViscosity( string const & consistencyTolerance )
{
  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
  setupProblemFromXML( state.getProblemManager(), generateBarXml( BarCase::thermalViscosity, consistencyTolerance ).c_str() );
  DomainPartition & domain = state.getProblemManager().getDomainPartition();
  SinglePhaseMixedMFD & solver = state.getProblemManager().getPhysicsSolverManager().getGroup< SinglePhaseMixedMFD >( "flow" );

  solver.execute( 0.0, 1.0e15, 0, 0, 0, domain );

  real64 const referenceViscosity = 1.0e-3;
  real64 const mobilityCoefficient = 1000.0 * 1.0e-13;                     // rho k
  real64 const mu0 = referenceViscosity * ( 1.0 - barViscosityCoefficient * ( barTemperature - 350.0 ) );
  real64 const dMu = referenceViscosity * barViscosityCoefficient * barTemperatureDrop;
  real64 const massFlux = mobilityCoefficient * barPressureDrop / ( mu0 + 0.5 * dMu );
  real64 const h = 0.2;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();
  CellElementSubRegion const & subRegion = mesh.getElemManager().getRegion( 0 ).getSubRegion< CellElementSubRegion >( 0 );
  arrayView2d< real64 const > const center = subRegion.getElementCenter();
  arrayView1d< real64 const > const temp = subRegion.getField< fields::flow::temperature >();
  arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
  for( localIndex ei = 0; ei < subRegion.size(); ++ei )
  {
    real64 const x = center[ei][0];
    EXPECT_NEAR( temp[ei], barTemperature - barTemperatureDrop * x, 1.0e-7 * barTemperatureDrop ) << "cell at x = " << x;
    real64 const exactPressure = barPressure - massFlux / mobilityCoefficient * ( mu0 * x + 0.5 * dMu * x * x );
    real64 const midpointDefect = massFlux * dMu * h * h / ( 8.0 * mobilityCoefficient );
    EXPECT_NEAR( pres[ei], exactPressure - midpointDefect, 1.0e-7 * barPressureDrop ) << "cell at x = " << x;
  }

  FaceManager const & faces = mesh.getFaceManager();
  arrayView2d< real64 const > const faceNormal = faces.faceNormal();
  arrayView1d< real64 const > const flux = faces.getField< fields::mixedMimetic::faceMassFlux >();
  for( localIndex kf = 0; kf < faces.size(); ++kf )
  {
    real64 const exact = LvArray::math::abs( faceNormal[kf][0] ) > 0.5 ? massFlux * barArea : 0.0;
    EXPECT_NEAR( LvArray::math::abs( flux[kf] ), exact, 1.0e-7 * massFlux * barArea ) << "face " << kf;
  }
}

// outward heat flux of the Robin face: alpha ( x_f - g ) with x_f = ( c x_0 + alpha g ) / ( c + alpha )
real64 constexpr barRobinFlux = barAlpha * ( ( barConductivity * barTemperature + barAlpha * barAmbient ) / ( barConductivity + barAlpha ) - barAmbient );

} // namespace

TEST( SinglePhaseMixedMFDExactBar, temperatureDependentViscosity_tpfa )
{
  testTemperatureDependentViscosity( "1.0e+20" );
}

TEST( SinglePhaseMixedMFDExactBar, temperatureDependentViscosity_consistent )
{
  testTemperatureDependentViscosity( "1.0e-20" );
}

TEST( SinglePhaseMixedMFDExactBar, flowNeumann_tpfa )
{
  testBar( BarCase::flowNeumann, "1.0e+20", barPressure, barMassFlux, barMobility );
}

TEST( SinglePhaseMixedMFDExactBar, flowNeumann_consistent )
{
  testBar( BarCase::flowNeumann, "1.0e-20", barPressure, barMassFlux, barMobility );
}

TEST( SinglePhaseMixedMFDExactBar, heatNeumann_tpfa )
{
  testBar( BarCase::heatNeumann, "1.0e+20", barTemperature, barHeatFlux, barConductivity );
}

TEST( SinglePhaseMixedMFDExactBar, heatNeumann_consistent )
{
  testBar( BarCase::heatNeumann, "1.0e-20", barTemperature, barHeatFlux, barConductivity );
}

TEST( SinglePhaseMixedMFDExactBar, heatRobin_tpfa )
{
  testBar( BarCase::heatRobin, "1.0e+20", barTemperature, barRobinFlux, barConductivity );
}

TEST( SinglePhaseMixedMFDExactBar, heatRobin_consistent )
{
  testBar( BarCase::heatRobin, "1.0e-20", barTemperature, barRobinFlux, barConductivity );
}

//======================== hydrostatic equilibrium ===================================
// Fluid at rest in the unit cube: p = p_top on z = 1, every other boundary face impervious, constant density.
// The exact solution is m_f = 0 and p = p_top + rho |g| ( 1 - z ). With m = 0 the Darcy law of each cell reduces to
// pi_f = p_K - rho ( gamma_K - gamma_f ) for every inner product and every cell shape: the discrete pressure is exact
// at the cell centres, on planar and on non-planar faces alike.

namespace
{

real64 constexpr hydroTopPressure = 1.0e7;
real64 constexpr hydroDensity = 1000.0;
real64 constexpr hydroGravity = 9.81;

/**
 * @brief Write a mesh of n^3 hexahedra of the unit cube whose nodes are displaced by at most amplitude / n.
 *        A node of a boundary plane stays on that plane: the domain is unchanged and the interior faces are non-planar.
 */
string writeDistortedHexMesh( integer const n, real64 const amplitude )
{
  string const fileName = GEOS_FMT( "{}/hydrostatic_hexahedra_{}_{}.vtu", TEST_BINARY_DIR, n, integer( 100 * amplitude ) );
  integer const numNodes1d = n + 1;
  real64 const h = 1.0 / n;
  std::uint64_t seed = 12345;
  auto const random = [&seed]()
  {
    seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
    return 2.0 * ( real64( seed >> 11 ) / real64( 1ULL << 53 ) ) - 1.0;
  };

  std::ofstream os( fileName );
  os << std::setprecision( 17 );
  os << "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
     << "<UnstructuredGrid>\n<Piece NumberOfPoints=\"" << numNodes1d * numNodes1d * numNodes1d
     << "\" NumberOfCells=\"" << n * n * n << "\">\n"
     << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for( integer k = 0; k < numNodes1d; ++k )
  {
    for( integer j = 0; j < numNodes1d; ++j )
    {
      for( integer i = 0; i < numNodes1d; ++i )
      {
        integer const index[3] = { i, j, k };
        for( integer d = 0; d < 3; ++d )
        {
          bool const onBoundaryPlane = index[d] == 0 || index[d] == n;
          real64 const shift = amplitude * h * random();
          os << index[d] * h + ( onBoundaryPlane ? 0.0 : shift ) << " ";
        }
        os << "\n";
      }
    }
  }
  os << "</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
  auto const node = [numNodes1d]( integer const i, integer const j, integer const k )
  {
    return i + numNodes1d * ( j + numNodes1d * k );
  };
  for( integer k = 0; k < n; ++k )
  {
    for( integer j = 0; j < n; ++j )
    {
      for( integer i = 0; i < n; ++i )
      {
        os << node( i, j, k ) << " " << node( i + 1, j, k ) << " " << node( i + 1, j + 1, k ) << " " << node( i, j + 1, k ) << " "
           << node( i, j, k + 1 ) << " " << node( i + 1, j, k + 1 ) << " " << node( i + 1, j + 1, k + 1 ) << " " << node( i, j + 1, k + 1 ) << "\n";
      }
    }
  }
  os << "</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
  for( integer c = 1; c <= n * n * n; ++c )
  {
    os << 8 * c << " ";
  }
  os << "\n</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for( integer c = 0; c < n * n * n; ++c )
  {
    os << "12 ";
  }
  os << "\n</DataArray>\n</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
  return fileName;
}

string generateHydrostaticXml( string const & meshFile, string const & consistencyTolerance )
{
  std::ostringstream os;
  os <<
    R"xml(
  <Problem>
    <Solvers gravityVector="{ 0.0, 0.0, -9.81 }">
      <SinglePhaseMixedMFD name="flow" logLevel="0" discretization="mixedMFD" targetRegions="{ domain }">
        <NonlinearSolverParameters newtonTol="1.0e-10" newtonMaxIter="10"/>
        <LinearSolverParameters directParallel="0"/>
      </SinglePhaseMixedMFD>
    </Solvers>
    <NumericalMethods>
      <MixedMimetic>
        <MixedMimeticDiscretization name="mixedMFD" innerProductType="RT" adaptiveConsistency="1"
                                    consistencyTolerance=")xml" << consistencyTolerance <<
    R"xml(" nominalGradient="{ 1.0, 1.0, 1.0 }"/>
      </MixedMimetic>
    </NumericalMethods>
    <Mesh>
      <VTKMesh name="mesh" partitionRefinement="0" useGlobalIds="0" file=")xml" << meshFile <<
    R"xml("/>
    </Mesh>
    <Geometry>
      <Box name="top" xMin="{ -0.01, -0.01, 0.999 }" xMax="{ 1.01, 1.01, 1.001 }"/>
    </Geometry>
    <ElementRegions>
      <CellElementRegion name="domain" cellBlocks="{ * }" materialList="{ rock, fluid }"/>
    </ElementRegions>
    <Constitutive>
      <NullModel name="nullSolid"/>
      <PressurePorosity name="rockPorosity" defaultReferencePorosity="0.2" referencePressure="1.0e7" compressibility="0.0"/>
      <ConstantPermeability name="rockPerm" permeabilityComponents="{ 1.0e-13, 2.0e-13, 0.5e-13 }"/>
      <CompressibleSolidConstantPermeability name="rock" solidModelName="nullSolid" porosityModelName="rockPorosity"
                                             permeabilityModelName="rockPerm"/>
      <CompressibleSinglePhaseFluid name="fluid" defaultDensity="1000" defaultViscosity="0.001" referencePressure="1.0e7"
                                    referenceDensity="1000" referenceViscosity="0.001" compressibility="0.0" viscosibility="0.0"/>
    </Constitutive>
    <FieldSpecifications>
      <FieldSpecification name="initialPressure" initialCondition="1" setNames="{ all }" objectPath="ElementRegions"
                          fieldName="pressure" scale="1.0e7"/>
      <FieldSpecification name="topPressure" setNames="{ top }" objectPath="faceManager" fieldName="bcPressure" scale="1.0e7"/>
    </FieldSpecifications>
  </Problem>
  )xml";
  return os.str();
}

void testHydrostaticEquilibrium( string const & meshFile, string const & consistencyTolerance )
{
  SCOPED_TRACE( GEOS_FMT( "mesh {}, consistencyTolerance {}", meshFile, consistencyTolerance ) );
  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
  setupProblemFromXML( state.getProblemManager(), generateHydrostaticXml( meshFile, consistencyTolerance ).c_str() );
  DomainPartition & domain = state.getProblemManager().getDomainPartition();
  SinglePhaseMixedMFD & solver = state.getProblemManager().getPhysicsSolverManager().getGroup< SinglePhaseMixedMFD >( "flow" );

  solver.execute( 0.0, 1.0e6, 0, 0, 0, domain );

  real64 const gradient = hydroDensity * hydroGravity;
  real64 const tolerance = 1.0e-12 * hydroTopPressure;
  MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();

  // least-squares line p = a + b z through the cell values
  real64 sum1 = 0.0, sumZ = 0.0, sumZZ = 0.0, sumP = 0.0, sumZP = 0.0;
  mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
  {
    arrayView2d< real64 const > const center = subRegion.getElementCenter();
    arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
    for( localIndex ei = 0; ei < subRegion.size(); ++ei )
    {
      // exact pressure at the cell centre
      real64 const z = center[ei][2];
      EXPECT_NEAR( pres[ei], hydroTopPressure + gradient * ( 1.0 - z ), tolerance ) << "cell at z = " << z;

      real64 const dp = pres[ei] - hydroTopPressure;
      sum1 += 1.0; sumZ += z; sumZZ += z * z; sumP += dp; sumZP += z * dp;
    }
  } );

  // the pressure is linear in z with the slope - rho |g| and the value p_top at z = 1
  real64 const slope = ( sum1 * sumZP - sumZ * sumP ) / ( sum1 * sumZZ - sumZ * sumZ );
  real64 const intercept = ( sumP - slope * sumZ ) / sum1;
  EXPECT_NEAR( slope, -gradient, 1.0e-9 * gradient );
  EXPECT_NEAR( intercept + slope, 0.0, 1.0e-9 * gradient );
  mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
  {
    arrayView2d< real64 const > const center = subRegion.getElementCenter();
    arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
    for( localIndex ei = 0; ei < subRegion.size(); ++ei )
    {
      EXPECT_NEAR( pres[ei] - hydroTopPressure, intercept + slope * center[ei][2], tolerance );
    }
  } );

  // the fluid is at rest: the scale is the flux driven by the hydrostatic pressure difference
  arrayView1d< real64 const > const flux = mesh.getFaceManager().getField< fields::mixedMimetic::faceMassFlux >();
  real64 const fluxScale = hydroDensity * 1.0e-13 / 1.0e-3 * gradient;
  for( localIndex kf = 0; kf < flux.size(); ++kf )
  {
    EXPECT_NEAR( flux[kf], 0.0, 1.0e-10 * fluxScale ) << "face " << kf;
  }
}

} // namespace

TEST( SinglePhaseMixedMFDHydrostatic, hexahedraPyramidsTetrahedra )
{
  string const meshFile = string( TEST_BINARY_DIR ) + "/hex_pyr_tet_nested_mixed.vtu";
  for( string const tolerance : { "1.0e+20", "0.1", "1.0e-20" } )
  {
    testHydrostaticEquilibrium( meshFile, tolerance );
  }
}

TEST( SinglePhaseMixedMFDHydrostatic, distortedHexahedra )
{
  for( real64 const amplitude : { 0.0, 0.15, 0.3 } )
  {
    string const meshFile = writeDistortedHexMesh( 4, amplitude );
    for( string const tolerance : { "1.0e+20", "0.1", "1.0e-20" } )
    {
      testHydrostaticEquilibrium( meshFile, tolerance );
    }
  }
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
