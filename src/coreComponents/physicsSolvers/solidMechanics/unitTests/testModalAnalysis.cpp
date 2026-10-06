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
 * @file testModalAnalysis.cpp
 *
 * Tests of the Modal time integration option of SolidMechanicsLagrangianFEM on a bar made of one column of
 * trilinear hexahedra, with Poisson ratio zero and lumped mass.
 *
 * When the lateral displacements are constrained, the axial motion that is uniform over the cross-section is
 * exactly a one-dimensional chain of linear bar elements with lumped mass. Its eigenvalues are computed here
 * with a dense solver and used as reference. Modes that are not uniform over the cross-section also have shear
 * energy, have no closed-form reference, and are recognized by their zero participation factor.
 */

#include "denseLinearAlgebra/interfaces/blaslapack/BlasLapackLA.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/initialization.hpp"
#include "mesh/DomainPartition.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

using namespace geos;

CommandLineOptions g_commandLineOptions;

namespace
{

/// Bar length and number of elements. E = rho = A = 1, so the element size is 1.
integer constexpr numElements = 8;
real64 constexpr barLength = numElements;

/// Frequency of the shift sigma = -1 (in eigenvalue units): -1 / ( 2 pi )
char const shiftFrequency[] = "-0.15915494309189535";

string replaceAll( string text, string const & token, string const & value )
{
  for( size_t pos = text.find( token ); pos != string::npos; pos = text.find( token, pos + value.size() ) )
  {
    text.replace( pos, token.size(), value );
  }
  return text;
}

string const lateralConstraints =
  R"xml(
    <FieldSpecification name="fixYneg" objectPath="nodeManager" fieldName="totalDisplacement" component="1" scale="0.0" setNames="{ yneg }"/>
    <FieldSpecification name="fixYpos" objectPath="nodeManager" fieldName="totalDisplacement" component="1" scale="0.0" setNames="{ ypos }"/>
    <FieldSpecification name="fixZneg" objectPath="nodeManager" fieldName="totalDisplacement" component="2" scale="0.0" setNames="{ zneg }"/>
    <FieldSpecification name="fixZpos" objectPath="nodeManager" fieldName="totalDisplacement" component="2" scale="0.0" setNames="{ zpos }"/>
)xml";

string const clampedEnd =
  R"xml(
    <FieldSpecification name="fixXneg" objectPath="nodeManager" fieldName="totalDisplacement" component="0" scale="0.0" setNames="{ xneg }"/>
)xml";

string makeInput( string const & constraints, integer const numModes, integer const blockSize )
{
  string xml =
    R"xml(
<Problem>
  <Solvers>
    <SolidMechanicsLagrangianFEM name="solid"
                                 timeIntegrationOption="Modal"
                                 discretization="FE1"
                                 targetRegions="{ Region }"
                                 modalNumModes="@MODES@"
                                 modalShiftFrequency="@SHIFT@"
                                 modalBlockSize="@BLOCK@"
                                 modalTolerance="1e-9">
      <LinearSolverParameters solverType="cg"
                              preconditionerType="jacobi"
                              krylovTol="1e-13"
                              krylovMaxIter="5000"/>
    </SolidMechanicsLagrangianFEM>
  </Solvers>

  <Mesh>
    <InternalMesh name="mesh"
                  elementTypes="{ C3D8 }"
                  xCoords="{ 0, @LENGTH@ }"
                  yCoords="{ 0, 1 }"
                  zCoords="{ 0, 1 }"
                  nx="{ @NX@ }"
                  ny="{ 1 }"
                  nz="{ 1 }"
                  cellBlockNames="{ cb }"/>
  </Mesh>

  <Events maxTime="1.0">
    <PeriodicEvent name="modalAnalysis" forceDt="1.0" target="/Solvers/solid"/>
  </Events>

  <NumericalMethods>
    <FiniteElements>
      <FiniteElementSpace name="FE1" order="1"/>
    </FiniteElements>
  </NumericalMethods>

  <ElementRegions>
    <CellElementRegion name="Region" cellBlocks="{ cb }" materialList="{ rock }"/>
  </ElementRegions>

  <Constitutive>
    <ElasticIsotropic name="rock" defaultDensity="1" defaultYoungModulus="1" defaultPoissonRatio="0"/>
  </Constitutive>

  <FieldSpecifications>
@CONSTRAINTS@
  </FieldSpecifications>
</Problem>
)xml";
  xml = replaceAll( xml, "@MODES@", std::to_string( numModes ) );
  xml = replaceAll( xml, "@SHIFT@", shiftFrequency );
  xml = replaceAll( xml, "@BLOCK@", std::to_string( blockSize ) );
  xml = replaceAll( xml, "@LENGTH@", std::to_string( barLength ) );
  xml = replaceAll( xml, "@NX@", std::to_string( numElements ) );
  xml = replaceAll( xml, "@CONSTRAINTS@", constraints );
  return xml;
}

struct ModalResult
{
  std::vector< real64 > eigenvalues;
  std::vector< real64 > residuals;
  std::vector< std::array< real64, 3 > > participation;
  std::vector< real64 > frequencies;
  bool shapesRegistered = false;
};

ModalResult runModalAnalysis( string const & xml )
{
  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
  ProblemManager & problem = state.getProblemManager();

  // Partition the bar along its axis when running with several MPI ranks
  dataRepository::Group & commandLine = problem.getGroup< dataRepository::Group >( problem.groupKeys.commandLine );
  commandLine.getReference< integer >( problem.viewKeys.xPartitionsOverride ) = MpiWrapper::commSize( MPI_COMM_GEOS );

  problem.parseInputString( xml );
  problem.problemSetup();
  problem.applyInitialConditions();
  EXPECT_FALSE( problem.runSimulation() ) << "Simulation exited early.";

  SolidMechanicsLagrangianFEM & solver = problem.getGroupByPath< SolidMechanicsLagrangianFEM >( "/Solvers/solid" );

  ModalResult result;
  arrayView1d< real64 const > const lambda = solver.modalEigenvalues();
  arrayView1d< real64 const > const residual = solver.modalResiduals();
  arrayView1d< real64 const > const frequency = solver.modalFrequencies();
  arrayView2d< real64 const > const gamma = solver.modalParticipationFactors();
  for( localIndex k = 0; k < lambda.size(); ++k )
  {
    result.eigenvalues.push_back( lambda[k] );
    result.residuals.push_back( residual[k] );
    result.frequencies.push_back( frequency[k] );
    result.participation.push_back( { gamma( k, 0 ), gamma( k, 1 ), gamma( k, 2 ) } );
  }

  NodeManager const & nodes = problem.getDomainPartition().getMeshBody( 0 ).getBaseDiscretization().getNodeManager();
  result.shapesRegistered = nodes.hasWrapper( SolidMechanicsLagrangianFEM::modeShapeFieldName( 1 ) ) &&
                            nodes.hasWrapper( SolidMechanicsLagrangianFEM::modeShapeFieldName( LvArray::integerConversion< integer >( lambda.size() ) ) );
  return result;
}

/**
 * @brief Eigenvalues of a chain of linear bar elements (E = A = rho = 1) with lumped mass.
 * @param clamped if true, the first node is fixed
 * @return the eigenvalues of K x = lambda M x in ascending order
 */
std::vector< real64 > chainEigenvalues( bool const clamped )
{
  integer const numNodes = numElements + 1;
  real64 const h = barLength / numElements;

  std::vector< std::vector< real64 > > K( numNodes, std::vector< real64 >( numNodes, 0.0 ) );
  std::vector< real64 > m( numNodes, h );
  m.front() = 0.5 * h;
  m.back() = 0.5 * h;
  for( integer e = 0; e < numElements; ++e )
  {
    K[e][e] += 1.0 / h;
    K[e + 1][e + 1] += 1.0 / h;
    K[e][e + 1] -= 1.0 / h;
    K[e + 1][e] -= 1.0 / h;
  }

  integer const first = clamped ? 1 : 0;
  integer const n = numNodes - first;
  array2d< real64, MatrixLayout::COL_MAJOR_PERM > S( n, n );
  array2d< real64, MatrixLayout::COL_MAJOR_PERM > V( n, n );
  array1d< real64 > lambda( n );
  for( integer i = 0; i < n; ++i )
  {
    for( integer j = 0; j < n; ++j )
    {
      S( i, j ) = K[first + i][first + j] / std::sqrt( m[first + i] * m[first + j] );
    }
  }
  BlasLapackLA::matrixSymmetricEigen( S.toSliceConst(), lambda.toSlice(), V.toSlice() );
  return std::vector< real64 >( lambda.begin(), lambda.end() );
}

bool nearlyEqual( real64 const a, real64 const b )
{
  return std::fabs( a - b ) <= 1.0e-8 + 1.0e-6 * std::fabs( b );
}

/**
 * @brief Check the computed spectrum against the chain eigenvalues.
 * @param result the modal analysis result
 * @param reference the chain eigenvalues
 * @param netParticipation if true, the axial modes are recognized by a non-zero x participation factor (true when
 *        the end is clamped); if false (free end: elastic modes are M-orthogonal to the rigid translation), the
 *        reference eigenvalues are only required to be contained in the computed spectrum
 */
void checkAxialModes( ModalResult const & result, std::vector< real64 > const & reference, bool const netParticipation )
{
  ASSERT_FALSE( result.eigenvalues.empty() );
  EXPECT_TRUE( result.shapesRegistered );
  EXPECT_TRUE( std::is_sorted( result.eigenvalues.begin(), result.eigenvalues.end() ) );

  real64 const maxFound = result.eigenvalues.back();
  real64 const participationThreshold = 1.0e-6;

  for( size_t k = 0; k < result.eigenvalues.size(); ++k )
  {
    EXPECT_LT( result.residuals[k], 1.0e-6 ) << "mode " << k + 1;

    // The signed frequency is consistent with the eigenvalue
    real64 const f = result.frequencies[k];
    EXPECT_NEAR( ( f < 0.0 ? -1.0 : 1.0 ) * f * f * 4.0 * M_PI * M_PI, result.eigenvalues[k], 1.0e-8 + 1.0e-10 * std::fabs( result.eigenvalues[k] ) );

    if( netParticipation && std::fabs( result.participation[k][0] ) > participationThreshold )
    {
      bool found = false;
      for( real64 const ref : reference )
      {
        found = found || nearlyEqual( result.eigenvalues[k], ref );
      }
      EXPECT_TRUE( found ) << "mode " << k + 1 << " has eigenvalue " << result.eigenvalues[k] << " not in the chain spectrum";
    }
  }

  // Every chain eigenvalue below the largest computed eigenvalue must have been found: this detects missed modes
  for( real64 const ref : reference )
  {
    if( ref < maxFound - 1.0e-6 * std::fabs( maxFound ) )
    {
      bool found = false;
      for( size_t k = 0; k < result.eigenvalues.size(); ++k )
      {
        found = found || ( nearlyEqual( result.eigenvalues[k], ref ) &&
                           ( !netParticipation || std::fabs( result.participation[k][0] ) > participationThreshold ) );
      }
      EXPECT_TRUE( found ) << "chain eigenvalue " << ref << " was not found";
    }
  }
}

} // namespace

TEST( SolidMechanicsModal, clampedBarAxialModes )
{
  ModalResult const result = runModalAnalysis( makeInput( lateralConstraints + clampedEnd, 16, 1 ) );
  ASSERT_EQ( result.eigenvalues.size(), 16u );
  checkAxialModes( result, chainEigenvalues( true ), true );
}

TEST( SolidMechanicsModal, freeBarAxialModes )
{
  // Only the x translation is a rigid-body mode: K is singular and the shifted operator K + M is not
  ModalResult const result = runModalAnalysis( makeInput( lateralConstraints, 16, 1 ) );
  ASSERT_EQ( result.eigenvalues.size(), 16u );
  EXPECT_NEAR( result.eigenvalues[0], 0.0, 1.0e-8 );
  EXPECT_GT( std::fabs( result.participation[0][0] ), 1.0 );
  checkAxialModes( result, chainEigenvalues( false ), false );
}

TEST( SolidMechanicsModal, freeBodyRigidModes )
{
  // No constraint at all: six rigid-body modes at zero frequency, followed by elastic modes
  integer const numModes = 14;
  ModalResult const single = runModalAnalysis( makeInput( "", numModes, 1 ) );
  ASSERT_EQ( single.eigenvalues.size(), static_cast< size_t >( numModes ) );

  for( integer k = 0; k < 6; ++k )
  {
    EXPECT_NEAR( single.eigenvalues[k], 0.0, 1.0e-8 ) << "rigid mode " << k + 1;
  }
  EXPECT_GT( single.eigenvalues[6], 1.0e-3 );

  // The six rigid modes span the rigid-body space: for each direction, the squared participation factors
  // add up to the total mass rho * A * L
  for( integer d = 0; d < 3; ++d )
  {
    real64 sum = 0.0;
    for( integer k = 0; k < 6; ++k )
    {
      sum += single.participation[k][d] * single.participation[k][d];
    }
    EXPECT_NEAR( sum, barLength, 1.0e-6 * barLength ) << "direction " << d;
  }
  // Elastic modes are M-orthogonal to the rigid translations
  for( integer k = 6; k < numModes; ++k )
  {
    for( integer d = 0; d < 3; ++d )
    {
      EXPECT_NEAR( single.participation[k][d], 0.0, 1.0e-6 ) << "mode " << k + 1;
    }
    EXPECT_LT( single.residuals[k], 1.0e-6 );
  }

  // A block size equal to the multiplicity of the rigid modes gives the same spectrum
  ModalResult const block = runModalAnalysis( makeInput( "", numModes, 6 ) );
  ASSERT_EQ( block.eigenvalues.size(), single.eigenvalues.size() );
  for( size_t k = 6; k < single.eigenvalues.size(); ++k )
  {
    EXPECT_NEAR( block.eigenvalues[k], single.eigenvalues[k], 1.0e-6 * single.eigenvalues[k] ) << "mode " << k + 1;
  }
  for( size_t k = 0; k < 6; ++k )
  {
    EXPECT_NEAR( block.eigenvalues[k], 0.0, 1.0e-8 );
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
