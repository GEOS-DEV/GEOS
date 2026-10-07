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
 * Tests of the SolidMechanicsModalAnalysis solver on a bar made of one column of
 * trilinear hexahedra, with Poisson ratio zero and lumped mass.
 *
 * When the lateral displacements are constrained, the axial motion that is uniform over the cross-section is
 * exactly a one-dimensional chain of linear bar elements with lumped mass. Its eigenvalues are computed here
 * with a dense solver and used as reference. Modes that are not uniform over the cross-section also have shear
 * energy, have no closed-form reference, and are recognized by their zero participation factor.
 */

#include "denseLinearAlgebra/interfaces/blaslapack/BlasLapackLA.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/initialization.hpp"
#include "mesh/DomainPartition.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsModalAnalysis.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"

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

string const fullyClampedEnd =
  R"xml(
    <FieldSpecification name="clampX" objectPath="nodeManager" fieldName="totalDisplacement" component="0" scale="0.0" setNames="{ xneg }"/>
    <FieldSpecification name="clampY" objectPath="nodeManager" fieldName="totalDisplacement" component="1" scale="0.0" setNames="{ xneg }"/>
    <FieldSpecification name="clampZ" objectPath="nodeManager" fieldName="totalDisplacement" component="2" scale="0.0" setNames="{ xneg }"/>
)xml";

string makeInput( string const & constraints, integer const numModes, integer const blockSize,
                  string const & solverType = "arnoldi", integer const deflateRigidBodyModes = 0 )
{
  string xml =
    R"xml(
<Problem>
  <Solvers>
    <SolidMechanicsModalAnalysis name="solid"
                                 discretization="FE1"
                                 targetRegions="{ Region }"
                                 modalNumModes="@MODES@"
                                 modalShiftFrequency="@SHIFT@"
                                 modalBlockSize="@BLOCK@"
                                 modalSolverType="@EIGENSOLVER@"
                                 modalSubspaceSize="@SUBSPACE@"
                                 modalDeflateRigidBodyModes="@DEFLATE@"
                                 modalTolerance="1e-9">
      <LinearSolverParameters solverType="cg"
                              preconditionerType="jacobi"
                              krylovTol="1e-13"
                              krylovMaxIter="5000"/>
    </SolidMechanicsModalAnalysis>
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
  xml = replaceAll( xml, "@EIGENSOLVER@", solverType );
  xml = replaceAll( xml, "@DEFLATE@", std::to_string( deflateRigidBodyModes ) );
  // LOBPCG: guard vectors keep it from cutting a cluster of repeated eigenvalues at the last requested mode
  xml = replaceAll( xml, "@SUBSPACE@", solverType == "lobpcg" ? std::to_string( numModes + 6 ) : "0" );
  xml = replaceAll( xml, "@LENGTH@", std::to_string( barLength ) );
  xml = replaceAll( xml, "@NX@", std::to_string( numElements ) );
  xml = replaceAll( xml, "@CONSTRAINTS@", constraints );
  return xml;
}

struct ModalResult
{
  stdVector< real64 > eigenvalues;
  stdVector< real64 > residuals;
  stdVector< std::array< real64, 3 > > participation;
  stdVector< real64 > frequencies;
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

  SolidMechanicsModalAnalysis & solver = problem.getGroupByPath< SolidMechanicsModalAnalysis >( "/Solvers/solid" );

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
  result.shapesRegistered = nodes.hasWrapper( SolidMechanicsModalAnalysis::modeShapeFieldName( 1 ) ) &&
                            nodes.hasWrapper( SolidMechanicsModalAnalysis::modeShapeFieldName( LvArray::integerConversion< integer >( lambda.size() ) ) );

  if( MpiWrapper::commSize( MPI_COMM_GEOS ) > 1 )
  {
    // An independent all-reduce of owned values checks every received mode-field
    // component. This also exercises host/device validity after the halo exchange.
    arrayView1d< globalIndex const > const globalIds = nodes.localToGlobalMap();
    arrayView1d< integer const > const ghostRanks = nodes.ghostRank();
    globalIds.move( hostMemorySpace, false );
    ghostRanks.move( hostMemorySpace, false );
    globalIndex maxId = -1;
    for( localIndex a = 0; a < nodes.size(); ++a )
      maxId = std::max( maxId, globalIds[a] );
    globalIndex const globalNodes = MpiWrapper::max( maxId ) + 1;
    stdVector< real64 > ownedValues( static_cast< size_t >( globalNodes * lambda.size() * 3 ), 0.0 );
    for( localIndex k = 0; k < lambda.size(); ++k )
    {
      auto const field = nodes.getReference< fields::solidMechanics::array2dLayoutTotalDisplacement >(
        SolidMechanicsModalAnalysis::modeShapeFieldName( k + 1 ) ).toViewConst();
      field.move( hostMemorySpace, false );
      for( localIndex a = 0; a < nodes.size(); ++a )
        if( ghostRanks[a] < 0 )
          for( integer d = 0; d < 3; ++d )
            ownedValues[( k * globalNodes + globalIds[a] ) * 3 + d] = field( a, d );
    }
    stdVector< real64 > ownerValues( ownedValues.size(), 0.0 );
    MpiWrapper::allReduce( ownedValues, ownerValues, MpiWrapper::Reduction::Sum );
    for( localIndex k = 0; k < lambda.size(); ++k )
    {
      auto const field = nodes.getReference< fields::solidMechanics::array2dLayoutTotalDisplacement >(
        SolidMechanicsModalAnalysis::modeShapeFieldName( k + 1 ) ).toViewConst();
      array2d< real64 > readback( nodes.size(), 3 );
      arrayView2d< real64 > const readbackView = readback.toView();
      geos::forAll< geos::parallelDevicePolicy<> >( nodes.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        for( integer d = 0; d < 3; ++d )
          readbackView( a, d ) = field( a, d );
      } );
      readback.move( hostMemorySpace, false );
      for( localIndex a = 0; a < nodes.size(); ++a )
        if( ghostRanks[a] >= 0 )
          for( integer d = 0; d < 3; ++d )
            EXPECT_DOUBLE_EQ( readback( a, d ), ownerValues[( k * globalNodes + globalIds[a] ) * 3 + d] )
              << "mode " << k + 1 << ", global node " << globalIds[a] << ", component " << d;
    }
  }
  return result;
}

/**
 * @brief Eigenvalues of a chain of linear bar elements (E = A = rho = 1) with lumped mass.
 * @param clamped if true, the first node is fixed
 * @return the eigenvalues of K x = lambda M x in ascending order
 */
stdVector< real64 > chainEigenvalues( bool const clamped )
{
  integer const numNodes = numElements + 1;
  real64 const h = barLength / numElements;

  stdVector< stdVector< real64 > > K( numNodes, stdVector< real64 >( numNodes, 0.0 ) );
  stdVector< real64 > m( numNodes, h );
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
  return stdVector< real64 >( lambda.begin(), lambda.end() );
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
void checkAxialModes( ModalResult const & result, stdVector< real64 > const & reference, bool const netParticipation )
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


/// The bar of makeInput extended by a second cell region that the solver does not target
string makeInputWithOutsideRegion( string const & constraints, integer const numModes )
{
  string xml = makeInput( constraints, numModes, 1 );
  xml = replaceAll( xml, "xCoords=\"{ 0, " + std::to_string( barLength ) + " }\"",
                    "xCoords=\"{ 0, " + std::to_string( barLength ) + ", " + std::to_string( 2 * barLength ) + " }\"" );
  xml = replaceAll( xml, "nx=\"{ " + std::to_string( numElements ) + " }\"",
                    "nx=\"{ " + std::to_string( numElements ) + ", " + std::to_string( numElements ) + " }\"" );
  xml = replaceAll( xml, "cellBlockNames=\"{ cb }\"", "cellBlockNames=\"{ cb, cbOutside }\"" );
  xml = replaceAll( xml, "<CellElementRegion name=\"Region\" cellBlocks=\"{ cb }\" materialList=\"{ rock }\"/>",
                    "<CellElementRegion name=\"Region\" cellBlocks=\"{ cb }\" materialList=\"{ rock }\"/>\n"
                    "    <CellElementRegion name=\"Outside\" cellBlocks=\"{ cbOutside }\" materialList=\"{ rock }\"/>" );
  return xml;
}

} // namespace

TEST( SolidMechanicsModal, clampedBarAxialModes )
{
  for( string const solverType : { "arnoldi", "lobpcg" } )
  {
    SCOPED_TRACE( solverType );
    ModalResult const result = runModalAnalysis( makeInput( lateralConstraints + clampedEnd, 16, 1, solverType ) );
    ASSERT_EQ( result.eigenvalues.size(), 16u );
    checkAxialModes( result, chainEigenvalues( true ), true );
  }
}

TEST( SolidMechanicsModal, freeBarAxialModes )
{
  // Only the x translation is a rigid-body mode: K is singular and the shifted operator K + M is not
  for( string const solverType : { "arnoldi", "lobpcg" } )
  {
    SCOPED_TRACE( solverType );
    ModalResult const result = runModalAnalysis( makeInput( lateralConstraints, 16, 1, solverType ) );
    ASSERT_EQ( result.eigenvalues.size(), 16u );
    EXPECT_NEAR( result.eigenvalues[0], 0.0, 1.0e-8 );
    EXPECT_GT( std::fabs( result.participation[0][0] ), 1.0 );
    checkAxialModes( result, chainEigenvalues( false ), false );
  }
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

  // A block size equal to the multiplicity of the rigid modes, and the LOBPCG solver, give the same spectrum
  ModalResult const block = runModalAnalysis( makeInput( "", numModes, 6 ) );
  ModalResult const lobpcg = runModalAnalysis( makeInput( "", numModes, 1, "lobpcg" ) );
  ASSERT_EQ( block.eigenvalues.size(), single.eigenvalues.size() );
  ASSERT_EQ( lobpcg.eigenvalues.size(), single.eigenvalues.size() );
  for( size_t k = 6; k < single.eigenvalues.size(); ++k )
  {
    EXPECT_NEAR( block.eigenvalues[k], single.eigenvalues[k], 1.0e-6 * single.eigenvalues[k] ) << "mode " << k + 1;
    EXPECT_NEAR( lobpcg.eigenvalues[k], single.eigenvalues[k], 1.0e-6 * single.eigenvalues[k] ) << "mode " << k + 1;
  }
  for( size_t k = 0; k < 6; ++k )
  {
    EXPECT_NEAR( block.eigenvalues[k], 0.0, 1.0e-8 );
    EXPECT_NEAR( lobpcg.eigenvalues[k], 0.0, 1.0e-8 );
  }
  for( integer d = 0; d < 3; ++d )
  {
    real64 sum = 0.0;
    for( integer k = 0; k < 6; ++k )
    {
      sum += lobpcg.participation[k][d] * lobpcg.participation[k][d];
    }
    EXPECT_NEAR( sum, barLength, 1.0e-6 * barLength ) << "LOBPCG direction " << d;
  }
}

TEST( SolidMechanicsModal, deflatedRigidBodyModes )
{
  // The rigid-body modes are computed analytically and deflated: the spectrum is the one of the free body
  integer const numModes = 14;
  ModalResult const reference = runModalAnalysis( makeInput( "", numModes, 1 ) );
  ASSERT_EQ( reference.eigenvalues.size(), static_cast< size_t >( numModes ) );

  for( string const solverType : { "arnoldi", "lobpcg" } )
  {
    SCOPED_TRACE( solverType );
    ModalResult const result = runModalAnalysis( makeInput( "", numModes, 1, solverType, 1 ) );
    ASSERT_EQ( result.eigenvalues.size(), static_cast< size_t >( numModes ) );
    for( size_t k = 0; k < 6; ++k )
    {
      EXPECT_NEAR( result.eigenvalues[k], 0.0, 1.0e-8 ) << "rigid mode " << k + 1;
    }
    for( size_t k = 6; k < result.eigenvalues.size(); ++k )
    {
      EXPECT_NEAR( result.eigenvalues[k], reference.eigenvalues[k], 1.0e-6 * reference.eigenvalues[k] ) << "mode " << k + 1;
      EXPECT_LT( result.residuals[k], 1.0e-6 ) << "mode " << k + 1;
    }
    for( integer d = 0; d < 3; ++d )
    {
      real64 sum = 0.0;
      for( integer k = 0; k < 6; ++k )
      {
        sum += result.participation[k][d] * result.participation[k][d];
      }
      EXPECT_NEAR( sum, barLength, 1.0e-10 * barLength ) << "direction " << d;
    }
  }
}

TEST( SolidMechanicsModal, constraintsOutsideTargetRegionsAreIgnored )
{
  // A displacement condition on nodes that no target region owns has no degree of freedom to remove. These nodes
  // have the degree of freedom number -1, and a component offset must not turn it into the number of another node.
  string const outsideConstraint =
    R"xml(
    <FieldSpecification name="fixOutside" objectPath="nodeManager" fieldName="totalDisplacement" component="1" scale="0.0" setNames="{ xpos }"/>
)xml";
  integer const numModes = 14;
  ModalResult const reference = runModalAnalysis( makeInputWithOutsideRegion( "", numModes ) );
  ModalResult const constrained = runModalAnalysis( makeInputWithOutsideRegion( outsideConstraint, numModes ) );
  ASSERT_EQ( reference.eigenvalues.size(), static_cast< size_t >( numModes ) );
  ASSERT_EQ( constrained.eigenvalues.size(), static_cast< size_t >( numModes ) );

  // The target region is a free body: six rigid-body modes
  for( size_t k = 0; k < 6; ++k )
  {
    EXPECT_NEAR( reference.eigenvalues[k], 0.0, 1.0e-8 ) << "rigid mode " << k + 1;
    EXPECT_NEAR( constrained.eigenvalues[k], 0.0, 1.0e-8 ) << "rigid mode " << k + 1;
  }
  for( size_t k = 6; k < reference.eigenvalues.size(); ++k )
  {
    EXPECT_NEAR( constrained.eigenvalues[k], reference.eigenvalues[k], 1.0e-6 * reference.eigenvalues[k] ) << "mode " << k + 1;
  }
}

TEST( SolidMechanicsModal, clampedBarDoubleBendingModes )
{
  // The bar has a square cross-section, so every bending mode is double. The spectrum is: bending (double),
  // torsion (single), second bending (double). The shift is far below the spectrum: the Ritz values of the
  // completeness check are then clustered, and a check that is too weak misses the second copy of a double mode.
  for( integer const numModes : { 3, 5 } )
  {
    SCOPED_TRACE( numModes );
    string const xml = replaceAll( makeInput( fullyClampedEnd, numModes, 1 ), shiftFrequency, "-1.0" );
    ModalResult const result = runModalAnalysis( xml );
    ASSERT_EQ( result.eigenvalues.size(), static_cast< size_t >( numModes ) );
    stdVector< real64 > const & lambda = result.eigenvalues;
    EXPECT_NEAR( lambda[1], lambda[0], 1.0e-6 * lambda[0] ) << "first bending pair";
    EXPECT_GT( lambda[2], 1.01 * lambda[1] ) << "torsion mode";
    if( numModes >= 5 )
    {
      EXPECT_GT( lambda[3], 1.01 * lambda[2] ) << "second bending pair";
      EXPECT_NEAR( lambda[4], lambda[3], 1.0e-6 * lambda[3] ) << "second bending pair";
    }
  }
}

TEST( SolidMechanicsModal, fewModesWithDefaultBasis )
{
  // One or two modes, with the default Krylov basis size, must converge to the same values as a larger solve
  ModalResult const reference = runModalAnalysis( makeInput( fullyClampedEnd, 6, 1 ) );
  ASSERT_EQ( reference.eigenvalues.size(), 6u );
  for( integer const numModes : { 1, 2 } )
  {
    SCOPED_TRACE( numModes );
    ModalResult const result = runModalAnalysis( makeInput( fullyClampedEnd, numModes, 1 ) );
    ASSERT_EQ( result.eigenvalues.size(), static_cast< size_t >( numModes ) );
    for( integer k = 0; k < numModes; ++k )
    {
      EXPECT_NEAR( result.eigenvalues[k], reference.eigenvalues[k], 1.0e-7 * reference.eigenvalues[k] ) << "mode " << k + 1;
      EXPECT_LT( result.residuals[k], 1.0e-5 ) << "mode " << k + 1;
    }
  }
}

namespace
{

/**
 * @brief Modal analysis of the free unit tetrahedron (vertices 0, e_x, e_y, e_z), E = 1, nu = 0.3, rho = 1.
 * @param massType the modalMassType
 * @return the result with eight eigenvalues, the six rigid-body modes being deflated
 */
ModalResult runFreeTetrahedron( string const & massType )
{
  // Leave room for the Arnoldi search basis in the six-dimensional elastic complement.
  string xml = makeInput( "", 8, 1, "arnoldi", 1 );
  size_t const meshStart = xml.find( "  <Mesh>" );
  size_t const meshEnd = xml.find( "  </Mesh>", meshStart ) + string( "  </Mesh>" ).size();
  xml.replace( meshStart, meshEnd - meshStart,
               "<Mesh><VTKMesh name=\"mesh\" file=\"" GEOS_MODAL_TEST_DATA_DIR "/free_tetra.vtu\"/></Mesh>" );
  xml = replaceAll( xml, "cellBlocks=\"{ cb }\"", "cellBlocks=\"{ tetrahedra }\"" );
  xml = replaceAll( xml, "defaultPoissonRatio=\"0\"", "defaultPoissonRatio=\"0.3\"" );
  xml = replaceAll( xml, "modalTolerance=\"1e-9\"", "modalTolerance=\"1e-9\" modalMassType=\"" + massType + "\" modalVerifyFreeBody=\"1\"" );
  return runModalAnalysis( xml );
}

/**
 * @brief Dense reference for the free unit tetrahedron.
 * @param consistent if true, the consistent mass M = ( I + J ) / 120 (x) I_3, otherwise the lumped mass I / 24
 * @return the twelve eigenvalues of K x = lambda M x in ascending order
 */
stdVector< real64 > freeTetrahedronEigenvalues( bool const consistent )
{
  // For the unit simplex, grad N = (-1,-1,-1), e_x, e_y, e_z.
  real64 const gradients[4][3] = { {-1, -1, -1}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };
  real64 B[6][12] = {};
  for( integer a = 0; a < 4; ++a )
  {
    real64 const x = gradients[a][0], y = gradients[a][1], z = gradients[a][2];
    B[0][3*a] = x; B[1][3*a+1] = y; B[2][3*a+2] = z;
    B[3][3*a] = y; B[3][3*a+1] = x;
    B[4][3*a+1] = z; B[4][3*a+2] = y;
    B[5][3*a] = z; B[5][3*a+2] = x;
  }
  real64 const mu = 1.0 / 2.6, lambdaLame = 0.3 / ( 1.3 * 0.4 );
  real64 K[12][12] = {}, W[12][12] = {};
  // W = M^{-1/2}. Consistent mass: analytically, using J^2 = 4 J. Lumped mass: M = I / 24.
  for( integer i = 0; i < 12; ++i )
    for( integer j = 0; j < 12; ++j )
    {
      if( consistent )
        W[i][j] = i%3 == j%3 ? std::sqrt( 120.0 ) * ( ( i == j ? 1.0 : 0.0 ) - ( 1.0 - 1.0/std::sqrt( 5.0 ) )/4.0 ) : 0.0;
      else
        W[i][j] = i == j ? std::sqrt( 24.0 ) : 0.0;
      for( integer p = 0; p < 6; ++p )
        for( integer q = 0; q < 6; ++q )
        {
          real64 const D = ( p == q ? ( p < 3 ? 2*mu : mu ) : 0.0 ) + ( p < 3 && q < 3 ? lambdaLame : 0.0 );
          K[i][j] += B[p][i] * D * B[q][j] / 6.0;
        }
    }
  array2d< real64, MatrixLayout::COL_MAJOR_PERM > S( 12, 12 ), V( 12, 12 );
  array1d< real64 > reference( 12 );
  S.zero();
  for( integer i = 0; i < 12; ++i )
    for( integer j = 0; j < 12; ++j )
      for( integer p = 0; p < 12; ++p )
        for( integer q = 0; q < 12; ++q )
          S( i, j ) += W[i][p] * K[p][q] * W[q][j];
  BlasLapackLA::matrixSymmetricEigen( S.toSliceConst(), reference.toSlice(), V.toSlice() );
  return stdVector< real64 >( reference.begin(), reference.end() );
}

} // namespace

TEST( SolidMechanicsModal, consistentFreeTetrahedron )
{
  if( MpiWrapper::commSize( MPI_COMM_GEOS ) != 1 )
    GTEST_SKIP() << "A one-cell body is a serial assembly oracle.";
  ModalResult const result = runFreeTetrahedron( "consistent" );
  stdVector< real64 > const reference = freeTetrahedronEigenvalues( true );
  ASSERT_EQ( result.eigenvalues.size(), 8u );
  for( integer k = 0; k < 8; ++k )
    EXPECT_NEAR( result.eigenvalues[k], reference[k], 1e-8 + 1e-8 * std::fabs( reference[k] ) ) << "mode " << k + 1;
}

TEST( SolidMechanicsModal, lumpedFreeTetrahedron )
{
  // The lumped mass of a Tet4 is a quarter of the element mass per node. The nodal mass of the base solver used to
  // be six times too large on tetrahedra, which scaled every eigenvalue by 1/6.
  if( MpiWrapper::commSize( MPI_COMM_GEOS ) != 1 )
    GTEST_SKIP() << "A one-cell body is a serial assembly oracle.";
  ModalResult const result = runFreeTetrahedron( "lumped" );
  stdVector< real64 > const reference = freeTetrahedronEigenvalues( false );
  ASSERT_EQ( result.eigenvalues.size(), 8u );
  for( integer k = 0; k < 8; ++k )
    EXPECT_NEAR( result.eigenvalues[k], reference[k], 1e-8 + 1e-8 * std::fabs( reference[k] ) ) << "mode " << k + 1;
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
