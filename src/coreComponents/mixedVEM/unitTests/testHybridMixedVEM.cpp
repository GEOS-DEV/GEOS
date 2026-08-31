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

#include "denseLinearAlgebra/interfaces/blaslapack/BlasLapackLA.hpp"
#include "mixedVEM/HybridMixedVEM.hpp"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <vector>

using namespace geos;
using namespace geos::mixedVEM;

namespace
{

using Point = std::array< real64, 3 >;

/// A small conforming mesh: global nodes, global face node loops, cell to face maps.
struct Mesh
{
  std::vector< Point > nodes;
  std::vector< std::vector< int > > faces;
  std::vector< std::vector< int > > cells;
};

/// Coordinate accessor over one node loop.
struct LoopCoordinates
{
  std::vector< Point > const * nodes;
  std::vector< int > const * loop;

  real64 operator()( integer const k, integer const i ) const
  { return ( *nodes )[ static_cast< std::size_t >( ( *loop )[ static_cast< std::size_t >( k ) ] ) ][ static_cast< std::size_t >( i ) ]; }
};

void loopNormal( Mesh const & mesh,
                 std::vector< int > const & loop,
                 real64 (& normal)[3] )
{
  normal[0] = normal[1] = normal[2] = 0.0;

  std::size_t const n = loop.size();
  for( std::size_t k = 0; k < n; ++k )
  {
    Point const & p0 = mesh.nodes[ static_cast< std::size_t >( loop[k] ) ];
    Point const & p1 = mesh.nodes[ static_cast< std::size_t >( loop[ ( k + 1 ) % n ] ) ];

    normal[0] += ( p0[1] - p1[1] ) * ( p0[2] + p1[2] );
    normal[1] += ( p0[2] - p1[2] ) * ( p0[0] + p1[0] );
    normal[2] += ( p0[0] - p1[0] ) * ( p0[1] + p1[1] );
  }

  LvArray::tensorOps::normalize< 3 >( normal );
}

/// Two hexahedra sharing the face x = 1, the smallest mesh with an interior face.
Mesh twoCellMesh()
{
  Mesh mesh;

  mesh.nodes = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 },
    { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 },
    { 2, 0, 0 }, { 2, 1, 0 }, { 2, 0, 1 }, { 2, 1, 1 } };

  mesh.faces = { { 0, 3, 2, 1 },     // 0  z = 0 of the left cell
    { 4, 5, 6, 7 },                  // 1  z = 1 of the left cell
    { 0, 1, 5, 4 },                  // 2  y = 0 of the left cell
    { 1, 2, 6, 5 },                  // 3  x = 1, the interior face
    { 2, 3, 7, 6 },                  // 4  y = 1 of the left cell
    { 3, 0, 4, 7 },                  // 5  x = 0
    { 1, 2, 9, 8 },                  // 6  z = 0 of the right cell
    { 5, 6, 11, 10 },                // 7  z = 1 of the right cell
    { 1, 8, 10, 5 },                 // 8  y = 0 of the right cell
    { 2, 9, 11, 6 },                 // 9  y = 1 of the right cell
    { 8, 9, 11, 10 } };              // 10 x = 2

  mesh.cells = { { 0, 1, 2, 3, 4, 5 },
    { 3, 6, 7, 8, 9, 10 } };

  // slightly distort the interior face so the two cells are not mirror images
  mesh.nodes[6] = { 1.0, 1.0, 1.0 };

  return mesh;
}

/// Everything one cell needs, plus its operators and its condensed blocks.
struct Cell
{
  std::vector< FaceGeometry > faceGeom;
  ElementMoments moments;
  real64 center[3];
  real64 diameter;
  integer numFaces;
  integer numStressDof;

  array2d< real64 > divergence;
  array2d< real64 > divReconstruction;
  array2d< real64 > projection;
  array2d< real64 > workspace;
  array2d< real64 > stiffness;

  array2d< real64 > factorization;
  array2d< real64 > couplingTranspose;
  array2d< real64 > schur;
  real64 inverseDivGram[NUM_RM_DOF][NUM_RM_DOF] = { { 0.0 } };
};

void buildCell( Mesh const & mesh,
                integer const cellIndex,
                real64 const (&compliance)[NUM_SYM_COMP][NUM_SYM_COMP],
                Cell & cell )
{
  std::vector< int > const & cellFaces = mesh.cells[ static_cast< std::size_t >( cellIndex ) ];

  cell.numFaces = LvArray::integerConversion< integer >( cellFaces.size() );
  cell.numStressDof = NUM_FACE_DOF * cell.numFaces;
  cell.faceGeom.resize( static_cast< std::size_t >( cell.numFaces ) );

  // x_E is the average of the cell nodes, gathered from its faces
  std::vector< int > cellNodes;
  for( int face : cellFaces )
  {
    for( int node : mesh.faces[ static_cast< std::size_t >( face ) ] )
    {
      if( std::find( cellNodes.begin(), cellNodes.end(), node ) == cellNodes.end() )
      {
        cellNodes.push_back( node );
      }
    }
  }

  cell.center[0] = cell.center[1] = cell.center[2] = 0.0;
  for( int node : cellNodes )
  {
    for( integer i = 0; i < 3; ++i )
    {
      cell.center[i] += mesh.nodes[ static_cast< std::size_t >( node ) ][ static_cast< std::size_t >( i ) ];
    }
  }
  for( integer i = 0; i < 3; ++i )
  {
    cell.center[i] /= cellNodes.size();
  }

  cell.diameter = 0.0;
  for( std::size_t a = 0; a < cellNodes.size(); ++a )
  {
    for( std::size_t b = a + 1; b < cellNodes.size(); ++b )
    {
      real64 distanceSquared = 0.0;
      for( integer i = 0; i < 3; ++i )
      {
        real64 const d = mesh.nodes[ static_cast< std::size_t >( cellNodes[a] ) ][ static_cast< std::size_t >( i ) ]
                         - mesh.nodes[ static_cast< std::size_t >( cellNodes[b] ) ][ static_cast< std::size_t >( i ) ];
        distanceSquared += d * d;
      }
      cell.diameter = std::max( cell.diameter, std::sqrt( distanceSquared ) );
    }
  }

  initialize( cell.moments );

  for( integer lf = 0; lf < cell.numFaces; ++lf )
  {
    std::vector< int > const & loop = mesh.faces[ static_cast< std::size_t >( cellFaces[ static_cast< std::size_t >( lf ) ] ) ];

    LoopCoordinates const X { &mesh.nodes, &loop };
    integer const numNodes = LvArray::integerConversion< integer >( loop.size() );

    real64 normal[3];
    loopNormal( mesh, loop, normal );

    FaceGeometry & geom = cell.faceGeom[ static_cast< std::size_t >( lf ) ];

    computeFaceGeometry( X, numNodes, normal, geom );
    orientFace( cell.center, geom );
    accumulateElementMoments( X, numNodes, geom, cell.center, cell.moments );
  }

  cell.divergence.resize( NUM_RM_DOF, cell.numStressDof );
  cell.divReconstruction.resize( NUM_RM_DOF, cell.numStressDof );
  cell.projection.resize( NUM_SYM_COMP, cell.numStressDof );
  cell.workspace.resize( NUM_SYM_COMP, cell.numStressDof );
  cell.stiffness.resize( cell.numStressDof, cell.numStressDof );

  computeElementOperators( cell.faceGeom.data(),
                           cell.numFaces,
                           cell.center,
                           cell.diameter,
                           cell.moments,
                           compliance,
                           cell.divergence.toSlice(),
                           cell.divReconstruction.toSlice(),
                           cell.projection.toSlice(),
                           cell.workspace.toSlice(),
                           cell.stiffness.toSlice() );

  cell.factorization.resize( cell.numStressDof, cell.numStressDof );
  cell.couplingTranspose.resize( NUM_RM_DOF, cell.numStressDof );
  cell.schur.resize( cell.numStressDof, cell.numStressDof );

  bool const success = computeLocalCondensation( cell.stiffness.toSliceConst(),
                                                 cell.divergence.toSliceConst(),
                                                 cell.numFaces,
                                                 cell.factorization.toSlice(),
                                                 cell.couplingTranspose.toSlice(),
                                                 cell.inverseDivGram,
                                                 cell.schur.toSlice() );
  ASSERT_TRUE( success );
}

/// A reproducible pseudo random sequence, so the tests do not depend on <random>.
real64 sample( int const k )
{
  return std::sin( 1.0 + 0.7 * k ) + 0.3 * std::cos( 2.3 * k );
}

void solveDense( array2d< real64 > & matrix,
                 std::vector< real64 > const & rhs,
                 std::vector< real64 > & solution )
{
  integer const n = LvArray::integerConversion< integer >( rhs.size() );

  array1d< real64 > b( n );
  array1d< real64 > x( n );
  for( integer i = 0; i < n; ++i )
  {
    b[i] = rhs[ static_cast< std::size_t >( i ) ];
  }

  arrayView2d< real64 > const matrixView = matrix.toView();
  arrayView1d< real64 const > const bView = b.toViewConst();
  arrayView1d< real64 > const xView = x.toView();

  BlasLapackLA::solveLinearSystem( matrixView.toSlice(), bView.toSliceConst(), xView.toSlice() );

  solution.resize( static_cast< std::size_t >( n ) );
  for( integer i = 0; i < n; ++i )
  {
    solution[ static_cast< std::size_t >( i ) ] = x[i];
  }
}

} // namespace

TEST( HybridMixedVEM, localCondensationInvertsTheSaddleBlock )
{
  real64 compliance[NUM_SYM_COMP][NUM_SYM_COMP];
  makeIsotropicCompliance( 1.9, 0.8, compliance );

  Mesh const mesh = twoCellMesh();

  Cell cell;
  buildCell( mesh, 0, compliance, cell );

  integer const numStressDof = cell.numStressDof;

  // S_E is symmetric
  for( integer i = 0; i < numStressDof; ++i )
  {
    for( integer j = 0; j < numStressDof; ++j )
    {
      EXPECT_NEAR( cell.schur( i, j ), cell.schur( j, i ), 1e-14 );
    }
  }

  // S_E annihilates range(B_E^T): S_E B_E^T = W - W G^{-1} G = 0. Unlike the same statement
  // written on W, this one is independent of the basis chosen for T_h(f).
  for( integer k = 0; k < NUM_RM_DOF; ++k )
  {
    for( integer i = 0; i < numStressDof; ++i )
    {
      real64 value = 0.0;
      for( integer j = 0; j < numStressDof; ++j )
      {
        value += cell.schur( i, j ) * cell.divergence( k, j );
      }
      EXPECT_NEAR( value, 0.0, 1e-10 ) << "kernel mode " << k << " row " << i;
    }
  }

  // S_E is positive semidefinite
  for( int trial = 0; trial < 5; ++trial )
  {
    std::vector< real64 > x( static_cast< std::size_t >( numStressDof ) );
    for( integer i = 0; i < numStressDof; ++i )
    {
      x[ static_cast< std::size_t >( i ) ] = sample( 17 * trial + i );
    }

    real64 quadratic = 0.0;
    for( integer i = 0; i < numStressDof; ++i )
    {
      for( integer j = 0; j < numStressDof; ++j )
      {
        quadratic += x[ static_cast< std::size_t >( i ) ] * cell.schur( i, j ) * x[ static_cast< std::size_t >( j ) ];
      }
    }
    EXPECT_GT( quadratic, 0.0 );
  }

  // the recovered pair solves the element saddle system exactly
  std::vector< real64 > multiplier( static_cast< std::size_t >( numStressDof ) );
  std::vector< real64 > stressRhs( static_cast< std::size_t >( numStressDof ) );
  for( integer i = 0; i < numStressDof; ++i )
  {
    multiplier[ static_cast< std::size_t >( i ) ] = sample( 3 + i );
    stressRhs[ static_cast< std::size_t >( i ) ] = 0.1 * sample( 101 + i );
  }

  real64 bodyForce[NUM_RM_DOF];
  for( integer k = 0; k < NUM_RM_DOF; ++k )
  {
    bodyForce[k] = 0.2 * sample( 57 + k );
  }

  std::vector< real64 > stress( static_cast< std::size_t >( numStressDof ) );
  real64 displacement[NUM_RM_DOF];

  recoverElementSolution( cell.schur.toSliceConst(),
                          cell.couplingTranspose.toSliceConst(),
                          cell.inverseDivGram,
                          cell.faceGeom.data(),
                          cell.numFaces,
                          multiplier.data(),
                          stressRhs.data(),
                          bodyForce,
                          stress.data(),
                          displacement );

  // K_E sigma + B_E^T u = g_E + C_E^T lambda
  for( integer i = 0; i < numStressDof; ++i )
  {
    real64 value = 0.0;
    for( integer j = 0; j < numStressDof; ++j )
    {
      value += cell.stiffness( i, j ) * stress[ static_cast< std::size_t >( j ) ];
    }
    for( integer k = 0; k < NUM_RM_DOF; ++k )
    {
      value += cell.divergence( k, i ) * displacement[k];
    }

    real64 const expected = stressRhs[ static_cast< std::size_t >( i ) ]
                            + cell.faceGeom[ static_cast< std::size_t >( i / NUM_FACE_DOF ) ].outwardSign
                            * multiplier[ static_cast< std::size_t >( i ) ];

    EXPECT_NEAR( value, expected, 1e-9 ) << "stress row " << i;
  }

  // B_E sigma = -f_E
  for( integer k = 0; k < NUM_RM_DOF; ++k )
  {
    real64 value = 0.0;
    for( integer i = 0; i < numStressDof; ++i )
    {
      value += cell.divergence( k, i ) * stress[ static_cast< std::size_t >( i ) ];
    }
    EXPECT_NEAR( value, -bodyForce[k], 1e-9 ) << "displacement row " << k;
  }
}

TEST( HybridMixedVEM, hybridReproducesTheMixedSolution )
{
  real64 compliance[NUM_SYM_COMP][NUM_SYM_COMP];
  makeIsotropicCompliance( 1.4, 1.1, compliance );

  Mesh const mesh = twoCellMesh();

  integer const numCells = LvArray::integerConversion< integer >( mesh.cells.size() );
  integer const numFaces = LvArray::integerConversion< integer >( mesh.faces.size() );

  // face 3 is the only interior face of this mesh
  std::vector< int > faceCellCount( static_cast< std::size_t >( numFaces ), 0 );
  for( auto const & cellFaces : mesh.cells )
  {
    for( int face : cellFaces )
    {
      ++faceCellCount[ static_cast< std::size_t >( face ) ];
    }
  }

  std::vector< Cell > cells( static_cast< std::size_t >( numCells ) );
  for( integer e = 0; e < numCells; ++e )
  {
    buildCell( mesh, e, compliance, cells[ static_cast< std::size_t >( e ) ] );
  }

  integer const numGlobalStressDof = NUM_FACE_DOF * numFaces;
  integer const numGlobalDispDof = NUM_RM_DOF * numCells;
  integer const numMixedDof = numGlobalStressDof + numGlobalDispDof;

  // arbitrary Dirichlet data on the boundary faces, arbitrary load on the cells
  std::vector< real64 > dirichlet( static_cast< std::size_t >( numGlobalStressDof ), 0.0 );
  for( integer f = 0; f < numFaces; ++f )
  {
    if( faceCellCount[ static_cast< std::size_t >( f ) ] == 1 )
    {
      for( integer j = 0; j < NUM_FACE_DOF; ++j )
      {
        dirichlet[ static_cast< std::size_t >( NUM_FACE_DOF * f + j ) ] = 0.05 * sample( 7 * f + j );
      }
    }
  }

  std::vector< real64 > load( static_cast< std::size_t >( numGlobalDispDof ) );
  for( integer i = 0; i < numGlobalDispDof; ++i )
  {
    load[ static_cast< std::size_t >( i ) ] = 0.1 * sample( 211 + i );
  }

  // the mixed saddle point system
  array2d< real64 > mixedMatrix( numMixedDof, numMixedDof );
  mixedMatrix.zero();

  std::vector< real64 > mixedRhs( static_cast< std::size_t >( numMixedDof ), 0.0 );

  for( integer i = 0; i < numGlobalStressDof; ++i )
  {
    mixedRhs[ static_cast< std::size_t >( i ) ] = dirichlet[ static_cast< std::size_t >( i ) ];
  }
  for( integer i = 0; i < numGlobalDispDof; ++i )
  {
    mixedRhs[ static_cast< std::size_t >( numGlobalStressDof + i ) ] = -load[ static_cast< std::size_t >( i ) ];
  }

  for( integer e = 0; e < numCells; ++e )
  {
    Cell const & cell = cells[ static_cast< std::size_t >( e ) ];
    std::vector< int > const & cellFaces = mesh.cells[ static_cast< std::size_t >( e ) ];

    std::vector< integer > stressDof( static_cast< std::size_t >( cell.numStressDof ) );
    for( integer lf = 0; lf < cell.numFaces; ++lf )
    {
      for( integer j = 0; j < NUM_FACE_DOF; ++j )
      {
        stressDof[ static_cast< std::size_t >( NUM_FACE_DOF * lf + j ) ] =
          NUM_FACE_DOF * cellFaces[ static_cast< std::size_t >( lf ) ] + j;
      }
    }

    for( integer i = 0; i < cell.numStressDof; ++i )
    {
      integer const row = stressDof[ static_cast< std::size_t >( i ) ];

      for( integer j = 0; j < cell.numStressDof; ++j )
      {
        mixedMatrix( row, stressDof[ static_cast< std::size_t >( j ) ] ) += cell.stiffness( i, j );
      }

      for( integer k = 0; k < NUM_RM_DOF; ++k )
      {
        integer const column = numGlobalStressDof + NUM_RM_DOF * e + k;
        mixedMatrix( row, column ) += cell.divergence( k, i );
        mixedMatrix( column, row ) += cell.divergence( k, i );
      }
    }
  }

  std::vector< real64 > mixedSolution;
  solveDense( mixedMatrix, mixedRhs, mixedSolution );

  // the hybridized system: one multiplier block per interior face
  std::vector< integer > multiplierDof( static_cast< std::size_t >( numFaces ), -1 );
  integer numMultiplierDof = 0;
  for( integer f = 0; f < numFaces; ++f )
  {
    if( faceCellCount[ static_cast< std::size_t >( f ) ] == 2 )
    {
      multiplierDof[ static_cast< std::size_t >( f ) ] = numMultiplierDof;
      numMultiplierDof += NUM_FACE_DOF;
    }
  }
  ASSERT_EQ( numMultiplierDof, NUM_FACE_DOF );

  array2d< real64 > interfaceMatrix( numMultiplierDof, numMultiplierDof );
  interfaceMatrix.zero();

  std::vector< real64 > interfaceRhs( static_cast< std::size_t >( numMultiplierDof ), 0.0 );

  std::vector< std::vector< real64 > > stressRhs( static_cast< std::size_t >( numCells ) );
  std::vector< std::array< real64, NUM_RM_DOF > > bodyForce( static_cast< std::size_t >( numCells ) );

  for( integer e = 0; e < numCells; ++e )
  {
    Cell & cell = cells[ static_cast< std::size_t >( e ) ];
    std::vector< int > const & cellFaces = mesh.cells[ static_cast< std::size_t >( e ) ];

    // g_E carries the Dirichlet data, lambda carries the interior faces
    stressRhs[ static_cast< std::size_t >( e ) ].assign( static_cast< std::size_t >( cell.numStressDof ), 0.0 );
    for( integer lf = 0; lf < cell.numFaces; ++lf )
    {
      int const face = cellFaces[ static_cast< std::size_t >( lf ) ];
      if( faceCellCount[ static_cast< std::size_t >( face ) ] == 1 )
      {
        for( integer j = 0; j < NUM_FACE_DOF; ++j )
        {
          stressRhs[ static_cast< std::size_t >( e ) ][ static_cast< std::size_t >( NUM_FACE_DOF * lf + j ) ] =
            dirichlet[ static_cast< std::size_t >( NUM_FACE_DOF * face + j ) ];
        }
      }
    }

    for( integer k = 0; k < NUM_RM_DOF; ++k )
    {
      bodyForce[ static_cast< std::size_t >( e ) ][ static_cast< std::size_t >( k ) ] =
        load[ static_cast< std::size_t >( NUM_RM_DOF * e + k ) ];
    }

    // sigma_E at lambda = 0, whose jump is the interface right hand side
    std::vector< real64 > zeroMultiplier( static_cast< std::size_t >( cell.numStressDof ), 0.0 );
    std::vector< real64 > stress( static_cast< std::size_t >( cell.numStressDof ) );
    real64 displacement[NUM_RM_DOF];

    recoverElementSolution( cell.schur.toSliceConst(),
                            cell.couplingTranspose.toSliceConst(),
                            cell.inverseDivGram,
                            cell.faceGeom.data(),
                            cell.numFaces,
                            zeroMultiplier.data(),
                            stressRhs[ static_cast< std::size_t >( e ) ].data(),
                            bodyForce[ static_cast< std::size_t >( e ) ].data(),
                            stress.data(),
                            displacement );

    // H = sum_E C_E S_E C_E^T, restricted to the multiplier unknowns
    array2d< real64 > contribution( cell.numStressDof, cell.numStressDof );
    for( integer i = 0; i < cell.numStressDof; ++i )
    {
      for( integer j = 0; j < cell.numStressDof; ++j )
      {
        contribution( i, j ) = cell.schur( i, j );
      }
    }
    applyContinuityOperator( cell.faceGeom.data(), cell.numFaces, contribution.toSlice() );

    for( integer i = 0; i < cell.numStressDof; ++i )
    {
      integer const rowFace = multiplierDof[ static_cast< std::size_t >( cellFaces[ static_cast< std::size_t >( i / NUM_FACE_DOF ) ] ) ];
      if( rowFace < 0 )
      {
        continue;
      }
      integer const row = rowFace + i % NUM_FACE_DOF;

      interfaceRhs[ static_cast< std::size_t >( row ) ] -=
        cell.faceGeom[ static_cast< std::size_t >( i / NUM_FACE_DOF ) ].outwardSign * stress[ static_cast< std::size_t >( i ) ];

      for( integer j = 0; j < cell.numStressDof; ++j )
      {
        integer const columnFace = multiplierDof[ static_cast< std::size_t >( cellFaces[ static_cast< std::size_t >( j / NUM_FACE_DOF ) ] ) ];
        if( columnFace < 0 )
        {
          continue;
        }
        interfaceMatrix( row, columnFace + j % NUM_FACE_DOF ) += contribution( i, j );
      }
    }
  }

  // the interface operator is symmetric positive definite
  for( integer i = 0; i < numMultiplierDof; ++i )
  {
    for( integer j = 0; j < numMultiplierDof; ++j )
    {
      EXPECT_NEAR( interfaceMatrix( i, j ), interfaceMatrix( j, i ), 1e-12 );
    }
  }

  {
    array2d< real64 > copy( numMultiplierDof, numMultiplierDof );
    for( integer i = 0; i < numMultiplierDof; ++i )
    {
      for( integer j = 0; j < numMultiplierDof; ++j )
      {
        copy( i, j ) = interfaceMatrix( i, j );
      }
    }
    EXPECT_TRUE( choleskyFactorize( copy.toSlice(), numMultiplierDof ) );
  }

  std::vector< real64 > multiplierSolution;
  solveDense( interfaceMatrix, interfaceRhs, multiplierSolution );

  // recover (sigma_E, u_E) elementwise and compare with the mixed solution
  for( integer e = 0; e < numCells; ++e )
  {
    Cell const & cell = cells[ static_cast< std::size_t >( e ) ];
    std::vector< int > const & cellFaces = mesh.cells[ static_cast< std::size_t >( e ) ];

    std::vector< real64 > multiplier( static_cast< std::size_t >( cell.numStressDof ), 0.0 );
    for( integer i = 0; i < cell.numStressDof; ++i )
    {
      integer const rowFace = multiplierDof[ static_cast< std::size_t >( cellFaces[ static_cast< std::size_t >( i / NUM_FACE_DOF ) ] ) ];
      if( rowFace >= 0 )
      {
        multiplier[ static_cast< std::size_t >( i ) ] =
          multiplierSolution[ static_cast< std::size_t >( rowFace + i % NUM_FACE_DOF ) ];
      }
    }

    std::vector< real64 > stress( static_cast< std::size_t >( cell.numStressDof ) );
    real64 displacement[NUM_RM_DOF];

    recoverElementSolution( cell.schur.toSliceConst(),
                            cell.couplingTranspose.toSliceConst(),
                            cell.inverseDivGram,
                            cell.faceGeom.data(),
                            cell.numFaces,
                            multiplier.data(),
                            stressRhs[ static_cast< std::size_t >( e ) ].data(),
                            bodyForce[ static_cast< std::size_t >( e ) ].data(),
                            stress.data(),
                            displacement );

    for( integer i = 0; i < cell.numStressDof; ++i )
    {
      integer const globalDof = NUM_FACE_DOF * cellFaces[ static_cast< std::size_t >( i / NUM_FACE_DOF ) ] + i % NUM_FACE_DOF;

      EXPECT_NEAR( stress[ static_cast< std::size_t >( i ) ],
                   mixedSolution[ static_cast< std::size_t >( globalDof ) ], 1e-8 )
        << "cell " << e << " stress " << i;
    }

    for( integer k = 0; k < NUM_RM_DOF; ++k )
    {
      EXPECT_NEAR( displacement[k],
                   mixedSolution[ static_cast< std::size_t >( numGlobalStressDof + NUM_RM_DOF * e + k ) ], 1e-8 )
        << "cell " << e << " displacement " << k;
    }
  }
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
