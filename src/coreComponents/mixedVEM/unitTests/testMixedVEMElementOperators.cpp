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

#include "mixedVEM/MixedVEMAssembly.hpp"
#include "mixedVEM/MixedVEMCellOutput.hpp"
#include "mixedVEM/MixedVEMCompliance.hpp"
#include "mixedVEM/MixedVEMElementOperators.hpp"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <vector>

using namespace geos;
using namespace geos::mixedVEM;

namespace
{

using Point = std::array< real64, 3 >;

/// A polyhedron given by its nodes and by the node loop of each face.
struct Polyhedron
{
  std::vector< Point > nodes;
  std::vector< std::vector< int > > faces;
};

/// Coordinate accessor over one node loop.
struct LoopCoordinates
{
  std::vector< Point > const * nodes;
  std::vector< int > const * loop;

  real64 operator()( integer const k, integer const i ) const
  { return ( *nodes )[ static_cast< std::size_t >( ( *loop )[ static_cast< std::size_t >( k ) ] ) ][ static_cast< std::size_t >( i ) ]; }
};

/// Coordinate accessor over every node of the polyhedron.
struct AllCoordinates
{
  std::vector< Point > const * nodes;

  real64 operator()( integer const k, integer const i ) const
  { return ( *nodes )[ static_cast< std::size_t >( k ) ][ static_cast< std::size_t >( i ) ]; }
};

/// Newell normal of a node loop, the face-intrinsic reference normal of the tests.
void loopNormal( Polyhedron const & poly,
                 std::vector< int > const & loop,
                 real64 (& normal)[3] )
{
  normal[0] = normal[1] = normal[2] = 0.0;

  std::size_t const n = loop.size();
  for( std::size_t k = 0; k < n; ++k )
  {
    Point const & p0 = poly.nodes[ static_cast< std::size_t >( loop[k] ) ];
    Point const & p1 = poly.nodes[ static_cast< std::size_t >( loop[ ( k + 1 ) % n ] ) ];

    normal[0] += ( p0[1] - p1[1] ) * ( p0[2] + p1[2] );
    normal[1] += ( p0[2] - p1[2] ) * ( p0[0] + p1[0] );
    normal[2] += ( p0[0] - p1[0] ) * ( p0[1] + p1[1] );
  }

  LvArray::tensorOps::normalize< 3 >( normal );
}

/// Everything the element operators need for one polyhedron.
struct ElementData
{
  std::vector< FaceGeometry > faceGeom;
  ElementMoments moments;
  real64 center[3];
  real64 diameter;
  integer numFaces;
  integer numStressDof;
};

ElementData buildElement( Polyhedron const & poly,
                          real64 const (&offset)[3] = { 0.0, 0.0, 0.0 },
                          bool const flipNormals = false )
{
  ElementData data;

  data.numFaces = LvArray::integerConversion< integer >( poly.faces.size() );
  data.numStressDof = NUM_FACE_DOF * data.numFaces;
  data.faceGeom.resize( static_cast< std::size_t >( data.numFaces ) );

  // x_E is the node average, displaced by offset to exercise a non barycentric center
  data.center[0] = data.center[1] = data.center[2] = 0.0;
  for( Point const & p : poly.nodes )
  {
    for( integer i = 0; i < 3; ++i )
    {
      data.center[i] += p[ static_cast< std::size_t >( i ) ];
    }
  }
  for( integer i = 0; i < 3; ++i )
  {
    data.center[i] = data.center[i] / poly.nodes.size() + offset[i];
  }

  initialize( data.moments );

  for( integer lf = 0; lf < data.numFaces; ++lf )
  {
    std::vector< int > const & loop = poly.faces[ static_cast< std::size_t >( lf ) ];

    LoopCoordinates const X { &poly.nodes, &loop };
    integer const numNodes = LvArray::integerConversion< integer >( loop.size() );

    real64 normal[3];
    loopNormal( poly, loop, normal );

    // a reference normal opposed to the node loop, the branch where loopOrientation is -1
    if( flipNormals )
    {
      LvArray::tensorOps::scale< 3 >( normal, -1.0 );
    }

    FaceGeometry & geom = data.faceGeom[ static_cast< std::size_t >( lf ) ];

    computeFaceGeometry( X, numNodes, normal, geom );
    orientFace( data.center, geom );
    accumulateElementMoments( X, numNodes, geom, data.center, data.moments );
  }

  AllCoordinates const allX { &poly.nodes };
  data.diameter = computeElementDiameter( allX, LvArray::integerConversion< integer >( poly.nodes.size() ) );

  return data;
}

Polyhedron unitCube()
{
  return { { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 },
    { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 } },
    { { 0, 3, 2, 1 }, { 4, 5, 6, 7 },
      { 0, 1, 5, 4 }, { 1, 2, 6, 5 },
      { 2, 3, 7, 6 }, { 3, 0, 4, 7 } } };
}

Polyhedron distortedHex()
{
  Polyhedron poly = unitCube();

  // T_h(f) is built on planar faces, so every face touching the moved node is triangulated
  poly.nodes[6] = { 1.35, 1.2, 1.4 };

  poly.faces = { { 0, 3, 2, 1 },
    { 4, 5, 6 }, { 4, 6, 7 },
    { 0, 1, 5, 4 },
    { 1, 2, 6 }, { 1, 6, 5 },
    { 2, 3, 7 }, { 2, 7, 6 },
    { 3, 0, 4, 7 } };

  return poly;
}

Polyhedron shearedHex()
{
  Polyhedron poly = unitCube();

  // an affine shear keeps every quadrilateral face planar
  for( Point & p : poly.nodes )
  {
    p = { p[0] + 0.4 * p[2] + 0.15 * p[1], 1.3 * p[1] - 0.2 * p[2], 0.8 * p[2] };
  }

  return poly;
}

Polyhedron tetrahedron()
{
  return { { { 0.0, 0.0, 0.0 }, { 1.3, 0.1, 0.0 }, { 0.2, 1.1, 0.1 }, { 0.3, 0.2, 0.9 } },
    { { 0, 2, 1 }, { 0, 1, 3 }, { 1, 2, 3 }, { 2, 0, 3 } } };
}

Polyhedron hexagonalPrism()
{
  Polyhedron poly;

  std::vector< int > bottom, top;
  for( int k = 0; k < 6; ++k )
  {
    real64 const angle = 2.0 * M_PI * k / 6.0;
    real64 const radius = 1.0 + 0.15 * k;

    poly.nodes.push_back( { radius * std::cos( angle ), radius * std::sin( angle ), 0.0 } );
    poly.nodes.push_back( { radius * std::cos( angle ), radius * std::sin( angle ), 0.7 } );

    bottom.push_back( 2 * k );
    top.push_back( 2 * k + 1 );
  }

  std::vector< int > reversedBottom( bottom.rbegin(), bottom.rend() );
  poly.faces.push_back( reversedBottom );
  poly.faces.push_back( top );

  for( int k = 0; k < 6; ++k )
  {
    int const kp1 = ( k + 1 ) % 6;
    poly.faces.push_back( { 2 * k, 2 * kp1, 2 * kp1 + 1, 2 * k + 1 } );
  }

  return poly;
}

/// Degrees of freedom of a constant symmetric stress: (pi n)|_f expands on the constant modes.
std::vector< real64 > constantStressDofs( ElementData const & data,
                                          real64 const (&sigma)[NUM_SYM_COMP] )
{
  std::vector< real64 > dofs( static_cast< std::size_t >( data.numStressDof ), 0.0 );

  for( integer lf = 0; lf < data.numFaces; ++lf )
  {
    FaceGeometry const & geom = data.faceGeom[ static_cast< std::size_t >( lf ) ];

    real64 traction[3];
    symBasisTraction( sigma, geom.normal, traction );

    std::size_t const offset = static_cast< std::size_t >( NUM_FACE_DOF * lf );

    for( integer i = 0; i < 3; ++i )
    {
      dofs[offset + static_cast< std::size_t >( i )] = geom.area * traction[i];
    }
  }

  return dofs;
}

/// Rigid body motion R_i evaluated at a point.
void evaluateRigidMotion( integer const i,
                          real64 const (&r)[3],
                          real64 (& value)[3] )
{
  if( i < 3 )
  {
    value[0] = value[1] = value[2] = 0.0;
    value[i] = 1.0;
  }
  else
  {
    real64 axis[3] = { 0.0, 0.0, 0.0 };
    axis[i - 3] = 1.0;
    LvArray::tensorOps::crossProduct( value, axis, r );
  }
}

/// Element operators plus their storage, so that a test can hold on to all of them.
struct Operators
{
  array2d< real64 > divergence;
  array2d< real64 > divReconstruction;
  array2d< real64 > projection;
  array2d< real64 > workspace;
  array2d< real64 > stiffness;
};

Operators computeOperators( ElementData const & data,
                            real64 const (&compliance)[NUM_SYM_COMP][NUM_SYM_COMP] )
{
  Operators ops;

  ops.divergence.resize( NUM_RM_DOF, data.numStressDof );
  ops.divReconstruction.resize( NUM_RM_DOF, data.numStressDof );
  ops.projection.resize( NUM_SYM_COMP, data.numStressDof );
  ops.workspace.resize( NUM_SYM_COMP, data.numStressDof );
  ops.stiffness.resize( data.numStressDof, data.numStressDof );

  computeElementOperators( data.faceGeom.data(),
                           data.numFaces,
                           data.center,
                           data.diameter,
                           data.moments,
                           compliance,
                           ops.divergence.toSlice(),
                           ops.divReconstruction.toSlice(),
                           ops.projection.toSlice(),
                           ops.workspace.toSlice(),
                           ops.stiffness.toSlice() );

  return ops;
}

/// Second order rule on the tetrahedral fan, exact for every volume integrand below.
template< typename FUNC >
void forEachVolumeQuadraturePoint( Polyhedron const & poly,
                                   ElementData const & data,
                                   FUNC && func )
{
  constexpr real64 alpha = 0.5854101966249685;
  constexpr real64 beta = 0.1381966011250153;

  real64 const barycentric[4][4] = { { alpha, beta, beta, beta },
    { beta, alpha, beta, beta },
    { beta, beta, alpha, beta },
    { beta, beta, beta, alpha } };

  for( integer lf = 0; lf < data.numFaces; ++lf )
  {
    FaceGeometry const & geom = data.faceGeom[ static_cast< std::size_t >( lf ) ];
    std::vector< int > const & loop = poly.faces[ static_cast< std::size_t >( lf ) ];

    real64 const orientation = geom.loopOrientation * geom.outwardSign;

    std::size_t const n = loop.size();
    for( std::size_t k = 0; k < n; ++k )
    {
      Point const & p2 = poly.nodes[ static_cast< std::size_t >( loop[k] ) ];
      Point const & p3 = poly.nodes[ static_cast< std::size_t >( loop[ ( k + 1 ) % n ] ) ];

      real64 const vertices[4][3] = { { data.center[0], data.center[1], data.center[2] },
        { geom.center[0], geom.center[1], geom.center[2] },
        { p2[0], p2[1], p2[2] },
        { p3[0], p3[1], p3[2] } };

      real64 a[3], b[3], c[3], cross[3];
      for( integer i = 0; i < 3; ++i )
      {
        a[i] = vertices[1][i] - vertices[0][i];
        b[i] = vertices[2][i] - vertices[0][i];
        c[i] = vertices[3][i] - vertices[0][i];
      }
      LvArray::tensorOps::crossProduct( cross, b, c );

      real64 const volume = orientation * LvArray::tensorOps::AiBi< 3 >( a, cross ) / 6.0;

      for( integer q = 0; q < 4; ++q )
      {
        real64 point[3] = { 0.0, 0.0, 0.0 };
        for( integer v = 0; v < 4; ++v )
        {
          for( integer i = 0; i < 3; ++i )
          {
            point[i] += barycentric[q][v] * vertices[v][i];
          }
        }
        func( point, 0.25 * volume );
      }
    }
  }
}

/// Second order rule on the triangular fan of one face.
template< typename FUNC >
void forEachFaceQuadraturePoint( Polyhedron const & poly,
                                 ElementData const & data,
                                 integer const lf,
                                 FUNC && func )
{
  constexpr real64 alpha = 2.0 / 3.0;
  constexpr real64 beta = 1.0 / 6.0;

  real64 const barycentric[3][3] = { { alpha, beta, beta },
    { beta, alpha, beta },
    { beta, beta, alpha } };

  FaceGeometry const & geom = data.faceGeom[ static_cast< std::size_t >( lf ) ];
  std::vector< int > const & loop = poly.faces[ static_cast< std::size_t >( lf ) ];

  std::size_t const n = loop.size();
  for( std::size_t k = 0; k < n; ++k )
  {
    Point const & p1 = poly.nodes[ static_cast< std::size_t >( loop[k] ) ];
    Point const & p2 = poly.nodes[ static_cast< std::size_t >( loop[ ( k + 1 ) % n ] ) ];

    real64 const vertices[3][3] = { { geom.center[0], geom.center[1], geom.center[2] },
      { p1[0], p1[1], p1[2] },
      { p2[0], p2[1], p2[2] } };

    real64 e0[3], e1[3], cross[3];
    for( integer i = 0; i < 3; ++i )
    {
      e0[i] = vertices[1][i] - vertices[0][i];
      e1[i] = vertices[2][i] - vertices[0][i];
    }
    LvArray::tensorOps::crossProduct( cross, e0, e1 );

    real64 const area = 0.5 * geom.loopOrientation * LvArray::tensorOps::AiBi< 3 >( cross, geom.normal );

    for( integer q = 0; q < 3; ++q )
    {
      real64 point[3] = { 0.0, 0.0, 0.0 };
      for( integer v = 0; v < 3; ++v )
      {
        for( integer i = 0; i < 3; ++i )
        {
          point[i] += barycentric[q][v] * vertices[v][i];
        }
      }
      func( point, area / 3.0 );
    }
  }
}

/// p_a(x) = pi_a (x - x_E), the linear field whose symmetric gradient is pi_a.
void evaluateProjectionField( integer const a,
                              real64 const (&r)[3],
                              real64 (& value)[3] )
{
  real64 coefficients[NUM_SYM_COMP] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  coefficients[a] = 1.0;

  symBasisTraction( coefficients, r, value );
}

bool isPositiveDefinite( array2d< real64 > const & matrix )
{
  integer const n = LvArray::integerConversion< integer >( matrix.size( 0 ) );

  std::vector< real64 > factor( static_cast< std::size_t >( n * n ), 0.0 );

  for( integer i = 0; i < n; ++i )
  {
    for( integer j = 0; j <= i; ++j )
    {
      real64 sum = matrix( i, j );
      for( integer k = 0; k < j; ++k )
      {
        sum -= factor[ static_cast< std::size_t >( i * n + k ) ] * factor[ static_cast< std::size_t >( j * n + k ) ];
      }

      if( i == j )
      {
        if( sum <= 0.0 )
        {
          return false;
        }
        factor[ static_cast< std::size_t >( i * n + j ) ] = std::sqrt( sum );
      }
      else
      {
        factor[ static_cast< std::size_t >( i * n + j ) ] = sum / factor[ static_cast< std::size_t >( j * n + j ) ];
      }
    }
  }

  return true;
}

std::vector< std::pair< std::string, Polyhedron > > testElements()
{
  return { { "unitCube", unitCube() },
    { "shearedHex", shearedHex() },
    { "distortedHex", distortedHex() },
    { "tetrahedron", tetrahedron() },
    { "hexagonalPrism", hexagonalPrism() } };
}

} // namespace

TEST( MixedVEMGeometry, unitSquareFace )
{
  Polyhedron const poly = unitCube();
  std::vector< int > const loop = { 0, 1, 2, 3 };

  LoopCoordinates const X { &poly.nodes, &loop };
  real64 const normal[3] = { 0.0, 0.0, 1.0 };

  FaceGeometry geom;
  computeFaceGeometry( X, 4, normal, geom );

  EXPECT_NEAR( geom.area, 1.0, 1e-14 );
  EXPECT_NEAR( geom.center[0], 0.5, 1e-14 );
  EXPECT_NEAR( geom.center[1], 0.5, 1e-14 );
  EXPECT_NEAR( geom.center[2], 0.0, 1e-14 );

  EXPECT_NEAR( geom.m20, 1.0 / 12.0, 1e-14 );
  EXPECT_NEAR( geom.m02, 1.0 / 12.0, 1e-14 );
  EXPECT_NEAR( geom.m11, 0.0, 1e-14 );
}

TEST( MixedVEMGeometry, unitCubeMoments )
{
  ElementData const data = buildElement( unitCube() );

  EXPECT_NEAR( data.moments.volume, 1.0, 1e-13 );

  for( integer i = 0; i < 3; ++i )
  {
    EXPECT_NEAR( data.moments.firstMoment[i], 0.0, 1e-13 );

    for( integer j = 0; j < 3; ++j )
    {
      real64 const expected = ( i == j ) ? 1.0 / 12.0 : 0.0;
      EXPECT_NEAR( data.moments.secondMoment[i][j], expected, 1e-13 );
    }
  }

  EXPECT_NEAR( data.diameter, std::sqrt( 3.0 ), 1e-13 );
}

TEST( MixedVEMGeometry, volumeMatchesQuadrature )
{
  for( auto const & entry : testElements() )
  {
    ElementData const data = buildElement( entry.second );

    real64 quadratureVolume = 0.0;
    forEachVolumeQuadraturePoint( entry.second, data,
                                  [&]( real64 const (&)[3], real64 const weight )
    {
      quadratureVolume += weight;
    } );

    EXPECT_NEAR( quadratureVolume, data.moments.volume, 1e-12 * data.moments.volume ) << entry.first;
    EXPECT_GT( data.moments.volume, 0.0 ) << entry.first;
  }
}

TEST( MixedVEMElementOperators, constantStressPatchTest )
{
  real64 compliance[NUM_SYM_COMP][NUM_SYM_COMP];
  makeIsotropicCompliance( 1.7, 0.9, compliance );

  real64 const stresses[3][NUM_SYM_COMP] = { { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 },
    { 0.4, -1.1, 0.7, 0.3, -0.6, 0.9 } };

  // a deliberately non barycentric x_E, so the reconstruction cannot rely on m1 = 0
  real64 const offsets[2][3] = { { 0.0, 0.0, 0.0 }, { 0.07, -0.05, 0.03 } };

  for( auto const & entry : testElements() )
  {
    for( auto const & offset : offsets )
    {
      for( bool const flipNormals : { false, true } )
      {
        ElementData const data = buildElement( entry.second, offset, flipNormals );
        Operators const ops = computeOperators( data, compliance );

        real64 const scale = data.moments.volume + data.diameter;

        for( auto const & sigma : stresses )
        {
          std::vector< real64 > const dofs = constantStressDofs( data, sigma );

          // Pi_E reproduces constant symmetric stresses exactly
          for( integer a = 0; a < NUM_SYM_COMP; ++a )
          {
            real64 value = 0.0;
            for( integer j = 0; j < data.numStressDof; ++j )
            {
              value += ops.projection( a, j ) * dofs[ static_cast< std::size_t >( j ) ];
            }
            EXPECT_NEAR( value, sigma[a], 1e-11 ) << entry.first << " component " << a;
          }

          // div of a constant stress vanishes, hence B_E sigma = 0 and D_E sigma = 0
          for( integer i = 0; i < NUM_RM_DOF; ++i )
          {
            real64 divergenceValue = 0.0;
            real64 reconstructionValue = 0.0;
            for( integer j = 0; j < data.numStressDof; ++j )
            {
              divergenceValue += ops.divergence( i, j ) * dofs[ static_cast< std::size_t >( j ) ];
              reconstructionValue += ops.divReconstruction( i, j ) * dofs[ static_cast< std::size_t >( j ) ];
            }
            EXPECT_NEAR( divergenceValue, 0.0, 1e-11 * scale ) << entry.first << " mode " << i;
            EXPECT_NEAR( reconstructionValue, 0.0, 1e-11 * scale ) << entry.first << " mode " << i;
          }

          // the stabilization annihilates sigma, so K_E sigma = |E| P_E^T D sigma
          for( integer i = 0; i < data.numStressDof; ++i )
          {
            real64 value = 0.0;
            for( integer j = 0; j < data.numStressDof; ++j )
            {
              value += ops.stiffness( i, j ) * dofs[ static_cast< std::size_t >( j ) ];
            }

            real64 expected = 0.0;
            for( integer a = 0; a < NUM_SYM_COMP; ++a )
            {
              for( integer b = 0; b < NUM_SYM_COMP; ++b )
              {
                expected += ops.projection( a, i ) * compliance[a][b] * sigma[b];
              }
            }
            expected *= data.moments.volume;

            EXPECT_NEAR( value, expected, 1e-10 * scale ) << entry.first << " row " << i;
          }
        }
      }
    }
  }
}

TEST( MixedVEMElementOperators, stiffnessIsSymmetricPositiveDefinite )
{
  real64 compliance[NUM_SYM_COMP][NUM_SYM_COMP];
  makeIsotropicCompliance( 3.0, 1.4, compliance );

  for( auto const & entry : testElements() )
  {
    ElementData const data = buildElement( entry.second );
    Operators const ops = computeOperators( data, compliance );

    real64 maxEntry = 0.0;
    for( integer i = 0; i < data.numStressDof; ++i )
    {
      for( integer j = 0; j < data.numStressDof; ++j )
      {
        maxEntry = std::max( maxEntry, std::abs( ops.stiffness( i, j ) ) );
      }
    }

    for( integer i = 0; i < data.numStressDof; ++i )
    {
      for( integer j = 0; j < data.numStressDof; ++j )
      {
        EXPECT_NEAR( ops.stiffness( i, j ), ops.stiffness( j, i ), 1e-12 * maxEntry ) << entry.first;
      }
    }

    EXPECT_TRUE( isPositiveDefinite( ops.stiffness ) ) << entry.first;
  }
}

TEST( MixedVEMElementOperators, matchesDirectQuadrature )
{
  real64 compliance[NUM_SYM_COMP][NUM_SYM_COMP];
  makeIsotropicCompliance( 2.2, 1.1, compliance );

  for( auto const & entry : testElements() )
  {
    Polyhedron const & poly = entry.second;

    ElementData const data = buildElement( poly, { 0.05, 0.02, -0.04 }, true );
    Operators const ops = computeOperators( data, compliance );

    integer const numStressDof = data.numStressDof;
    real64 const scale = data.moments.volume + data.diameter;

    // B_E from face quadrature of (phi_j n_E) . R_i
    array2d< real64 > referenceDivergence( NUM_RM_DOF, numStressDof );
    referenceDivergence.zero();

    // face part of |E| P_E, and the stabilization block, from the same quadrature
    array2d< real64 > referenceProjection( NUM_SYM_COMP, numStressDof );
    referenceProjection.zero();

    array2d< real64 > referenceStabilization( numStressDof, numStressDof );
    referenceStabilization.zero();

    for( integer lf = 0; lf < data.numFaces; ++lf )
    {
      FaceGeometry const & geom = data.faceGeom[ static_cast< std::size_t >( lf ) ];
      integer const offset = NUM_FACE_DOF * lf;

      forEachFaceQuadraturePoint( poly, data, lf,
                                  [&]( real64 const (&point)[3], real64 const weight )
      {
        real64 r[3];
        for( integer i = 0; i < 3; ++i )
        {
          r[i] = point[i] - data.center[i];
        }

        real64 phi[NUM_FACE_DOF][3];
        evaluateFaceBasis( geom, point, phi );

        for( integer j = 0; j < NUM_FACE_DOF; ++j )
        {
          for( integer i = 0; i < NUM_RM_DOF; ++i )
          {
            real64 rigidMotion[3];
            evaluateRigidMotion( i, r, rigidMotion );

            referenceDivergence( i, offset + j ) +=
              weight * geom.outwardSign * LvArray::tensorOps::AiBi< 3 >( phi[j], rigidMotion );
          }

          for( integer a = 0; a < NUM_SYM_COMP; ++a )
          {
            real64 field[3];
            evaluateProjectionField( a, r, field );

            referenceProjection( a, offset + j ) +=
              weight * geom.outwardSign * LvArray::tensorOps::AiBi< 3 >( phi[j], field );
          }
        }

        // residual (I - Pi_E) phi_i n on this face, for every element degree of freedom
        std::vector< std::array< real64, 3 > > residual( static_cast< std::size_t >( numStressDof ) );

        for( integer i = 0; i < numStressDof; ++i )
        {
          real64 coefficients[NUM_SYM_COMP];
          for( integer a = 0; a < NUM_SYM_COMP; ++a )
          {
            coefficients[a] = ops.projection( a, i );
          }

          real64 traction[3];
          symBasisTraction( coefficients, geom.normal, traction );

          for( integer p = 0; p < 3; ++p )
          {
            residual[ static_cast< std::size_t >( i ) ][ static_cast< std::size_t >( p ) ] = traction[p];
          }
        }

        for( integer j = 0; j < NUM_FACE_DOF; ++j )
        {
          for( integer p = 0; p < 3; ++p )
          {
            residual[ static_cast< std::size_t >( offset + j ) ][ static_cast< std::size_t >( p ) ] -= phi[j][p];
          }
        }

        for( integer i = 0; i < numStressDof; ++i )
        {
          for( integer j = 0; j < numStressDof; ++j )
          {
            real64 product = 0.0;
            for( integer p = 0; p < 3; ++p )
            {
              product += residual[ static_cast< std::size_t >( i ) ][ static_cast< std::size_t >( p ) ]
                         * residual[ static_cast< std::size_t >( j ) ][ static_cast< std::size_t >( p ) ];
            }
            referenceStabilization( i, j ) += weight * product;
          }
        }
      } );
    }

    // volume part of |E| P_E: -int_E div phi_j . p_a dE, with div phi_j read from D_E
    forEachVolumeQuadraturePoint( poly, data,
                                  [&]( real64 const (&point)[3], real64 const weight )
    {
      real64 r[3];
      for( integer i = 0; i < 3; ++i )
      {
        r[i] = point[i] - data.center[i];
      }

      for( integer j = 0; j < numStressDof; ++j )
      {
        real64 omega[3], divergenceValue[3];
        for( integer i = 0; i < 3; ++i )
        {
          omega[i] = ops.divReconstruction( 3 + i, j );
        }
        LvArray::tensorOps::crossProduct( divergenceValue, omega, r );

        for( integer i = 0; i < 3; ++i )
        {
          divergenceValue[i] += ops.divReconstruction( i, j );
        }

        for( integer a = 0; a < NUM_SYM_COMP; ++a )
        {
          real64 field[3];
          evaluateProjectionField( a, r, field );

          referenceProjection( a, j ) -= weight * LvArray::tensorOps::AiBi< 3 >( divergenceValue, field );
        }
      }
    } );

    for( integer i = 0; i < NUM_RM_DOF; ++i )
    {
      for( integer j = 0; j < numStressDof; ++j )
      {
        EXPECT_NEAR( ops.divergence( i, j ), referenceDivergence( i, j ), 1e-11 * scale )
          << entry.first << " B(" << i << "," << j << ")";
      }
    }

    for( integer a = 0; a < NUM_SYM_COMP; ++a )
    {
      for( integer j = 0; j < numStressDof; ++j )
      {
        EXPECT_NEAR( ops.projection( a, j ), referenceProjection( a, j ) / data.moments.volume, 1e-10 * scale )
          << entry.first << " P(" << a << "," << j << ")";
      }
    }

    real64 const stabScale = compliance[3][3] * data.diameter;

    for( integer i = 0; i < numStressDof; ++i )
    {
      for( integer j = 0; j < numStressDof; ++j )
      {
        real64 consistency = 0.0;
        for( integer a = 0; a < NUM_SYM_COMP; ++a )
        {
          for( integer b = 0; b < NUM_SYM_COMP; ++b )
          {
            consistency += ops.projection( a, i ) * compliance[a][b] * ops.projection( b, j );
          }
        }
        consistency *= data.moments.volume;

        real64 const expected = consistency + stabScale * referenceStabilization( i, j );

        EXPECT_NEAR( ops.stiffness( i, j ), expected, 1e-9 * scale )
          << entry.first << " K(" << i << "," << j << ")";
      }
    }
  }
}

TEST( MixedVEMCompliance, voigtRoundTrip )
{
  real64 const lambda = 1.3;
  real64 const mu = 0.8;

  real64 voigtStiffness[NUM_SYM_COMP][NUM_SYM_COMP] = { { 0.0 } };
  for( integer i = 0; i < 3; ++i )
  {
    for( integer j = 0; j < 3; ++j )
    {
      voigtStiffness[i][j] = lambda + ( ( i == j ) ? 2.0 * mu : 0.0 );
    }
    voigtStiffness[3 + i][3 + i] = mu;
  }

  real64 stiffness[NUM_SYM_COMP][NUM_SYM_COMP];
  convertVoigtStiffness( voigtStiffness, stiffness );

  real64 compliance[NUM_SYM_COMP][NUM_SYM_COMP];
  invertStiffness( stiffness, compliance );

  real64 expected[NUM_SYM_COMP][NUM_SYM_COMP];
  makeIsotropicCompliance( lambda, mu, expected );

  for( integer i = 0; i < NUM_SYM_COMP; ++i )
  {
    for( integer j = 0; j < NUM_SYM_COMP; ++j )
    {
      EXPECT_NEAR( compliance[i][j], expected[i][j], 1e-12 );
    }
  }
}

TEST( MixedVEMAssembly, scatterReproducesDenseBlock )
{
  real64 compliance[NUM_SYM_COMP][NUM_SYM_COMP];
  makeIsotropicCompliance( 1.0, 1.0, compliance );

  Polyhedron const poly = unitCube();
  ElementData const data = buildElement( poly );
  Operators const ops = computeOperators( data, compliance );

  integer const numStressDof = data.numStressDof;
  integer const numDof = numStressDof + NUM_RM_DOF;

  // one face per element face, numbered in order
  array1d< localIndex > elemToFaces( data.numFaces );
  array1d< globalIndex > faceDofNumber( data.numFaces );
  for( integer lf = 0; lf < data.numFaces; ++lf )
  {
    elemToFaces[lf] = lf;
    faceDofNumber[lf] = NUM_FACE_DOF * lf;
  }

  std::vector< globalIndex > stressDofIndices( static_cast< std::size_t >( numStressDof ) );
  arrayView1d< localIndex const > const elemToFacesView = elemToFaces.toViewConst();
  gatherStressDofIndices( elemToFacesView.toSliceConst(),
                          data.numFaces,
                          faceDofNumber.toViewConst(),
                          stressDofIndices.data() );

  for( integer j = 0; j < numStressDof; ++j )
  {
    EXPECT_EQ( stressDofIndices[ static_cast< std::size_t >( j ) ], j );
  }

  globalIndex dispDofIndices[NUM_RM_DOF];
  for( integer k = 0; k < NUM_RM_DOF; ++k )
  {
    dispDofIndices[k] = numStressDof + k;
  }

  CRSMatrix< real64, globalIndex > matrix;
  matrix.resize( numDof, numDof, numDof );

  std::vector< globalIndex > allColumns( static_cast< std::size_t >( numDof ) );
  std::vector< real64 > zeros( static_cast< std::size_t >( numDof ), 0.0 );
  for( integer j = 0; j < numDof; ++j )
  {
    allColumns[ static_cast< std::size_t >( j ) ] = j;
  }
  for( integer i = 0; i < numDof; ++i )
  {
    matrix.insertNonZeros( i, allColumns.data(), zeros.data(), numDof );
  }

  std::vector< real64 > rowBuffer( static_cast< std::size_t >( numStressDof ), 0.0 );

  addElementToMatrix< serialAtomic >( matrix.toViewConstSizes(),
                                      0,
                                      stressDofIndices.data(),
                                      numStressDof,
                                      dispDofIndices,
                                      ops.stiffness.toSliceConst(),
                                      ops.divergence.toSliceConst(),
                                      rowBuffer.data() );

  for( integer i = 0; i < numDof; ++i )
  {
    arraySlice1d< globalIndex const > const columns = matrix.getColumns( i );
    arraySlice1d< real64 const > const values = matrix.getEntries( i );

    for( integer k = 0; k < numDof; ++k )
    {
      integer const j = LvArray::integerConversion< integer >( columns[k] );

      real64 expected = 0.0;
      if( i < numStressDof && j < numStressDof )
      {
        expected = ops.stiffness( i, j );
      }
      else if( i < numStressDof )
      {
        expected = ops.divergence( j - numStressDof, i );
      }
      else if( j < numStressDof )
      {
        // the rigid motion rows carry -B_E
        expected = -ops.divergence( i - numStressDof, j );
      }

      EXPECT_NEAR( values[k], expected, 1e-14 ) << "entry (" << i << "," << j << ")";
    }
  }
}

TEST( MixedVEMAssembly, meshAdaptorMatchesDirectGeometry )
{
  Polyhedron const poly = unitCube();
  ElementData const reference = buildElement( poly );

  integer const numNodes = LvArray::integerConversion< integer >( poly.nodes.size() );
  integer const numFaces = reference.numFaces;

  array2d< real64, nodes::REFERENCE_POSITION_PERM > nodePositions( numNodes, 3 );
  for( integer k = 0; k < numNodes; ++k )
  {
    for( integer i = 0; i < 3; ++i )
    {
      nodePositions( k, i ) = poly.nodes[ static_cast< std::size_t >( k ) ][ static_cast< std::size_t >( i ) ];
    }
  }

  ArrayOfArrays< localIndex > faceToNodes( numFaces, 8 );
  array2d< real64 > faceNormals( numFaces, 3 );
  array1d< localIndex > elemToFaces( numFaces );

  for( integer lf = 0; lf < numFaces; ++lf )
  {
    std::vector< int > const & loop = poly.faces[ static_cast< std::size_t >( lf ) ];
    for( int node : loop )
    {
      faceToNodes.emplaceBack( lf, node );
    }

    real64 normal[3];
    loopNormal( poly, loop, normal );
    for( integer i = 0; i < 3; ++i )
    {
      faceNormals( lf, i ) = normal[i];
    }

    elemToFaces[lf] = lf;
  }

  std::vector< FaceGeometry > faceGeom( static_cast< std::size_t >( numFaces ) );
  ElementMoments moments;

  arrayView1d< localIndex const > const elemToFacesView = elemToFaces.toViewConst();

  buildElementGeometry( nodePositions.toViewConst(),
                        faceToNodes.toViewConst(),
                        faceNormals.toViewConst(),
                        elemToFacesView.toSliceConst(),
                        numFaces,
                        reference.center,
                        faceGeom.data(),
                        moments );

  EXPECT_NEAR( moments.volume, reference.moments.volume, 1e-13 );

  for( integer lf = 0; lf < numFaces; ++lf )
  {
    FaceGeometry const & expected = reference.faceGeom[ static_cast< std::size_t >( lf ) ];
    FaceGeometry const & actual = faceGeom[ static_cast< std::size_t >( lf ) ];

    EXPECT_NEAR( actual.area, expected.area, 1e-13 );
    EXPECT_NEAR( actual.m20, expected.m20, 1e-13 );
    EXPECT_NEAR( actual.m02, expected.m02, 1e-13 );
    EXPECT_NEAR( actual.outwardSign, expected.outwardSign, 1e-13 );

    for( integer i = 0; i < 3; ++i )
    {
      EXPECT_NEAR( actual.center[i], expected.center[i], 1e-13 );
      EXPECT_NEAR( actual.normal[i], expected.normal[i], 1e-13 );
      EXPECT_NEAR( actual.t1[i], expected.t1[i], 1e-13 );
    }
  }
}

TEST( MixedVEMCellOutput, cellFieldsAreWritable )
{
  real64 compliance[NUM_SYM_COMP][NUM_SYM_COMP];
  makeIsotropicCompliance( 1.0, 1.0, compliance );

  real64 const sigma[NUM_SYM_COMP] = { 0.4, -1.1, 0.7, 0.3, -0.6, 0.9 };

  for( auto const & entry : testElements() )
  {
    ElementData const data = buildElement( entry.second );
    Operators const ops = computeOperators( data, compliance );

    std::vector< real64 > const dofs = constantStressDofs( data, sigma );

    // the cell stress is the plain tensor, the orthonormal sqrt(2) removed
    real64 stress[NUM_SYM_COMP];
    computeCellStress( ops.projection.toSliceConst(), data.numFaces, dofs.data(), stress );

    for( integer a = 0; a < NUM_SYM_COMP; ++a )
    {
      real64 const expected = ( a < 3 ) ? sigma[a] : INV_SQRT_2 * sigma[a];
      EXPECT_NEAR( stress[a], expected, 1e-11 ) << entry.first << " component " << a;
    }

    // the cell displacement is u_h(x_E), the rotation is omega
    real64 const rigidMotion[NUM_RM_DOF] = { 0.3, -0.2, 0.5, 0.1, 0.4, -0.3 };

    real64 displacement[3], rotation[3];
    computeCellDisplacement( rigidMotion, displacement, rotation );

    real64 atCenter[3];
    evaluateCellDisplacement( rigidMotion, data.center, data.center, atCenter );

    for( integer i = 0; i < 3; ++i )
    {
      EXPECT_NEAR( displacement[i], rigidMotion[i], 1e-14 );
      EXPECT_NEAR( rotation[i], rigidMotion[3 + i], 1e-14 );
      EXPECT_NEAR( atCenter[i], displacement[i], 1e-14 );
    }

    // away from x_E the rigid motion adds omega ^ (x - x_E)
    real64 const point[3] = { data.center[0] + 0.3, data.center[1] - 0.2, data.center[2] + 0.1 };
    real64 offset[3] = { 0.3, -0.2, 0.1 };
    real64 expectedValue[3];
    LvArray::tensorOps::crossProduct( expectedValue, rotation, offset );

    real64 value[3];
    evaluateCellDisplacement( rigidMotion, data.center, point, value );

    for( integer i = 0; i < 3; ++i )
    {
      EXPECT_NEAR( value[i], displacement[i] + expectedValue[i], 1e-13 );
    }
  }
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
