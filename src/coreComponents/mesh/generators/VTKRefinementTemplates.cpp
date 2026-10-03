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
 * @file VTKRefinementTemplates.cpp
 */

#include "VTKRefinementTemplates.hpp"
#include "LvArray/src/tensorOps.hpp"

#include <vtkCell.h>
#include <vtkCellType.h>
#include <vtkPoints.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace geos::vtk::refinement
{
namespace
{
Coordinates subtract( Coordinates const & a, Coordinates const & b )
{
  return { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
}
double determinant( Coordinates const & a, Coordinates const & b, Coordinates const & c )
{
  // The mesh stores coordinates in std::array for VTK interchange. Use the
  // existing LvArray tensor operations through fixed-size stack arrays.
  double const matrix[3][3] = { { a[0], b[0], c[0] }, { a[1], b[1], c[1] }, { a[2], b[2], c[2] } };
  return LvArray::tensorOps::determinant< 3 >( matrix );
}
Coordinates crossProduct( Coordinates const & a, Coordinates const & b )
{
  double const lhs[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( a );
  double const rhs[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( b );
  double result[3];
  LvArray::tensorOps::crossProduct( result, lhs, rhs );
  return { result[0], result[1], result[2] };
}
double normSquared( Coordinates const & coordinates )
{
  double const vector[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( coordinates );
  return LvArray::tensorOps::l2NormSquared< 3 >( vector );
}
double tetraMeasure( Connectivity const & p, PointRegistry const & registry )
{
  Coordinates const & a = registry.position( p[0] );
  return determinant( subtract( registry.position( p[1] ), a ), subtract( registry.position( p[2] ), a ),
                      subtract( registry.position( p[3] ), a ) ) /
         6;
}
std::uint64_t add( std::uint64_t a, std::uint64_t b )
{
  if( b > std::numeric_limits< std::uint64_t >::max() - a )
    throw std::overflow_error( "Uniform refinement cell-count addition overflow" );
  return a + b;
}
std::uint64_t multiply( std::uint64_t a, std::uint64_t b )
{
  if( b && a > std::numeric_limits< std::uint64_t >::max() / b )
    throw std::overflow_error( "Uniform refinement cell-count multiplication overflow" );
  return a * b;
}
Connectivity ids( vtkCell & cell )
{
  Connectivity result;
  for( vtkIdType i = 0; i < cell.GetNumberOfPoints(); ++i )
    result.push_back( cell.GetPointId( i ) );
  return result;
}

// Derivatives of first-order parent maps. Pyramid derivatives omit the positive
// (1-z)^2 collapse factor, so validity can be checked without rejecting its
// apex.
double jacobian( Cell const & cell, PointRegistry const & registry, double x, double y, double z )
{
  std::array< Coordinates, 3 > derivatives{};
  auto accumulate = [&]( int i, double dx, double dy, double dz )
  {
    Coordinates const delta = subtract( registry.position( cell.points[i] ), registry.position( cell.points[0] ) );
    for( int d = 0; d < 3; ++d )
    {
      derivatives[0][d] += dx * delta[d];
      derivatives[1][d] += dy * delta[d];
      derivatives[2][d] += dz * delta[d];
    }
  };
  constexpr int cube[8][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 } };
  if( cell.vtkType == VTK_HEXAHEDRON )
    for( int i = 0; i < 8; ++i )
    {
      double const a = cube[i][0] ? x : 1 - x;
      double const b = cube[i][1] ? y : 1 - y;
      double const c = cube[i][2] ? z : 1 - z;
      accumulate( i, ( cube[i][0] ? 1 : -1 ) * b * c, ( cube[i][1] ? 1 : -1 ) * a * c, ( cube[i][2] ? 1 : -1 ) * a * b );
    }
  else if( cell.vtkType == VTK_WEDGE )
  {
    double const triangle[3] = { 1 - x - y, x, y };
    constexpr int dx[3] = { -1, 1, 0 }, dy[3] = { -1, 0, 1 };
    for( int i = 0; i < 6; ++i )
      accumulate( i, dx[i % 3] * ( i < 3 ? 1 - z : z ), dy[i % 3] * ( i < 3 ? 1 - z : z ), triangle[i % 3] * ( i < 3 ? -1 : 1 ) );
  }
  else if( cell.vtkType == VTK_PYRAMID )
  {
    for( int i = 0; i < 4; ++i )
    {
      double const a = cube[i][0] ? x : 1 - x, b = cube[i][1] ? y : 1 - y;
      accumulate( i, ( cube[i][0] ? 1 : -1 ) * b, ( cube[i][1] ? 1 : -1 ) * a, -a * b );
    }
    accumulate( 4, 0, 0, 1 );
  }
  else
    throw std::invalid_argument( "No refinement Jacobian for this cell type" );
  // GEOS' existing import permutation and the VTK 9.4 fixtures use clockwise
  // lower-cap wedge order. Keep children in that convention, including its
  // sign.
  double const result = determinant( derivatives[0], derivatives[1], derivatives[2] );
  return cell.vtkType == VTK_WEDGE ? -result : result;
}

// A triquadratic determinant is positive everywhere when all its Bernstein
// coefficients are positive. Subdivision tightens the bound for warped cells.
// Unlike sampling alone, a successful return certifies the entire parent map.
void certifyBox( Cell const & cell, PointRegistry const & registry, Coordinates const & low, Coordinates const & high, double tolerance,
                 int depth )
{
  double b[3][3][3];
  for( int i = 0; i < 3; ++i )
    for( int j = 0; j < 3; ++j )
      for( int k = 0; k < 3; ++k )
      {
        b[i][j][k] = jacobian( cell, registry, low[0] + ( high[0] - low[0] ) * i / 2, low[1] + ( high[1] - low[1] ) * j / 2,
                               low[2] + ( high[2] - low[2] ) * k / 2 );
        if( !std::isfinite( b[i][j][k] ) || b[i][j][k] <= tolerance )
          throw std::invalid_argument( "Degenerate or inverted uniform refinement cell map" );
      }
  for( int j = 0; j < 3; ++j )
    for( int k = 0; k < 3; ++k )
      b[1][j][k] = 2 * b[1][j][k] - ( b[0][j][k] + b[2][j][k] ) / 2;
  for( int i = 0; i < 3; ++i )
    for( int k = 0; k < 3; ++k )
      b[i][1][k] = 2 * b[i][1][k] - ( b[i][0][k] + b[i][2][k] ) / 2;
  for( int i = 0; i < 3; ++i )
    for( int j = 0; j < 3; ++j )
      b[i][j][1] = 2 * b[i][j][1] - ( b[i][j][0] + b[i][j][2] ) / 2;
  bool positive = true;
  for( auto const & plane : b )
    for( auto const & row : plane )
      for( double value : row )
        positive = positive && value > tolerance;
  if( positive )
    return;
  if( depth == 8 )
    throw std::invalid_argument( "Cannot certify positivity of warped refinement cell map" );
  Coordinates middle;
  for( int d = 0; d < 3; ++d )
    middle[d] = ( low[d] + high[d] ) / 2;
  for( int octant = 0; octant < 8; ++octant )
  {
    Coordinates a, c;
    for( int d = 0; d < 3; ++d )
    {
      a[d] = octant & ( 1 << d ) ? middle[d] : low[d];
      c[d] = octant & ( 1 << d ) ? high[d] : middle[d];
    }
    certifyBox( cell, registry, a, c, tolerance, depth + 1 );
  }
}

char const * cellTypeName( Cell const & cell )
{
  if( cell.prismSides )
    return "polygonal prism";
  switch( cell.vtkType )
  {
    case VTK_TETRA: return "tetrahedron";
    case VTK_PYRAMID: return "pyramid";
    case VTK_WEDGE: return "wedge";
    case VTK_HEXAHEDRON: return "hexahedron";
    default: return "cell";
  }
}

// Smallest sampled Jacobian, divided by the cube of the cell size. It only
// explains a rejected cell; validateGeometry is the validity proof.
std::string sampledScaledJacobian( Cell const & cell, PointRegistry const & registry )
{
  double scale = 0;
  for( vtkIdType p : cell.points )
    for( double x : subtract( registry.position( p ), registry.position( cell.points[0] ) ) )
      scale = std::max( scale, std::abs( x ) );
  double lowest = std::numeric_limits< double >::infinity();
  try
  {
    int constexpr samples = 8;
    for( int i = 0; i <= samples; ++i )
      for( int j = 0; j <= samples; ++j )
        for( int k = 0; k <= samples; ++k )
        {
          double const x = double( i ) / samples, y = double( j ) / samples, z = double( k ) / samples;
          if( cell.vtkType == VTK_WEDGE && x + y > 1 )
            continue;
          lowest = std::min( lowest, jacobian( cell, registry, x, y, z ) );
        }
  }
  catch( std::exception const & )
  {
    return "not available";
  }
  std::ostringstream text;
  text.precision( 2 );
  text << std::scientific << lowest / ( scale * scale * scale );
  return text.str();
}
} // namespace

std::vector< Connectivity > cellFaces( Cell const & cell )
{
  std::vector< Connectivity > faces;
  if( cell.prismSides )
  {
    int const n = cell.prismSides;
    Connectivity lower, upper;
    for( int i = 0; i < n; ++i )
    {
      lower.push_back( n - 1 - i );
      upper.push_back( n + i );
    }
    faces = { lower, upper };
    for( int i = 0; i < n; ++i )
      faces.push_back( { i, ( i + 1 ) % n, n + ( i + 1 ) % n, n + i } );
  }
  else
    switch( cell.vtkType )
    {
      case VTK_TETRA:
        faces = { { 0, 2, 1 }, { 0, 1, 3 }, { 1, 2, 3 }, { 2, 0, 3 } };
        break;
      case VTK_HEXAHEDRON:
        faces = { { 0, 3, 2, 1 }, { 4, 5, 6, 7 }, { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 0, 4, 7 } };
        break;
      case VTK_WEDGE:
        faces = { { 0, 1, 2 }, { 3, 5, 4 }, { 0, 3, 4, 1 }, { 1, 4, 5, 2 }, { 2, 5, 3, 0 } };
        break;
      case VTK_PYRAMID:
        faces = { { 0, 3, 2, 1 }, { 0, 1, 4 }, { 1, 2, 4 }, { 2, 3, 4 }, { 3, 0, 4 } };
        break;
      default:
        throw std::invalid_argument( "Unsupported uniform refinement volume cell type" );
    }
  for( auto & face : faces )
    for( vtkIdType & i : face )
      i = cell.points.at( i );
  return faces;
}

std::vector< EntitySupport > cellEntitySupports( Cell const & cell, vtkIdType globalCellId, Connectivity const & globalPointIds,
                                                 int localRank, std::uint64_t meshNamespace )
{
  std::map< EntityKey, EntitySupport > entities;
  auto append = [&]( EntityKind kind, Connectivity const & corners )
  {
    Connectivity ids;
    for( vtkIdType corner : corners )
      ids.push_back( globalPointIds.at( corner ) );
    EntityKey key = entityKey( kind, std::move( ids ), meshNamespace );
    entities.emplace( key, EntitySupport{ key, corners, { localRank } } );
  };
  for( vtkIdType point : cell.points )
    append( EntityKind::vertex, { point } );
  for( Connectivity const & face : cellFaces( cell ) )
  {
    append( EntityKind::face, face );
    for( std::size_t i = 0; i < face.size(); ++i )
      append( EntityKind::edge, { face[i], face[( i + 1 ) % face.size()] } );
  }
  EntityKey const key = entityKey( EntityKind::cell, { globalCellId }, meshNamespace );
  entities.emplace( key, EntitySupport{ key, cell.points, { localRank } } );
  std::vector< EntitySupport > result;
  result.reserve( entities.size() );
  for( auto & entry : entities )
    result.push_back( std::move( entry.second ) );
  return result;
}

std::vector< EntitySupport > meshEntitySupports( std::vector< Cell > const & cells, Connectivity const & globalCellIds,
                                                 Connectivity const & globalPointIds, int localRank, std::uint64_t meshNamespace )
{
  if( cells.size() != globalCellIds.size() )
    throw std::invalid_argument( "Refinement cell ID count does not match local volumes" );
  std::unordered_map< EntityKey, std::size_t, EntityKeyHash > indices;
  indices.reserve( cells.size() );
  std::vector< EntitySupport > result;
  for( std::size_t i = 0; i < cells.size(); ++i )
    for( EntitySupport & entity : cellEntitySupports( cells[i], globalCellIds[i], globalPointIds, localRank, meshNamespace ) )
    {
      auto const [previous, inserted] = indices.emplace( entity.key, result.size() );
      if( inserted )
        result.push_back( std::move( entity ) );
      else if( entity.key.kind == EntityKind::cell )
        throw std::invalid_argument( "Duplicate local refinement volume ID" );
      else
      {
        Connectivity corners = entity.localCorners, expected = result[previous->second].localCorners;
        std::sort( corners.begin(), corners.end() );
        std::sort( expected.begin(), expected.end() );
        if( corners != expected )
          throw std::invalid_argument( "Conflicting local refinement entity incidence" );
      }
    }
  return result;
}

Cell normalizeCell( vtkCell & cell )
{
  Cell result{ cell.GetCellType(), ids( cell ), 0 };
  std::set< vtkIdType > unique( result.points.begin(), result.points.end() );
  if( unique.size() != result.points.size() )
    throw std::invalid_argument( "Repeated refinement cell vertex" );
  switch( result.vtkType )
  {
    case VTK_TETRA:
    case VTK_PYRAMID:
    case VTK_WEDGE:
    case VTK_HEXAHEDRON:
      return result;
    case VTK_VOXEL:
      std::swap( result.points[2], result.points[3] );
      std::swap( result.points[6], result.points[7] );
      result.vtkType = VTK_HEXAHEDRON;
      return result;
    case VTK_PENTAGONAL_PRISM:
      result.prismSides = 5;
      return result;
    case VTK_HEXAGONAL_PRISM:
      result.prismSides = 6;
      return result;
    case VTK_POLYHEDRON:
      break;
    default:
      throw std::invalid_argument( "Unsupported uniform refinement cell (including high-order cells)" );
  }

  std::vector< Connectivity > faces;
  std::map< Connectivity, int > edgeCounts;
  std::set< Connectivity > faceSets;
  for( int f = 0; f < cell.GetNumberOfFaces(); ++f )
  {
    Connectivity face = ids( *cell.GetFace( f ) );
    Connectivity const key = canonicalCycle( face );
    if( !faceSets.insert( key ).second )
      throw std::invalid_argument( "Duplicate polyhedron face" );
    for( std::size_t i = 0; i < face.size(); ++i )
    {
      if( !unique.count( face[i] ) )
        throw std::invalid_argument( "Polyhedron face vertex not in cell" );
      Connectivity edge{ face[i], face[( i + 1 ) % face.size()] };
      std::sort( edge.begin(), edge.end() );
      ++edgeCounts[edge];
    }
    faces.push_back( std::move( face ) );
  }
  for( auto const & edge : edgeCounts )
    if( edge.second != 2 )
      throw std::invalid_argument( "Polyhedron is not a closed manifold incidence graph" );
  if( unique.size() + faces.size() != edgeCounts.size() + 2 )
    throw std::invalid_argument( "Unsupported polyhedron topology" );

  auto orient = [&]( Cell & candidate )
  {
    std::map< vtkIdType, Coordinates > xyz;
    for( vtkIdType i = 0; i < cell.GetNumberOfPoints(); ++i )
      std::copy_n( cell.GetPoints()->GetPoint( i ), 3, xyz[cell.GetPointId( i )].begin() );
    // Avoid geometry arrays indexed by mesh-wide local IDs: remap the
    // candidate.
    Cell local = candidate;
    std::vector< Coordinates > compact;
    for( vtkIdType & p : local.points )
    {
      compact.push_back( xyz[p] );
      p = compact.size() - 1;
    }
    Connectivity gids( compact.size() );
    for( std::size_t i = 0; i < gids.size(); ++i )
      gids[i] = i;
    PointRegistry registry( std::move( compact ), std::move( gids ) );
    if( signedMeasure( local, registry ) < 0 )
    {
      if( candidate.vtkType == VTK_TETRA )
        std::swap( candidate.points[1], candidate.points[2] );
      else
      {
        int const n = candidate.vtkType == VTK_PYRAMID ? 4 : static_cast< int >( candidate.points.size() / 2 );
        std::reverse( candidate.points.begin() + 1, candidate.points.begin() + n );
        if( candidate.vtkType != VTK_PYRAMID )
          std::reverse( candidate.points.begin() + n + 1, candidate.points.end() );
      }
    }
  };
  if( unique.size() == 4 && faces.size() == 4 && edgeCounts.size() == 6 )
  {
    for( auto const & face : faces )
      if( face.size() != 3 )
        throw std::invalid_argument( "Malformed polyhedral tetrahedron" );
    result.vtkType = VTK_TETRA;
    orient( result );
    return result;
  }
  if( unique.size() == 5 && faces.size() == 5 && edgeCounts.size() == 8 )
  {
    auto base = std::find_if( faces.begin(), faces.end(), []( auto const & f ) { return f.size() == 4; } );
    if( base == faces.end() )
      throw std::invalid_argument( "Malformed polyhedral pyramid" );
    result.points = *base;
    for( vtkIdType p : unique )
      if( std::find( base->begin(), base->end(), p ) == base->end() )
        result.points.push_back( p );
    result.vtkType = VTK_PYRAMID;
    auto expected = cellFaces( result );
    for( auto const & f : expected )
      if( !faceSets.count( canonicalCycle( f ) ) )
        throw std::invalid_argument( "Malformed polyhedral pyramid incidence" );
    orient( result );
    return result;
  }
  int const n = static_cast< int >( unique.size() / 2 );
  if( unique.size() % 2 || n < 3 || n > 11 || faces.size() != static_cast< std::size_t >( n + 2 ) )
    throw std::invalid_argument( "Unrecognized polyhedron; uniform refinement "
                                 "supports canonical shapes only" );
  for( auto const & base : faces )
  {
    if( base.size() != static_cast< std::size_t >( n ) )
      continue;
    Connectivity top;
    bool valid = true;
    for( vtkIdType p : base )
    {
      Connectivity partners;
      for( auto const & edge : edgeCounts )
        if( edge.first[0] == p || edge.first[1] == p )
        {
          vtkIdType const q = edge.first[0] == p ? edge.first[1] : edge.first[0];
          if( std::find( base.begin(), base.end(), q ) == base.end() )
            partners.push_back( q );
        }
      if( partners.size() != 1 )
      {
        valid = false;
        break;
      }
      top.push_back( partners[0] );
    }
    if( !valid || std::set< vtkIdType >( top.begin(), top.end() ).size() != top.size() )
      continue;
    result.points = base;
    result.points.insert( result.points.end(), top.begin(), top.end() );
    result.prismSides = n;
    auto expected = cellFaces( result );
    for( auto const & f : expected )
      valid = valid && faceSets.count( canonicalCycle( f ) );
    if( !valid )
      continue;
    result.prismSides = n > 4 ? n : 0;
    result.vtkType = n == 3 ? VTK_WEDGE : n == 4 ? VTK_HEXAHEDRON : VTK_POLYHEDRON;
    orient( result );
    return result;
  }
  throw std::invalid_argument( "Invalid polyhedral prism cap correspondence" );
}

Cell normalizeCoarseCell( vtkCell & input, PointRegistry const & registry )
{
  Cell original = normalizeCell( input );
  auto id = [&]( vtkIdType point )
  {
    if( point < 0 || point >= registry.originalSize() )
      throw std::invalid_argument( "Coarse reference frame requires original input vertices" );
    return registry.points()[point].key.corners.front();
  };
  for( vtkIdType point : original.points )
    id( point );
  if( original.vtkType == VTK_TETRA )
  {
    Cell canonical = original;
    std::sort( canonical.points.begin(), canonical.points.end(), [&]( vtkIdType a, vtkIdType b ) { return id( a ) < id( b ); } );
    // Restrict to even permutations of the input frame. This preserves an
    // inverted native parent's sign so normal geometry validation rejects it.
    int inversions = 0;
    for( std::size_t i = 0; i < canonical.points.size(); ++i )
      for( std::size_t j = i + 1; j < canonical.points.size(); ++j )
        inversions += std::find( original.points.begin(), original.points.end(), canonical.points[i] ) >
                      std::find( original.points.begin(), original.points.end(), canonical.points[j] );
    if( inversions % 2 )
      std::swap( canonical.points[2], canonical.points[3] );
    return canonical;
  }
  int const n = original.prismSides ? original.prismSides : ( original.vtkType == VTK_WEDGE ? 3 : 4 );
  auto const faces = cellFaces( original );
  std::map< vtkIdType, std::set< vtkIdType > > adjacent;
  for( auto const & face : faces )
    for( std::size_t i = 0; i < face.size(); ++i )
    {
      vtkIdType const a = face[i], b = face[( i + 1 ) % face.size()];
      adjacent[a].insert( b );
      adjacent[b].insert( a );
    }
  Cell best;
  bool found = false;
  for( auto const & face : faces )
  {
    if( face.size() != static_cast< std::size_t >( n ) )
      continue;
    for( int direction : { -1, 1 } )
      for( int start = 0; start < n; ++start )
      {
        Cell candidate = original;
        candidate.points.clear();
        for( int i = 0; i < n; ++i )
          candidate.points.push_back( face[( start + ( direction == 1 ? i : n - i ) ) % n] );
        if( original.vtkType == VTK_PYRAMID )
          candidate.points.push_back( original.points[4] );
        else
          for( int i = 0; i < n; ++i )
          {
            Connectivity partners;
            for( vtkIdType neighbor : adjacent.at( candidate.points[i] ) )
              if( std::find( face.begin(), face.end(), neighbor ) == face.end() )
                partners.push_back( neighbor );
            if( partners.size() != 1 )
              throw std::invalid_argument( "Invalid coarse reference-frame cap correspondence" );
            candidate.points.push_back( partners.front() );
          }
        // Matching outward cycles selects the orientation-preserving shape
        // automorphisms, without repeated Jacobian integrations per candidate.
        auto const lower = cellFaces( candidate ).front();
        auto const begin = std::find( face.begin(), face.end(), lower.front() );
        bool sameOrientation = true;
        for( int i = 0; i < n; ++i )
          sameOrientation = sameOrientation && lower[i] == face[( std::distance( face.begin(), begin ) + i ) % n];
        if( !sameOrientation )
          continue;
        if( !found || std::lexicographical_compare( candidate.points.begin(), candidate.points.end(), best.points.begin(),
                                                    best.points.end(), [&]( vtkIdType a, vtkIdType b ) { return id( a ) < id( b ); } ) )
        {
          best = std::move( candidate );
          found = true;
        }
      }
  }
  if( !found )
    throw std::invalid_argument( "No orientation-preserving coarse reference frame" );
  return best;
}

std::vector< Connectivity > subdivideFace( Connectivity const & p, PointRegistry & registry )
{
  canonicalCycle( p ); // Check local face incidence even for triangles.
  Connectivity e;
  for( std::size_t i = 0; i < p.size(); ++i )
    e.push_back( registry.edge( p[i], p[( i + 1 ) % p.size()] ) );
  if( p.size() == 3 )
    return { { p[0], e[0], e[2] }, { e[0], p[1], e[1] }, { e[2], e[1], p[2] }, { e[0], e[1], e[2] } };
  vtkIdType const c = registry.face( p );
  std::vector< Connectivity > children;
  for( std::size_t i = 0; i < p.size(); ++i )
    children.push_back( { p[i], e[i], c, e[( i + p.size() - 1 ) % p.size()] } );
  return children;
}

Subdivision subdivideCell( Cell const & cell, vtkIdType globalCellId, PointRegistry & registry, bool includeFaceChildren )
{
  try
  {
    validateGeometry( cell, registry );
  }
  catch( std::exception const & error )
  {
    throw std::invalid_argument( std::string( "the coarse " ) + cellTypeName( cell ) + " is degenerate or inverted (" + error.what() +
                                 "; minimum sampled scaled Jacobian " + sampledScaledJacobian( cell, registry ) + ")" );
  }
  Subdivision result;
  result.children.reserve( refinedCellCount( cell, 1 ) );
  if( includeFaceChildren )
    for( auto const & face : cellFaces( cell ) )
      result.faceChildren.push_back( subdivideFace( face, registry ) );
  auto const & p = cell.points;
  auto edge = [&]( int a, int b ) { return registry.edge( p[a], p[b] ); };
  auto child = [&]( int type, Connectivity corners ) { result.children.push_back( { type, std::move( corners ), 0 } ); };
  if( cell.prismSides )
  {
    int const n = cell.prismSides;
    Connectivity bottom( p.begin(), p.begin() + n ), top( p.begin() + n, p.end() );
    vtkIdType const B = registry.face( bottom ), T = registry.face( top ), C = registry.cell( globalCellId, p );
    Connectivity eb, et, m, q;
    for( int i = 0; i < n; ++i )
    {
      eb.push_back( edge( i, ( i + 1 ) % n ) );
      et.push_back( edge( n + i, n + ( i + 1 ) % n ) );
      m.push_back( edge( i, n + i ) );
      q.push_back( registry.face( { p[i], p[( i + 1 ) % n], p[n + ( i + 1 ) % n], p[n + i] } ) );
    }
    for( int i = 0; i < n; ++i )
    {
      int const j = ( i + n - 1 ) % n;
      child( VTK_HEXAHEDRON, { p[i], eb[i], B, eb[j], m[i], q[i], C, q[j] } );
      child( VTK_HEXAHEDRON, { m[i], q[i], C, q[j], p[n + i], et[i], T, et[j] } );
    }
  }
  else if( cell.vtkType == VTK_HEXAHEDRON )
  {
    constexpr int corners[8][3] = { { 0, 0, 0 }, { 2, 0, 0 }, { 2, 2, 0 }, { 0, 2, 0 },
      { 0, 0, 2 }, { 2, 0, 2 }, { 2, 2, 2 }, { 0, 2, 2 } };
    vtkIdType lattice[3][3][3];
    for( int x = 0; x < 3; ++x )
      for( int y = 0; y < 3; ++y )
        for( int z = 0; z < 3; ++z )
        {
          Connectivity support;
          for( int i = 0; i < 8; ++i )
            if( ( x == 1 || x == corners[i][0] ) && ( y == 1 || y == corners[i][1] ) && ( z == 1 || z == corners[i][2] ) )
              support.push_back( p[i] );
          if( support.size() == 1 )
            lattice[x][y][z] = support[0];
          else if( support.size() == 2 )
            lattice[x][y][z] = registry.edge( support[0], support[1] );
          else if( support.size() == 4 )
          {
            // Input VTK vertex order is not a cyclic order on every face.
            for( auto const & f : cellFaces( cell ) )
              if( std::all_of( f.begin(), f.end(),
                               [&]( vtkIdType a ) { return std::find( support.begin(), support.end(), a ) != support.end(); } ) )
              {
                lattice[x][y][z] = registry.face( f );
                break;
              }
          }
          else
            lattice[x][y][z] = registry.cell( globalCellId, p );
        }
    for( int x = 0; x < 2; ++x )
      for( int y = 0; y < 2; ++y )
        for( int z = 0; z < 2; ++z )
          child( VTK_HEXAHEDRON,
                 { lattice[x][y][z], lattice[x + 1][y][z], lattice[x + 1][y + 1][z], lattice[x][y + 1][z], lattice[x][y][z + 1],
                   lattice[x + 1][y][z + 1], lattice[x + 1][y + 1][z + 1], lattice[x][y + 1][z + 1] } );
  }
  else if( cell.vtkType == VTK_PYRAMID )
  {
    vtkIdType e[4], s[4];
    for( int i = 0; i < 4; ++i )
    {
      e[i] = edge( i, ( i + 1 ) % 4 );
      s[i] = edge( i, 4 );
    }
    vtkIdType const c = registry.face( { p[0], p[1], p[2], p[3] } );
    for( int i = 0; i < 4; ++i )
      child( VTK_PYRAMID, { p[i], e[i], c, e[( i + 3 ) % 4], s[i] } );
    child( VTK_PYRAMID, { s[0], s[1], s[2], s[3], p[4] } );
    child( VTK_PYRAMID, { s[3], s[2], s[1], s[0], c } );
    for( int i = 0; i < 4; ++i )
      child( VTK_TETRA, { e[i], s[i], s[( i + 1 ) % 4], c } );
  }
  else if( cell.vtkType == VTK_WEDGE )
  {
    vtkIdType layers[3][6];
    for( int z = 0; z < 3; ++z )
    {
      for( int i = 0; i < 3; ++i )
        layers[z][i] = z == 1 ? edge( i, i + 3 ) : p[i + 3 * ( z / 2 )];
      for( int i = 0; i < 3; ++i )
        layers[z][3 + i] = z == 1 ? registry.face( { p[i], p[( i + 1 ) % 3], p[3 + ( i + 1 ) % 3], p[3 + i] } )
                                  : edge( i + 3 * ( z / 2 ), ( i + 1 ) % 3 + 3 * ( z / 2 ) );
    }
    constexpr int triangles[4][3] = { { 0, 3, 5 }, { 3, 1, 4 }, { 5, 4, 2 }, { 3, 4, 5 } };
    for( int z = 0; z < 2; ++z )
      for( auto const & t : triangles )
        child( VTK_WEDGE,
               { layers[z][t[0]], layers[z][t[1]], layers[z][t[2]], layers[z + 1][t[0]], layers[z + 1][t[1]], layers[z + 1][t[2]] } );
  }
  else if( cell.vtkType == VTK_TETRA )
  {
    vtkIdType const a = edge( 0, 1 ), b = edge( 0, 2 ), c = edge( 0, 3 ), d = edge( 1, 2 ), e = edge( 1, 3 ), f = edge( 2, 3 );
    child( VTK_TETRA, { p[0], a, b, c } );
    child( VTK_TETRA, { a, p[1], d, e } );
    child( VTK_TETRA, { b, d, p[2], f } );
    child( VTK_TETRA, { c, e, f, p[3] } );
    std::array< Connectivity, 3 > diagonals{ Connectivity{ a, f, b, c, e, d }, Connectivity{ b, e, a, c, f, d },
                                             Connectivity{ c, d, a, b, f, e } };
    double best = -1;
    std::vector< Cell > interior;
    for( auto const & diag : diagonals )
    {
      std::vector< Cell > candidate;
      double quality = std::numeric_limits< double >::max();
      for( int i = 0; i < 4; ++i )
      {
        Cell t{ VTK_TETRA, { diag[0], diag[1], diag[2 + i], diag[2 + ( i + 1 ) % 4] }, 0 };
        if( tetraMeasure( t.points, registry ) < 0 )
          std::swap( t.points[0], t.points[1] );
        double longest = 0;
        for( int j = 0; j < 4; ++j )
          for( int k = 0; k < j; ++k )
          {
            Coordinates delta = subtract( registry.position( t.points[j] ), registry.position( t.points[k] ) );
            longest = std::max( longest, normSquared( delta ) );
          }
        quality = std::min( quality, tetraMeasure( t.points, registry ) / std::pow( longest, 1.5 ) );
        candidate.push_back( std::move( t ) );
      }
      // Tie-breaking uses inherited VTK template ordering, never allocated IDs.
      if( quality > best + 64 * std::numeric_limits< double >::epsilon() )
      {
        best = quality;
        interior = std::move( candidate );
      }
    }
    result.children.insert( result.children.end(), interior.begin(), interior.end() );
  }
  else
    throw std::invalid_argument( "Unsupported volume refinement template" );
  for( std::size_t k = 0; k < result.children.size(); ++k )
  {
    Cell const & c = result.children[k];
    try
    {
      validateGeometry( c, registry );
    }
    catch( std::exception const & error )
    {
      std::string message = std::string( "the parent " ) + cellTypeName( cell ) + " is valid, but child " + std::to_string( k ) + " (a " +
                            cellTypeName( c ) + ") of its refinement template is degenerate or inverted (" + error.what() +
                            "). Minimum sampled scaled Jacobian: parent " + sampledScaledJacobian( cell, registry ) + ", child " +
                            sampledScaledJacobian( c, registry ) + ".";
      // The pyramid template's inverted center child has its apex at the base
      // center and its base at the midpoints of the lateral edges.
      if( cell.vtkType == VTK_PYRAMID && c.vtkType == VTK_PYRAMID && k == 5 )
        message += " The pyramid is too flat for its warped base: the center of its base lies above the midpoints of its lateral edges.";
      throw std::invalid_argument( message );
    }
  }
  double sum = 0;
  for( Cell const & c : result.children )
    sum += signedMeasure( c, registry );
  double const parent = signedMeasure( cell, registry );
  // Midpoint coordinates are averaged at the mesh's physical scale. On large
  // translated meshes, summing child tetrahedra can accumulate roundoff above
  // the previous 1e-10 relative threshold while preserving the parent volume.
  if( std::abs( sum - parent ) > 1e-6 * parent )
    throw std::invalid_argument( "Refinement does not preserve signed parent measure" );
  return result;
}

double signedMeasure( Cell const & cell, PointRegistry const & registry )
{
  if( cell.vtkType == VTK_TETRA )
    return tetraMeasure( cell.points, registry );
  if( cell.prismSides )
  {
    // The polygon-cap geometry is its common center/quad fan, including
    // mildly warped caps from legacy files. Integrate exactly those bilinear
    // patches, rather than choosing a triangle diagonal on each rank.
    Coordinates const origin = registry.position( cell.points[0] );
    double volume = 0;
    constexpr double g[2] = { 0.21132486540518713, 0.78867513459481287 };
    auto integrateQuad = [&]( std::array< Coordinates, 4 > const & p )
    {
      for( double x : g )
        for( double y : g )
        {
          Coordinates a{}, dx{}, dy{};
          double const w[4] = { ( 1 - x ) * ( 1 - y ), x * ( 1 - y ), x * y, ( 1 - x ) * y };
          double const wx[4] = { -( 1 - y ), 1 - y, y, -y }, wy[4] = { -( 1 - x ), -x, x, 1 - x };
          for( int i = 0; i < 4; ++i )
            for( int d = 0; d < 3; ++d )
            {
              a[d] += w[i] * p[i][d];
              dx[d] += wx[i] * p[i][d];
              dy[d] += wy[i] * p[i][d];
            }
          volume += determinant( a, dx, dy ) / 12;
        }
    };
    for( auto const & face : cellFaces( cell ) )
    {
      std::vector< Coordinates > p;
      Coordinates center{};
      for( vtkIdType id : face )
      {
        p.push_back( subtract( registry.position( id ), origin ) );
        for( int d = 0; d < 3; ++d )
          center[d] += p.back()[d] / face.size();
      }
      if( face.size() == 4 )
        integrateQuad( { p[0], p[1], p[2], p[3] } );
      else
        for( std::size_t i = 0; i < face.size(); ++i )
        {
          Coordinates next{}, previous{};
          for( int d = 0; d < 3; ++d )
          {
            next[d] = ( p[i][d] + p[( i + 1 ) % p.size()][d] ) / 2;
            previous[d] = ( p[i][d] + p[( i + p.size() - 1 ) % p.size()][d] ) / 2;
          }
          integrateQuad( { p[i], next, center, previous } );
        }
    }
    return volume;
  }
  constexpr double g[2] = { 0.21132486540518713, 0.78867513459481287 };
  double volume = 0;
  if( cell.vtkType == VTK_WEDGE )
    for( double z : g )
      volume += jacobian( cell, registry, 1. / 3, 1. / 3, z ) / 4;
  else if( cell.vtkType == VTK_PYRAMID )
    for( double x : g )
      for( double y : g )
        volume += jacobian( cell, registry, x, y, 0 ) / 12;
  else if( cell.vtkType == VTK_HEXAHEDRON )
    for( double x : g )
      for( double y : g )
        for( double z : g )
          volume += jacobian( cell, registry, x, y, z ) / 8;
  else
    throw std::invalid_argument( "Unsupported signed refinement measure" );
  return volume;
}

void validateGeometry( Cell const & cell, PointRegistry const & registry )
{
  double scale = 0;
  for( vtkIdType p : cell.points )
    for( double x : subtract( registry.position( p ), registry.position( cell.points[0] ) ) )
      scale = std::max( scale, std::abs( x ) );
  double const tolerance = 128 * std::numeric_limits< double >::epsilon() * scale * scale * scale;
  double const measure = signedMeasure( cell, registry );
  if( !std::isfinite( measure ) || measure <= tolerance )
    throw std::invalid_argument( "Degenerate or inverted refinement parent/child" );
  if( cell.vtkType == VTK_HEXAHEDRON || cell.vtkType == VTK_PYRAMID )
    certifyBox( cell, registry, { 0, 0, 0 }, { 1, 1, 1 }, tolerance, 0 );
  if( cell.vtkType == VTK_WEDGE )
  {
    // The determinant is affine on the triangle and quadratic along extrusion.
    for( auto const & xy : std::array< std::array< double, 2 >, 3 >{ { { 0, 0 }, { 1, 0 }, { 0, 1 } } } )
      certifyBox( cell, registry, { xy[0], xy[1], 0 }, { xy[0], xy[1], 1 }, tolerance, 0 );
  }
  if( cell.prismSides )
  {
    int const n = cell.prismSides;
    for( int cap = 0; cap < 2; ++cap )
    {
      Coordinates center{};
      Coordinates const origin = registry.position( cell.points[n * cap] );
      for( int i = 0; i < n; ++i )
        for( int d = 0; d < 3; ++d )
          center[d] += ( registry.position( cell.points[n * cap + i] )[d] - origin[d] ) / n;
      Coordinates normal{};
      for( int i = 0; i < n; ++i )
      {
        Coordinates a = subtract( registry.position( cell.points[n * cap + i] ), origin );
        Coordinates b = subtract( registry.position( cell.points[n * cap + ( i + 1 ) % n] ), origin );
        Coordinates const faceNormal = crossProduct( a, b );
        for( int d = 0; d < 3; ++d )
          normal[d] += faceNormal[d];
      }
      double const norm = std::sqrt( normSquared( normal ) );
      if( norm <= 128 * std::numeric_limits< double >::epsilon() * scale * scale )
        throw std::invalid_argument( "Degenerate prism cap" );
      double winding = 0;
      for( int i = 0; i < n; ++i )
      {
        Coordinates a = subtract( registry.position( cell.points[n * cap + i] ), origin );
        Coordinates b = subtract( registry.position( cell.points[n * cap + ( i + 1 ) % n] ), origin );
        Coordinates const ca = subtract( a, center ), cb = subtract( b, center );
        double const cross = determinant( ca, cb, normal );
        if( cross <= tolerance * norm / scale )
          throw std::invalid_argument( "Prism cap center is outside its admissible oriented fan" );
        double dot = 0, na = 0, nb = 0;
        for( int d = 0; d < 3; ++d )
        {
          dot += ca[d] * cb[d];
          na += ca[d] * normal[d];
          nb += cb[d] * normal[d];
        }
        dot -= na * nb / ( norm * norm );
        winding += std::atan2( cross / norm, dot );
      }
      if( std::abs( winding - 2 * std::acos( -1. ) ) > 1e-8 )
        throw std::invalid_argument( "Self-intersecting or overlapping polygonal prism cap fan" );
    }
  }
}

void CellCounts::include( Cell const & cell )
{
  if( cell.prismSides )
  {
    if( cell.prismSides < 5 || cell.prismSides > 11 )
      throw std::invalid_argument( "Unsupported prism arity in refinement growth counts" );
    prisms[cell.prismSides] = add( prisms[cell.prismSides], 1 );
    return;
  }
  switch( cell.vtkType )
  {
    case VTK_HEXAHEDRON: hexahedra = add( hexahedra, 1 ); break;
    case VTK_TETRA: tetrahedra = add( tetrahedra, 1 ); break;
    case VTK_WEDGE: wedges = add( wedges, 1 ); break;
    case VTK_PYRAMID: pyramids = add( pyramids, 1 ); break;
    default: throw std::invalid_argument( "Unsupported uniform refinement count type" );
  }
}

std::uint64_t refinedCellCount( Cell const & cell, int levels )
{
  if( levels < 0 )
    throw std::invalid_argument( "Negative uniform refinement level" );
  if( levels == 0 )
    return 1;
  if( cell.vtkType == VTK_TRIANGLE || cell.vtkType == VTK_QUAD || cell.vtkType == VTK_POLYGON )
  {
    if( cell.vtkType == VTK_POLYGON && cell.points.size() < 3 )
      throw std::invalid_argument( "Refinement polygon requires at least three corners" );
    std::uint64_t count = cell.vtkType == VTK_POLYGON && cell.points.size() > 4 ? cell.points.size() : 4;
    for( int l = 1; l < levels; ++l )
      count = multiply( count, 4 );
    return count;
  }
  CellCounts counts;
  counts.include( cell );
  for( int l = 0; l < levels; ++l )
    counts = counts.next();
  return counts.total();
}

std::uint64_t CellCounts::total() const
{
  std::uint64_t result = add( add( hexahedra, tetrahedra ), add( wedges, pyramids ) );
  for( auto count : prisms )
    result = add( result, count );
  return result;
}
CellCounts CellCounts::next() const
{
  for( int n = 0; n < 5; ++n )
    if( prisms[n] )
      throw std::invalid_argument( "Unsupported prism arity in refinement growth counts" );
  CellCounts result;
  result.hexahedra = multiply( hexahedra, 8 );
  result.wedges = multiply( wedges, 8 );
  result.tetrahedra = add( multiply( tetrahedra, 8 ), multiply( pyramids, 4 ) );
  result.pyramids = multiply( pyramids, 6 );
  for( int n = 5; n <= 11; ++n )
    result.hexahedra = add( result.hexahedra, multiply( prisms[n], 2 * n ) );
  result.total();
  return result;
}
} // namespace geos::vtk::refinement
