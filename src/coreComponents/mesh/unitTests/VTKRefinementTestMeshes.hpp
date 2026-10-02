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

/** @file VTKRefinementTestMeshes.hpp */
#ifndef GEOS_VTK_REFINEMENT_TEST_MESHES_HPP
#define GEOS_VTK_REFINEMENT_TEST_MESHES_HPP

#include "../generators/VTKRefinementTemplates.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vtkCellType.h>

namespace geos::vtk::refinement::testMeshes
{
inline Cell regularPrism( int n, std::vector< Coordinates > & xyz )
{
  Cell cell{ VTK_POLYHEDRON, {}, n };
  for( int z = 0; z < 2; ++z )
    for( int i = 0; i < n; ++i )
    {
      double const angle = 2 * std::acos( -1. ) * i / n;
      xyz.push_back( { std::cos( angle ) + 0.3 * z, std::sin( angle ) - 0.15 * z, static_cast< double >( z ) } );
      cell.points.push_back( xyz.size() - 1 );
    }
  return cell;
}

struct ReferenceCell
{
  Cell cell;
  std::vector< Coordinates > xyz;
  int face;
};
inline ReferenceCell referenceCell( int type, int face )
{
  ReferenceCell r{ { type, {}, 0 }, {}, face };
  switch( type )
  {
    case VTK_TETRA:
      r.xyz = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
      break;
    case VTK_PYRAMID:
      r.xyz = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { .5, .5, 1 } };
      break;
    case VTK_WEDGE:
      r.xyz = { { 0, 0, 0 }, { 0, 1, 0 }, { 1, 0, 0 }, { 0, 0, 1 }, { 0, 1, 1 }, { 1, 0, 1 } };
      break;
    case VTK_HEXAHEDRON:
      r.xyz = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 } };
      break;
    default:
      r.cell = regularPrism( type, r.xyz );
      return r;
  }
  r.cell.points.resize( r.xyz.size() );
  std::iota( r.cell.points.begin(), r.cell.points.end(), 0 );
  return r;
}

// Attach a canonical reference face to a common triangle/square by an affine
// map. Only the explicitly corresponding face corners share IDs.
inline Cell attach( ReferenceCell const & r, int sign, std::vector< Coordinates > & xyz )
{
  auto const face = cellFaces( r.cell )[r.face];
  auto const origin = r.xyz[face[0]];
  Coordinates u{}, v{}, normal{};
  for( int d = 0; d < 3; ++d )
  {
    u[d] = r.xyz[face[1]][d] - origin[d];
    v[d] = r.xyz[face.back()][d] - origin[d];
  }
  normal = { u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0] };
  double uu = 0, vv = 0, uv = 0, nn = 0;
  for( int d = 0; d < 3; ++d )
  {
    uu += u[d] * u[d];
    vv += v[d] * v[d];
    uv += u[d] * v[d];
    nn += normal[d] * normal[d];
  }
  Cell result = r.cell;
  for( std::size_t i = 0; i < r.xyz.size(); ++i )
  {
    auto const found = std::find( face.begin(), face.end(), static_cast< vtkIdType >( i ) );
    if( found != face.end() )
    {
      auto const index = std::distance( face.begin(), found );
      result.points[i] = sign > 0 ? index : ( face.size() - index ) % face.size();
    }
    else
    {
      double a = 0, b = 0, z = 0;
      for( int d = 0; d < 3; ++d )
      {
        double const delta = r.xyz[i][d] - origin[d];
        a += delta * u[d];
        b += delta * v[d];
        z += delta * normal[d] / std::sqrt( nn );
      }
      double const x = ( a * vv - b * uv ) / ( uu * vv - uv * uv ), y = ( b * uu - a * uv ) / ( uu * vv - uv * uv );
      result.points[i] = xyz.size();
      xyz.push_back( sign > 0 ? Coordinates{ x, y, z } : Coordinates{ y, x, -z } );
    }
  }
  return result;
}
} // namespace geos::vtk::refinement::testMeshes
#endif
