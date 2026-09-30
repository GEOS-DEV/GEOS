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

/** @file VTKRefinementAssociations.cpp */
#include "VTKRefinementAssociations.hpp"

#include <algorithm>
#include <set>
#include <stdexcept>

namespace geos::vtk::refinement
{
SurfaceAssociations::SurfaceAssociations( std::vector< MainFace > faces, std::uint64_t mainNamespace ) : m_namespace( mainNamespace )
{
  std::unordered_map< EntityKey, std::size_t, EntityKeyHash > unique;
  std::unordered_map< Connectivity, Connectivity, ConnectivityHash > cycles;
  for( auto & face : faces )
  {
    EntityKey const key = entityKey( EntityKind::face, face.globalCorners, m_namespace );
    if( face.globalCorners.size() > 11 || face.owners.empty() )
    {
      throw std::invalid_argument( "Invalid main face association arity or owners" );
    }
    std::sort( face.owners.begin(), face.owners.end() );
    face.owners.erase( std::unique( face.owners.begin(), face.owners.end() ), face.owners.end() );
    if( face.owners.front() < 0 )
    {
      throw std::invalid_argument( "Negative volume-face owner" );
    }
    auto set = face.globalCorners;
    std::sort( set.begin(), set.end() );
    auto const [cycle, insertedCycle] = cycles.emplace( std::move( set ), key.corners );
    if( !insertedCycle && cycle->second != key.corners )
    {
      throw std::invalid_argument( "Crossed main face association cycles" );
    }
    auto const [previous, inserted] = unique.emplace( key, m_faces.size() );
    if( !inserted )
    {
      auto & owners = m_faces[previous->second].owners;
      owners.insert( owners.end(), face.owners.begin(), face.owners.end() );
      std::sort( owners.begin(), owners.end() );
      owners.erase( std::unique( owners.begin(), owners.end() ), owners.end() );
      continue;
    }
    for( vtkIdType id : face.globalCorners )
    {
      m_incident[id].push_back( m_faces.size() );
    }
    m_faces.push_back( std::move( face ) );
  }
}

std::vector< SurfaceSide > SurfaceAssociations::match( Connectivity const & surfaceCorners,
                                                       std::vector< Connectivity > const & collocationBuckets ) const
{
  if( surfaceCorners.size() < 3 || surfaceCorners.size() > 11 )
  {
    throw std::invalid_argument( "Unsupported associated surface polygon" );
  }
  auto unique = surfaceCorners;
  std::sort( unique.begin(), unique.end() );
  if( std::adjacent_find( unique.begin(), unique.end() ) != unique.end() )
  {
    throw std::invalid_argument( "Repeated associated surface corner" );
  }
  for( vtkIdType corner : surfaceCorners )
  {
    if( corner < 0 || static_cast< std::size_t >( corner ) >= collocationBuckets.size() || collocationBuckets[corner].empty() )
    {
      throw std::invalid_argument( "Missing surface collocation bucket" );
    }
    for( vtkIdType id : collocationBuckets[corner] )
    {
      if( id < 0 )
      {
        throw std::invalid_argument( "Unnormalized negative surface collocation entry" );
      }
    }
  }
  std::set< std::size_t > candidates;
  for( vtkIdType id : collocationBuckets[surfaceCorners.front()] )
  {
    auto const found = m_incident.find( id );
    if( found != m_incident.end() )
    {
      candidates.insert( found->second.begin(), found->second.end() );
    }
  }
  std::vector< SurfaceSide > result;
  for( std::size_t index : candidates )
  {
    MainFace const & face = m_faces[index];
    if( face.globalCorners.size() != surfaceCorners.size() )
    {
      continue;
    }
    Connectivity mapped( surfaceCorners.size(), -1 );
    bool complete = true, ambiguous = false;
    for( vtkIdType id : face.globalCorners )
    {
      int matches = 0;
      for( std::size_t i = 0; i < surfaceCorners.size(); ++i )
      {
        auto const & bucket = collocationBuckets[surfaceCorners[i]];
        if( std::find( bucket.begin(), bucket.end(), id ) != bucket.end() )
        {
          ++matches;
          ambiguous = ambiguous || mapped[i] != -1;
          mapped[i] = id;
        }
      }
      complete = complete && matches > 0;
      ambiguous = ambiguous || matches > 1;
    }
    if( !complete )
    {
      continue;
    }
    if( ambiguous || std::find( mapped.begin(), mapped.end(), -1 ) != mapped.end() )
    {
      throw std::invalid_argument( "Ambiguous surface-to-side vertex association" );
    }
    auto const key = entityKey( EntityKind::face, face.globalCorners, m_namespace );
    if( !( entityKey( EntityKind::face, mapped, m_namespace ) == key ) )
    {
      throw std::invalid_argument( "Surface and main face have incompatible edge cycles" );
    }
    result.push_back( { key, std::move( mapped ), face.owners } );
  }
  return result;
}

EntityKey SurfaceAssociations::supportKey( EntityKind kind, Connectivity const & surfaceSupport, Connectivity const & surfaceFace,
                                           SurfaceSide const & side )
{
  if( surfaceFace.size() != side.mainCornersBySurface.size() || kind == EntityKind::cell )
  {
    throw std::invalid_argument( "Invalid surface support projection" );
  }
  Connectivity indices, ids;
  for( vtkIdType corner : surfaceSupport )
  {
    auto const found = std::find( surfaceFace.begin(), surfaceFace.end(), corner );
    if( found == surfaceFace.end() )
    {
      throw std::invalid_argument( "Surface support spans unrelated parent faces" );
    }
    auto const index = static_cast< std::size_t >( found - surfaceFace.begin() );
    indices.push_back( index );
    ids.push_back( side.mainCornersBySurface[index] );
  }
  if( kind == EntityKind::edge )
  {
    if( indices.size() != 2 || ( ( indices[0] + 1 ) % surfaceFace.size() != static_cast< std::size_t >( indices[1] ) &&
                                 ( indices[1] + 1 ) % surfaceFace.size() != static_cast< std::size_t >( indices[0] ) ) )
    {
      throw std::invalid_argument( "Surface support is not an actual parent edge" );
    }
  }
  if( kind == EntityKind::face )
  {
    auto actual = surfaceSupport, expected = surfaceFace;
    std::sort( actual.begin(), actual.end() );
    std::sort( expected.begin(), expected.end() );
    if( actual != expected )
    {
      throw std::invalid_argument( "Incomplete surface face-center support" );
    }
    // Point recipes list affine support in sorted ID order, which need not be
    // an edge cycle. The validated actual side face supplies the cycle.
    return side.mainFace;
  }
  return entityKey( kind, std::move( ids ), side.mainFace.meshNamespace );
}
} // namespace geos::vtk::refinement
