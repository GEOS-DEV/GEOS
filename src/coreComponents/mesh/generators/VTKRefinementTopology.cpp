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
 * @file VTKRefinementTopology.cpp
 */

#include "VTKRefinementTopology.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>

namespace geos::vtk::refinement
{

bool EntityKey::operator<( EntityKey const & other ) const
{
  return std::tie( meshNamespace, kind, corners ) < std::tie( other.meshNamespace, other.kind, other.corners );
}

bool EntityKey::operator==( EntityKey const & other ) const
{
  return meshNamespace == other.meshNamespace && kind == other.kind && corners == other.corners;
}

Connectivity canonicalCycle( Connectivity const & corners )
{
  Connectivity sorted = corners;
  std::sort( sorted.begin(), sorted.end() );
  if( corners.size() < 3 || std::adjacent_find( sorted.begin(), sorted.end() ) != sorted.end() )
    throw std::invalid_argument( "Refinement face requires at least three distinct corners" );
  Connectivity best = corners;
  for( int direction : { -1, 1 } )
    for( std::size_t start = 0; start < corners.size(); ++start )
    {
      Connectivity candidate;
      for( std::size_t i = 0; i < corners.size(); ++i )
        candidate.push_back( corners[( start + ( direction == 1 ? i : corners.size() - i ) ) % corners.size()] );
      best = std::min( best, candidate );
    }
  return best;
}

std::uint64_t stableHash( EntityKey const & key )
{
  std::uint64_t hash = UINT64_C( 14695981039346656037 );
  auto append = [&]( std::uint64_t value )
  {
    for( int i = 0; i < 8; ++i )
    {
      hash ^= ( value >> ( 8 * i ) ) & 255;
      hash *= UINT64_C( 1099511628211 );
    }
  };
  append( key.meshNamespace );
  append( static_cast< std::uint64_t >( key.kind ) );
  append( key.corners.size() );
  for( vtkIdType id : key.corners )
    append( static_cast< std::uint64_t >( id ) );
  return hash;
}

std::size_t ConnectivityHash::operator()( Connectivity const & corners ) const
{
  std::uint64_t hash = UINT64_C( 14695981039346656037 );
  for( vtkIdType id : corners )
    for( int i = 0; i < 8; ++i )
    {
      hash ^= ( static_cast< std::uint64_t >( id ) >> ( 8 * i ) ) & 255;
      hash *= UINT64_C( 1099511628211 );
    }
  return hash;
}

EntityKey entityKey( EntityKind kind, Connectivity corners, std::uint64_t meshNamespace )
{
  if( kind != EntityKind::vertex && kind != EntityKind::edge && kind != EntityKind::face && kind != EntityKind::cell )
    throw std::invalid_argument( "Unknown refinement entity kind" );
  for( vtkIdType corner : corners )
    if( corner < 0 )
      throw std::invalid_argument( "Negative refinement entity corner ID" );
  if( kind == EntityKind::face )
    corners = canonicalCycle( corners );
  else
  {
    if( ( kind == EntityKind::vertex && corners.size() != 1 ) ||
        ( kind == EntityKind::edge && ( corners.size() != 2 || corners[0] == corners[1] ) ) ||
        ( kind == EntityKind::cell && corners.size() != 1 ) )
      throw std::invalid_argument( "Invalid refinement entity key arity" );
    std::sort( corners.begin(), corners.end() );
  }
  return { meshNamespace, kind, std::move( corners ) };
}

PointRegistry::PointRegistry( std::vector< Coordinates > coordinates, Connectivity globalIds, std::uint64_t meshNamespace )
  : m_namespace( meshNamespace ), m_globalIds( std::move( globalIds ) )
{
  if( coordinates.size() != m_globalIds.size() )
    throw std::invalid_argument( "Refinement coordinates/global ID tuple counts differ" );
  if( coordinates.size() > static_cast< std::size_t >( std::numeric_limits< vtkIdType >::max() ) )
    throw std::overflow_error( "Refinement input point count exceeds vtkIdType" );
  m_points.reserve( coordinates.size() );
  m_indices.reserve( coordinates.size() );
  for( std::size_t i = 0; i < coordinates.size(); ++i )
  {
    if( m_globalIds[i] < 0 )
      throw std::invalid_argument( "Negative refinement point global ID" );
    for( double x : coordinates[i] )
      if( !std::isfinite( x ) )
        throw std::invalid_argument( "Nonfinite refinement point coordinates" );
    EntityKey key{ m_namespace, EntityKind::vertex, { m_globalIds[i] } };
    if( !m_indices.emplace( key, i ).second )
      throw std::invalid_argument( "Duplicate local refinement point global ID" );
    m_points.push_back( { std::move( key ), { static_cast< vtkIdType >( i ) }, coordinates[i] } );
  }
}

vtkIdType PointRegistry::insert( EntityKey key, Connectivity corners )
{
  auto found = m_indices.find( key );
  if( found != m_indices.end() )
    return found->second;
  if( corners.empty() )
    throw std::invalid_argument( "Empty refinement point support" );
  std::sort( corners.begin(), corners.end(), [&]( vtkIdType a, vtkIdType b ) { return m_globalIds.at( a ) < m_globalIds.at( b ); } );
  Coordinates position{};
  // Shift by the first corner before averaging to avoid overflow/cancellation
  // in translated meshes. Every participant evaluates the same canonical order.
  Coordinates const & origin = m_points.at( corners.front() ).position;
  for( vtkIdType corner : corners )
    for( int d = 0; d < 3; ++d )
      position[d] += ( m_points.at( corner ).position[d] - origin[d] ) / corners.size();
  for( int d = 0; d < 3; ++d )
  {
    position[d] += origin[d];
    if( !std::isfinite( position[d] ) )
      throw std::overflow_error( "Nonfinite refinement point recipe" );
  }
  if( m_points.size() >= static_cast< std::size_t >( std::numeric_limits< vtkIdType >::max() ) )
    throw std::overflow_error( "Refinement point count exceeds vtkIdType" );
  vtkIdType const result = m_points.size();
  m_points.push_back( { key, std::move( corners ), position } );
  m_indices.emplace( std::move( key ), result );
  return result;
}

vtkIdType PointRegistry::edge( vtkIdType a, vtkIdType b )
{
  if( a == b )
    throw std::invalid_argument( "Repeated refinement edge endpoint" );
  Connectivity ids{ m_globalIds.at( a ), m_globalIds.at( b ) };
  std::sort( ids.begin(), ids.end() );
  return insert( { m_namespace, EntityKind::edge, std::move( ids ) }, { a, b } );
}

vtkIdType PointRegistry::face( Connectivity const & corners )
{
  Connectivity ids;
  for( vtkIdType corner : corners )
    ids.push_back( m_globalIds.at( corner ) );
  Connectivity cycle = canonicalCycle( ids );
  std::sort( ids.begin(), ids.end() );
  auto const [it, inserted] = m_faceCycles.emplace( ids, cycle );
  if( !inserted && it->second != cycle )
    throw std::invalid_argument( "Incompatible edge cycles on a shared refinement face" );
  return insert( { m_namespace, EntityKind::face, std::move( cycle ) }, corners );
}

vtkIdType PointRegistry::cell( vtkIdType globalCellId, Connectivity const & corners )
{
  if( globalCellId < 0 )
    throw std::invalid_argument( "Negative refinement cell global ID" );
  return insert( { m_namespace, EntityKind::cell, { globalCellId } }, corners );
}

vtkIdType PointRegistry::pointForKey( EntityKey const & key ) const
{
  auto const found = m_indices.find( key );
  if( found == m_indices.end() )
    throw std::invalid_argument( "Refinement interface point was not planned" );
  return found->second;
}

vtkIdType PointRegistry::findPoint( EntityKey const & key ) const
{
  auto const found = m_indices.find( key );
  return found == m_indices.end() ? -1 : found->second;
}

vtkIdType PointRegistry::findEdge( vtkIdType a, vtkIdType b ) const
{
  Connectivity ids{ m_globalIds.at( a ), m_globalIds.at( b ) };
  std::sort( ids.begin(), ids.end() );
  return findPoint( { m_namespace, EntityKind::edge, std::move( ids ) } );
}

vtkIdType PointRegistry::findFace( Connectivity const & corners ) const
{
  Connectivity ids;
  for( vtkIdType corner : corners )
    ids.push_back( m_globalIds.at( corner ) );
  return findPoint( { m_namespace, EntityKind::face, canonicalCycle( ids ) } );
}

SharingInheritance::SharingInheritance( PointRegistry const & points, std::vector< EntitySupport > entities, int localRank )
  : m_entities( std::move( entities ) ), m_incident( points.originalSize() )
{
  if( localRank < 0 )
    throw std::invalid_argument( "Negative refinement participant rank" );
  std::unordered_map< EntityKey, std::size_t, EntityKeyHash > unique;
  for( std::size_t i = 0; i < m_entities.size(); ++i )
  {
    EntitySupport & entity = m_entities[i];
    if( !( entity.key == entityKey( entity.key.kind, entity.key.corners, entity.key.meshNamespace ) ) )
      throw std::invalid_argument( "Noncanonical refinement support key" );
    std::sort( entity.localCorners.begin(), entity.localCorners.end() );
    if( entity.localCorners.empty() ||
        std::adjacent_find( entity.localCorners.begin(), entity.localCorners.end() ) != entity.localCorners.end() )
      throw std::invalid_argument( "Empty or repeated refinement support corner" );
    if( entity.participants.empty() || !std::is_sorted( entity.participants.begin(), entity.participants.end() ) ||
        entity.participants.front() < 0 ||
        std::adjacent_find( entity.participants.begin(), entity.participants.end() ) != entity.participants.end() ||
        !std::binary_search( entity.participants.begin(), entity.participants.end(), localRank ) )
      throw std::invalid_argument( "Invalid inherited refinement participants" );
    if( entity.key.kind == EntityKind::cell && entity.participants.size() != 1 )
      throw std::invalid_argument( "A refinement volume cell must have exactly one owner" );
    Connectivity globalCorners;
    for( vtkIdType corner : entity.localCorners )
    {
      if( corner < 0 || corner >= points.originalSize() )
        throw std::invalid_argument( "Refinement support corner outside previous level" );
      globalCorners.push_back( points.points()[corner].key.corners.front() );
    }
    if( entity.key.kind != EntityKind::cell )
    {
      std::sort( globalCorners.begin(), globalCorners.end() );
      Connectivity expected = entity.key.corners;
      std::sort( expected.begin(), expected.end() );
      if( globalCorners != expected )
        throw std::invalid_argument( "Refinement support key/corner mismatch" );
    }
    auto const [previous, inserted] = unique.emplace( entity.key, i );
    if( !inserted )
    {
      EntitySupport const & prior = m_entities[previous->second];
      if( prior.localCorners != entity.localCorners || prior.participants != entity.participants )
        throw std::invalid_argument( "Conflicting inherited refinement support" );
      continue;
    }
    for( vtkIdType corner : entity.localCorners )
      m_incident[corner].push_back( i );
  }
  m_pointSupports.reserve( points.points().size() );
  for( PointRecipe const & recipe : points.points() )
  {
    Connectivity corners = recipe.support;
    std::sort( corners.begin(), corners.end() );
    auto const entity = unique.find( recipe.key );
    if( entity == unique.end() || m_entities[entity->second].localCorners != corners )
      throw std::invalid_argument( "Missing exact previous-level point support" );
    m_pointSupports.push_back( std::move( corners ) );
  }
}

EntitySupport const & SharingInheritance::support( Connectivity const & fineCorners ) const
{
  Connectivity oldCorners;
  for( vtkIdType corner : fineCorners )
  {
    if( corner < 0 || static_cast< std::size_t >( corner ) >= m_pointSupports.size() )
      throw std::invalid_argument( "Fine refinement corner outside point registry" );
    Connectivity const & defining = m_pointSupports[corner];
    oldCorners.insert( oldCorners.end(), defining.begin(), defining.end() );
  }
  if( oldCorners.empty() )
    throw std::invalid_argument( "Empty fine refinement entity" );
  std::sort( oldCorners.begin(), oldCorners.end() );
  oldCorners.erase( std::unique( oldCorners.begin(), oldCorners.end() ), oldCorners.end() );
  // Start from the least incident corner: lookup cost follows local valence,
  // rather than scanning every old cell for each fine entity.
  auto const least = std::min_element( oldCorners.begin(), oldCorners.end(),
                                       [&]( vtkIdType a, vtkIdType b ) { return m_incident[a].size() < m_incident[b].size(); } );
  EntitySupport const * best = nullptr;
  bool ambiguous = false;
  for( std::size_t index : m_incident[*least] )
  {
    EntitySupport const & candidate = m_entities[index];
    if( !std::includes( candidate.localCorners.begin(), candidate.localCorners.end(), oldCorners.begin(), oldCorners.end() ) )
      continue;
    if( best == nullptr || candidate.key.kind < best->key.kind )
    {
      best = &candidate;
      ambiguous = false;
    }
    else if( candidate.key.kind == best->key.kind && !( candidate.key == best->key ) )
      ambiguous = true;
  }
  if( best == nullptr )
    throw std::invalid_argument( "Fine entity spans unrelated previous-level supports" );
  if( ambiguous )
    throw std::invalid_argument( "Ambiguous previous-level refinement support" );
  return *best;
}

} // namespace geos::vtk::refinement
