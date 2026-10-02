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

/** @file VTKRefinementSharing.cpp */
#include "VTKRefinementSharing.hpp"

#include <algorithm>
#include <stdexcept>
#include <unordered_set>

namespace geos::vtk::refinement
{
namespace
{
Connectivity globalCorners( Connectivity const & corners, Connectivity const & globalIds )
{
  Connectivity result;
  result.reserve( corners.size() );
  for( vtkIdType corner : corners )
  {
    if( corner < 0 || static_cast< std::size_t >( corner ) >= globalIds.size() )
      throw std::invalid_argument( "Refinement interface corner outside local points" );
    result.push_back( globalIds[corner] );
  }
  return result;
}

bool sameCorners( Connectivity a, Connectivity b )
{
  std::sort( a.begin(), a.end() );
  std::sort( b.begin(), b.end() );
  return a == b;
}
} // namespace

CoarseBoundary coarseBoundary( std::vector< Cell > const & cells, Connectivity const & globalCellIds, Connectivity const & globalPointIds,
                               int localRank, std::uint64_t meshNamespace )
{
  if( localRank < 0 || cells.size() != globalCellIds.size() )
    throw std::invalid_argument( "Invalid coarse refinement rank or cell ID count" );
  std::unordered_set< vtkIdType > points;
  points.reserve( globalPointIds.size() );
  for( vtkIdType id : globalPointIds )
    if( id < 0 || !points.insert( id ).second )
      throw std::invalid_argument( "Negative or repeated local coarse point ID" );
  struct Face
  {
    EntitySupport entity;
    int count = 1;
  };
  std::vector< Face > faces;
  std::unordered_map< Connectivity, std::size_t, ConnectivityHash > faceIndices;
  std::unordered_set< vtkIdType > volumes;
  CoarseBoundary result;
  result.volumeIds.reserve( cells.size() );
  for( std::size_t i = 0; i < cells.size(); ++i )
  {
    vtkIdType const id = globalCellIds[i];
    if( id < 0 || !volumes.insert( id ).second )
      throw std::invalid_argument( "Negative or duplicate local coarse volume ID" );
    result.volumeIds.push_back( entityKey( EntityKind::cell, { id }, meshNamespace ) );
    for( Connectivity const & corners : cellFaces( cells[i] ) )
    {
      Connectivity cycle = globalCorners( corners, globalPointIds );
      EntityKey const key = entityKey( EntityKind::face, cycle, meshNamespace );
      std::sort( cycle.begin(), cycle.end() );
      auto const [previous, inserted] = faceIndices.emplace( cycle, faces.size() );
      if( inserted )
        faces.push_back( { { key, corners, { localRank } }, 1 } );
      else
      {
        Face & face = faces[previous->second];
        if( !( face.entity.key == key ) )
          throw std::invalid_argument( "Incompatible local volume face edge cycles" );
        if( ++face.count > 2 )
          throw std::invalid_argument( "Nonmanifold local volume face incidence" );
      }
    }
  }
  std::unordered_map< EntityKey, std::size_t, EntityKeyHash > indices;
  auto append = [&]( EntityKind kind, Connectivity const & corners )
  {
    EntityKey key = entityKey( kind, globalCorners( corners, globalPointIds ), meshNamespace );
    if( indices.emplace( key, result.entities.size() ).second )
      result.entities.push_back( { std::move( key ), corners, { localRank } } );
  };
  for( Face const & face : faces )
    if( face.count == 1 )
    {
      append( EntityKind::face, face.entity.localCorners );
      auto const & corners = face.entity.localCorners;
      for( std::size_t i = 0; i < corners.size(); ++i )
      {
        append( EntityKind::vertex, { corners[i] } );
        append( EntityKind::edge, { corners[i], corners[( i + 1 ) % corners.size()] } );
      }
    }
  return result;
}

InterfaceSharing::InterfaceSharing( PointRegistry & points, std::vector< EntitySupport > sharedEntities, int localRank )
  : m_points( points ), m_local{ localRank }, m_entities( std::move( sharedEntities ) )
{
  if( localRank < 0 )
    throw std::invalid_argument( "Negative refinement interface rank" );
  m_indices.reserve( m_entities.size() );
  for( std::size_t i = 0; i < m_entities.size(); ++i )
  {
    auto const & entity = m_entities[i];
    auto const & ranks = entity.participants;
    if( entity.key.kind == EntityKind::cell || ranks.size() < 2 || ranks.front() < 0 || !std::is_sorted( ranks.begin(), ranks.end() ) ||
        std::adjacent_find( ranks.begin(), ranks.end() ) != ranks.end() || !std::binary_search( ranks.begin(), ranks.end(), localRank ) )
      throw std::invalid_argument( "Invalid shared refinement interface participants" );
    Connectivity ids;
    for( vtkIdType corner : entity.localCorners )
    {
      if( corner < 0 || corner >= points.originalSize() )
        throw std::invalid_argument( "Shared refinement corner outside previous level" );
      auto const & vertex = points.points()[corner].key;
      if( vertex.meshNamespace != entity.key.meshNamespace )
        throw std::invalid_argument( "Wrong shared refinement mesh namespace" );
      ids.push_back( vertex.corners.front() );
    }
    if( !( entityKey( entity.key.kind, std::move( ids ), entity.key.meshNamespace ) == entity.key ) )
      throw std::invalid_argument( "Shared refinement key/corner cycle mismatch" );
    auto const [previous, inserted] = m_indices.emplace( entity.key, i );
    if( !inserted )
    {
      auto const & existing = m_entities[previous->second];
      if( existing.participants != ranks || !sameCorners( existing.localCorners, entity.localCorners ) )
        throw std::invalid_argument( "Conflicting shared refinement interface" );
      continue;
    }
    for( vtkIdType corner : entity.localCorners )
      m_incident[corner].push_back( i );
  }
  // A missing shared edge cannot be approximated by its face participants.
  for( auto const & entity : m_entities )
  {
    auto required = [&]( EntityKind kind, Connectivity ids )
    {
      auto const found = m_indices.find( entityKey( kind, std::move( ids ), entity.key.meshNamespace ) );
      if( found == m_indices.end() )
        throw std::invalid_argument( "Shared refinement incidence is not closed" );
      auto const & ranks = m_entities[found->second].participants;
      if( !std::includes( ranks.begin(), ranks.end(), entity.participants.begin(), entity.participants.end() ) )
        throw std::invalid_argument( "Shared refinement subentity omits parent participants" );
    };
    for( vtkIdType id : entity.key.corners )
      required( EntityKind::vertex, { id } );
    if( entity.key.kind == EntityKind::face )
      for( std::size_t i = 0; i < entity.key.corners.size(); ++i )
        required( EntityKind::edge, { entity.key.corners[i], entity.key.corners[( i + 1 ) % entity.key.corners.size()] } );
  }
}

Participants const & InterfaceSharing::participants( Connectivity const & fineCorners ) const
{
  if( fineCorners.empty() )
    throw std::invalid_argument( "Empty fine refinement interface" );
  Connectivity support;
  for( vtkIdType corner : fineCorners )
  {
    if( corner < 0 || static_cast< std::size_t >( corner ) >= m_points.points().size() )
      throw std::invalid_argument( "Fine refinement corner outside point registry" );
    auto const & point = m_points.points()[corner];
    if( fineCorners.size() == 1 )
    {
      auto const found = m_indices.find( point.key );
      return found == m_indices.end() ? m_local : m_entities[found->second].participants;
    }
    support.insert( support.end(), point.support.begin(), point.support.end() );
  }
  std::sort( support.begin(), support.end() );
  support.erase( std::unique( support.begin(), support.end() ), support.end() );
  std::vector< std::size_t > const * candidates = nullptr;
  for( vtkIdType corner : support )
  {
    auto const found = m_incident.find( corner );
    if( found == m_incident.end() )
      return m_local;
    if( candidates == nullptr || found->second.size() < candidates->size() )
      candidates = &found->second;
  }
  EntitySupport const * best = nullptr;
  bool ambiguous = false;
  for( std::size_t index : *candidates )
  {
    auto const & entity = m_entities[index];
    if( !std::all_of(
          support.begin(), support.end(), [&]( vtkIdType corner )
    { return std::find( entity.localCorners.begin(), entity.localCorners.end(), corner ) != entity.localCorners.end(); } ) )
      continue;
    if( best == nullptr || entity.key.kind < best->key.kind )
    {
      best = &entity;
      ambiguous = false;
    }
    else if( entity.key.kind == best->key.kind && !( entity.key == best->key ) )
      ambiguous = true;
  }
  if( ambiguous )
    throw std::invalid_argument( "Ambiguous shared refinement support" );
  return best == nullptr ? m_local : best->participants;
}

std::vector< EntitySupport > InterfaceSharing::fineSupports( Connectivity const & fineGlobalIds )
{
  std::size_t const plannedPoints = m_points.points().size();
  if( fineGlobalIds.size() != plannedPoints )
    throw std::invalid_argument( "Fine interface ID count mismatch" );
  for( auto const & entity : m_entities )
    if( entity.key.kind == EntityKind::edge || ( entity.key.kind == EntityKind::face && entity.localCorners.size() >= 4 ) )
      m_points.pointForKey( entity.key );
  std::vector< EntitySupport > result;
  std::unordered_map< EntityKey, std::size_t, EntityKeyHash > indices;
  auto append = [&]( EntityKind kind, Connectivity const & corners, std::uint64_t meshNamespace )
  {
    auto const & ranks = participants( corners );
    if( ranks.size() < 2 )
      return;
    EntityKey key = entityKey( kind, globalCorners( corners, fineGlobalIds ), meshNamespace );
    auto const [previous, inserted] = indices.emplace( key, result.size() );
    if( inserted )
      result.push_back( { std::move( key ), corners, ranks } );
    else if( result[previous->second].participants != ranks )
      throw std::invalid_argument( "Conflicting inherited fine interface participants" );
  };
  for( auto const & entity : m_entities )
  {
    auto const ns = entity.key.meshNamespace;
    if( entity.key.kind == EntityKind::vertex )
      append( EntityKind::vertex, entity.localCorners, ns );
    else if( entity.key.kind == EntityKind::edge )
    {
      vtkIdType const a = entity.localCorners[0], b = entity.localCorners[1];
      vtkIdType const middle = m_points.pointForKey( entity.key );
      if( m_points.points().size() != plannedPoints )
        throw std::invalid_argument( "Shared edge midpoint was not planned by a volume template" );
      append( EntityKind::vertex, { middle }, ns );
      append( EntityKind::edge, { a, middle }, ns );
      append( EntityKind::edge, { middle, b }, ns );
    }
    else
    {
      auto const children = subdivideFace( entity.localCorners, m_points );
      if( m_points.points().size() != plannedPoints )
        throw std::invalid_argument( "Shared face points were not planned by a volume template" );
      for( auto const & face : children )
      {
        append( EntityKind::face, face, ns );
        for( std::size_t i = 0; i < face.size(); ++i )
        {
          append( EntityKind::vertex, { face[i] }, ns );
          append( EntityKind::edge, { face[i], face[( i + 1 ) % face.size()] }, ns );
        }
      }
    }
  }
  return result;
}

} // namespace geos::vtk::refinement
