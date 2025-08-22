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
 * @file ScatterDataHistoryCollection.cpp
 */

#include "ScatterDataHistoryCollection.hpp"

#include "mesh/DomainPartition.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "common/MpiWrapper.hpp"
#include <cstring>

namespace geos
{

using namespace dataRepository;

ScatterDataHistoryCollection::ScatterDataHistoryCollection( string const & name, Group * const parent )
  : HistoryCollectionBase( name, parent ),
  m_scatterDataProvider( nullptr )
{
  registerWrapper( viewKeyStruct::solverNameString(), &m_solverName )
    .setInputFlag( InputFlags::REQUIRED )
    .setDescription( "Name of the physics solver that provides scatter data" );

  registerWrapper( viewKeyStruct::includeCoordinatesString(), &m_includeCoordinates )
    .setInputFlag( InputFlags::OPTIONAL )
    .setApplyDefaultValue( 1 )
    .setDescription( "Flag to include coordinates in output (1=yes, 0=no)" );

  registerWrapper( viewKeyStruct::includeMetadataString(), &m_includeMetadata )
    .setInputFlag( InputFlags::OPTIONAL )
    .setApplyDefaultValue( 0 )
    .setDescription( "Flag to include metadata in output (1=yes, 0=no)" );
}

void ScatterDataHistoryCollection::initializePostInitialConditionsPreSubGroups()
{
  HistoryCollectionBase::initializePostInitialConditionsPreSubGroups();

  // Find the domain partition to locate the physics solver
  findScatterDataProvider();

  // Set up the collection count
  m_collectionCount = 1; // Main scatter data

  if( m_includeCoordinates && m_scatterDataProvider != nullptr )
  {
    array2d< real64 > const & coordinates = m_scatterDataProvider->getScatterCoordinates();
    if( coordinates.size( 0 ) > 0 && coordinates.size( 1 ) > 0 )
    {
      m_collectionCount += coordinates.size( 1 ); // Add X, Y, Z coordinates
    }
  }
}

void ScatterDataHistoryCollection::findScatterDataProvider()
{
  // Get the problem manager from the hierarchy
  Group const & problemManager = this->getGroupByPath( "/Problem" );

  // Get the physics solver manager
  Group const & physicsSolverManager = problemManager.getGroup( "Solvers" );

  // Try to find the solver by name
  Group const * solver = physicsSolverManager.getGroupPointer( m_solverName );
  GEOS_THROW_IF( solver == nullptr,
                 GEOS_FMT( "Could not find physics solver named '{}'", m_solverName ),
                 InputError );

  // Cast to ScatterDataProvider
  m_scatterDataProvider = dynamic_cast< ScatterDataProvider const * >( solver );
  GEOS_THROW_IF( m_scatterDataProvider == nullptr,
                 GEOS_FMT( "Physics solver '{}' does not implement ScatterDataProvider interface", m_solverName ),
                 InputError );

  GEOS_LOG_RANK_0( GEOS_FMT( "Found scatter data provider '{}' with {} points", m_solverName, m_scatterDataProvider->getNumScatterPoints()) );
}

void ScatterDataHistoryCollection::collect( DomainPartition const & domain )
{
  GEOS_UNUSED_VAR( domain );

  if( m_scatterDataProvider == nullptr )
  {
    GEOS_THROW( "Scatter data provider not initialized", std::runtime_error );
  }

  // Get the current scatter data values
  array1d< real64 > const & scatterData = m_scatterDataProvider->getScatterData();
  localIndex const numPoints = m_scatterDataProvider->getNumScatterPoints();

  GEOS_THROW_IF( scatterData.size() != numPoints,
                 GEOS_FMT( "Scatter data size ({}) does not match number of points ({})", scatterData.size(), numPoints ),
                 std::runtime_error );
}

void ScatterDataHistoryCollection::collect( DomainPartition const & domain,
                                            localIndex const collectionIdx,
                                            buffer_unit_type * & buffer )
{
  GEOS_UNUSED_VAR( domain );
  GEOS_MARK_FUNCTION;

  if( m_scatterDataProvider == nullptr )
  {
    return;
  }

  localIndex const numPoints = m_scatterDataProvider->getNumScatterPoints();

  // Collection index 0 is the main scatter data
  if( collectionIdx == 0 )
  {
    array1d< real64 > const & scatterData = m_scatterDataProvider->getScatterData();

    // Copy data to buffer
    std::memcpy( buffer, scatterData.data(), numPoints * sizeof( real64 ) );
    buffer += numPoints * sizeof( real64 );
  }
  // Subsequent indices are coordinates if requested
  else if( m_includeCoordinates )
  {
    array2d< real64 > const & coordinates = m_scatterDataProvider->getScatterCoordinates();
    localIndex const numDim = coordinates.size( 1 );
    localIndex const coordIdx = collectionIdx - 1;

    if( coordIdx < numDim && coordinates.size( 0 ) == numPoints )
    {
      // Extract the specific coordinate component
      real64 * coordData = reinterpret_cast< real64 * >( buffer );
      for( localIndex i = 0; i < numPoints; ++i )
      {
        coordData[i] = coordinates[i][coordIdx];
      }
      buffer += numPoints * sizeof( real64 );
    }
  }
}

localIndex ScatterDataHistoryCollection::numCollectors() const
{
  localIndex count = 1; // Main scatter data

  if( m_includeCoordinates && m_scatterDataProvider != nullptr )
  {
    array2d< real64 > const & coordinates = m_scatterDataProvider->getScatterCoordinates();
    if( coordinates.size( 0 ) > 0 && coordinates.size( 1 ) > 0 )
    {
      count += coordinates.size( 1 ); // Add X, Y, Z coordinates
    }
  }

  return count;
}

const string & ScatterDataHistoryCollection::getTargetName() const
{
  return m_solverName;
}

HistoryMetadata ScatterDataHistoryCollection::getMetaData( DomainPartition const & domain, localIndex collectionIdx ) const
{
  GEOS_UNUSED_VAR( domain );

  if( m_scatterDataProvider == nullptr )
  {
    return HistoryMetadata();
  }

  localIndex const numPoints = m_scatterDataProvider->getNumScatterPoints();

  // Collection index 0 is always the main scatter data
  if( collectionIdx == 0 )
  {
    return HistoryMetadata( "scatterData", numPoints, std::type_index( typeid( real64 ) ) );
  }

  // Subsequent indices are coordinates if requested
  if( m_includeCoordinates )
  {
    array2d< real64 > const & coordinates = m_scatterDataProvider->getScatterCoordinates();
    localIndex const numDim = coordinates.size( 1 );

    if( collectionIdx <= numDim )
    {
      string coordName;
      switch( collectionIdx - 1 )
      {
        case 0: coordName = "coordinateX"; break;
        case 1: coordName = "coordinateY"; break;
        case 2: coordName = "coordinateZ"; break;
        default: coordName = "coordinate" + std::to_string( collectionIdx - 1 ); break;
      }
      return HistoryMetadata( coordName, numPoints, std::type_index( typeid( real64 ) ) );
    }
  }

  return HistoryMetadata();
}

void ScatterDataHistoryCollection::updateSetsIndices( DomainPartition const & domain )
{
  GEOS_UNUSED_VAR( domain );
}

localIndex ScatterDataHistoryCollection::numMetaDataCollectors() const
{
  return 0;
}

HistoryCollection & ScatterDataHistoryCollection::getMetaDataCollector( localIndex metaIdx )
{
  GEOS_UNUSED_VAR( metaIdx );
  GEOS_THROW( "Scatter data collections do not have metadata collectors",
              std::runtime_error );
}

REGISTER_CATALOG_ENTRY( TaskBase, ScatterDataHistoryCollection, string const &, Group * const )

} // namespace geos
