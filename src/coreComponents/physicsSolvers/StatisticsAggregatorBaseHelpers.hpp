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

#include "mesh/MeshBody.hpp"
#include "physicsSolvers/StatisticsAggregatorBase.hpp"

#include "LvArray/src/math.hpp"
#include "common/DataTypes.hpp"
#include "common/format/Format.hpp"
#include "common/format/StringUtilities.hpp"
#include "common/logger/Logger.hpp"
#include "dataRepository/DataContext.hpp"
#include "dataRepository/Group.hpp"
#include "mesh/CellElementRegion.hpp"
#include "mesh/CellElementSubRegion.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/MeshLevel.hpp"

namespace geos
{

inline RegionStatisticsBase::RegionStatisticsBase( string const & targetName,
                                                   dataRepository::Group * const parent,
                                                   bool const statsOutputEnabled ):
  dataRepository::Group( targetName, parent ),
  m_time( std::numeric_limits< double >::lowest() )
{
  Group::setRestartFlags( statsOutputEnabled ?
                          dataRepository::RestartFlags::WRITE_AND_READ :
                          dataRepository::RestartFlags::NO_WRITE );

  // TODO : registerWrappers to store results in HDF5 (but need repairing of 1D HDF5 outputs)
}

template< typename Impl >
StatsAggregatorBase< Impl >::StatsAggregatorBase( dataRepository::DataContext const & ownerDataContext,
                                                  dataRepository::Group & meshBodies,
                                                  bool const statsOutputEnabled ):
  m_ownerDataContext( ownerDataContext ),
  m_statsOutputEnabled( statsOutputEnabled ),
  m_meshBodies( meshBodies )
{}

template< typename Impl >
void
StatsAggregatorBase< Impl >::initStatisticsAggregation( SolverType & solver )
{
  solver.forDiscretizationOnMeshTargets( m_meshBodies, [&] ( string const & meshBodyName,
                                                             MeshLevel & mesh,
                                                             string_array const & regionNames )
  {
    // getting the container of all requesters statistics groups (can be already initialized)
    dataRepository::Group * meshStatsGroup = mesh.getGroupPointer( ViewKeys::statisticsString() );
    if( meshStatsGroup == nullptr )
      meshStatsGroup = &mesh.registerGroup( ViewKeys::statisticsString() );

    // registering the container of instance statistics groups (must be unique for this instance)
    string const & ownerName = getOwnerName();
    GEOS_ERROR_IF_NE_MSG( meshStatsGroup->hasGroup( ownerName ), false,
                          GEOS_FMT( "A statistics aggregator have already been requested for '{}'.",
                                    ownerName ),
                          m_ownerDataContext );
    meshStatsGroup->registerGroup( ownerName );

    // remembering the path of this discretization
    MeshBody const & body = m_meshBodies.getGroup< MeshBody >( meshBodyName );
    DiscretizationGroupPath const path {
      /* .m_meshBody = */ body.getIndexInParent(),
      /* .m_meshLevel = */ mesh.getIndexInParent(),
      /* .m_regionNames = */ regionNames,
    };
    m_discretizationsPaths.push_back( path );
  } );
}

template< typename Impl >
void
StatsAggregatorBase< Impl >::enableRegionStatisticsAggregation( RegionStatsRegisterFunc && registerStatsFunc )
{
  integer regionCount = 0;
  integer subRegionCount = 0;

  for( auto const & path : m_discretizationsPaths )
  {
    MeshLevel & mesh = getMeshLevel( path );
    ElementRegionManager & elemManager = mesh.getElemManager();
    dataRepository::Group & statisticsGroup = getInstanceStatisticsGroup( mesh );
    StatsGroupType & meshRegionsStats = registerStatsFunc( statisticsGroup,
                                                           ViewKeys::regionsStatisticsString() );

    for( size_t i = 0; i < path.m_regionNames.size(); ++i )
    {
      CellElementRegion & region = elemManager.getRegion< CellElementRegion >( path.m_regionNames[i] );
      StatsGroupType & regionStats = registerStatsFunc( meshRegionsStats,
                                                        region.getName() );

      region.forElementSubRegions< CellElementSubRegion >( [&] ( CellElementSubRegion & subRegion )
      {
        registerStatsFunc( regionStats,
                           subRegion.getName() );
        ++subRegionCount;
      } );
      ++regionCount;
    }
  }

  GEOS_ERROR_IF( regionCount == 0 || subRegionCount == 0,
                 GEOS_FMT( "Missing region for computing statistics:\n- {} regions,\n- {} sub-regions.",
                           getOwnerName(), regionCount, subRegionCount ),
                 m_ownerDataContext );

  m_regionStatsState.m_isEnabled = true;
  m_regionStatsState.m_isDirty = true;
}

template< typename Impl >
void
StatsAggregatorBase< Impl >::forRegionStatistics( RegionStatsFunc< MeshLevel > const & func ) const
{
  for( auto const & path : m_discretizationsPaths )
  {
    MeshLevel & mesh = getMeshLevel( path );
    StatsGroupType & meshRegionsStats = getRegionsStatistics( mesh );

    func( mesh, meshRegionsStats );
  }
}

template< typename Impl >
void
StatsAggregatorBase< Impl >::forRegionStatistics( MeshLevel & mesh,
                                                  StatsGroupType & meshRegionsStatistics,
                                                  RegionStatsFunc< CellElementRegion > const & func ) const
{
  ElementRegionManager & elemManager = mesh.getElemManager();
  meshRegionsStatistics.template forSubGroups< StatsGroupType >( [&] ( StatsGroupType & regionStatistics )
  {
    string_view targetName = regionStatistics.getTargetName();
    CellElementRegion & region = elemManager.getRegion< CellElementRegion >( string( targetName ) );

    func( region, regionStatistics );
  } );
}

template< typename Impl >
void
StatsAggregatorBase< Impl >::forRegionStatistics( CellElementRegion & region,
                                                  StatsGroupType & regionStatistics,
                                                  RegionStatsFunc< CellElementSubRegion > const & func ) const
{
  regionStatistics.template forSubGroups< StatsGroupType >( [&] ( StatsGroupType & subRegionStatistics )
  {
    string_view targetName = subRegionStatistics.getTargetName();
    CellElementSubRegion & subRegion = region.getSubRegion< CellElementSubRegion >( string( targetName ) );
    func( subRegion, subRegionStatistics );
  } );
}

template< typename Impl >
bool
StatsAggregatorBase< Impl >::isComputed( real64 const timeRequest, StatsGroupType const & stats )
{
  real64 const timePrecisionScale = LvArray::math::max( LvArray::math::abs( timeRequest ),
                                                        LvArray::math::abs( stats.m_time ) );
  static constexpr real64 timeRelTol = 1.0e-12;

  return
    !m_regionStatsState.m_isDirty &&
    LvArray::math::abs( timeRequest - stats.m_time ) < timeRelTol * timePrecisionScale;
}

template< typename Impl >
void
StatsAggregatorBase< Impl >::setDirty()
{
  m_regionStatsState.m_isDirty = true;
}

template< typename Impl >
bool
StatsAggregatorBase< Impl >::computeRegionsStatistics( real64 const timeRequest )
{
  GEOS_MARK_FUNCTION;

  m_warnings.clear();

  // computation of sub region stats
  forRegionStatistics( [&, timeRequest] ( MeshLevel & mesh, StatsGroupType & meshRegionsStats )
  {
    forRegionStatistics( mesh,
                         meshRegionsStats,
                         [&, timeRequest] ( CellElementRegion & region, StatsGroupType & regionStats )
    {
      forRegionStatistics( region,
                           regionStats,
                           [&, timeRequest] ( CellElementSubRegion & subRegion, StatsGroupType & subRegionStats )
      {
        initStats( subRegionStats, timeRequest );
        computeSubRegionRankStats( subRegion, subRegionStats );
      } );
    } );
  } );

  // aggregation of computations from the sub regions
  forRegionStatistics( [&, timeRequest] ( MeshLevel & mesh, StatsGroupType & meshRegionsStats )
  {
    initStats( meshRegionsStats, timeRequest );

    forRegionStatistics( mesh,
                         meshRegionsStats,
                         [&, timeRequest] ( CellElementRegion & region, StatsGroupType & regionStats )
    {
      initStats( regionStats, timeRequest );

      forRegionStatistics( region,
                           regionStats,
                           [&] ( CellElementSubRegion &, StatsGroupType & subRegionStats )
      {
        aggregateStats( regionStats, subRegionStats );

        mpiAggregateStats( subRegionStats );
        postAggregateStats( subRegionStats );
      } );

      aggregateStats( meshRegionsStats, regionStats );

      mpiAggregateStats( regionStats );
      postAggregateStats( regionStats );
    } );

    mpiAggregateStats( meshRegionsStats );
    postAggregateStats( meshRegionsStats );
  } );

  m_regionStatsState.m_isDirty = false;

  return true;
}

template< typename Impl >
MeshLevel &
StatsAggregatorBase< Impl >::getMeshLevel( DiscretizationGroupPath const & path ) const
{
  MeshBody & body = m_meshBodies.getGroup< MeshBody >( path.m_meshBody );
  MeshLevel & mesh = body.getMeshLevel( path.m_meshLevel );
  return mesh;
}

template< typename Impl >
dataRepository::Group &
StatsAggregatorBase< Impl >::getInstanceStatisticsGroup( MeshLevel & mesh ) const
{
  // considering everything is initialized, or else, crash gracefully
  dataRepository::Group & meshStatsGroup = mesh.getGroup( ViewKeys::statisticsString() );
  dataRepository::Group & instanceStatisticsGroup = meshStatsGroup.getGroup( getOwnerName() );
  return instanceStatisticsGroup;
}

template< typename Impl >
typename StatsAggregatorBase< Impl >::StatsGroupType &
StatsAggregatorBase< Impl >::
getRegionsStatistics( MeshLevel & mesh ) const
{
  // considering everything is initialized, or else, crash gracefully
  dataRepository::Group & instanceStatisticsGroup = getInstanceStatisticsGroup( mesh );
  return instanceStatisticsGroup.getGroup< StatsGroupType >( ViewKeys::regionsStatisticsString() );
}

template< typename Impl >
typename StatsAggregatorBase< Impl >::StatsGroupType &
StatsAggregatorBase< Impl >::getRegionStatistics( MeshLevel & mesh,
                                                  string_view regionName ) const
{
  StatsGroupType & meshRegionsStats = getRegionsStatistics( mesh );
  StatsGroupType * const stats = meshRegionsStats.template getGroupPointer< StatsGroupType >( string( regionName ) );
  GEOS_THROW_IF( stats == nullptr,
                 GEOS_FMT( "Region '{}' not found to get region statistics, is it a target of the reservoir solver?\n"
                           "Available target regions:\n- {}",
                           regionName, stringutilities::join( meshRegionsStats.getSubGroupsNames(), "\n- " ) ),
                 InputError, m_ownerDataContext );
  return *stats;
}

} /* namespace geos */
