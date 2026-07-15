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
                                                   string_view setName,
                                                   bool const dataOutputEnabled ):
  dataRepository::Group( targetName, parent ),
  m_time( std::numeric_limits< double >::lowest() ),
  m_setName( setName )
{
  Group::setRestartFlags( dataOutputEnabled ?
                          dataRepository::RestartFlags::WRITE_AND_READ :
                          dataRepository::RestartFlags::NO_WRITE );

  // TODO : registerWrappers to store results in HDF5 (but need repairing of 1D HDF5 outputs)
}

template< typename Impl >
StatsAggregatorBase< Impl >::StatsAggregatorBase( dataRepository::DataContext const & ownerDataContext,
                                                  dataRepository::Group & meshBodies,
                                                  bool const dataOutputEnabled ):
  m_ownerDataContext( ownerDataContext ),
  m_dataOutputEnabled( dataOutputEnabled ),
  m_meshBodies( meshBodies ),
  m_isAnySetsIntersecting( false )
{}

template< typename Impl >
std::set< string >
StatsAggregatorBase< Impl >::getMeshLevelPartialSetNames( MeshLevel const & mesh,
                                                          string_array const & regionNames )
{
  std::set< string > foundSetNames;

  mesh.getElemManager().forElementRegions< CellElementRegion >( regionNames,
                                                                [&] ( localIndex &,
                                                                      CellElementRegion const & region )
  {
    region.forElementSubRegions< CellElementSubRegion >( [&] ( CellElementSubRegion const & subRegion ) {
      dataRepository::Group const & sets = subRegion.sets();
      sets.forWrappers< SetType >( [&] ( dataRepository::Wrapper< SetType > const & set )
      {
        foundSetNames.emplace( set.getName() );
      } );
    } );
  } );

  return foundSetNames;
}

template< typename Impl >
void
StatsAggregatorBase< Impl >::initStatisticsAggregation( SolverType & solver )
{
  solver.forDiscretizationOnMeshTargets( m_meshBodies, [&] ( string const & meshBodyName,
                                                             MeshLevel & mesh,
                                                             string_array const & regionNames )
  {
    // getting the container of all requesters statistics groups (can be already initialized)
    dataRepository::Group * meshStatsGroup = mesh.getGroupPointer( string( ViewKeys::statisticsString()) );
    if( meshStatsGroup == nullptr )
      meshStatsGroup = &mesh.registerGroup( string( ViewKeys::statisticsString() ) );

    // registering the container of instance statistics groups (must be unique for this instance)
    string const & ownerName = getOwnerName();
    GEOS_ERROR_IF_NE_MSG( meshStatsGroup->hasGroup( ownerName ), false,
                          GEOS_FMT( "A statistics aggregator have already been requested for '{}'.",
                                    ownerName ),
                          m_ownerDataContext );
    meshStatsGroup->registerGroup( ownerName );

    MeshBody const & body = m_meshBodies.getGroup< MeshBody >( meshBodyName );

    /// finding all sets in this mesh level
    std::set< string > foundSetNames = getMeshLevelPartialSetNames( mesh, regionNames );

    { // remembering the path of this discretization
      DiscretizationSetPath const path {
        /* .m_meshBody = */ body.getIndexInParent(),
        /* .m_meshLevel = */ mesh.getIndexInParent(),
        /* .m_setNames = */ string_array( foundSetNames.begin(), foundSetNames.end()),
        /* .m_regionNames = */ regionNames,
        /* .m_setsCompound = */ {},
      };
      m_discretizationsPaths.emplace_back( path );
    }
  } );
}

template< typename Impl >
void
StatsAggregatorBase< Impl >::enableRegionStatisticsAggregation( RegionStatsRegisterFunc && registerStatsFunc,
                                                                string_array const & setNames )
{
  m_setNames.insert( setNames.begin(), setNames.end() );
  if( m_setNames.empty())
    m_setNames.emplace( "all" );

  bool regionFound = false;
  bool subRegionFound = false;
  std::set< string > confirmedSets;

  for( auto & path : m_discretizationsPaths )
  {
    MeshLevel & mesh = getMeshLevel( path );
    ElementRegionManager & elemManager = mesh.getElemManager();
    dataRepository::Group & statisticsGroup = getInstanceStatisticsGroup( mesh );
    StatsGroupType & meshStats = registerStatsFunc( statisticsGroup, string( ViewKeys::setsStatisticsString() ),
                                                    "", m_dataOutputEnabled );
    SetsCompound meshSetsCompounds;
    bool isAnySetIntersects = false;

    // registering all stats Group in the dataRepository for each set
    for( string const & setName : m_setNames )
    {
      StatsGroupType & meshSetStats = registerStatsFunc( meshStats, setName,
                                                         setName, m_dataOutputEnabled );

      for( size_t regionId = 0; regionId < path.m_regionNames.size(); ++regionId )
      {
        ElementRegionBase & region = elemManager.getRegion< ElementRegionBase >( path.m_regionNames[regionId] );
        if( dynamic_cast< CellElementRegion * >(&region) != nullptr )
        { // we only process cell regions
          StatsGroupType & setRegionStats = registerStatsFunc( meshSetStats, region.getName(),
                                                               setName, m_dataOutputEnabled );
          stdMap< string, SetType > & regionSetsCompounds = meshSetsCompounds.get_inserted( region.getName() );
          regionFound = true;

          region.forElementSubRegions< CellElementSubRegion >( [&] ( CellElementSubRegion & subRegion )
          {
            registerStatsFunc( setRegionStats, subRegion.getName(),
                               setName, false );
            subRegionFound = true;

            auto const * const subRegionSetWrapper = subRegion.sets().getWrapperPointer< SetType >( setName );
            if( subRegionSetWrapper != nullptr )
            {
              confirmedSets.emplace( setName );

              // Insert the set elements ids in the mesh-level compound. If doubles are found, sets intersect.
              SetViewType const & subRegionSet = subRegionSetWrapper->reference();
              arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();
              SetType & subRegionSetsCompound = regionSetsCompounds.get_inserted( subRegion.getName() );
              for( localIndex setIter = 0; setIter < subRegionSet.size(); ++setIter )
              {
                localIndex elemId = subRegionSet[setIter];
                if( elemGhostRank[elemId] < 0 )
                {
                  if( !subRegionSetsCompound.insert( elemId ) )
                    isAnySetIntersects = true;
                }
              }
            }
          } );
        }
      }
    }

    // if needed, registering stats Group in the dataRepository for the compound set
    isAnySetIntersects = (bool)( MpiWrapper::max( (unsigned char)( isAnySetIntersects ) ) );
    if( isAnySetIntersects )
    {
      string_view setName = ViewKeys::compoundSetNameString();
      StatsGroupType & meshSetStats = registerStatsFunc( meshStats, string( setName ),
                                                         setName, false );
      for( size_t regionId = 0; regionId < path.m_regionNames.size(); ++regionId )
      {
        ElementRegionBase & region = elemManager.getRegion< ElementRegionBase >( path.m_regionNames[regionId] );
        if( dynamic_cast< CellElementRegion * >(&region) != nullptr )
        { // we only process cell regions
          StatsGroupType & setRegionStats = registerStatsFunc( meshSetStats, region.getName(),
                                                               setName, false );

          region.forElementSubRegions< CellElementSubRegion >( [&] ( CellElementSubRegion & subRegion )
          {
            registerStatsFunc( setRegionStats, subRegion.getName(),
                               setName, false );

            path.m_setsCompound
              .get_inserted( region.getName())
              .insert( meshSetsCompounds
                         .at( region.getName() )
                         .extract( subRegion.getName() ) );
          } );
        }
      }
    }
  }

  m_regionStatsState.m_isEnabled = true;
  m_regionStatsState.m_isDirty = true;

  GEOS_ERROR_IF( regionFound == 0 || subRegionFound == 0,
                 GEOS_FMT( "Missing region for computing statistics:\n- {} regions,\n- {} sub-regions.",
                           getOwnerName(),
                           regionFound ? "found" : "no",
                           subRegionFound ? "found" : "no" ),
                 m_ownerDataContext );

  // at this point, m_setNames content is element sets user *requests*: they are the same on every ranks,
  // but we don't know yet if any of these set is missing or empty (=> confirmedSets).
  std::set< string > notFoundSets;
  for( string const & requestedSetName : m_setNames )
  {
    if( MpiWrapper::max( confirmedSets.count( requestedSetName ) ) == 0 )
      notFoundSets.emplace( requestedSetName );
  }
  if( MpiWrapper::commRank() == 0 )
    GEOS_ERROR_IF( !notFoundSets.empty(),
                   GEOS_FMT( "During retrieval of statistics element sets, the following set(s) could not be "
                             "found or does not contain any element:\n- {}",
                             stringutilities::join( notFoundSets, "\n- " ) ),
                   m_ownerDataContext );
}

template< typename Impl >
MeshLevel &
StatsAggregatorBase< Impl >::getMeshLevel( DiscretizationSetPath const & path ) const
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
  dataRepository::Group & meshStatsGroup = mesh.getGroup( string( ViewKeys::statisticsString() ) );
  dataRepository::Group & instanceStatisticsGroup = meshStatsGroup.getGroup( getOwnerName() );
  return instanceStatisticsGroup;
}

template< typename Impl >
typename StatsAggregatorBase< Impl >::StatsGroupType &
StatsAggregatorBase< Impl >::
getMeshLevelStatistics( MeshLevel & mesh ) const
{
  // considering everything is initialized, or else, crash gracefully
  dataRepository::Group & instanceStatisticsGroup = getInstanceStatisticsGroup( mesh );
  return instanceStatisticsGroup.getGroup< StatsGroupType >( string( ViewKeys::setsStatisticsString()) );
}

template< typename Impl >
typename StatsAggregatorBase< Impl >::StatsGroupType &
StatsAggregatorBase< Impl >::getSetStatistics( MeshLevel & mesh,
                                               string_view requestedSetName ) const
{
  // considering mesh-level stats structure is initialized, or else, crash gracefully
  StatsGroupType & meshStats = getMeshLevelStatistics( mesh );
  string const setName = string( requestedSetName.empty() ? "all" : requestedSetName );
  StatsGroupType * const stats = meshStats.template getGroupPointer< StatsGroupType >( setName );
  GEOS_THROW_IF( stats == nullptr,
                 GEOS_FMT( "Element set '{}' statistics data structure not found, has the set been requested?\nRequested element sets:\n- {}",
                           setName, stringutilities::join( m_setNames, "\n- " ) ),
                 InputError, m_ownerDataContext );
  return *stats;
}

template< typename Impl >
typename StatsAggregatorBase< Impl >::StatsGroupType &
StatsAggregatorBase< Impl >::getCompoundSetStatistics( MeshLevel & mesh ) const
{
  // considering mesh-level stats structure is initialized, or else, crash gracefully
  StatsGroupType & meshStats = getMeshLevelStatistics( mesh );
  string const setName = string( ViewKeys::compoundSetNameString() );
  StatsGroupType * const stats = meshStats.template getGroupPointer< StatsGroupType >( setName );
  GEOS_THROW_IF( stats == nullptr,
                 GEOS_FMT( "Compound element set data structure not found ({}), is the aggregator initialized? Was the compound needed?\nRequested element sets:\n- {}",
                           setName, stringutilities::join( m_setNames, "\n- " ) ),
                 InputError, m_ownerDataContext );
  return *stats;
}

template< typename Impl >
bool
StatsAggregatorBase< Impl >::isCompoundSetStatistics( StatsGroupType const & setStats ) const
{
  return setStats.getSetName() == ViewKeys::compoundSetNameString();
}

template< typename Impl >
typename StatsAggregatorBase< Impl >::StatsGroupType &
StatsAggregatorBase< Impl >::getRegionStatistics( MeshLevel & mesh,
                                                  string_view regionName,
                                                  string_view setName ) const
{
  StatsGroupType & setStats = getSetStatistics( mesh, setName );
  StatsGroupType * const stats = setStats.template getGroupPointer< StatsGroupType >( string( regionName ) );
  GEOS_THROW_IF( stats == nullptr,
                 GEOS_FMT( "Region '{}' not found to get region statistics, is it a target of the reservoir solver?\n"
                           "Available target regions:\n- {}",
                           regionName, stringutilities::join( setStats.getSubGroupsNames(), "\n- " ) ),
                 InputError, m_ownerDataContext );
  return *stats;
}

template< typename Impl >
void
StatsAggregatorBase< Impl >::forRegionStatistics( RegionStatsFunc< MeshLevel & > const & func ) const
{
  for( auto const & path : m_discretizationsPaths )
  {
    MeshLevel & mesh = getMeshLevel( path );
    StatsGroupType & meshLevelStats = getMeshLevelStatistics( mesh );

    func( mesh, meshLevelStats );
  }
}

template< typename Impl >
void
StatsAggregatorBase< Impl >::forRegionStatistics( MeshLevel & mesh,
                                                  StatsGroupType & meshLevelStats,
                                                  bool enableCompoundSet,
                                                  RegionStatsFunc< MeshLevelSet > const & func ) const
{
  meshLevelStats.template forSubGroups< StatsGroupType >( [&] ( StatsGroupType & setStats )
  {
    if( enableCompoundSet || setStats.getSetName() != ViewKeys::compoundSetNameString() )
    {
      MeshLevelSet meshSet {
        /* .mesh = */ mesh,
        /* .setName = */ setStats.getTargetName(),
      };

      func( meshSet, setStats );
    }
  } );
}

template< typename Impl >
void
StatsAggregatorBase< Impl >::forRegionStatistics( MeshLevelSet const meshSet,
                                                  StatsGroupType & setStats,
                                                  RegionStatsFunc< CellElementRegion & > const & func ) const
{
  ElementRegionManager & elemManager = meshSet.mesh.getElemManager();
  setStats.template forSubGroups< StatsGroupType >( [&] ( StatsGroupType & setRegionStats )
  {
    string_view regionName = setRegionStats.getTargetName();
    ElementRegionBase & region = elemManager.getRegion< ElementRegionBase >( string( regionName ) );
    CellElementRegion * cellRegion = dynamic_cast< CellElementRegion * >(&region); // we only process cell regions

    if( cellRegion != nullptr )
      func( *cellRegion, setRegionStats );
  } );
}

template< typename Impl >
void
StatsAggregatorBase< Impl >::forRegionStatistics( CellElementRegion & setRegion,
                                                  StatsGroupType & setRegionStats,
                                                  RegionStatsFunc< CellElementSubRegion & > const & func ) const
{
  setRegionStats.template forSubGroups< StatsGroupType >( [&] ( StatsGroupType & subRegionStatistics )
  {
    string_view subRegionName = subRegionStatistics.getTargetName();
    CellElementSubRegion & subRegion = setRegion.getSubRegion< CellElementSubRegion >( string( subRegionName ) );

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

  for( auto const & path : m_discretizationsPaths )
  {
    MeshLevel & mesh = getMeshLevel( path );
    StatsGroupType & meshLevelStats = getMeshLevelStatistics( mesh );
    bool const isSetsCompoundEnabled = !path.m_setsCompound.empty();

    // computation of sub region stats for each selected set
    initStats( meshLevelStats, timeRequest );
    forRegionStatistics( mesh, meshLevelStats, true,
                         [&, timeRequest] ( MeshLevelSet meshSet, StatsGroupType & setStats )
    {
      initStats( setStats, timeRequest );
      forRegionStatistics( meshSet, setStats,
                           [&, timeRequest] ( CellElementRegion & region, StatsGroupType & regionStats )
      {
        initStats( regionStats, timeRequest );
        forRegionStatistics( region, regionStats,
                             [&, timeRequest] ( CellElementSubRegion & subRegion, StatsGroupType & subRegionStats )
        {
          initStats( subRegionStats, timeRequest );

          SetType const & targetSet = meshSet.setName == ViewKeys::compoundSetNameString() ?
                                      path.m_setsCompound.at( region.getName() ).at( subRegion.getName() ) :
                                      subRegion.sets().getReference< SetType >( string( meshSet.setName ) );
          computeSubRegionRankStats( subRegion, subRegionStats, targetSet );
        } );
      } );
    } );

    // aggregation of computations from the sub regions
    forRegionStatistics( mesh, meshLevelStats, true,
                         [&] ( MeshLevelSet meshSet, StatsGroupType & setStats )
    {
      forRegionStatistics( meshSet, setStats,
                           [&] ( CellElementRegion & region, StatsGroupType & regionStats )
      {
        forRegionStatistics( region, regionStats,
                             [&] ( CellElementSubRegion &, StatsGroupType & subRegionStats )
        {
          aggregateStats( regionStats, subRegionStats );
          // sub-region stats finalization is disabled, it does not seem useful
          // mpiAggregateStats( subRegionStats );
          // postAggregateStats( subRegionStats );
        } );

        aggregateStats( setStats, regionStats );

        mpiAggregateStats( regionStats );
        postAggregateStats( regionStats );
      } );

      // if there is no compound set, we can simply aggregate each set statistics to the mesh-level statistics.
      // if there is one, we will only take into account the statistics of the compound set.
      if( !isSetsCompoundEnabled || isCompoundSetStatistics( setStats ) )
      {
        aggregateStats( meshLevelStats, setStats );
      }

      mpiAggregateStats( setStats );
      postAggregateStats( setStats );
    } );

    mpiAggregateStats( meshLevelStats );
    postAggregateStats( meshLevelStats );
  }

  m_regionStatsState.m_isDirty = false;

  return true;
}

} /* namespace geos */
