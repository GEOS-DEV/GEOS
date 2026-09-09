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
 * @file StatisticsAggregatorBase.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_STATISTICSAGGREGATOR_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_STATISTICSAGGREGATOR_HPP_

#include "common/DataTypes.hpp"
#include "dataRepository/DataContext.hpp"
#include "dataRepository/Group.hpp"
#include "mesh/CellElementRegion.hpp"
#include "mesh/MeshBody.hpp"
#include "mesh/MeshLevel.hpp"

namespace geos
{

/**
 * @brief Output data group to contain the result of a given stat aggregator on the dataRepository.
 *        Attributes are public since the class is a POD.
 * @todo repair 1D HDF5 outputs to enable stats HDF5 outputs
 */
class RegionStatisticsBase : public dataRepository::Group
{
public:

  /// Time of statistics computation (std::numeric_limits< double >::lowest() if not initialized)
  real64 m_time;

  /**
   * @brief Construct a new Region Statistics Base object
   * @param targetName name of the data-repository object that is targeted by the statistics
   *                   (mesh level / region / sub-region).
   * @param parent the instance parent in data-repository
   * @param statsOutputEnabled If true, the stats are saved in the output HDF5
   *                           (through dataRepository::RestartFlags, but not functional for this output for now).
   */
  RegionStatisticsBase( string const & targetName,
                        dataRepository::Group * const parent,
                        bool statsOutputEnabled );

  /**
   * @return the name of the data-repository object that is targeted by the statistics
   *         (mesh level / region / sub-region).
   */
  string_view getTargetName() const
  { return getName(); }

};

template< typename T >
struct StatsAggregatorTraits;

/**
 * @brief Reponsible of computing physical statistics over the grid, registering the result in the
 *        data repository, but not storing / outputing it by itself. It does not have mutable state
 *        except the encountered issues.
 * @todo repair 1D HDF5 outputs to enable stats HDF5 outputs
 * @tparam Impl the derived type of the statistics aggregator which contains all necessary implementations (CRTP)
 */
template< typename Impl >
class StatsAggregatorBase
{
public:

  using SolverType = typename StatsAggregatorTraits< Impl >::SolverType;

  using StatsGroupType = typename StatsAggregatorTraits< Impl >::StatsGroupType;

  /**
   * @brief Standard function signature for any functor that applies on statistics group instances (StatsGroupType)
   *        - param 0: OwnerType &, the group instance containing the data for which we want to aggregate the statistics (MeshLevel,
   * CellElementRegion...)
   *        - param 1: StatsAggregateGroupType &, the statistics aggregate Group where to store the data
   * @tparam OwnerType the concrete type of the OwnerType param
   */
  template< typename OwnerType >
  using RegionStatsFunc = std::function< void ( OwnerType &,
                                                StatsGroupType & ) >;

  /**
   * @brief A functor that can be used to register a statistics Group instance. Parameters:
   *        - dataRepository::Group &, the parent Group,
   *        - string const &, name of the statistics target (i.e. region name).
   */
  using RegionStatsRegisterFunc = std::function< StatsGroupType & ( dataRepository::Group &,
                                                                    string const & ) >;

  /**
   * @brief the associated view keys
   */
  struct ViewKeys
  {
    /// String for the discretization statistics group
    constexpr static char const * statisticsString() { return "statistics"; }
    /// String for the region statistics group
    constexpr static char const * regionsStatisticsString() { return "regionsStatistics"; }
  };

  /**
   * @brief Construct a new Stats Aggregator object
   * @param ownerName the unique name of the entity requesting the statistics.
   *                  An error is thrown if not unique in this context.
   * @param meshBodies The Group containing the MeshBody objects
   * @param statsOutputEnabled If true, the stats are saved in the output HDF5
   *                           (through dataRepository::RestartFlags, but not functional for this output for now).
   */
  StatsAggregatorBase( dataRepository::DataContext const & ownerDataContext,
                       dataRepository::Group & meshBodies,
                       bool statsOutputEnabled );

  /**
   * @brief Enable the computation of any statistics, initialize data structure to collect them.
   *        Register the resulting data wrappers so they will be targeted by TimeHistory output
   * @param solver flow solver object to retrieve:
                   - the simulated regions,
                   - fields for statistics computation.
   */
  void initStatisticsAggregation( SolverType & solver );

  void forRegionStatistics( RegionStatsFunc< MeshLevel > const & functor ) const;

  void forRegionStatistics( MeshLevel & mesh,
                            StatsGroupType & meshRegionsStatistics,
                            RegionStatsFunc< CellElementRegion > const & functor ) const;

  void forRegionStatistics( CellElementRegion & region,
                            StatsGroupType & regionStatistics,
                            RegionStatsFunc< CellElementSubRegion > const & functor ) const;

  /**
   * @param[in] timeRequest The time for which we want to know if the statistics are computed.
   * @param[in] stats the statistics data structure we want to know if it has been computed
   * @return true if the statistics have been computed.
   */
  bool isComputed( real64 const timeRequest, StatsGroupType const & stats );

  /**
   * @brief set the statistics as dirty, ensuring isComputed() will be false until the next computation.
   */
  void setDirty();

  /**
   * @brief Compute statistics on the mesh discretizations (average field pressure, etc)
   *        Results are reduced on rank 0, and broadcasted over all ranks.
   * @param[in] timeRequest The time for which we want to compute the statistics.
   * @return false if there was a problem that prevented the statistics to be computed correctly.
   */
  bool computeRegionsStatistics( real64 const timeRequest );

  /**
   * @return the name of the entity that needs the statistics.
   */
  string const & getOwnerName() const
  { return m_ownerDataContext.getTargetName(); }

  /**
   * @return The encountered issues during the last computing method call.
   */
  stdVector< string > const & getWarnings() const
  { return m_warnings; }

  dataRepository::Group & getInstanceStatisticsGroup( MeshLevel & mesh ) const;

  StatsGroupType & getRegionsStatistics( MeshLevel & mesh ) const;

  /**
   * @return a specific statistics Group instance.
   * @param mesh The desired mesh-level
   * @param regionName The name of the desired region
   * @throw InputError if no statistics data is found for the given region name.
   */
  StatsGroupType & getRegionStatistics( MeshLevel & mesh, string_view regionName ) const;

protected:

  struct StatsState
  {
    bool m_isEnabled = false;
    bool m_isDirty = false;
  };

  struct DiscretizationGroupPath
  {
    localIndex m_meshBody;
    localIndex m_meshLevel;
    string_array m_regionNames;
  };

  /// @see getOwnerName()
  dataRepository::DataContext const & m_ownerDataContext;

  /// If true, the stats are to save in the output HDF5
  bool const m_statsOutputEnabled;

  dataRepository::Group & m_meshBodies;

  stdVector< DiscretizationGroupPath > m_discretizationsPaths;

  /// @see getWarnings()
  stdVector< string > m_warnings;

  /// The current state of the region statistics
  StatsState m_regionStatsState;

  /**
   * @brief Enable the computation of region statistics, initialize data structure to collect them.
   *        Register the resulting data wrappers so they will be targeted by TimeHistory output
   * @note Must be called in or after the "registerDataOnMesh" initialization phase
   * @param registerStatsFunc The functor which register each statistics group whithin the regions hierarchy
   */
  void enableRegionStatisticsAggregation( RegionStatsRegisterFunc && registerStatsFunc );

  /**
   * @param path the path of the discretization group in the data-repository.
   * @return MeshLevel& the MeshLevel Group for the given discretisation
   */
  MeshLevel & getMeshLevel( DiscretizationGroupPath const & path ) const;

  /**
   * @brief Initialize all statistics values to aggregable default values,
   *        before any computation / reduction for the current timestep.
   * @param stats the statistics instance
   * @param time start time of the current timestep (s)
   * @note Must be implemented for each type that implements this template (CRTP).
   */
  void initStats( StatsGroupType & stats, real64 time ) const
  { static_cast< Impl const * >(this)->initStats( stats, time ); }

  /**
   * @brief Compute the rank-local stats for the given sub-region and store the results in the given stats group.
   * @param subRegion
   * @param subRegionStats the stats group instance for the subregion
   * @note Must be implemented for each type that implements this template (CRTP).
   */
  void computeSubRegionRankStats( CellElementSubRegion & subRegion, StatsGroupType & subRegionStats ) const
  { static_cast< Impl const * >(this)->computeSubRegionRankStats( subRegion, subRegionStats ); }

  /**
   * @brief Aggregate all instance statistics with those of another instance on the current rank.
   * @param stats the statistics instance
   * @param other the other instance to aggregate with.
   * @note Must be implemented for each type that implements this template (CRTP).
   */
  void aggregateStats( StatsGroupType & stats, StatsGroupType const & other ) const
  { static_cast< Impl const * >(this)->aggregateStats( stats, other ); }

  /**
   * @brief Aggregate all instance statistics with those of other ranks.
   * @param stats the statistics instance
   * @note Must be implemented for each type that implements this template (CRTP).
   */
  void mpiAggregateStats( StatsGroupType & stats ) const
  { static_cast< Impl const * >(this)->mpiAggregateStats( stats ); }

  /**
   * @brief Do the final computations for the statistics. Must be called after computations & aggregations.
   * @param stats the statistics instance
   * @note Must be implemented for each type that implements this template (CRTP).
   */
  void postAggregateStats( StatsGroupType & stats )
  { static_cast< Impl * >(this)->postAggregateStats( stats ); }

};

} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_STATISTICSAGGREGATOR_HPP_ */
