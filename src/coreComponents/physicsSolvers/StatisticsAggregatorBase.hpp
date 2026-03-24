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
 * @file CompositionalMultiphaseStatisticsAg/gregator.hpp
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
   * @param setName Optional name of a set, restricting the statistics on given mesh element.
   *                If none is specified (empty), "all" elements set is targeted.
   * @param dataOutputEnabled If true, the stats are saved in the output HDF5 (dataRepository::RestartFlags, not functional for this output
   * for now).
   */
  RegionStatisticsBase( string const & targetName,
                        dataRepository::Group * const parent,
                        string_view setName,
                        bool dataOutputEnabled );

  /**
   * @return the name of the data-repository object that is targeted by the statistics
   *         (mesh level / region / sub-region).
   */
  string_view getTargetName() const
  { return getName(); }

  /**
   * @return Optional name of a set, restricting the statistics on given mesh element.
   *         If none is specified (empty), "all" elements set is targeted.
   */
  string_view getSetName() const
  { return m_setName; }

private:
  /// @see getSetName()
  string const m_setName;
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

  using SetType = SortedArray< localIndex >;

  using SetViewType = SortedArrayView< localIndex const >;

  /**
   * @brief Standard function signature for any functor that applies on statistics group instances (StatsGroupType)
   *        - param 0: OwnerType &, the group instance containing the data for which we want to aggregate the statistics (MeshLevel,
   * CellElementRegion...)
   *        - param 1: StatsAggregateGroupType &, the statistics aggregate Group where to store the data
   * @tparam OwnerType the concrete type of the OwnerType param
   */
  template< typename OwnerType >
  using RegionStatsFunc = std::function< void ( OwnerType,
                                                StatsGroupType & ) >;

  /**
   * @brief TODO Document
   *        - dataRepository::Group & targetName: the parent Group,
   *        - string const & parent: name of the statistics target (i.e. region name).
   *        - string_view setName: name of the target element set. If none is specified (empty), "all" elements set is targeted.
   *        - bool dataOutputEnabled: enable Group output features
   */
  using RegionStatsRegisterFunc = std::function< StatsGroupType & ( dataRepository::Group &,
                                                                    string const &,
                                                                    string_view,
                                                                    bool ) >;

  /**
   * @brief the associated view keys
   */
  struct ViewKeys
  {
    /// String for the mesh-level statistics group
    constexpr static char const * statisticsString() { return "statistics"; }
    /// String for the region statistics group
    constexpr static char const * setsStatisticsString() { return "setsStatistics"; }
    /// String for the region statistics group
    constexpr static char const * regionsStatisticsString() { return "regionsStatistics"; }
    /// on-purpose generated compound set name
    constexpr static char const * compoundSetNameString() { return "__compound"; }
  };

  /**
   * @brief Allow to reference a set in a given MeshLevel structure.
   */
  struct MeshLevelSet
  {
    MeshLevel & mesh;
    string_view setName;
  };

  /**
   * @brief Construct a new Stats Aggregator object
   * @param ownerName the unique name of the entity requesting the statistics.
   *                  An error is thrown if not unique in this context.
   */
  StatsAggregatorBase( dataRepository::DataContext const & ownerDataContext,
                       bool statsOutputEnabled );

  /**
   * @brief Enable the computation of any statistics, explores the data structure to collect them.
   * @param solver flow solver object to retrieve:
                   - the simulated regions,
                   - fields for statistics computation.
   * @param meshBodies The Group containing the MeshBody objects
   */
  void initStatisticsAggregation( dataRepository::Group & meshBodies,
                                  SolverType & solver );

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
   * @name Data Structure Traversal Strategy
   */
  ///@{

  /**
   * @brief Execute the given functor on each mesh-level (discretization) available for this instance.
   * @param functor the function to execute.
   */
  void forRegionStatistics( RegionStatsFunc< MeshLevel & > const & functor ) const;

  /**
   * @brief Execute the given functor on each set present in the given mesh-level (discretization).
   * @param mesh the target mesh-level.
   * @param meshSetsStats the instance mesh-level statistics data structure.
   * @param enableCompoundSet if true, the functor will also be applied on the compound set (which
   *                          is an optional internal computation data structure).
   * @param functor the function to execute.
   */
  void forRegionStatistics( MeshLevel & mesh,
                            StatsGroupType & meshSetsStats,
                            bool enableCompoundSet,
                            RegionStatsFunc< MeshLevelSet > const & functor ) const;

  /**
   * @brief Execute the given functor on each region available in the given mesh-level for the given set.
   * @param meshSet the target mesh-level set.
   * @param setStats the mesh-level statistics data structure for the given set.
   * @param functor the function to execute.
   */
  void forRegionStatistics( MeshLevelSet meshSet,
                            StatsGroupType & setStats,
                            RegionStatsFunc< CellElementRegion & > const & functor ) const;

  /**
   * @brief Execute the given functor on each sub-region available in the given region for the given set.
   * @param meshSet the target region set.
   * @param setRegionStats the region statistics data structure for the given set.
   * @param functor the function to execute.
   */
  void forRegionStatistics( CellElementRegion & region,
                            StatsGroupType & setRegionStats,
                            RegionStatsFunc< CellElementSubRegion & > const & functor ) const;

  /**
   * @brief Get the Instance Statistics Group for the given mesh-level.
   * @param mesh the mesh-level (discretization).
   * @return A Group instance hierarchically holding all statistics data structures.
   */
  dataRepository::Group & getInstanceStatisticsGroup( MeshLevel & mesh ) const;

  /**
   * @return A statistics data structure (see StatsGroupType), aggregated from all targeted regions
   *         and sets in the given mesh-level.
   * @param mesh the mesh-level (discretization).
   */
  StatsGroupType & getMeshLevelStatistics( MeshLevel & mesh ) const;

  /**
   * @return A statistics data structure (see StatsGroupType), aggregated from all targeted regions
   *         for the given set in the given mesh-level.
   * @param mesh the mesh-level (discretization).
   * @param setName the name of element set. If none is specified (empty), "all" elements set is targeted.
   * @throw InputError if "all" is targeted, but the instance does not target "all" element set.
   */
  StatsGroupType & getSetStatistics( MeshLevel & mesh, string_view setName = "" ) const;

  /**
   * @return A statistics data structure (see StatsGroupType), for the given set, in the given region,
   *         in the given mesh-level.
   * @param mesh The mesh-level (discretization)
   * @param setName the name of element set. If none is specified (empty), "all" elements set is targeted.
   * @param regionName The name of the desired region
   * @throw InputError - if no statistics data is found for the given region name.
   *                   - if "all" is targeted, but the instance does not target "all" element set.
   */
  StatsGroupType & getRegionStatistics( MeshLevel & mesh, string_view regionName, string_view setName = "" ) const;

  ///@}

  /**
   * @param[in] timeRequest The time for which we want to know if the statistics are computed.
   * @param[in] stats the statistics data structure we want to know if it has been computed
   * @return true if the statistics have been computed.
   */
  bool isComputed( real64 const timeRequest, StatsGroupType const & stats );

  /**
   * @return true if the region statistics target given sets.
   *         false if the statistics are targeted on only one mesh-level-wide set
   */
  bool isRestrictedToSets()
  {return m_setNames.empty() || ( m_setNames.size() == 1 && m_setNames.count( "all" ) > 0 ); }

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

protected:

  struct StatsState
  {
    bool m_isEnabled = false;
    bool m_isDirty = false;
  };

  struct DiscretizationSetPath
  {
    /// The id of the mesh body in the "meshBodies" Group
    localIndex m_meshBody;
    /// The id of the mesh level in the mesh body.
    localIndex m_meshLevel;
    /// The names of this mesh-level element sets. If none is specified (empty), "all" elements set is targeted.
    string_array m_setNames;
    /// the regions names in the mesh-level / mesh level
    string_array m_regionNames;
    /// if at least one set intersects with another, a compound set is needed to compute mesh-level statistics.
    SetType m_setsCompound;
  };

  /// @see getOwnerName()
  dataRepository::DataContext const & m_ownerDataContext;

  bool const m_dataOutputEnabled;

  dataRepository::Group * m_meshBodies = nullptr;

  /// Data repository paths of the element sets per mesh levels
  /// TODO: can it be generalized with the "all" set
  std::vector< DiscretizationSetPath > m_discretizationsPaths;

  /// The list of mesh element sets to restrict the region statistics to.
  /// Cannot be empty: if the whole mesh-level is to process, "all" must be targeted.
  std::set< string > m_setNames;

  /// if true, at least ont set intersects with another, threfore we must compute mesh-level statistics with a coupound.
  bool m_isAnySetsIntersecting;

  /// @see getWarnings()
  stdVector< string > m_warnings;

  /// The current state of the region statistics
  StatsState m_regionStatsState;

  /**
   * @brief Enable the computation of region statistics, initialize data structure to collect them.
   *        Register the resulting data wrappers so they will be targeted by TimeHistory output
   * @note Must be called in or after the "registerDataOnMesh" initialization phase
   * @param registerStatsFunc The functor which register each statistics group whithin the regions hierarchy
   * @param setNames The list of mesh element sets to restrict the statistics to.
   *                 If empty, the whole mesh-level is processed.
   */
  void enableRegionStatisticsAggregation( RegionStatsRegisterFunc && registerStatsFunc,
                                          string_array const & setNames );

  /**
   * @param path the path of the mesh-level group in the data-repository.
   * @return MeshLevel& the MeshLevel Group for the given discretisation
   */
  MeshLevel & getMeshLevel( DiscretizationSetPath const & path ) const;

  /**
   * @return An optional statistics data structure (see StatsGroupType), aggregated from all sets
   *         over all targeted regions, for the given set in the given mesh-level.
   * @param mesh the mesh-level (discretization).
   * @throw InputError if no compound set exists, for instance if isSetsCompoundEnabled == false.
   */
  StatsGroupType & getCompoundSetStatistics( MeshLevel & mesh ) const;

  /**
   * @brief TODO
   */
  static std::set< string > getMeshLevelPartialSetNames( MeshLevel const & mesh,
                                                         string_array const & regionNames );

  /**
   * @brief Set the statstics target sets to restrict the statistics to.
   * @param setNames the list of target sets. If "all" is included, or none is present, all the discretisation is targeted.
   */
  void setTargetSets( string_array const & setNames );

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
   * @param targetSet
   * @note Must be implemented for each type that implements this template (CRTP).
   */
  void computeSubRegionRankStats( CellElementSubRegion & subRegion,
                                  StatsGroupType & subRegionStats,
                                  SetType const & targetSet ) const
  { static_cast< Impl const * >(this)->computeSubRegionRankStats( subRegion, subRegionStats, targetSet ); }

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
