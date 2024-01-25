/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SourceFluxStatistics.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_SOURCEFLUXSTATISTICS_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_SOURCEFLUXSTATISTICS_HPP_

#include "physicsSolvers/FieldStatisticsBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"

class SourceFluxBoundaryCondition;

namespace geos
{

/**
 * @class SourceFluxStatsAggregator
 *
 * Task class allowing for the computation of aggregate statistics of SourceFluxBoundaryCondition
 */
class SourceFluxStatsAggregator : public FieldStatisticsBase< FlowSolverBase >
{
public:

  /**
   * @brief Potentially aggregated statistics flux data
   */
  struct StatData
  {
    /// fluid mass produced by the flux(es) (kg). Negative if injecting.
    real64 m_producedMass = 0.0;
    /// flux(es) production rate (kg/s). Negative if injecting.
    real64 m_productionRate = 0.0;
    /// Number of elements in which we are producing / injecting
    integer m_elementCount = 0;

    /**
     * @brief Aggregate the statistics of the instance with those of another one.
     * @param other the other stats structure.
     */
    void combine( StatData const & other );
    /**
     * @brief Aggregate the statistics of the instance with those from all instances from other MPI ranks.
     */
    void mpiReduce();
  };
  /**
   * @brief Class that aggregate statistics of a flux over multiple time-steps for a given
   * SourceFluxStatsAggregator and a for a given mesh part (i.e. a subregion, a region...).
   */
  class WrappedStats
  {
public:
    /**
     * @brief Set the subjects targeted by the stats.
     * @param aggregatorName the name of the targeted SourceFluxStatsAggregator
     * @param fluxName the name of the targeted SourceFluxBoundaryCondition
     */
    void setTarget( string_view aggregatorName, string_view fluxName );

    /**
     * @brief Set the current time step stats.
     * @param overwriteTimeStepStats false when this is the first time we're writing the current
     *                               timestep data, true otherwise (i.e. new solver iteration)
     * @param dt                     time delta of the current timestep
     * @param productedMass          time-step producted mass (see StatData::m_producedMass).
     * @param elementCount           number of cell elements concerned by this instance
     */
    void setTimeStepStats( real64 dt, real64 productedMass, integer elementCount,
                           bool overwriteTimeStepStats );

    /**
     * @brief Finalize the statistics period and render data.
     * @return the accumulated statistics of the period.
     */
    StatData finalizePeriod();

    /**
     * @return get the name of the SourceFluxStatsAggregator that want to collect data on this instance.
     */
    string_view getAggregatorName() const
    { return m_aggregatorName; }

    /**
     * @return get the name of the SourceFluxBoundaryCondition from which we are collecting data on this instance.
     */
    string_view getFluxName() const
    { return m_fluxName; }
private:
    /// producted mass of the current time-step.
    real64 m_currentTimeStepMass = 0.0;
    /// producted mass sum from all previous time-step of the current period.
    real64 m_pendingPeriodMass = 0.0;
    /// delta time the current period
    real64 m_periodDeltaTime = 0.0;
    /// number of cell elements concerned by this instance
    integer m_elementCount = 0;

    /// Name of the SourceFluxStatsAggregator that want to collect data on this instance.
    string m_aggregatorName;
    /// Name of the SourceFluxBoundaryCondition from which we are collecting data on this instance.
    string m_fluxName;
  };


  /**
   * @brief Constructor for the statistics class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  SourceFluxStatsAggregator( const string & name,
                                  Group * const parent );

  /// Accessor for the catalog name
  static string catalogName() { return "SourceFluxStatistics"; }

  /**
   * @defgroup Tasks Interface Functions
   *
   * This function implements the interface defined by the abstract TaskBase class
   */
  /**@{*/

  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /**@}*/

  /**
   * @return a WrappedStats struct that contains the statistics of the flux for the
   * SourceFluxStatsAggregator instance in the container.
   * @note To be retrieved, the WrappedStats struct must be registered on the container during the
   * registerDataOnMesh() call.
   * @param container the container from which we want the statistics.
   * @param fluxName  the name of the flux from which we want the statistics.
   * @throw           a GEOS_ERROR if the flux was not found.
   */
  WrappedStats & getFluxStatData( Group & container, string_view fluxName );

  /**
   * @brief Apply a functor to all WrappedStats of the given group that target a given flux (and
   * potentially multiple SourceFluxStatsAggregator).
   * The functor takes in parameter the reference to the currently processed WrappedStats.
   * @note To be retrieved, the WrappedStats structs must be registered on the container during the
   * registerDataOnMesh() call.
   * @tparam LAMBDA   the type of lambda function to call in the function
   * @param container the container from which we want the statistics.
   * @param fluxName  the name of the flux from which we want the statistics.
   * @param lambda    the functor that will be called for each WrappedStats. Takes in parameter the
   *                  reference to the currently processed WrappedStats.
   */
  template< typename LAMBDA >
  static void forAllFluxStatData( Group & container, string_view fluxName, LAMBDA && lambda );

  /**
   * @return a string used to name the wrapper that is added to each region that is simulated by the solver.
   * The string is unique within the region for the SourceFluxBoundaryCondition and the SourceFluxStatsAggregator.
   */
  inline string getRegionStatDataName( string_view fluxName ) const;


  /**
   * @brief View keys
   */
  struct viewKeyStruct
  {
    /// @return The key for setName
    constexpr static char const * fluxNamesString() { return "fluxNames"; }
  };

protected:

  /**
   * @copydoc Group::postProcessInput()
   */
  void postProcessInput() override;

private:
  using Base = FieldStatisticsBase< FlowSolverBase >;

  /// the names of the SourceFlux(s) for which we want the statistics
  string_array m_fluxNames;

  /**
   * @copydoc Group::registerDataOnMesh(Group &)
   */
  void registerDataOnMesh( Group & meshBodies ) override;

  /**
   * @brief Output in the log the given statistics.
   * @param regionName the name of the element group (a region, a sub-region...) from which we want
   * to output the data.
   * @param minLogLevel the min log level to output any line.
   * @param subSetName the region / sub-subregion name concerned by the statistics.
   * @param fluxName the flux name concerned by the statistics.
   * @param stats the statistics that must be output in the log.
   */
  void writeStatData( integer minLogLevel, string_view subSetName, string_view fluxName, StatData const & stats );

};


template< typename LAMBDA >
void SourceFluxStatsAggregator::forAllFluxStatData( Group & container,
                                                         string_view fluxName,
                                                         LAMBDA && lambda )
{
  container.forWrappers< WrappedStats >( [&]( dataRepository::Wrapper< WrappedStats > & statsWrapper )
  {
    if( statsWrapper.referenceAsView().getFluxName() == fluxName )
    {
      lambda( statsWrapper.reference() );
    }
  } );
}

inline string SourceFluxStatsAggregator::getRegionStatDataName( string_view fluxName ) const
{ return GEOS_FMT( "{}_region_stats_for_{}", fluxName, getName() ); }


/******************************** SourceFluxStatisticsKernel ********************************/
struct SourceFluxStatisticsKernel
{};


} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_SOURCEFLUXSTATISTICS_HPP_ */
