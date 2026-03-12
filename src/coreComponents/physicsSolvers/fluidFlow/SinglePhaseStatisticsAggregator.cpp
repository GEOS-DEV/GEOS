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
 * @file SinglePhaseStatistics.cpp
 */

#include "SinglePhaseStatisticsAggregator.hpp"

#include "mesh/DomainPartition.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "physicsSolvers/StatisticsAggregatorBaseHelpers.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/StatisticsKernel.hpp"
#include "common/format/table/TableData.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "common/format/table/TableLayout.hpp"

namespace geos
{

namespace singlePhaseStatistics
{

using namespace constitutive;
using namespace dataRepository;

RegionStatistics::RegionStatistics( string const & name,
                                    dataRepository::Group * const parent,
                                    bool const statsOutputEnabled ):
  RegionStatisticsBase( name, parent, statsOutputEnabled )
{}

StatsAggregator::StatsAggregator( DataContext const & ownerDataContext,
                                  bool const statsOutputEnabled ):
  Base( ownerDataContext, statsOutputEnabled )
{}

void StatsAggregator::enableRegionStatisticsAggregation()
{
  auto const registerStats = [=] ( Group & parent,
                                   string const & targetName ) -> RegionStatistics &
  {
    return parent.registerGroup( targetName,
                                 std::make_unique< RegionStatistics >( targetName,
                                                                       &parent,
                                                                       m_statsOutputEnabled ) );
  };

  Base::enableRegionStatisticsAggregation( registerStats );
}

void StatsAggregator::initStats( RegionStatistics & stats, real64 const time ) const
{
  stats.m_time = time;

  stats.m_averagePressure = 0.0;
  stats.m_maxPressure = 0.0;
  stats.m_minPressure = LvArray::NumericLimits< real64 >::max;

  stats.m_maxDeltaPressure = -LvArray::NumericLimits< real64 >::max;
  stats.m_minDeltaPressure = LvArray::NumericLimits< real64 >::max;

  stats.m_averageTemperature = 0.0;
  stats.m_maxTemperature = 0.0;
  stats.m_minTemperature = LvArray::NumericLimits< real64 >::max;

  stats.m_totalDynamicPoreVolume = 0.0;
  stats.m_totalUncompactedPoreVolume = 0.0;

  stats.m_totalMass = 0.0;
}

void StatsAggregator::computeSubRegionRankStats( CellElementSubRegion & subRegion,
                                                 RegionStatistics & subRegionStats ) const
{
  static constexpr string_view solidNamesVK = SinglePhaseBase::viewKeyStruct::solidNamesString();
  static constexpr string_view fluidNamesVK = FlowSolverBase::viewKeyStruct::fluidNamesString();
  static constexpr string_view modelsVK = ElementSubRegionBase::groupKeyStruct::constitutiveModelsString();

  arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();
  arrayView1d< real64 const > const volume = subRegion.getElementVolume();
  arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
  arrayView1d< real64 const > const deltaPres = subRegion.getField< fields::flow::deltaPressure >();
  arrayView1d< real64 const > const temp = subRegion.getField< fields::flow::temperature >();

  string const & solidName = subRegion.getReference< string >( string( solidNamesVK ) );
  Group const & constitutiveModels = subRegion.getGroup( string( modelsVK ) );
  CoupledSolidBase const & solid = constitutiveModels.getGroup< CoupledSolidBase >( solidName );
  arrayView1d< real64 const > const refPorosity = solid.getReferencePorosity();
  arrayView2d< real64 const > const porosity = solid.getPorosity();

  string const & fluidName = subRegion.template getReference< string >( string( fluidNamesVK ) );
  SingleFluidBase const & fluid = constitutiveModels.getGroup< SingleFluidBase >( fluidName );
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const densities = fluid.density();

  singlePhaseBaseKernels::StatisticsKernel::launch( subRegion.size(),
                                                    elemGhostRank,
                                                    volume,
                                                    pres,
                                                    deltaPres,
                                                    temp,
                                                    refPorosity,
                                                    porosity,
                                                    densities,
                                                    subRegionStats.m_minPressure,
                                                    subRegionStats.m_averagePressure,
                                                    subRegionStats.m_maxPressure,
                                                    subRegionStats.m_minDeltaPressure,
                                                    subRegionStats.m_maxDeltaPressure,
                                                    subRegionStats.m_minTemperature,
                                                    subRegionStats.m_averageTemperature,
                                                    subRegionStats.m_maxTemperature,
                                                    subRegionStats.m_totalUncompactedPoreVolume,
                                                    subRegionStats.m_totalDynamicPoreVolume,
                                                    subRegionStats.m_totalMass );
}


void StatsAggregator::aggregateStats( RegionStatistics & stats,
                                      RegionStatistics const & other ) const
{
  stats.m_averagePressure += other.m_averagePressure;
  stats.m_minPressure = LvArray::math::min( stats.m_minPressure, other.m_minPressure );
  stats.m_maxPressure = LvArray::math::max( stats.m_maxPressure, other.m_maxPressure );

  stats.m_minDeltaPressure = LvArray::math::min( stats.m_minDeltaPressure, other.m_minDeltaPressure );
  stats.m_maxDeltaPressure = LvArray::math::max( stats.m_maxDeltaPressure, other.m_maxDeltaPressure );

  stats.m_averageTemperature += other.m_averageTemperature;
  stats.m_minTemperature = LvArray::math::min( stats.m_minTemperature, other.m_minTemperature );
  stats.m_maxTemperature = LvArray::math::max( stats.m_maxTemperature, other.m_maxTemperature );

  stats.m_totalUncompactedPoreVolume += other.m_totalUncompactedPoreVolume;
  stats.m_totalDynamicPoreVolume += other.m_totalDynamicPoreVolume;

  stats.m_totalMass += other.m_totalMass;
}

void StatsAggregator::mpiAggregateStats( RegionStatistics & stats ) const
{
  stats.m_averagePressure = MpiWrapper::sum( stats.m_averagePressure );
  stats.m_minPressure = MpiWrapper::min( stats.m_minPressure );
  stats.m_maxPressure = MpiWrapper::max( stats.m_maxPressure );

  stats.m_minDeltaPressure = MpiWrapper::min( stats.m_minDeltaPressure );
  stats.m_maxDeltaPressure = MpiWrapper::max( stats.m_maxDeltaPressure );

  stats.m_averageTemperature = MpiWrapper::sum( stats.m_averageTemperature );
  stats.m_minTemperature = MpiWrapper::min( stats.m_minTemperature );
  stats.m_maxTemperature = MpiWrapper::max( stats.m_maxTemperature );

  stats.m_totalUncompactedPoreVolume = MpiWrapper::sum( stats.m_totalUncompactedPoreVolume );
  stats.m_totalDynamicPoreVolume = MpiWrapper::sum( stats.m_totalDynamicPoreVolume );

  stats.m_totalMass = MpiWrapper::sum( stats.m_totalMass );
}

void StatsAggregator::postAggregateStats( RegionStatistics & stats )
{
  if( stats.m_totalUncompactedPoreVolume > 0 )
  {
    float invTotalUncompactedPoreVolume = 1.0 / stats.m_totalUncompactedPoreVolume;
    stats.m_averagePressure *= invTotalUncompactedPoreVolume;
    stats.m_averageTemperature *= invTotalUncompactedPoreVolume;
  }
  else
  {
    stats.m_averagePressure = 0.0;
    stats.m_averageTemperature = 0.0;
    m_warnings.emplace_back( GEOS_FMT( "Cannot compute average pressure for '{}' because pore volume is zero in '{}'.",
                                       getOwnerName(), stats.getTargetName() ) );
  }
}

} /* namespace singlePhaseStatistics */

template class StatsAggregatorBase< singlePhaseStatistics::StatsAggregator >;

} /* namespace geos */
