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
 * @file CompositionalMultiphaseStatistics.cpp
 * @details Region statistics data is stored as follow:

 * Problem : ProblemManager
 * |-> domain : DomainPartition
 *     |-> MeshBodies : Group
 *         |-> cartesianMesh : MeshBody
 *             |-> meshLevels : Group
 *                 |-> Level0 : MeshLevel
 *                 |   |-> nodeManager : NodeManager
 *                 |   |   |-> sets : Group
 *                 |   |       | * all : Wrapper< index array >
 *                 |   |       | * xneg : Wrapper< index array >
 *                 |   |       [...] (other element sets)
 *                 |   |
 *                 |   |-> ElementRegions : ElementRegionManager
 *                 |   |   |-> Channel : CellElementRegion
 *                 |   |   |   |-> cb-0_0_0 : CellElementSubRegion
 *                 |   |   |   |   | * pressure : Wrapper< real64 array >
 *                 |   |   |   |   | * temperature : Wrapper< real64 array >
 *                 |   |   |   |   [...] (other fields)
 *                 |   |   |   |
 *                 |   |   |   |-> cb-0_0_1 : CellElementSubRegion
 *                 |   |   |   |   | * pressure : Wrapper< real64 array >
 *                 |   |   |   |   | * temperature : Wrapper< real64 array >
 *                 |   |   |   |   [...] (other fields)
 *                 |   |   |   |
 *                 |   |   |   [...] (other sub-regions)
 *                 |   |   |
 *                 |   |   |-> Barrier : CellElementRegion
 *                 |   |       |-> cb-1_0_0 : CellElementSubRegion
 *                 |   |       |-> cb-1_0_1 : CellElementSubRegion
 *                 |   |       [...] (other sub-regions)
 *                 |   |
 *                 |   [...] (other element managers)
 *          ____   |   |
 *          |      |   |-> statistics : Group (storage for all stats)
 *          |      |       |-> compFlowStats : Group (storage for this instance stats)
 *          |      |       |   |-> cflStatistics : CFLStatistics
 *          |      |       |   |-> regionsStatistics : RegionStatistics (aggregate)
 *          |      |       |       |-> Channel : RegionStatistics (aggregate, mpi reduced)
 *          |      |       |       |   |-> cb-0_0_0 : RegionStatistics (compute read-back)
 *  stats   |      |       |       |   |-> cb-0_0_1 : RegionStatistics (compute read-back)
 *  data -> |      |       |       |   [...] (other sub-regions stats)
 *          |      |       |       |
 *          |      |       |       |-> Barrier : RegionStatistics (aggregate, mpi reduced)
 *          |      |       |           |-> cb-1_0_0 : RegionStatistics (compute read-back)
 *          |      |       |           |-> cb-1_0_1 : RegionStatistics (compute read-back)
 *          |      |       |           [...] (other sub-regions stats)
 *          |      |       |
 *          |___   |       [...] (other stats storages)
 *                 |
 *                 [...] (other discretizations)
 */

#include "CompositionalMultiphaseStatisticsAggregator.hpp"

#include "physicsSolvers/StatisticsAggregatorBaseHelpers.hpp"
#include "mesh/DomainPartition.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/StatisticsKernel.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;

namespace compositionalMultiphaseStatistics
{

RegionStatistics::RegionStatistics( string const & name,
                                    dataRepository::Group * const parent,
                                    string_view setName,
                                    bool dataOutputEnabled,
                                    integer const numPhases,
                                    integer const numComponents ):
  RegionStatisticsBase( name, parent, setName, dataOutputEnabled ),
  m_phaseDynamicPoreVolume( numPhases ),
  m_phaseMass( numPhases ),
  m_trappedPhaseMass( numPhases ),
  m_nonTrappedPhaseMass( numPhases ),
  m_immobilePhaseMass( numPhases ),
  m_mobilePhaseMass( numPhases ),
  m_componentMass( numPhases, numComponents )
{
  // TODO : registerWrappers to store results in HDF5 (but need repairing of 1D HDF5 outputs)
}

CFLStatistics::CFLStatistics( string const & name,
                              dataRepository::Group * const parent,
                              bool const dataOutputEnabled ):
  RegionStatisticsBase( name,
                        parent,
                        "", // for now, CFL numbers are not for given sets
                        dataOutputEnabled )
{
  // TODO : registerWrappers to store results in HDF5 (but need repairing of 1D HDF5 outputs)
}

StatsAggregator::StatsAggregator( DataContext const & ownerDataContext,
                                  bool const dataOutputEnabled ):
  Base( ownerDataContext, dataOutputEnabled ),
  m_params()
{}

void StatsAggregator::initStatisticsAggregation( dataRepository::Group & meshBodies,
                                                 CompositionalMultiphaseBase & solver )
{
  m_solver = &solver;
  m_numPhases = solver.numFluidPhases();
  m_numComponents = solver.numFluidComponents();

  Base::initStatisticsAggregation( meshBodies, solver );
}

void StatsAggregator::enableRegionStatisticsAggregation( string_array const & setNames )
{
  auto const registerStats = [=] ( Group & parent,
                                   string const & targetName,
                                   string_view setName,
                                   bool const dataOutputEnabled ) -> RegionStatistics &
  {
    return parent.registerGroup( targetName,
                                 std::make_unique< RegionStatistics >( targetName,
                                                                       &parent,
                                                                       setName,
                                                                       dataOutputEnabled,
                                                                       m_numPhases,
                                                                       m_numComponents ) );
  };

  Base::enableRegionStatisticsAggregation( registerStats, setNames );
}

void StatsAggregator::enableCFLStatistics()
{
  if( m_solver == nullptr )
    return;

  m_solver->registerDataForCFL( *m_meshBodies );
  for( auto const & path : m_discretizationsPaths )
  {
    MeshLevel & mesh = getMeshLevel( path );
    Group & statisticsGroup = getInstanceStatisticsGroup( mesh );
    auto const cflStatsName = ViewKeys::cflStatisticsString();
    statisticsGroup.registerGroup< CFLStatistics >( cflStatsName,
                                                    std::make_unique< CFLStatistics >( cflStatsName,
                                                                                       &statisticsGroup,
                                                                                       m_dataOutputEnabled ) );
  }

  m_cflStatsState.m_isEnabled = true;
  m_cflStatsState.m_isDirty = true;
}

void StatsAggregator::setDirty()
{
  Base::setDirty();
  m_cflStatsState.m_isDirty = true;
}

bool StatsAggregator::computeCFLNumbers( real64 const time,
                                         real64 const dt,
                                         DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  real64 maxPhaseCFL, maxCompCFL;
  CFLStatistics * stats = getCFLStatistics( domain );

  m_warnings.clear();

  if( !m_cflStatsState.m_isEnabled )
  {
    m_warnings.emplace_back( "CFL numbers computation is not enabled." );
    return false;
  }
  if( stats == nullptr )
  {
    m_warnings.emplace_back( GEOS_FMT( "No statistics structure to compute CFL numbers for domain '{}'.", domain.getName() ));
    return false;
  }

  m_solver->computeCFLNumbers( domain, dt, maxPhaseCFL, maxCompCFL );

  stats->m_time = time;
  stats->m_maxPhaseCFL = maxPhaseCFL;
  stats->m_maxCompCFL = maxCompCFL;

  m_cflStatsState.m_isDirty = false;

  return true;
}

CFLStatistics * StatsAggregator::getCFLStatistics( DomainPartition & domain ) const
{
  return domain.getGroupPointer< CFLStatistics >( ViewKeys::cflStatisticsString() );
}

CFLStatistics & StatsAggregator::getCflStatistics( MeshLevel & mesh ) const
{
  // considering everything is initialized, or else, crash gracefully
  Group & statisticsGroup = getInstanceStatisticsGroup( mesh );
  return statisticsGroup.getGroup< CFLStatistics >( ViewKeys::cflStatisticsString() );
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

  stats.m_totalMass = 0.0;

  stats.m_totalPoreVolume = 0.0;
  stats.m_totalUncompactedPoreVolume = 0.0;
  stats.m_phaseDynamicPoreVolume.setValues< serialPolicy >( 0.0 );

  stats.m_phaseMass.setValues< serialPolicy >( 0.0 );
  stats.m_trappedPhaseMass.setValues< serialPolicy >( 0.0 );
  stats.m_nonTrappedPhaseMass.setValues< serialPolicy >( 0.0 );
  stats.m_immobilePhaseMass.setValues< serialPolicy >( 0.0 );
  stats.m_mobilePhaseMass.setValues< serialPolicy >( 0.0 );
  stats.m_componentMass.setValues< serialPolicy >( 0.0 );
}

void StatsAggregator::computeSubRegionRankStats( CellElementSubRegion & subRegion,
                                                 RegionStatistics & subRegionStats,
                                                 SetType const & targetSet ) const
{
  arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();
  arrayView1d< real64 const > const volume = subRegion.getElementVolume();
  arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
  arrayView1d< real64 const > const temp = subRegion.getField< fields::flow::temperature >();
  arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
    subRegion.getField< fields::flow::phaseVolumeFraction >();
  arrayView1d< real64 const > const deltaPres = subRegion.getField< fields::flow::deltaPressure >();

  Group const & constitutiveModels = subRegion.getGroup( ElementSubRegionBase::groupKeyStruct::constitutiveModelsString() );

  string const & solidName = subRegion.getReference< string >( CompositionalMultiphaseBase::viewKeyStruct::solidNamesString() );
  CoupledSolidBase const & solid = constitutiveModels.getGroup< CoupledSolidBase >( solidName );
  arrayView1d< real64 const > const refPorosity = solid.getReferencePorosity();
  arrayView2d< real64 const > const porosity = solid.getPorosity();

  string const & fluidName = subRegion.getReference< string >( CompositionalMultiphaseBase::viewKeyStruct::fluidNamesString() );
  MultiFluidBase const & fluid = constitutiveModels.getGroup< MultiFluidBase >( fluidName );
  arrayView3d< real64 const, multifluid::USD_PHASE > const phaseDensity = fluid.phaseDensity();
  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const phaseCompFraction = fluid.phaseCompFraction();

  //get min vol fraction for each phase to dispactche immobile/mobile mass
  string const & relpermName = subRegion.getReference< string >( CompositionalMultiphaseBase::viewKeyStruct::relPermNamesString() );
  RelativePermeabilityBase const & relperm = constitutiveModels.getGroup< RelativePermeabilityBase >( relpermName );
  arrayView3d< real64 const, relperm::USD_RELPERM > const phaseTrappedVolFrac = relperm.phaseTrappedVolFraction();
  arrayView3d< real64 const, relperm::USD_RELPERM > const phaseRelperm = relperm.phaseRelPerm();

  isothermalCompositionalMultiphaseBaseKernels::
    StatisticsKernel::
    launch< parallelDevicePolicy<> >( targetSet.toViewConst(),
                                      m_numComponents,
                                      m_numPhases,
                                      m_params.m_relpermThreshold,
                                      elemGhostRank,
                                      volume,
                                      pres,
                                      deltaPres,
                                      temp,
                                      refPorosity,
                                      porosity,
                                      phaseDensity,
                                      phaseCompFraction,
                                      phaseVolFrac,
                                      phaseTrappedVolFrac,
                                      phaseRelperm,
                                      subRegionStats.m_minPressure,
                                      subRegionStats.m_averagePressure,
                                      subRegionStats.m_maxPressure,
                                      subRegionStats.m_minDeltaPressure,
                                      subRegionStats.m_maxDeltaPressure,
                                      subRegionStats.m_minTemperature,
                                      subRegionStats.m_averageTemperature,
                                      subRegionStats.m_maxTemperature,
                                      subRegionStats.m_totalUncompactedPoreVolume,
                                      subRegionStats.m_phaseDynamicPoreVolume.toView(),
                                      subRegionStats.m_phaseMass.toView(),
                                      subRegionStats.m_trappedPhaseMass.toView(),
                                      subRegionStats.m_immobilePhaseMass.toView(),
                                      subRegionStats.m_componentMass.toView() );
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

  for( integer ip = 0; ip < m_numPhases; ++ip )
  {
    stats.m_phaseDynamicPoreVolume[ip] += other.m_phaseDynamicPoreVolume[ip];
    stats.m_phaseMass[ip] += other.m_phaseMass[ip];
    stats.m_trappedPhaseMass[ip] += other.m_trappedPhaseMass[ip];
    stats.m_immobilePhaseMass[ip] += other.m_immobilePhaseMass[ip];

    for( integer ic = 0; ic < m_numComponents; ++ic )
    {
      stats.m_componentMass[ip][ic] += other.m_componentMass[ip][ic];
    }
  }
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

  for( integer ip = 0; ip < m_numPhases; ++ip )
  {
    stats.m_phaseDynamicPoreVolume[ip] = MpiWrapper::sum( stats.m_phaseDynamicPoreVolume[ip] );
    stats.m_phaseMass[ip] = MpiWrapper::sum( stats.m_phaseMass[ip] );
    stats.m_trappedPhaseMass[ip] = MpiWrapper::sum( stats.m_trappedPhaseMass[ip] );
    stats.m_immobilePhaseMass[ip] = MpiWrapper::sum( stats.m_immobilePhaseMass[ip] );

    for( integer ic = 0; ic < m_numComponents; ++ic )
    {
      stats.m_componentMass[ip][ic] = MpiWrapper::sum( stats.m_componentMass[ip][ic] );
    }
  }
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

  for( integer ip = 0; ip < m_numPhases; ++ip )
  {
    stats.m_totalMass += stats.m_phaseMass[ip];
    stats.m_totalPoreVolume += stats.m_phaseDynamicPoreVolume[ip];
    stats.m_nonTrappedPhaseMass[ip] = stats.m_phaseMass[ip] - stats.m_trappedPhaseMass[ip];
    stats.m_mobilePhaseMass[ip] = stats.m_phaseMass[ip] - stats.m_immobilePhaseMass[ip];
  }
}

} /* namespace compositionalMultiphaseStatistics */

template class StatsAggregatorBase< compositionalMultiphaseStatistics::StatsAggregator >;

} /* namespace geos */
