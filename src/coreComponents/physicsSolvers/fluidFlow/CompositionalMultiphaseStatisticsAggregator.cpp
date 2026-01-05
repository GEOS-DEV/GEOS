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
 *
 * Problem : geos::ProblemManager
 * |-> domain : geos::DomainPartition
 *     |-> MeshBodies : geos::dataRepository::Group
 *         |-> cartesianMesh : geos::MeshBody
 *             |-> meshLevels : geos::dataRepository::Group
 *                 |-> Level0 : geos::MeshLevel
 *                 |   |-> nodeManager : geos::NodeManager
 *                 |   |   |-> sets : geos::dataRepository::Group
 *                 |   |       | * all : geos::dataRepository::Wrapper< index array >
 *                 |   |       | * xneg : geos::dataRepository::Wrapper< index array >
 *                 |   |       [...] (other element sets)
 *                 |   |
 *                 |   |-> ElementRegions : geos::ElementRegionManager
 *                 |   |   |-> Channel : geos::CellElementRegion
 *                 |   |   |   |-> cb-0_0_0 : geos::CellElementSubRegion
 *                 |   |   |   |   | * pressure : geos::dataRepository::Wrapper< real64 array >
 *                 |   |   |   |   | * temperature : geos::dataRepository::Wrapper< real64 array >
 *                 |   |   |   |   [...] (other fields)
 *                 |   |   |   |
 *                 |   |   |   |-> cb-0_0_1 : geos::CellElementSubRegion
 *                 |   |   |   |   | * pressure : geos::dataRepository::Wrapper< real64 array >
 *                 |   |   |   |   | * temperature : geos::dataRepository::Wrapper< real64 array >
 *                 |   |   |   |   [...] (other fields)
 *                 |   |   |   |
 *                 |   |   |   [...] (other sub-regions)
 *                 |   |   |
 *                 |   |   |-> Barrier : geos::CellElementRegion
 *                 |   |       |-> cb-1_0_0 : geos::CellElementSubRegion
 *                 |   |       |-> cb-1_0_1 : geos::CellElementSubRegion
 *                 |   |       [...] (other sub-regions)
 *                 |   |
 *                 |   [...] (other element managers)
 *          ____   |   |
 *          |      |   |-> statistics : geos::dataRepository::Group (storage for all stats)
 *          |      |       |-> compFlowStats : geos::dataRepository::Group (storage for this instance stats)
 *          |      |       |   |-> cflStatistics : geos::compositionalMultiphaseStatistics::CFLStatistics
 *          |      |       |   |-> regionsStats : geos::compositionalMultiphaseStatistics::RegionStatistics (aggregate)
 *          |      |       |       |-> Channel : geos::compositionalMultiphaseStatistics::RegionStatistics (aggregate, mpi reduced)
 *          |      |       |       |   |-> cb-0_0_0 : geos::compositionalMultiphaseStatistics::RegionStatistics (compute read-back)
 *  stats   |      |       |       |   |-> cb-0_0_1 : geos::compositionalMultiphaseStatistics::RegionStatistics (compute read-back)
 *  data -> |      |       |       |   [...] (other sub-regions stats)
 *          |      |       |       |
 *          |      |       |       |-> Barrier : geos::compositionalMultiphaseStatistics::RegionStatistics (aggregate, mpi reduced)
 *          |      |       |           |-> cb-1_0_0 : geos::compositionalMultiphaseStatistics::RegionStatistics (compute read-back)
 *          |      |       |           |-> cb-1_0_1 : geos::compositionalMultiphaseStatistics::RegionStatistics (compute read-back)
 *          |      |       |           [...] (other sub-regions stats)
 *          |      |       |
 *          |___   |       [...] (other stats storages)
 *                 |
 *                 [...] (other discretizations)
 *
 * not feasible because having a StatTask is not mandatory
 * Problem : geos::ProblemManager
 * |-> domain : geos::DomainPartition
 * |   |-> MeshBodies : geos::dataRepository::Group
 * |       |-> cartesianMesh : geos::MeshBody
 * |           |-> meshLevels : geos::dataRepository::Group
 * |               |-> Level0 : geos::MeshLevel
 * |               |   |-> nodeManager : geos::NodeManager
 * |               |   |   |-> sets : geos::dataRepository::Group
 * |               |   |       |-> all : geos::dataRepository::Wrapper< index array >
 * |               |   |       |-> xneg : geos::dataRepository::Wrapper< index array >
 * |               |   |       [...] (other element sets)
 * |               |   |
 * |               |   |-> ElementRegions : geos::ElementRegionManager
 * |               |   |   |-> Channel : geos::CellElementRegion
 * |               |   |   |   |-> cb-0_0_0 : geos::CellElementSubRegion
 * |               |   |   |   |   |-> pressure : geos::dataRepository::Wrapper< real64 array >
 * |               |   |   |   |   |-> temperature : geos::dataRepository::Wrapper< real64 array >
 * |               |   |   |   |   [...] (other fields)
 * |               |   |   |   |
 * |               |   |   |   |-> cb-0_0_1 : geos::CellElementSubRegion
 * |               |   |   |   |   |-> pressure : geos::dataRepository::Wrapper< real64 array >
 * |               |   |   |   |   |-> temperature : geos::dataRepository::Wrapper< real64 array >
 * |               |   |   |   |   [...] (other fields)
 * |               |   |   |   |
 * |               |   |   |   [...] (other sub-regions)
 * |               |   |   |
 * |               |   |   |-> Barrier : geos::CellElementRegion
 * |               |   |       |-> cb-1_0_0 : geos::CellElementSubRegion
 * |               |   |       |-> cb-1_0_1 : geos::CellElementSubRegion
 * |               |   |       [...] (other sub-regions)
 * |               |   |
 * |               |   [...] (other element managers)
 * |               |
 * |               [...] (other discretizations)
 * |               
 * |-> Tasks : geos::TaskManager
 *     |-> compFlowStats : geos::dataRepository::compositionalMultiphaseStatistics::StatTask (storage for this instance stats)
 *         |-> cartesianMesh : geos::dataRepository::Group
 *             |-> Level0 : geos::dataRepository::Group
 *          ____   |-> cartesianMesh_Level0 : geos::dataRepository::Group (storage for all stats)
 *          |      |   |   |-> cflStatistics : geos::compositionalMultiphaseStatistics::CFLStatistics
 *          |      |   |   |-> regionsStats : geos::compositionalMultiphaseStatistics::RegionStatistics (aggregate)
 *          |      |   |       |-> Channel : geos::compositionalMultiphaseStatistics::RegionStatistics (aggregate, mpi reduced)
 *          |      |   |       |   |-> cb-0_0_0 : geos::compositionalMultiphaseStatistics::RegionStatistics (compute read-back)
 *  stats   |      |   |       |   |-> cb-0_0_1 : geos::compositionalMultiphaseStatistics::RegionStatistics (compute read-back)
 *  data -> |      |   |       |   ... (other sub-regions stats)
 *          |      |   |       |
 *          |      |   |       |-> Barrier : geos::compositionalMultiphaseStatistics::RegionStatistics (aggregate, mpi reduced)
 *          |      |   |           |-> cb-1_0_0 : geos::compositionalMultiphaseStatistics::RegionStatistics (compute read-back)
 *          |      |   |           |-> cb-1_0_1 : geos::compositionalMultiphaseStatistics::RegionStatistics (compute read-back)
 *          |      |   |           ... (other sub-regions stats)
 *          |      |   |
 *          |___   |   ... (other stats storages)
 *                 |
 *                 ... (other discretizations)
 *
 * Stats data structures are much more scattered + the "statistics" Group is still needed
 * Problem : geos::ProblemManager
 * |-> domain : geos::DomainPartition
 *     |-> MeshBodies : geos::dataRepository::Group
 *         |-> cartesianMesh : geos::MeshBody
 *             |-> meshLevels : geos::dataRepository::Group
 *                 |-> Level0 : geos::MeshLevel
 *                 |   |-> nodeManager : geos::NodeManager
 *                 |   |   |-> sets : geos::dataRepository::Group
 *                 |   |       | * all : geos::dataRepository::Wrapper< index array >
 *                 |   |       | * xneg : geos::dataRepository::Wrapper< index array >
 *                 |   |       ... (other element sets)
 *                 |   |
 *                 |   |-> other element managers...
 *                 |   |
 *                 |   |-> ElementRegions : geos::ElementRegionManager
 *                 |       |-> Channel : geos::CellElementRegion
 *                 |       |   |-> cb-0_0_0 : geos::CellElementSubRegion
 *                 |       |   |   | * pressure : geos::dataRepository::Wrapper< real64 array >
 *                 |       |   |   | * temperature : geos::dataRepository::Wrapper< real64 array >
 *                 |       |   |   | [...] (other fields)
 *                 |       |   |   |-> compFlowStats : geos::compositionalMultiphaseStatistics::RegionStatistics (compute read-back)
 *                 |       |   |
 *                 |       |   |-> cb-0_0_1 : geos::CellElementSubRegion
 *                 |       |   |   | * pressure : geos::dataRepository::Wrapper< real64 array >
 *                 |       |   |   | * temperature : geos::dataRepository::Wrapper< real64 array >
 *                 |       |   |   | [...] (other fields)
 *                 |       |   |   |-> compFlowStats : geos::compositionalMultiphaseStatistics::RegionStatistics (compute read-back)
 *                 |       |   |
 *                 |       |   [...] (other sub-regions)
 *                 |       |   |
 *                 |       |   |-> compFlowStats : geos::compositionalMultiphaseStatistics::RegionStatistics (aggregate, mpi reduced)
 *                 |       |
 *                 |       |-> Barrier : geos::CellElementRegion
 *                 |       |   |-> cb-1_0_0 : geos::CellElementSubRegion
 *                 |       |   |-> cb-1_0_1 : geos::CellElementSubRegion
 *                 |       |   [...] (other sub-regions)
 *                 |       |
 *                 |       [...] (other regions)
 *                 |       |
 *          |      |       |-> statistics : geos::dataRepository::Group (storage for all stats)
 *                 |           |-> compFlowStats : geos::dataRepository::Group (storage for this instance stats)
 *                 |               |-> regionsStats : geos::compositionalMultiphaseStatistics::RegionStatistics (aggregate)
 *                 |               |-> cflStatistics : geos::compositionalMultiphaseStatistics::CFLStatistics
 *                 |
 *                 ... (other discretizations)
 */

#include "CompositionalMultiphaseStatisticsAggregator.hpp"

#include "common/DataTypes.hpp"
#include "common/format/Format.hpp"
#include "common/logger/Logger.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/StatisticsKernel.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include <ostream>


namespace geos
{

namespace compositionalMultiphaseStatistics
{

using namespace constitutive;
// using namespace fields;
using namespace dataRepository;

StatsAggregator::StatsAggregator():
  m_params(),
  m_isRegionStatsEnabled( false ),
  m_isCFLNumberEnabled( false )
{}

void StatsAggregator::enableRegionStatistics( dataRepository::Group & /*meshBodies*/ )
{
  GEOS_ERROR_IF_EQ_MSG( m_solver, nullptr, "Flow solver must be set." );

  // integer const numPhases = m_solver->numFluidPhases();
  // integer const numComps = m_solver->numFluidComponents();

  // m_solver->forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
  //                                                             MeshLevel & mesh,
  //                                                             string_array const & regionNames )
  // {
  //   ElementRegionManager & elemManager = mesh.getElemManager();

  //   for( size_t i = 0; i < regionNames.size(); ++i )
  //   {
  //     // ElementRegionBase & region = elemManager.getRegion( regionNames[i] );

  //     // RELOCALISER LES DONNEES region.registerWrapper< RegionStatistics >( VKStruct::regionStatisticsString() ).
  //     //   setRestartFlags( RestartFlags::NO_WRITE );
  //     // region.excludeWrappersFromPacking( { VKStruct::regionStatisticsString() } );
  //     // RegionStatistics & stats = region.getReference< RegionStatistics >( VKStruct::regionStatisticsString() );

  //     // stats.m_phasePoreVolume.resizeDimension< 0 >( numPhases );
  //     // stats.m_phaseMass.resizeDimension< 0 >( numPhases );
  //     // stats.m_trappedPhaseMass.resizeDimension< 0 >( numPhases );
  //     // stats.m_nonTrappedPhaseMass.resizeDimension< 0 >( numPhases );
  //     // stats.m_immobilePhaseMass.resizeDimension< 0 >( numPhases );
  //     // stats.m_mobilePhaseMass.resizeDimension< 0 >( numPhases );
  //     // stats.m_componentMass.resizeDimension< 0, 1 >( numPhases, numComps );
  //   }
  // } );
  // m_isRegionStatsEnabled = true;
}

void StatsAggregator::forRegionStatistics( MeshLevel const & mesh,
                                           StatsAggregator::RegionFunctor functor ) const
{
  GEOS_UNUSED_VAR( mesh );
  GEOS_UNUSED_VAR( functor );
}

void StatsAggregator::enableCFLStatistics( dataRepository::Group & meshBodies )
{
  GEOS_ERROR_IF_EQ_MSG( m_solver, nullptr, "Flow solver must be set." );

  m_solver->registerDataForCFL( meshBodies );
  m_isCFLNumberEnabled = true;
}

bool StatsAggregator::computeRegionsStatistics( real64 const time,
                                                MeshLevel & mesh,
                                                string_array const & regionNames )
{
  GEOS_MARK_FUNCTION;

  integer const numPhases = m_solver->numFluidPhases();
  integer const numComps = m_solver->numFluidComponents();

  // Step 1: initialize the average/min/max quantities
  ElementRegionManager & elemManager = mesh.getElemManager();
  for( size_t i = 0; i < regionNames.size(); ++i )
  {
    ElementRegionBase & region = elemManager.getRegion( regionNames[i] );
    RegionStatistics & stats = region.getReference< RegionStatistics >( VKStruct::regionStatisticsString() );

    stats.m_time = time;

    stats.m_averagePressure = 0.0;
    stats.m_maxPressure = 0.0;
    stats.m_minPressure = LvArray::NumericLimits< real64 >::max;

    stats.m_maxDeltaPressure = -LvArray::NumericLimits< real64 >::max;
    stats.m_minDeltaPressure = LvArray::NumericLimits< real64 >::max;

    stats.m_averageTemperature = 0.0;
    stats.m_maxTemperature = 0.0;
    stats.m_minTemperature = LvArray::NumericLimits< real64 >::max;

    stats.m_totalPoreVolume = 0.0;
    stats.m_totalUncompactedPoreVolume = 0.0;
    stats.m_phasePoreVolume.setValues< serialPolicy >( 0.0 );

    stats.m_phaseMass.setValues< serialPolicy >( 0.0 );
    stats.m_trappedPhaseMass.setValues< serialPolicy >( 0.0 );
    stats.m_nonTrappedPhaseMass.setValues< serialPolicy >( 0.0 );
    stats.m_immobilePhaseMass.setValues< serialPolicy >( 0.0 );
    stats.m_mobilePhaseMass.setValues< serialPolicy >( 0.0 );
    stats.m_componentMass.setValues< serialPolicy >( 0.0 );
  }

  // Step 2: increment the average/min/max quantities for all the subRegions
  elemManager.forElementSubRegions( regionNames, [&]( localIndex const,
                                                      ElementSubRegionBase & subRegion )
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

    real64 subRegionAvgPresNumerator = 0.0;
    real64 subRegionMinPres = 0.0;
    real64 subRegionMaxPres = 0.0;
    real64 subRegionMinDeltaPres = 0.0;
    real64 subRegionMaxDeltaPres = 0.0;
    real64 subRegionAvgTempNumerator = 0.0;
    real64 subRegionMinTemp = 0.0;
    real64 subRegionMaxTemp = 0.0;
    real64 subRegionTotalUncompactedPoreVol = 0.0;
    array1d< real64 > subRegionPhaseDynamicPoreVol( numPhases );
    array1d< real64 > subRegionPhaseMass( numPhases );
    array1d< real64 > subRegionTrappedPhaseMass( numPhases );
    array1d< real64 > subRegionImmobilePhaseMass( numPhases );
    array1d< real64 > subRegionRelpermPhaseMass( numPhases );
    array2d< real64 > subRegionComponentMass( numPhases, numComps );

    isothermalCompositionalMultiphaseBaseKernels::
      StatisticsKernel::
      launch< parallelDevicePolicy<> >( subRegion.size(),
                                        numComps,
                                        numPhases,
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
                                        subRegionMinPres,
                                        subRegionAvgPresNumerator,
                                        subRegionMaxPres,
                                        subRegionMinDeltaPres,
                                        subRegionMaxDeltaPres,
                                        subRegionMinTemp,
                                        subRegionAvgTempNumerator,
                                        subRegionMaxTemp,
                                        subRegionTotalUncompactedPoreVol,
                                        subRegionPhaseDynamicPoreVol.toView(),
                                        subRegionPhaseMass.toView(),
                                        subRegionTrappedPhaseMass.toView(),
                                        subRegionImmobilePhaseMass.toView(),
                                        subRegionComponentMass.toView() );

    ElementRegionBase & region = elemManager.getRegion( ElementRegionBase::getParentRegion( subRegion ).getName() );
    RegionStatistics & stats = region.getReference< RegionStatistics >( VKStruct::regionStatisticsString() );

    stats.m_averagePressure += subRegionAvgPresNumerator;
    if( subRegionMinPres < stats.m_minPressure )
    {
      stats.m_minPressure = subRegionMinPres;
    }
    if( subRegionMaxPres > stats.m_maxPressure )
    {
      stats.m_maxPressure = subRegionMaxPres;
    }

    if( subRegionMinDeltaPres < stats.m_minDeltaPressure )
    {
      stats.m_minDeltaPressure = subRegionMinDeltaPres;
    }
    if( subRegionMaxDeltaPres > stats.m_maxDeltaPressure )
    {
      stats.m_maxDeltaPressure = subRegionMaxDeltaPres;
    }

    stats.m_averageTemperature += subRegionAvgTempNumerator;
    if( subRegionMinTemp < stats.m_minTemperature )
    {
      stats.m_minTemperature = subRegionMinTemp;
    }
    if( subRegionMaxTemp > stats.m_maxTemperature )
    {
      stats.m_maxTemperature = subRegionMaxTemp;
    }

    stats.m_totalUncompactedPoreVolume += subRegionTotalUncompactedPoreVol;
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      stats.m_phasePoreVolume[ip] += subRegionPhaseDynamicPoreVol[ip];
      stats.m_phaseMass[ip] += subRegionPhaseMass[ip];
      stats.m_trappedPhaseMass[ip] += subRegionTrappedPhaseMass[ip];
      stats.m_immobilePhaseMass[ip] += subRegionImmobilePhaseMass[ip];

      for( integer ic = 0; ic < numComps; ++ic )
      {
        stats.m_componentMass[ip][ic] += subRegionComponentMass[ip][ic];
      }
    }

  } );

  // Step 3: synchronize the results over the MPI ranks
  for( size_t i = 0; i < regionNames.size(); ++i )
  {
    ElementRegionBase & region = elemManager.getRegion( regionNames[i] );
    RegionStatistics & stats = region.getReference< RegionStatistics >( VKStruct::regionStatisticsString() );

    stats.m_minPressure = MpiWrapper::min( stats.m_minPressure );
    stats.m_maxPressure = MpiWrapper::max( stats.m_maxPressure );
    stats.m_minDeltaPressure = MpiWrapper::min( stats.m_minDeltaPressure );
    stats.m_maxDeltaPressure = MpiWrapper::max( stats.m_maxDeltaPressure );
    stats.m_minTemperature = MpiWrapper::min( stats.m_minTemperature );
    stats.m_maxTemperature = MpiWrapper::max( stats.m_maxTemperature );
    stats.m_totalUncompactedPoreVolume = MpiWrapper::sum( stats.m_totalUncompactedPoreVolume );
    stats.m_totalPoreVolume = 0.0;
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      stats.m_phasePoreVolume[ip] = MpiWrapper::sum( stats.m_phasePoreVolume[ip] );
      stats.m_phaseMass[ip] = MpiWrapper::sum( stats.m_phaseMass[ip] );
      stats.m_trappedPhaseMass[ip] = MpiWrapper::sum( stats.m_trappedPhaseMass[ip] );
      stats.m_immobilePhaseMass[ip] = MpiWrapper::sum( stats.m_immobilePhaseMass[ip] );
      stats.m_totalPoreVolume += stats.m_phasePoreVolume[ip];
      for( integer ic = 0; ic < numComps; ++ic )
      {
        stats.m_componentMass[ip][ic] = MpiWrapper::sum( stats.m_componentMass[ip][ic] );
      }
    }
    stats.m_averagePressure = MpiWrapper::sum( stats.m_averagePressure );
    stats.m_averageTemperature = MpiWrapper::sum( stats.m_averageTemperature );
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
      m_issues.emplace_back( GEOS_FMT( "Cannot compute average pressure because region pore volume is zero in region '{}'.",
                                       regionNames[i] ) );
    }

    for( integer ip = 0; ip < numPhases; ++ip )
    {
      stats.m_nonTrappedPhaseMass[ip] = stats.m_phaseMass[ip] - stats.m_trappedPhaseMass[ip];
      stats.m_mobilePhaseMass[ip] = stats.m_phaseMass[ip] - stats.m_immobilePhaseMass[ip];
    }

    stdVector< string > phaseCompName;
    phaseCompName.reserve( numPhases*numComps );
    stdVector< string > massValues;
    phaseCompName.reserve( numPhases*numComps );

    ConstitutiveManager const & constitutiveManager = mesh.getGroupByPath< ConstitutiveManager >( "/Problem/domain/Constitutive" );
    MultiFluidBase const & fluid = constitutiveManager.getGroup< MultiFluidBase >( m_solver->referenceFluidModelName() );
    auto const phaseNames = fluid.phaseNames();
    auto const componentNames = fluid.componentNames();
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      for( integer ic = 0; ic < numComps; ++ic )
      {
        std::stringstream ss;
        ss << phaseNames[ip] << "/" << componentNames[ic];
        phaseCompName.push_back( ss.str() );
        massValues.push_back( GEOS_FMT( "{}", stats.m_componentMass[ip][ic] ) );
      }
    }
  }

  return true;
}

bool StatsAggregator::computeCFLNumbers( real64 const time,
                                         real64 const dt,
                                         DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  real64 maxPhaseCFL, maxCompCFL;
  CFLStatistics * stats = getCFLStatistics( domain );

  m_issues.clear();

  if( stats!=nullptr )
  {
    m_issues.emplace_back( GEOS_FMT( "No statistics structure to compute CFL numbers for domain '{}'.", domain.getName() ));
    return false;
  }

  m_solver->computeCFLNumbers( domain, dt, maxPhaseCFL, maxCompCFL );

  stats->m_time = time;
  stats->m_maxPhaseCFL = maxPhaseCFL;
  stats->m_maxCompCFL = maxCompCFL;

  return true;
}


} /* namespace compositionalMultiphaseStatistics */

} /* namespace geos */
