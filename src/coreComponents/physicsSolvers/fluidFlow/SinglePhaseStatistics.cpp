/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
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

#include "SinglePhaseStatistics.hpp"

#include "mainInterface/ProblemManager.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/LogLevelsInfo.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;

SinglePhaseStatistics::SinglePhaseStatistics( const string & name,
                                              Group * const parent ):
  Base( name, parent )
{
  addLogLevel< logInfo::Statistics >();
}

void SinglePhaseStatistics::registerDataOnMesh( Group & meshBodies )
{
  // the fields have to be registered in "registerDataOnMesh" (and not later)
  // otherwise they cannot be targeted by TimeHistory

  // for now, this guard is needed to avoid breaking the xml schema generation
  if( m_solver == nullptr )
  {
    return;
  }

  m_solver->forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                              MeshLevel & mesh,
                                                              arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    for( integer i = 0; i < regionNames.size(); ++i )
    {
      ElementRegionBase & region = elemManager.getRegion( regionNames[i] );
      region.registerGroup< RegionStatistics >( viewKeyStruct::regionStatisticsString() ).
        setRestartFlags( RestartFlags::NO_WRITE );

      // write output header
      if( m_writeCSV > 0 && MpiWrapper::commRank() == 0 )
      {
        std::ofstream outputFile( m_outputDir + "/" + regionNames[i] + ".csv" );
        outputFile <<
          "Time [s],Min pressure [Pa],Average pressure [Pa],Max pressure [Pa],Min delta pressure [Pa],Max delta pressure [Pa]," <<
          "Min temperature [Pa],Average temperature [Pa],Max temperature [Pa],Total dynamic pore volume [rm^3],Total fluid mass [kg]";
        outputFile << std::endl;
        outputFile.close();
      }
    }
  } );
}

bool SinglePhaseStatistics::execute( real64 const time_n,
                                     real64 const dt,
                                     integer const GEOS_UNUSED_PARAM( cycleNumber ),
                                     integer const GEOS_UNUSED_PARAM( eventCounter ),
                                     real64 const GEOS_UNUSED_PARAM( eventProgress ),
                                     DomainPartition & domain )
{
  m_solver->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                          MeshLevel & mesh,
                                                                          arrayView1d< string const > const & regionNames )
  {
    // current time is time_n + dt
    computeRegionStatistics( time_n + dt, mesh, regionNames );
  } );
  return false;
}

void SinglePhaseStatistics::computeRegionStatistics( real64 const time,
                                                     MeshLevel & mesh,
                                                     arrayView1d< string const > const & regionNames )
{
  GEOS_MARK_FUNCTION;

  // Step 1: initialize the average/min/max quantities
  ElementRegionManager & elemManager = mesh.getElemManager();
  for( integer i = 0; i < regionNames.size(); ++i )
  {
    ElementRegionBase & region = elemManager.getRegion( regionNames[i] );

    RegionStatistics & regionStatistics = region.getGroup< RegionStatistics >( viewKeyStruct::regionStatisticsString() );

    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::averagePressureString()).setApplyDefaultValue( 0.0 );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::maxPressureString()).setApplyDefaultValue( 0.0 );
    constexpr real64 max = LvArray::NumericLimits< real64 >::max;
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::minPressureString()).setApplyDefaultValue( max );

    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::maxDeltaPressureString()).setApplyDefaultValue( -max );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::minDeltaPressureString()).setApplyDefaultValue( max );

    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::totalMassString()).setApplyDefaultValue( 0.0 );

    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::averageTemperatureString()).setApplyDefaultValue( 0.0 );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::maxTemperatureString()).setApplyDefaultValue( 0.0 );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::minTemperatureString()).setApplyDefaultValue( max );

    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::totalPoreVolumeString()).setApplyDefaultValue( 0.0 );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::totalUncompactedPoreVolumeString()).setApplyDefaultValue( 0.0 );
  }

  // Step 2: increment the average/min/max quantities for all the subRegions
  elemManager.forElementSubRegions( regionNames, [&]( localIndex const,
                                                      ElementSubRegionBase & subRegion )
  {

    arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();
    arrayView1d< real64 const > const volume = subRegion.getElementVolume();
    arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
    arrayView1d< real64 const > const deltaPres = subRegion.getField< fields::flow::deltaPressure >();
    arrayView1d< real64 const > const temp = subRegion.getField< fields::flow::temperature >();

    string const & solidName = subRegion.getReference< string >( SinglePhaseBase::viewKeyStruct::solidNamesString() );
    Group const & constitutiveModels = subRegion.getGroup( ElementSubRegionBase::groupKeyStruct::constitutiveModelsString() );
    CoupledSolidBase const & solid = constitutiveModels.getGroup< CoupledSolidBase >( solidName );
    arrayView1d< real64 const > const refPorosity = solid.getReferencePorosity();
    arrayView2d< real64 const > const porosity = solid.getPorosity();

    string const & fluidName = subRegion.template getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() );
    SingleFluidBase const & fluid = constitutiveModels.getGroup< SingleFluidBase >( fluidName );
    arrayView2d< real64 const > const densities = fluid.density();

    real64 subRegionAvgPresNumerator = 0.0;
    real64 subRegionMinPres = 0.0;
    real64 subRegionMaxPres = 0.0;
    real64 subRegionMinDeltaPres = 0.0;
    real64 subRegionMaxDeltaPres = 0.0;
    real64 subRegionAvgTempNumerator = 0.0;
    real64 subRegionMinTemp = 0.0;
    real64 subRegionMaxTemp = 0.0;
    real64 subRegionTotalUncompactedPoreVol = 0.0;
    real64 subRegionTotalPoreVol = 0.0;
    real64 subRegionTotalMass = 0.0;

    singlePhaseBaseKernels::StatisticsKernel::
      launch( subRegion.size(),
              elemGhostRank,
              volume,
              pres,
              deltaPres,
              temp,
              refPorosity,
              porosity,
              densities,
              subRegionMinPres,
              subRegionAvgPresNumerator,
              subRegionMaxPres,
              subRegionMinDeltaPres,
              subRegionMaxDeltaPres,
              subRegionMinTemp,
              subRegionAvgTempNumerator,
              subRegionMaxTemp,
              subRegionTotalUncompactedPoreVol,
              subRegionTotalPoreVol,
              subRegionTotalMass );

    ElementRegionBase & region = elemManager.getRegion( ElementRegionBase::getParentRegion( subRegion ).getName() );
    RegionStatistics & regionStatistics = region.getGroup< RegionStatistics >( viewKeyStruct::regionStatisticsString() );

    real64 & averagePressure = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::averagePressureString());
    averagePressure += subRegionAvgPresNumerator;
    real64 & minPressure = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::minPressureString());
    if( subRegionMinPres < minPressure )
    {
      regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::minPressureString()).setApplyDefaultValue( subRegionMinPres );
    }
    real64 & maxPressure = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::maxPressureString());
    if( subRegionMaxPres > maxPressure )
    {
      regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::maxPressureString()).setApplyDefaultValue( subRegionMaxPres );
    }
    real64 & minDeltaPressure = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::minDeltaPressureString());
    if( subRegionMinDeltaPres < minDeltaPressure )
    {
      regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::minDeltaPressureString()).setApplyDefaultValue( subRegionMinDeltaPres );
    }
    real64 & maxDeltaPressure = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::maxDeltaPressureString());
    if( subRegionMaxDeltaPres > maxDeltaPressure )
    {
      regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::maxDeltaPressureString()).setApplyDefaultValue( subRegionMaxDeltaPres );
    }
    real64 & averageTemperature = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::averageTemperatureString());
    averageTemperature += subRegionAvgTempNumerator;
    real64 & minTemperature = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::minTemperatureString());
    if( subRegionMinTemp < minTemperature )
    {
      regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::minTemperatureString()).setApplyDefaultValue( subRegionMinTemp );
    }
    real64 & maxTemperature = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::maxTemperatureString());
    if( subRegionMaxTemp > maxTemperature )
    {
      regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::maxTemperatureString()).setApplyDefaultValue( subRegionMaxTemp );
    }

    real64 & totalUncompactedPoreVolume = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::totalPoreVolumeString());
    totalUncompactedPoreVolume += subRegionAvgPresNumerator;
    real64 & totalPoreVolume = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::totalUncompactedPoreVolumeString());
    totalPoreVolume += subRegionTotalPoreVol;
    real64 & totalMass = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::totalMassString());
    totalMass += subRegionTotalMass;
  } );

  // Step 3: synchronize the results over the MPI ranks
  for( integer i = 0; i < regionNames.size(); ++i )
  {
    ElementRegionBase & region = elemManager.getRegion( regionNames[i] );
    RegionStatistics & regionStatistics = region.getGroup< RegionStatistics >( viewKeyStruct::regionStatisticsString() );


    real64 minPressure = MpiWrapper::min( regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::minPressureString()) );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::minPressureString()).setApplyDefaultValue( minPressure );
    real64 maxPressure = MpiWrapper::max( regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::maxPressureString()) );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::maxPressureString()).setApplyDefaultValue( maxPressure );
    real64 averagePressure = MpiWrapper::sum( regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::averagePressureString()) );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::averagePressureString()).setApplyDefaultValue( averagePressure );

    real64 minDeltaPressure = MpiWrapper::min( regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::minDeltaPressureString()) );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::minDeltaPressureString()).setApplyDefaultValue( minDeltaPressure );
    real64 maxDeltaPressure = MpiWrapper::max( regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::maxDeltaPressureString()) );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::maxDeltaPressureString()).setApplyDefaultValue( maxDeltaPressure );

    real64 minTemperature = MpiWrapper::min( regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::minTemperatureString()) );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::minTemperatureString()).setApplyDefaultValue( minTemperature );
    real64 maxTemperature = MpiWrapper::max( regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::maxTemperatureString()) );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::maxTemperatureString()).setApplyDefaultValue( maxTemperature );
    real64 averageTemperature = MpiWrapper::sum( regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::averageTemperatureString()) );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::averageTemperatureString()).setApplyDefaultValue( averageTemperature );

    real64 totalUncompactedPoreVolume = MpiWrapper::sum( regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::totalUncompactedPoreVolumeString()) );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::totalUncompactedPoreVolumeString()).setApplyDefaultValue( totalUncompactedPoreVolume );
    real64 totalPoreVolume = MpiWrapper::sum( regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::totalPoreVolumeString()) );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::totalPoreVolumeString()).setApplyDefaultValue( totalPoreVolume );
    real64 totalMass = MpiWrapper::sum( regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::totalMassString()) );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::totalMassString()).setApplyDefaultValue( totalMass );

    if( totalUncompactedPoreVolume > 0 )
    {
      float invTotalUncompactedPoreVolume = 1.0 / totalUncompactedPoreVolume;
      averagePressure *= invTotalUncompactedPoreVolume;
      averageTemperature *= invTotalUncompactedPoreVolume;
    }
    else
    {
      averagePressure = 0.0;
      averageTemperature = 0.0;
      GEOS_WARNING( GEOS_FMT( "{}, {}: Cannot compute average pressure & temperature because region pore volume is zero.", getName(), regionNames[i] ) );
    }

    string statPrefix = GEOS_FMT( "{}, {} (time {} s):", getName(), regionNames[i], time );
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Statistics, GEOS_FMT( "{} Pressure (min, average, max): {}, {}, {} Pa",
                                                               statPrefix, minPressure, averagePressure, maxPressure ) );
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Statistics, GEOS_FMT( "{} Delta pressure (min, max): {}, {} Pa",
                                                               statPrefix, minDeltaPressure, maxDeltaPressure ) );
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Statistics, GEOS_FMT( "{} Temperature (min, average, max): {}, {}, {} K",
                                                               statPrefix, minTemperature, averageTemperature, maxTemperature ) );
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Statistics, GEOS_FMT( "{} Total dynamic pore volume: {} rm^3",
                                                               statPrefix, totalPoreVolume ) );
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Statistics, GEOS_FMT( "{} Total fluid mass: {} kg",
                                                               statPrefix, totalMass ) );
    if( m_writeCSV > 0 && MpiWrapper::commRank() == 0 )
    {
      std::ofstream outputFile( m_outputDir + "/" + regionNames[i] + ".csv", std::ios_base::app );
      outputFile << time << "," << minPressure << "," << averagePressure << "," << maxPressure << "," <<
        minDeltaPressure << "," << maxDeltaPressure << "," <<
        minTemperature << "," << averageTemperature << "," << maxTemperature << "," <<
        totalPoreVolume << "," << totalMass << std::endl;
      outputFile.close();
    }
  }
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        SinglePhaseStatistics,
                        string const &, dataRepository::Group * const )

} /* namespace geos */
