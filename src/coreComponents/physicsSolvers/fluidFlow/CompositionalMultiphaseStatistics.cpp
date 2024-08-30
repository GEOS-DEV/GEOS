/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CompositionalMultiphaseStatistics.cpp
 */

#include "CompositionalMultiphaseStatistics.hpp"

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseHybridFVM.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseFVMKernels.hpp"


namespace geos
{

using namespace constitutive;
using namespace dataRepository;

CompositionalMultiphaseStatistics::CompositionalMultiphaseStatistics( const string & name,
                                                                      Group * const parent ):
  Base( name, parent ),
  m_computeCFLNumbers( 0 ),
  m_computeRegionStatistics( 1 )
{
  registerWrapper( viewKeyStruct::computeCFLNumbersString(), &m_computeCFLNumbers ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to decide whether CFL numbers are computed or not" );

  registerWrapper( viewKeyStruct::computeRegionStatisticsString(), &m_computeRegionStatistics ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to decide whether region statistics are computed or not" );

  registerWrapper( viewKeyStruct::relpermThresholdString(), &m_relpermThreshold ).
    setApplyDefaultValue( 1e-6 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to decide whether a phase is considered mobile (when the relperm is above the threshold) or immobile (when the relperm is below the threshold) in metric 2" );

}

void CompositionalMultiphaseStatistics::postInputInitialization()
{
  Base::postInputInitialization();

  if( dynamicCast< CompositionalMultiphaseHybridFVM * >( m_solver ) && m_computeCFLNumbers != 0 )
  {
    GEOS_THROW( GEOS_FMT( "{} {}: the option to compute CFL numbers is incompatible with CompositionalMultiphaseHybridFVM",
                          catalogName(), getDataContext() ),
                InputError );
  }
}

void CompositionalMultiphaseStatistics::registerDataOnMesh( Group & meshBodies )
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

    integer const numPhases = m_solver->numFluidPhases();
    integer const numComps = m_solver->numFluidComponents();

    // if we have to report region statistics, we have to register them first here
    if( m_computeRegionStatistics )
    {

      for( integer i = 0; i < regionNames.size(); ++i )
      {
        ElementRegionBase & region = elemManager.getRegion( regionNames[i] );

        region.registerGroup< RegionStatistics >( viewKeyStruct::regionStatisticsString() ).
          setRestartFlags( RestartFlags::NO_WRITE );

        RegionStatistics & regionStatistics = region.getGroup< RegionStatistics >( viewKeyStruct::regionStatisticsString() );
        regionStatistics.init( numPhases, numComps );
        // write output header
        if( m_writeCSV > 0 && MpiWrapper::commRank() == 0 )
        {
          std::ofstream outputFile( m_outputDir + "/" + regionNames[i] + ".csv" );
          string_view massUnit = units::getSymbol( m_solver->getMassUnit() );
          outputFile <<
            "Time [s],Min pressure [Pa],Average pressure [Pa],Max pressure [Pa],Min delta pressure [Pa],Max delta pressure [Pa]," <<
            "Min temperature [Pa],Average temperature [Pa],Max temperature [Pa],Total dynamic pore volume [rm^3]";
          for( integer ip = 0; ip < numPhases; ++ip )
            outputFile << ",Phase " << ip << " dynamic pore volume [rm^3]";
          for( integer ip = 0; ip < numPhases; ++ip )
            outputFile << ",Phase " << ip << " mass [" << massUnit << "]";
          for( integer ip = 0; ip < numPhases; ++ip )
            outputFile << ",Trapped phase " << ip << " mass (metric 1) [" << massUnit << "]";
          for( integer ip = 0; ip < numPhases; ++ip )
            outputFile << ",Non-trapped phase " << ip << " mass (metric 1) [" << massUnit << "]";
          for( integer ip = 0; ip < numPhases; ++ip )
            outputFile << ",Immobile phase " << ip << " mass (metric 2) [" << massUnit << "]";
          for( integer ip = 0; ip < numPhases; ++ip )
            outputFile << ",Mobile phase " << ip << " mass (metric 2) [" << massUnit << "]";
          for( integer ip = 0; ip < numPhases; ++ip )
          {
            for( integer ic = 0; ic < numComps; ++ic )
              outputFile << ",Component " << ic << " (phase " << ip << ") mass [" << massUnit << "]";
          }
          outputFile << std::endl;
          outputFile.close();
        }
      }
    }

    // if we have to compute CFL numbers later, we need to register additional variables
    if( m_computeCFLNumbers )
    {
      elemManager.forElementSubRegions( regionNames, [&]( localIndex const,
                                                          ElementSubRegionBase & subRegion )
      {
        subRegion.registerField< fields::flow::phaseOutflux >( getName() ).
          reference().resizeDimension< 1 >( numPhases );
        subRegion.registerField< fields::flow::componentOutflux >( getName() ).
          reference().resizeDimension< 1 >( numComps );
        subRegion.registerField< fields::flow::phaseCFLNumber >( getName() );
        subRegion.registerField< fields::flow::componentCFLNumber >( getName() );
      } );
    }
  } );
}

bool CompositionalMultiphaseStatistics::execute( real64 const time_n,
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
    if( m_computeRegionStatistics )
    {
      // current time is time_n + dt
      computeRegionStatistics( time_n + dt, mesh, regionNames );
    }
  } );

  if( m_computeCFLNumbers )
  {
    // current time is time_n + dt
    computeCFLNumbers( time_n + dt, dt, domain );
  }

  return false;
}

void CompositionalMultiphaseStatistics::computeRegionStatistics( real64 const time,
                                                                 MeshLevel & mesh,
                                                                 arrayView1d< string const > const & regionNames )
{
  GEOS_MARK_FUNCTION;
  integer const numPhases = m_solver->numFluidPhases();
  integer const numComps = m_solver->numFluidComponents();

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

    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::averageTemperatureString()).setApplyDefaultValue( 0.0 );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::maxTemperatureString()).setApplyDefaultValue( 0.0 );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::minTemperatureString()).setApplyDefaultValue( max );

    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::totalPoreVolumeString()).setApplyDefaultValue( 0.0 );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::totalUncompactedPoreVolumeString()).setApplyDefaultValue( 0.0 );

    regionStatistics.getReference< array1d< real64 > >( RegionStatistics::viewKeyStruct::phasePoreVolumeString()).setValues< serialPolicy >( 0.0 );

    regionStatistics.getReference< array1d< real64 > >( RegionStatistics::viewKeyStruct::phaseMassString()).setValues< serialPolicy >( 0.0 );
    regionStatistics.getReference< array1d< real64 > >( RegionStatistics::viewKeyStruct::trappedPhaseMassString()).setValues< serialPolicy >( 0.0 );
    regionStatistics.getReference< array1d< real64 > >( RegionStatistics::viewKeyStruct::immobilePhaseMassString()).setValues< serialPolicy >( 0.0 );
    regionStatistics.getReference< array2d< real64 > >( RegionStatistics::viewKeyStruct::dissolvedComponentMassString()).setValues< serialPolicy >( 0.0 );
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
    array2d< real64 > subRegiondissolvedComponentMass( numPhases, numComps );

    isothermalCompositionalMultiphaseBaseKernels::
      StatisticsKernel::
      launch< parallelDevicePolicy<> >( subRegion.size(),
                                        numComps,
                                        numPhases,
                                        m_relpermThreshold,
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
                                        subRegiondissolvedComponentMass.toView() );

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
    real64 & totalUncompactedPoreVolume = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::totalUncompactedPoreVolumeString());
    totalUncompactedPoreVolume += subRegionTotalUncompactedPoreVol;
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      regionStatistics.getReference< array1d< real64 > >( RegionStatistics::viewKeyStruct::phasePoreVolumeString())[ip] += subRegionPhaseDynamicPoreVol[ip];
      regionStatistics.getReference< array1d< real64 > >( RegionStatistics::viewKeyStruct::phaseMassString())[ip] += subRegionPhaseMass[ip];
      regionStatistics.getReference< array1d< real64 > >( RegionStatistics::viewKeyStruct::trappedPhaseMassString())[ip] += subRegionTrappedPhaseMass[ip];
      regionStatistics.getReference< array1d< real64 > >( RegionStatistics::viewKeyStruct::immobilePhaseMassString())[ip] += subRegionImmobilePhaseMass[ip];

      for( integer ic = 0; ic < numComps; ++ic )
      {
        regionStatistics.getReference< array2d< real64 > >( RegionStatistics::viewKeyStruct::dissolvedComponentMassString())[ip][ic] += subRegiondissolvedComponentMass[ip][ic];
      }
    }

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

    real64 minDeltaPressure = MpiWrapper::min( regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::minDeltaPressureString()) );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::minDeltaPressureString()).setApplyDefaultValue( minDeltaPressure );

    real64 maxDeltaPressure = MpiWrapper::max( regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::maxDeltaPressureString()) );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::maxDeltaPressureString()).setApplyDefaultValue( maxDeltaPressure );

    real64 minTemperature = MpiWrapper::min( regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::minTemperatureString()) );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::minTemperatureString()).setApplyDefaultValue( minTemperature );

    real64 maxTemperature = MpiWrapper::max( regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::maxTemperatureString()) );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::maxTemperatureString()).setApplyDefaultValue( maxTemperature );

    real64 totalUncompactedPoreVolume = MpiWrapper::sum( regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::totalUncompactedPoreVolumeString()) );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::totalUncompactedPoreVolumeString()).setApplyDefaultValue( totalUncompactedPoreVolume );


    array1d< real64 > phasePoreVolume = regionStatistics.getReference< array1d< real64 > >( RegionStatistics::viewKeyStruct::phasePoreVolumeString());
    array1d< real64 > phaseMass = regionStatistics.getReference< array1d< real64 > >( RegionStatistics::viewKeyStruct::phaseMassString());
    array1d< real64 > trappedPhaseMass = regionStatistics.getReference< array1d< real64 > >( RegionStatistics::viewKeyStruct::trappedPhaseMassString());
    array1d< real64 > immobilePhaseMass = regionStatistics.getReference< array1d< real64 > >( RegionStatistics::viewKeyStruct::immobilePhaseMassString());
    array2d< real64 > dissolvedComponentMass = regionStatistics.getReference< array2d< real64 > >( RegionStatistics::viewKeyStruct::dissolvedComponentMassString());

    real64 totalPoreVolume = 0.0;
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      phasePoreVolume[ip] = MpiWrapper::sum( phasePoreVolume[ip] );
      phaseMass[ip] = MpiWrapper::sum( phaseMass[ip] );
      trappedPhaseMass[ip] = MpiWrapper::sum( trappedPhaseMass[ip] );
      immobilePhaseMass[ip] = MpiWrapper::sum( immobilePhaseMass[ip] );
      totalPoreVolume += phasePoreVolume[ip];
      for( integer ic = 0; ic < numComps; ++ic )
      {
        dissolvedComponentMass[ip][ic] = MpiWrapper::sum( dissolvedComponentMass[ip][ic] );
      }
    }
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::totalPoreVolumeString()).setApplyDefaultValue( totalPoreVolume );

    real64 & averagePressure = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::averagePressureString() );
    averagePressure = MpiWrapper::sum( averagePressure );

    real64 & averageTemperature = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::averageTemperatureString() );
    averageTemperature = MpiWrapper::sum( averageTemperature );

    if( totalUncompactedPoreVolume > 0 )
    {
      float invTotalUncompactedPoreVolume = 1.0 / totalUncompactedPoreVolume;
      averagePressure *= invTotalUncompactedPoreVolume;
      averageTemperature *= invTotalUncompactedPoreVolume;
    }
    else
    {
      averagePressure =  0.0;
      averageTemperature = 0.0;
      GEOS_LOG_LEVEL_RANK_0( 1, getName() << ", " << regionNames[i]
                                          << ": Cannot compute average pressure because region pore volume is zero." );
    }


    // helpers to report statistics
    array1d< real64 > nonTrappedPhaseMass( numPhases );
    array1d< real64 > mobilePhaseMass( numPhases );
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      nonTrappedPhaseMass[ip] = phaseMass[ip] - trappedPhaseMass[ip];
      mobilePhaseMass[ip] = phaseMass[ip] - immobilePhaseMass[ip];
    }

    string_view massUnit = units::getSymbol( m_solver->getMassUnit() );

    GEOS_LOG_LEVEL_RANK_0( 1, getName() << ", " << regionNames[i]
                                        << ": Pressure (min, average, max): "
                                        << minPressure << ", " << averagePressure << ", " << maxPressure << " Pa" );
    GEOS_LOG_LEVEL_RANK_0( 1, getName() << ", " << regionNames[i]
                                        << ": Delta pressure (min, max): "
                                        << minDeltaPressure << ", " << maxDeltaPressure << " Pa" );
    GEOS_LOG_LEVEL_RANK_0( 1, getName() << ", " << regionNames[i]
                                        << ": Temperature (min, average, max): "
                                        << minTemperature << ", " << averageTemperature << ", " << maxTemperature << " K" );
    GEOS_LOG_LEVEL_RANK_0( 1, getName() << ", " << regionNames[i]
                                        << ": Total dynamic pore volume: " << totalPoreVolume << " rm^3" );
    GEOS_LOG_LEVEL_RANK_0( 1, getName() << ", " << regionNames[i]
                                        << ": Phase dynamic pore volumes: " << phasePoreVolume << " rm^3" );
    GEOS_LOG_LEVEL_RANK_0( 1, getName() << ", " << regionNames[i]
                                        << ": Phase mass: " << phaseMass << " " << massUnit );

    // metric 1: trapping computed with the Land trapping coefficient (similar to Eclipse)
    GEOS_LOG_LEVEL_RANK_0( 1, getName() << ", " << regionNames[i]
                                        << ": Trapped phase mass (metric 1): " << trappedPhaseMass << " " << massUnit
                           );
    GEOS_LOG_LEVEL_RANK_0( 1, getName() << ", " << regionNames[i]
                                        << ": Non-trapped phase mass (metric 1): " << nonTrappedPhaseMass << " " << massUnit );

    // metric 2: immobile phase mass computed with a threshold on relative permeability
    GEOS_LOG_LEVEL_RANK_0( 1, getName() << ", " << regionNames[i]
                                        << ": Immobile phase mass (metric 2): " << immobilePhaseMass << " " << massUnit
                           );
    GEOS_LOG_LEVEL_RANK_0( 1, getName() << ", " << regionNames[i]
                                        << ": Mobile phase mass (metric 2): " << mobilePhaseMass << " " << massUnit );


    GEOS_LOG_LEVEL_RANK_0( 1, getName() << ", " << regionNames[i]
                                        << ": Dissolved component mass: " << dissolvedComponentMass << " " << massUnit
                           );

    if( m_writeCSV > 0 && MpiWrapper::commRank() == 0 )
    {
      std::ofstream outputFile( m_outputDir + "/" + regionNames[i] + ".csv", std::ios_base::app );
      outputFile << time << "," << minPressure << "," << averagePressure << "," << maxPressure << "," <<
        minDeltaPressure << "," << maxDeltaPressure << "," << minTemperature << "," <<
        averageTemperature << "," << maxTemperature << "," << totalPoreVolume;
      for( integer ip = 0; ip < numPhases; ++ip )
        outputFile << "," <<phasePoreVolume[ip];
      for( integer ip = 0; ip < numPhases; ++ip )
        outputFile << "," << phaseMass[ip];
      for( integer ip = 0; ip < numPhases; ++ip )
        outputFile << "," << trappedPhaseMass[ip];
      for( integer ip = 0; ip < numPhases; ++ip )
        outputFile << "," << nonTrappedPhaseMass[ip];
      for( integer ip = 0; ip < numPhases; ++ip )
        outputFile << "," << immobilePhaseMass[ip];
      for( integer ip = 0; ip < numPhases; ++ip )
        outputFile << "," << mobilePhaseMass[ip];
      for( integer ip = 0; ip < numPhases; ++ip )
      {
        for( integer ic = 0; ic < numComps; ++ic )
          outputFile << "," << dissolvedComponentMass[ip][ic];
      }
      outputFile << std::endl;
      outputFile.close();
    }
  }
}

void CompositionalMultiphaseStatistics::computeCFLNumbers( real64 const time,
                                                           real64 const dt,
                                                           DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  real64 maxPhaseCFL, maxCompCFL;
  m_solver->computeCFLNumbers( domain, dt, maxPhaseCFL, maxCompCFL );

  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{} (time {} s): Max phase CFL number: {}", getName(), time, maxPhaseCFL ) );
  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{} (time {} s): Max component CFL number: {}", getName(), time, maxCompCFL ) );
}


REGISTER_CATALOG_ENTRY( TaskBase,
                        CompositionalMultiphaseStatistics,
                        string const &, dataRepository::Group * const )

} /* namespace geos */