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

void CompositionalMultiphaseStatistics::postProcessInput()
{
  Base::postProcessInput();

  if( dynamicCast< CompositionalMultiphaseHybridFVM * >( m_solver ) && m_computeCFLNumbers != 0 )
  {
    GEOS_THROW( GEOS_FMT( "{} {}: the option to compute CFL numbers is incompatible with CompositionalMultiphaseHybridFVM",
                          catalogName(), getName() ),
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

bool CompositionalMultiphaseStatistics::execute( real64 const GEOS_UNUSED_PARAM( time_n ),
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
      computeRegionStatistics( mesh, regionNames );
    }
  } );

  if( m_computeCFLNumbers )
  {
    computeCFLNumbers( dt, domain );
  }

  return false;
}

void CompositionalMultiphaseStatistics::computeRegionStatistics( MeshLevel & mesh,
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

    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::averagePressureString()).setApplyDefaultValue(0.0);
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::maxPressureString()).setApplyDefaultValue(0.0);
    constexpr real64 max = LvArray::NumericLimits< real64 >::max;
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::minPressureString()).setApplyDefaultValue(max);

    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::maxDeltaPressureString()).setApplyDefaultValue(-max);
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::minDeltaPressureString()).setApplyDefaultValue(max);
    
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::averageTemperatureString()).setApplyDefaultValue(0.0);
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::maxTemperatureString()).setApplyDefaultValue(0.0);
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::minTemperatureString()).setApplyDefaultValue(max);

    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::totalPoreVolumeString()).setApplyDefaultValue(0.0);
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::totalUncompactedPoreVolumeString()).setApplyDefaultValue(0.0);
    
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
    arrayView1d< real64 > const deltaPres = subRegion.getField< fields::flow::deltaPressure >();
    arrayView1d< real64 const > const temp = subRegion.getField< fields::flow::temperature >();
    arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
      subRegion.getField< fields::flow::phaseVolumeFraction >();

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
    array2d< real64 > subRegionDissolvedComponentMass( numPhases, numComps );

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
                                        subRegionDissolvedComponentMass.toView() );

    ElementRegionBase & region = elemManager.getRegion( subRegion.getParent().getParent().getName() );
    RegionStatistics & regionStatistics = region.getGroup< RegionStatistics >( viewKeyStruct::regionStatisticsString() );

    real64& averagePressure = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::averagePressureString());    
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::averagePressureString()).setApplyDefaultValue(averagePressure + subRegionAvgPresNumerator);

    real64& minPressure = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::minPressureString());    
    if( subRegionMinPres < minPressure )
    {      
      regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::minPressureString()).setApplyDefaultValue(subRegionMinPres);
      real64& test = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::minPressureString());
      std::cout << "hahahahaah " << subRegionMinPres << " " << test << std::endl;
    }
    real64& maxPressure = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::maxPressureString());
    if( subRegionMaxPres > maxPressure )
    {
      regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::maxPressureString()).setApplyDefaultValue(subRegionMaxPres);
    }

    real64& minDeltaPressure = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::minDeltaPressureString());
    if( subRegionMinDeltaPres < minDeltaPressure )
    {
      regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::minDeltaPressureString()).setApplyDefaultValue(subRegionMinDeltaPres);
    }
    real64& maxDeltaPressure = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::maxDeltaPressureString());
    if( subRegionMaxDeltaPres > maxDeltaPressure )
    {
      regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::maxDeltaPressureString()).setApplyDefaultValue(subRegionMaxDeltaPres);
    }

    real64& averageTemperature = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::averageTemperatureString());
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::averageTemperatureString()).setApplyDefaultValue(averageTemperature + subRegionAvgTempNumerator);
    
    real64& minTemperature = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::minTemperatureString());
    if( subRegionMinTemp < minTemperature )
    {
      regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::minTemperatureString()).setApplyDefaultValue(subRegionMinTemp);
    }
    real64& maxTemperature = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::maxTemperatureString());
    if( subRegionMaxTemp > maxTemperature )
    {
      regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::maxTemperatureString()).setApplyDefaultValue(subRegionMaxTemp);
    }

    real64& totalUncompactedPoreVolume = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::totalUncompactedPoreVolumeString());
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::totalUncompactedPoreVolumeString()).setApplyDefaultValue(totalUncompactedPoreVolume + subRegionTotalUncompactedPoreVol);
     
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      regionStatistics.getReference< array1d< real64 > >( RegionStatistics::viewKeyStruct::phasePoreVolumeString())[ip] += subRegionPhaseDynamicPoreVol[ip];
      regionStatistics.getReference< array1d< real64 > >( RegionStatistics::viewKeyStruct::phaseMassString())[ip] += subRegionPhaseMass[ip];
      regionStatistics.getReference< array1d< real64 > >( RegionStatistics::viewKeyStruct::trappedPhaseMassString())[ip] += subRegionTrappedPhaseMass[ip];
      regionStatistics.getReference< array1d< real64 > >( RegionStatistics::viewKeyStruct::immobilePhaseMassString())[ip] += subRegionImmobilePhaseMass[ip];

      for( integer ic = 0; ic < numComps; ++ic )
      {
        regionStatistics.getReference< array2d< real64 > >( RegionStatistics::viewKeyStruct::dissolvedComponentMassString())[ip][ic] += subRegionDissolvedComponentMass[ip][ic];
      }
    }

  } );

  // Step 3: synchronize the results over the MPI ranks
  for( integer i = 0; i < regionNames.size(); ++i )
  {
    ElementRegionBase & region = elemManager.getRegion( regionNames[i] );
    RegionStatistics & regionStatistics = region.getGroup< RegionStatistics >( viewKeyStruct::regionStatisticsString() );
    
    real64 minPressure = MpiWrapper::min( regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::minPressureString()) );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::minPressureString()).setApplyDefaultValue(minPressure);
    
    real64 maxPressure = MpiWrapper::max( regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::maxPressureString()) );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::maxPressureString()).setApplyDefaultValue(maxPressure);

    real64 minDeltaPressure = MpiWrapper::min( regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::minDeltaPressureString()) );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::minDeltaPressureString()).setApplyDefaultValue(minDeltaPressure);

    real64 maxDeltaPressure = MpiWrapper::max( regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::maxDeltaPressureString()) );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::maxDeltaPressureString()).setApplyDefaultValue(maxDeltaPressure);
    
    real64 minTemperature = MpiWrapper::min( regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::minTemperatureString()) );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::minTemperatureString()).setApplyDefaultValue(minTemperature);

    real64 maxTemperature = MpiWrapper::max( regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::maxTemperatureString()) );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::maxTemperatureString()).setApplyDefaultValue(maxTemperature);

    real64 totalUncompactedPoreVolume = MpiWrapper::sum( regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::totalUncompactedPoreVolumeString()) );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::totalUncompactedPoreVolumeString()).setApplyDefaultValue(totalUncompactedPoreVolume);

    
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
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::totalPoreVolumeString()).setApplyDefaultValue(totalPoreVolume);

    real64 averagePressure = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::averagePressureString() );    
    averagePressure = MpiWrapper::sum( averagePressure );
    averagePressure /= regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::totalUncompactedPoreVolumeString() );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::averagePressureString() ).setApplyDefaultValue(averagePressure);

    real64 averageTemperature = regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::averageTemperatureString() );    
    averageTemperature = MpiWrapper::sum( averageTemperature );
    averageTemperature /= regionStatistics.getReference< real64 >( RegionStatistics::viewKeyStruct::totalUncompactedPoreVolumeString() );
    regionStatistics.getWrapper< real64 >( RegionStatistics::viewKeyStruct::averageTemperatureString() ).setApplyDefaultValue(averageTemperature);

    //helpers to report statistics
    array1d< real64 > nonTrappedPhaseMass( numPhases );
    array1d< real64 > mobilePhaseMass( numPhases );
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      nonTrappedPhaseMass[ip] = phaseMass[ip] - trappedPhaseMass[ip];
      mobilePhaseMass[ip] = phaseMass[ip] - immobilePhaseMass[ip];
    }

    integer const useMass = m_solver->getReference< integer >( CompositionalMultiphaseBase::viewKeyStruct::useMassFlagString() );
    string const massUnit = useMass ? "kg" : "mol";

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
  }
}

void CompositionalMultiphaseStatistics::computeCFLNumbers( real64 const & dt,
                                                           DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  integer const numPhases = m_solver->numFluidPhases();
  integer const numComps = m_solver->numFluidComponents();

  // Step 1: reset the arrays involved in the computation of CFL numbers
  m_solver->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                         MeshLevel & mesh,
                                                                         arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView2d< real64, compflow::USD_PHASE > const & phaseOutflux =
        subRegion.getField< fields::flow::phaseOutflux >();
      arrayView2d< real64, compflow::USD_COMP > const & compOutflux =
        subRegion.getField< fields::flow::componentOutflux >();
      phaseOutflux.zero();
      compOutflux.zero();
    } );

    // Step 2: compute the total volumetric outflux for each reservoir cell by looping over faces
    NumericalMethodsManager & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase & fluxApprox = fvManager.getFluxApproximation( m_solver->getDiscretizationName() );

    isothermalCompositionalMultiphaseFVMKernels::
      CFLFluxKernel::CompFlowAccessors compFlowAccessors( mesh.getElemManager(), getName() );
    isothermalCompositionalMultiphaseFVMKernels::
      CFLFluxKernel::MultiFluidAccessors multiFluidAccessors( mesh.getElemManager(), getName() );
    isothermalCompositionalMultiphaseFVMKernels::
      CFLFluxKernel::PermeabilityAccessors permeabilityAccessors( mesh.getElemManager(), getName() );
    isothermalCompositionalMultiphaseFVMKernels::
      CFLFluxKernel::RelPermAccessors relPermAccessors( mesh.getElemManager(), getName() );

    // TODO: find a way to compile with this modifiable accessors in CompFlowAccessors, and remove them from here
    ElementRegionManager::ElementViewAccessor< arrayView2d< real64, compflow::USD_PHASE > > const phaseOutfluxAccessor =
      mesh.getElemManager().constructViewAccessor< array2d< real64, compflow::LAYOUT_PHASE >,
                                                   arrayView2d< real64, compflow::USD_PHASE > >( fields::flow::phaseOutflux::key() );

    ElementRegionManager::ElementViewAccessor< arrayView2d< real64, compflow::USD_COMP > > const compOutfluxAccessor =
      mesh.getElemManager().constructViewAccessor< array2d< real64, compflow::LAYOUT_COMP >,
                                                   arrayView2d< real64, compflow::USD_COMP > >( fields::flow::componentOutflux::key() );


    fluxApprox.forAllStencils( mesh, [&] ( auto & stencil )
    {

      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

      // While this kernel is waiting for a factory class, pass all the accessors here
      isothermalCompositionalMultiphaseBaseKernels::KernelLaunchSelector1
      < isothermalCompositionalMultiphaseFVMKernels::CFLFluxKernel >( numComps,
                                                                      numPhases,
                                                                      dt,
                                                                      stencilWrapper,
                                                                      compFlowAccessors.get( fields::flow::pressure{} ),
                                                                      compFlowAccessors.get( fields::flow::gravityCoefficient{} ),
                                                                      compFlowAccessors.get( fields::flow::phaseVolumeFraction{} ),
                                                                      permeabilityAccessors.get( fields::permeability::permeability{} ),
                                                                      permeabilityAccessors.get( fields::permeability::dPerm_dPressure{} ),
                                                                      relPermAccessors.get( fields::relperm::phaseRelPerm{} ),
                                                                      multiFluidAccessors.get( fields::multifluid::phaseViscosity{} ),
                                                                      multiFluidAccessors.get( fields::multifluid::phaseDensity{} ),
                                                                      multiFluidAccessors.get( fields::multifluid::phaseMassDensity{} ),
                                                                      multiFluidAccessors.get( fields::multifluid::phaseCompFraction{} ),
                                                                      phaseOutfluxAccessor.toNestedView(),
                                                                      compOutfluxAccessor.toNestedView() );
    } );
  } );

  // Step 3: finalize the (cell-based) computation of the CFL numbers
  real64 localMaxPhaseCFLNumber = 0.0;
  real64 localMaxCompCFLNumber = 0.0;

  m_solver->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                         MeshLevel & mesh,
                                                                         arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView2d< real64 const, compflow::USD_PHASE > const & phaseOutflux =
        subRegion.getField< fields::flow::phaseOutflux >();
      arrayView2d< real64 const, compflow::USD_COMP > const & compOutflux =
        subRegion.getField< fields::flow::componentOutflux >();

      arrayView1d< real64 > const & phaseCFLNumber = subRegion.getField< fields::flow::phaseCFLNumber >();
      arrayView1d< real64 > const & compCFLNumber = subRegion.getField< fields::flow::componentCFLNumber >();

      arrayView1d< real64 const > const & volume = subRegion.getElementVolume();

      arrayView2d< real64 const, compflow::USD_COMP > const & compDens =
        subRegion.getField< fields::flow::globalCompDensity >();
      arrayView2d< real64 const, compflow::USD_COMP > const compFrac =
        subRegion.getField< fields::flow::globalCompFraction >();
      arrayView2d< real64, compflow::USD_PHASE > const phaseVolFrac =
        subRegion.getField< fields::flow::phaseVolumeFraction >();

      Group const & constitutiveModels = subRegion.getGroup( ElementSubRegionBase::groupKeyStruct::constitutiveModelsString() );

      string const & fluidName = subRegion.getReference< string >( CompositionalMultiphaseBase::viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = constitutiveModels.getGroup< MultiFluidBase >( fluidName );
      arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseVisc = fluid.phaseViscosity();

      string const & relpermName = subRegion.getReference< string >( CompositionalMultiphaseBase::viewKeyStruct::relPermNamesString() );
      RelativePermeabilityBase const & relperm = constitutiveModels.getGroup< RelativePermeabilityBase >( relpermName );
      arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseRelPerm = relperm.phaseRelPerm();
      arrayView4d< real64 const, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac = relperm.dPhaseRelPerm_dPhaseVolFraction();

      string const & solidName = subRegion.getReference< string >( CompositionalMultiphaseBase::viewKeyStruct::solidNamesString() );
      CoupledSolidBase const & solid = constitutiveModels.getGroup< CoupledSolidBase >( solidName );
      arrayView2d< real64 const > const & porosity = solid.getPorosity();

      real64 subRegionMaxPhaseCFLNumber = 0.0;
      real64 subRegionMaxCompCFLNumber = 0.0;

      isothermalCompositionalMultiphaseBaseKernels::KernelLaunchSelector2
      < isothermalCompositionalMultiphaseFVMKernels::CFLKernel >( numComps, numPhases,
                                                                  subRegion.size(),
                                                                  volume,
                                                                  porosity,
                                                                  compDens,
                                                                  compFrac,
                                                                  phaseVolFrac,
                                                                  phaseRelPerm,
                                                                  dPhaseRelPerm_dPhaseVolFrac,
                                                                  phaseVisc,
                                                                  phaseOutflux,
                                                                  compOutflux,
                                                                  phaseCFLNumber,
                                                                  compCFLNumber,
                                                                  subRegionMaxPhaseCFLNumber,
                                                                  subRegionMaxCompCFLNumber );

      localMaxPhaseCFLNumber = LvArray::math::max( localMaxPhaseCFLNumber, subRegionMaxPhaseCFLNumber );
      localMaxCompCFLNumber = LvArray::math::max( localMaxCompCFLNumber, subRegionMaxCompCFLNumber );

    } );
  } );

  real64 const globalMaxPhaseCFLNumber = MpiWrapper::max( localMaxPhaseCFLNumber );
  real64 const globalMaxCompCFLNumber = MpiWrapper::max( localMaxCompCFLNumber );

  GEOS_LOG_LEVEL_RANK_0( 1, getName() << ": Max phase CFL number: " << globalMaxPhaseCFLNumber );
  GEOS_LOG_LEVEL_RANK_0( 1, getName() << ": Max component CFL number: " << globalMaxCompCFLNumber );
}


REGISTER_CATALOG_ENTRY( TaskBase,
                        CompositionalMultiphaseStatistics,
                        string const &, dataRepository::Group * const )

} /* namespace geos */
