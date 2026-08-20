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
 * @file SinglePhaseWell.cpp
 */

#include "SinglePhaseWell.hpp"

#include "common/DataTypes.hpp"
#include "common/FieldSpecificationOps.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidSelector.hpp"
#include "dataRepository/Group.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "mesh/PerforationFields.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/wells/LogLevelsInfo.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/SolutionCheckHelpers.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWellFields.hpp"

#include "physicsSolvers/fluidFlow/wells/WellInjectionConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/WellProductionConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/WellBHPConstraints.hpp"
#include "physicsSolvers/fluidFlow/wells/WellVolumeRateConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/WellPhaseVolumeRateConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/WellMassRateConstraint.hpp"

#include "physicsSolvers/fluidFlow/wells/kernels/SinglePhaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/kernels/ThermalSinglePhaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/kernels/SinglePhasePerforationFluxKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/FluidUpdateKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/SolutionCheckKernel.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseStatisticsAggregator.hpp"
#include "physicsSolvers/fluidFlow/wells/kernels/SinglePhaseWellConstraintKernels.hpp"
#include "physicsSolvers/multiphysics/CoupledReservoirAndSinglePhaseWellKernels.hpp"
namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace fields;
using namespace singlePhaseWellKernels;
using namespace singlePhaseStatistics;

SinglePhaseBase & getFlowSolver( SinglePhaseWell & wellSolver )
{
  // TODO: change the way we access the flowSolver here
  return wellSolver.getParent().getGroup< SinglePhaseBase >( wellSolver.getFlowSolverName() );
}

SinglePhaseBase const & getFlowSolver( SinglePhaseWell const & wellSolver )
{
  // TODO: change the way we access the flowSolver here
  return wellSolver.getParent().getGroup< SinglePhaseBase >( wellSolver.getFlowSolverName() );
}

real64 getBHPReferenceGravityCoef( SinglePhaseWell const & wellSolver,
                                   ConstraintSourceId const source )
{
  real64 refGravCoef = 0.0;
  bool foundConstraint = false;

  if( wellSolver.isProducer() )
  {
    wellSolver.forSubGroups< MinimumBHPConstraint >( [&]( WellConstraintBase const & constraint )
    {
      if( constraint.isConstraintActive() && constraint.getConstraintSource() == source )
      {
        refGravCoef = static_cast< MinimumBHPConstraint const & >( constraint ).getReferenceGravityCoef();
        foundConstraint = true;
      }
    } );
  }
  else
  {
    wellSolver.forSubGroups< MaximumBHPConstraint >( [&]( WellConstraintBase const & constraint )
    {
      if( constraint.isConstraintActive() && constraint.getConstraintSource() == source )
      {
        refGravCoef = static_cast< MaximumBHPConstraint const & >( constraint ).getReferenceGravityCoef();
        foundConstraint = true;
      }
    } );
  }

  GEOS_THROW_IF( !foundConstraint,
                 GEOS_FMT( "Could not find active BHP constraint for well {}", wellSolver.getName() ),
                 InputError );

  return refGravCoef;
}

SinglePhaseWell::SinglePhaseWell( const string & name,
                                  Group * const parent ):
  WellControls( name, parent )
{
  m_numPhases = 1;
  m_numComponents = 1;

  this->registerWrapper( FlowSolverBase::viewKeyStruct::allowNegativePressureString(), &m_allowNegativePressure ).
    setApplyDefaultValue( 1 ). // negative pressure is allowed by default
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag indicating if negative pressure is allowed" );
}

SinglePhaseWell::~SinglePhaseWell() = default;

void SinglePhaseWell::registerWellDataOnMesh( WellElementSubRegion & subRegion )
{
  updateNumDofPerElement();
  WellControls::registerDataOnMesh( subRegion );
  setConstitutiveNames ( subRegion );

  if( m_referenceFluidModelName.empty() )
  {
    m_referenceFluidModelName = getConstitutiveName< SingleFluidBase >( subRegion );
  }
  subRegion.registerField< well::pressure >( getName() );
  subRegion.registerField< well::pressure_n >( getName() );

  subRegion.registerField< well::temperature >( getName() );
  if( isThermal() )
  {
    subRegion.registerField< well::temperature_n >( getName() );
  }
  subRegion.registerField< well::connectionRate_n >( getName() );
  subRegion.registerField< well::connectionRate >( getName() );
  subRegion.registerField< well::gravityCoefficient >( getName() );
  PerforationData & perforationData = *subRegion.getPerforationData();
  perforationData.registerField< well::gravityCoefficient >( getName() );

  perforationData.registerField< well::perforationRate >( getName() );
  perforationData.registerField< well::dPerforationRate >( getName() ).
    reference().resizeDimension< 1, 2 >( 2, 2 );
  if( isThermal() )
  {
    perforationData.registerField< well::energyPerforationFlux >( getName() );
    perforationData.registerField< well::dEnergyPerforationFlux >( getName() ).
      reference().resizeDimension< 1, 2 >( 2, 2 );
    perforationData.registerField< well::gravityCoefficient >( getName() );
  }

  registerWrapper< real64 >( viewKeyStruct::currentBHPString() );
  registerWrapper< real64 >( viewKeyStruct::currentVolRateString() );

  // write rates output header
  if( m_writeCSV > 0 && subRegion.isLocallyOwned())
  {
    string const fileName = GEOS_FMT( "{}/{}.csv", m_ratesOutputDir, getName() );
    string const conditionKey = useSurfaceConditions() ? "surface" : "reservoir";
    string const unitKey = useSurfaceConditions() ? "s" : "r";
    // format: time,bhp,total_rate,total_vol_rate
    makeDirsForPath( m_ratesOutputDir );
    GEOS_LOG( GEOS_FMT( "{}: Rates CSV generated at {}", getName(), fileName ) );
    std::ofstream outputFile( fileName );
    outputFile << "Time [s],BHP [Pa],Total rate [kg/s],Total " << conditionKey << " volumetric rate ["<<unitKey<<"m3/s]" << std::endl;
    outputFile.close();
  }

}

void SinglePhaseWell::setConstitutiveNames( ElementSubRegionBase & subRegion ) const
{
  setConstitutiveName< SingleFluidBase >( subRegion, viewKeyStruct::fluidNamesString(), "singlephase fluid" );
}

string SinglePhaseWell::resElementDofName() const
{
  return SinglePhaseBase::viewKeyStruct::elemDofFieldString();
}

void SinglePhaseWell::validateWellConstraints( real64 const & time_n,
                                               real64 const & GEOS_UNUSED_PARAM( dt ),
                                               WellElementSubRegion const & subRegion )
{
  GEOS_UNUSED_VAR( time_n );
  GEOS_UNUSED_VAR( subRegion );

  if( useSurfaceConditions() )
  {
    bool const useSegmentValues = referenceReservoirRegion().empty();

    static bool firstNoRefRegionMsg = true;
    if( useSegmentValues && firstNoRefRegionMsg )
    {
      GEOS_WARNING( WellControls::viewKeyStruct::referenceReservoirRegionString() <<
                    " not set: well constraint fluid property calculations will use top segement pressure and temp ",
                    getDataContext() );
      firstNoRefRegionMsg = false;
    }
  }
}

void SinglePhaseWell::initializeWellPostInitialConditionsPreSubGroups( WellElementSubRegion & subRegion )
{

  // set gravity coefficient
  setGravCoef( subRegion, getParent().getParent().getReference< R1Tensor >( PhysicsSolverManager::viewKeyStruct::gravityVectorString() ));

  // setup fluid model
  createSeparator( subRegion );
}
void SinglePhaseWell::initializePostInitialConditionsPreSubGroups()
{
  WellControls::initializePostInitialConditionsPreSubGroups();
}
void SinglePhaseWell::postRestartInitialization( )
{
  // setup fluid separator
  constitutive::SingleFluidBase & fluidSeparator =   getSingleFluidSeparator();
  fluidSeparator.allocateConstitutiveData( *this, 1 );
  fluidSeparator.resize( 1 );

}
void SinglePhaseWell::createSeparator( WellElementSubRegion & subRegion )
{

  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  SingleFluidBase & fluid = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );
  // setup fluid separator

  string const fluidSeparatorName =  getName() + "Separator";
  std::unique_ptr< constitutive::ConstitutiveBase >  fluidSeparatorPtr  = fluid.deliverClone( fluidSeparatorName, this );
  fluidSeparatorPtr->allocateConstitutiveData( *this, 1 );
  fluidSeparatorPtr->resize( 1 );
  setFluidSeparator( std::move( fluidSeparatorPtr ));

}

void SinglePhaseWell::updateBHPForConstraint( WellElementSubRegion & subRegion )
{
  GEOS_MARK_FUNCTION;

  // the rank that owns the reference well element is responsible for the calculations below.
  if( !subRegion.isLocallyOwned() )
  {
    return;
  }

  localIndex const iwelemRef = subRegion.getTopWellElementIndex();

  // subRegion data

  arrayView1d< real64 const > const pres =
    subRegion.getField< well::pressure >();

  arrayView1d< real64 const > const wellElemGravCoef =
    subRegion.getField< well::gravityCoefficient >();

  // fluid data
  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  SingleFluidBase & fluid = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & dens = fluid.density();
  arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const & dDens = fluid.dDensity();

  real64 const refGravCoef = getBHPReferenceGravityCoef( *this, ConstraintSourceId::USER );
  real64 & currentBHP = getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentBHPString() );

  // bring everything back to host, capture the scalars by reference
  forAll< serialPolicy >( 1, [pres,
                              dens,
                              dDens,
                              wellElemGravCoef,
                              &currentBHP,
                              &iwelemRef,
                              &refGravCoef] ( localIndex const )
  {
    real64 const diffGravCoef = refGravCoef - wellElemGravCoef[iwelemRef];
    currentBHP = pres[iwelemRef] + dens[iwelemRef][0] * diffGravCoef;
  } );

  GEOS_LOG_LEVEL_BY_RANK( logInfo::WellControl,
                          GEOS_FMT( "{}: The BHP (at the specified reference elevation) = {} Pa",
                                    getName(), currentBHP ) );

}

void SinglePhaseWell::calculateReferenceElementRates( WellElementSubRegion & subRegion )
{
  GEOS_MARK_FUNCTION;

  // the rank that owns the reference well element is responsible for the calculations below.
  if( !subRegion.isLocallyOwned() )
  {
    return;
  }

  localIndex const iwelemRef = subRegion.getTopWellElementIndex();

  // subRegion data

  arrayView1d< real64 const > const & connRate =
    subRegion.getField< well::connectionRate >();

  // fluid data
  constitutive::SingleFluidBase & fluidSeparator =   getSingleFluidSeparator();
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & dens = fluidSeparator.density();

  real64 & currentVolRate =
    getReference< real64 >( WellControls::viewKeyStruct::currentVolRateString() );

  // bring everything back to host, capture the scalars by reference
  forAll< serialPolicy >( 1, [connRate,
                              dens,
                              &currentVolRate,
                              &iwelemRef] ( localIndex const )
  {
    real64 const densInv = 1.0 / dens[iwelemRef][0];
    currentVolRate = connRate[iwelemRef] * densInv;
    // compute mass rate
  } );


}

void SinglePhaseWell::updateFluidModel( WellElementSubRegion & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  arrayView1d< real64 const > const pres = subRegion.getField< well::pressure >();
  arrayView1d< real64 const > const temp = subRegion.getField< well::temperature >();

  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  SingleFluidBase & fluid = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );

  constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    singlePhaseBaseKernels::FluidUpdateKernel::launch( fluidWrapper, pres, temp );
  } );
}
void SinglePhaseWell::updateSeparator( ElementRegionManager const & elemManager, WellElementSubRegion & subRegion )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( elemManager );

  // the rank that owns the reference well element is responsible for the calculations below.
  if( !subRegion.isLocallyOwned() )
  {
    return;
  }

  localIndex const iwelemRef = subRegion.getTopWellElementIndex();

  constitutive::SingleFluidBase & fluidSeparator = getSingleFluidSeparator();
  string const wellControlsName =  getName();
  bool const logSurfaceCondition = isLogLevelActive< logInfo::WellControl >( getLogLevel());
  integer const useSurfCond = useSurfaceConditions();
  ReferenceConditions const refConditions = getReferenceConditions( subRegion );

  constitutiveUpdatePassThru( fluidSeparator, [&]( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidSeparatorWrapper = castedFluid.createKernelWrapper();

    geos::internal::kernelLaunchSelectorThermalSwitch( isThermal(), [&] ( auto ISTHERMAL )
    {
      integer constexpr IS_THERMAL = ISTHERMAL();
      GEOS_UNUSED_VAR( IS_THERMAL );
      // bring everything back to host, capture the scalars by reference
      forAll< serialPolicy >( 1, [fluidSeparatorWrapper,
                                  logSurfaceCondition,
                                  &useSurfCond,
                                  refConditions,
                                  &iwelemRef,
                                  &wellControlsName] ( localIndex const )
      {
        // Refresh separator properties at the selected control/reference conditions.
        if constexpr ( IS_THERMAL )
        {
          fluidSeparatorWrapper.update( iwelemRef, 0, refConditions.pressure, refConditions.temperature );
        }
        else
        {
          fluidSeparatorWrapper.update( iwelemRef, 0, refConditions.pressure );
        }

        if( useSurfCond && logSurfaceCondition )
        {
          GEOS_LOG_RANK( GEOS_FMT( "{}: surface density computed with P_surface = {} Pa and T_surface = {} K",
                                   wellControlsName, refConditions.pressure, refConditions.temperature ) );
        }

#ifdef GEOS_USE_HIP
        GEOS_UNUSED_VAR( wellControlsName );
#endif
      } );
    } );
  } );
}

void SinglePhaseWell::precomputeReferenceConditions( real64 const time_n,
                                                     Group & meshBodies,
                                                     MeshBody & meshBody,
                                                     WellElementSubRegion const & subRegion )
{
  GEOS_UNUSED_VAR( subRegion );
  if( !useSurfaceConditions() )
  {
    string_view refRegionName = referenceReservoirRegion();
    bool const useSegmentValues = refRegionName.empty();
    if( useSegmentValues )
    {
      setRegionAveragePressure( -1 );
      setRegionAverageTemperature( -1 );
    }
    else
    {
      auto & flowSolver = getParent().getGroup< SinglePhaseBase >( getFlowSolverName() );
      MeshLevel & flowMeshLevel = meshBody.getMeshLevel( flowSolver.getDiscretizationName() );

      if( !m_reservoirStatsAggregator )
      { // lazily initialize the region statistics aggregator
        m_reservoirStatsAggregator = std::make_unique< StatsAggregator >( getDataContext(),
                                                                          meshBodies,
                                                                          false );
        m_reservoirStatsAggregator->initStatisticsAggregation( flowSolver );
        m_reservoirStatsAggregator->enableRegionStatisticsAggregation();
      }

      RegionStatistics & stats = m_reservoirStatsAggregator->getRegionStatistics( flowMeshLevel, refRegionName );

      // compute region stats only if needed (could have already been done for another subRegion)
      if( !m_reservoirStatsAggregator->isComputed( time_n, stats ) )
        m_reservoirStatsAggregator->computeRegionsStatistics( time_n );

      GEOS_WARNING_IF( stats.m_averagePressure <= 0.0,
                       GEOS_FMT( "No region average quantities computed in reference region '{}'.",
                                 referenceReservoirRegion() ),
                       getWrapperDataContext( WellControls::viewKeyStruct::referenceReservoirRegionString() ),
                       getDataContext() );

      setRegionAveragePressure( stats.m_averagePressure );
      setRegionAverageTemperature( stats.m_averageTemperature );
    }
  }
}

SinglePhaseWell::ReferenceConditions
SinglePhaseWell::getReferenceConditions( WellElementSubRegion const & subRegion )
{
  if( useSurfaceConditions() )
  {
    // use surface conditions
    return {
      /* .pressure = */ getSurfacePressure(),
      /* .temperature = */ getSurfaceTemperature(),
    };
  }
  else
  {
    if( getRegionAveragePressure() > 0.0 && getRegionAverageTemperature() > 0.0 )
    { // reference region condition properly computed, we can return them
      return {
        /* .pressure = */ getRegionAveragePressure(),
        /* .temperature = */ getRegionAverageTemperature(),
      };
    }
    else
    { // region average stats not initialized or initialized, fallback to top segment values
      static bool firstNoRefRegionMsg = true;
      if( firstNoRefRegionMsg )
      {
        GEOS_WARNING( "SinglePhaseWell: region average statsistics of reference region not initialized,"
                      " fallback to top segment values.",
                      getDataContext() );
        firstNoRefRegionMsg=false;
      }


      arrayView1d< real64 const > const & pres = subRegion.getField< well::pressure >();
      arrayView1d< real64 const > const & temp = subRegion.getField< well::temperature >();
      localIndex const iwelemRef = subRegion.getTopWellElementIndex();
      return {
        /* .pressure = */ pres[iwelemRef],
        /* .temperature = */ temp[iwelemRef],
      };
    }
  }
}

real64 SinglePhaseWell::updateSubRegionState( real64 const time_n,
                                              MeshBody const & meshBody,
                                              ElementRegionManager const & elemManager,
                                              WellElementSubRegion & subRegion )
{
  GEOS_UNUSED_VAR( time_n );
  GEOS_UNUSED_VAR( meshBody );

  if( getWellState())
  {

    // update volumetric rates for the well constraints
    // Warning! This must be called before updating the fluid model
    //calculateReferenceElementRates( subRegion );

    // update density in the well elements
    updateFluidModel( subRegion );
    updateSeparator( elemManager, subRegion );   //  Calculate fluid properties at control conditions

    // Calculate the reference element rates
    calculateReferenceElementRates( subRegion );
    // update the current BHP
    updateBHPForConstraint( subRegion );

    // Broad case the updated well state to other ranks
    // TODO: add the missing getters on SinglePhaseWell & WellElementSubRegion because look-up is not useful here.
    real64 & currentBHP =
      getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentBHPString() );
    real64 & currentTotalVolRate =
      getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentVolRateString() );
    integer topRank =
      subRegion.getReference< integer >( WellElementSubRegion::viewKeyStruct::topRankString() );
    MpiWrapper::broadcast( currentBHP, topRank );
    MpiWrapper::broadcast( currentTotalVolRate, topRank );
    WellConstraintBase * constraint = getCurrentConstraint();
    if( constraint != nullptr )
    {
      constraint->setBHP ( getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentBHPString() ));
      constraint->setTotalVolumeRate ( getReference< real64 >(
                                         SinglePhaseWell::viewKeyStruct::currentVolRateString() ));
    }
  }
  return 0.0; // change in phasevolume fraction doesnt apply
}
void SinglePhaseWell::initializeWell( DomainPartition & domain, Group & meshBodies, string const & meshBodyName, MeshLevel & mesh, WellElementSubRegion & subRegion, real64 const & time_n )
{
  GEOS_UNUSED_VAR( domain );

  PerforationData const & perforationData = *subRegion.getPerforationData();
  ElementRegionManager const & elemManager = mesh.getElemManager();
  // get the info stored on well elements
  arrayView1d< real64 const > const wellElemGravCoef =
    subRegion.getField< well::gravityCoefficient >();

  // get well primary variables on well elements
  arrayView1d< real64 > const wellElemPressure =
    subRegion.getField< well::pressure >();
  arrayView1d< real64 > const connRate =
    subRegion.getField< well::connectionRate >();
  arrayView1d< real64 > const wellElemTemperature =
    subRegion.getField< well::temperature >();
  // get the element region, subregion, index
  arrayView1d< localIndex const > const resElementRegion =
    perforationData.getField< perforation::reservoirElementRegion >();
  arrayView1d< localIndex const > const resElementSubRegion =
    perforationData.getField< perforation::reservoirElementSubRegion >();
  arrayView1d< localIndex const > const resElementIndex =
    perforationData.getField< perforation::reservoirElementIndex >();

  arrayView1d< real64 const > const & perfGravCoef =
    perforationData.getField< well::gravityCoefficient >();

  bool const hasNonZeroRate = MpiWrapper::max< integer >( hasNonZero( connRate ));

  if( time_n <= 0.0  || ( isWellOpen() && !hasNonZeroRate ) )
  {
    setWellState( true );
    if( getCurrentConstraint() == nullptr )
    {
      if( isProducer() )
      {
        forSubGroups< MinimumBHPConstraint, ProductionConstraint< VolumeRateConstraint >, ProductionConstraint< MassRateConstraint >,
                      ProductionConstraint< PhaseVolumeRateConstraint > >( [&]( auto & constraint )
        {
          if( ConstraintTypeId( getControl()) == constraint.getControl() )
          {
            setCurrentConstraint( &constraint );
          }
        } );
      }
      else
      {
        forSubGroups< MaximumBHPConstraint, InjectionConstraint< VolumeRateConstraint >, InjectionConstraint< MassRateConstraint >,
                      InjectionConstraint< PhaseVolumeRateConstraint > >( [&]( auto & constraint )
        {
          if( ConstraintTypeId( getControl()) == constraint.getControl() )
          {
            setCurrentConstraint( &constraint );
          }
        } );
      }
    }

    PresTempInitializationKernel::SinglePhaseFlowAccessors resSinglePhaseFlowAccessors( elemManager, getFlowSolverName());
    PresTempInitializationKernel::SingleFluidAccessors resSingleFluidAccessors( elemManager, getFlowSolverName() );
    real64 const refWellElemGravCoef = getBHPReferenceGravityCoef( *this, ConstraintSourceId::USER );
    // 1) Loop over all perforations to compute an average density
    // 2) Initialize the reference pressure
    // 3) Estimate the pressures in the well elements using the average density
    PresTempInitializationKernel::
      launch( isThermal(),
              perforationData.size(),
              subRegion.size(),
              perforationData.getNumPerforationsGlobal(),
              *this,
              refWellElemGravCoef,
              0.0,           // initialization done at t = 0
              resSinglePhaseFlowAccessors.get( flow::pressure{} ),
              resSinglePhaseFlowAccessors.get( flow::temperature{} ),
              resSingleFluidAccessors.get( fields::singlefluid::density{} ),
              resElementRegion,
              resElementSubRegion,
              resElementIndex,
              perfGravCoef,
              wellElemGravCoef,
              wellElemPressure,
              wellElemTemperature );

    // 4) Recompute the pressure-dependent properties
    // Note: I am leaving that here because I would like to use the perforationRates (computed in UpdateState)
    //       to better initialize the rates
    MeshBody & meshBody = domain.getMeshBody( meshBodyName );
    precomputeReferenceConditions( time_n, meshBodies, meshBody, subRegion );
    updateSubRegionState( time_n, meshBody, elemManager, subRegion );

    string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
    SingleFluidBase & fluid = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );
    arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & wellElemDens = fluid.density();

    // 5) Estimate the well rates
    RateInitializationKernel::launch( subRegion.size(),
                                      *this,
                                      0.0,     // initialization done at t = 0
                                      wellElemDens,
                                      connRate );

    calculateReferenceElementRates( subRegion );
    WellConstraintBase * constraint =  getCurrentConstraint();
    constraint->setBHP ( getReference< real64 >( WellControls::viewKeyStruct::currentBHPString() ));
    constraint->setTotalVolumeRate ( getReference< real64 >(
                                       SinglePhaseWell::viewKeyStruct::currentVolRateString() ));
    // 7) Copy well / fluid dofs to "prop"_n variables
    saveState( subRegion );
  }
  else if( !hasNonZeroRate )
  {
    setWellState( false );
  }
  else
  {
    setWellState( true );
    // setup for restart
    if( getCurrentConstraint() == nullptr )
    {
      if( isProducer() )
      {
        forSubGroups< MinimumBHPConstraint, ProductionConstraint< VolumeRateConstraint >, ProductionConstraint< MassRateConstraint >,
                      ProductionConstraint< PhaseVolumeRateConstraint > >( [&](
                                                                             auto
                                                                             & constraint )
        {
          if( ConstraintTypeId( getControl()) == constraint.getControl()  )
          {
            setCurrentConstraint( &constraint );
          }
        } );
      }
      else
      {
        forSubGroups< MaximumBHPConstraint, InjectionConstraint< VolumeRateConstraint >, InjectionConstraint< MassRateConstraint >, InjectionConstraint< PhaseVolumeRateConstraint > >( [&](
                                                                                                                                                                                          auto
                                                                                                                                                                                          &
                                                                                                                                                                                          constraint )
        {
          if( ConstraintTypeId( getControl()) == constraint.getControl()  )
          {
            setCurrentConstraint( &constraint );
          }
        } );
      }
      updateSubRegionState( time_n, domain.getMeshBody( meshBodyName ), elemManager, subRegion );
    }

  }

}

real64 SinglePhaseWell::updateWellState( MeshBody const & meshBody,
                                         ElementRegionManager const & elemManager,
                                         WellElementSubRegion & subRegion )
{
  GEOS_MARK_FUNCTION;

  updateSubRegionState( -1.0, meshBody, elemManager, subRegion );
  return 0.0;
}

void SinglePhaseWell::assembleWellFluxTerms( real64 const & time,
                                             real64 const & dt,
                                             WellElementSubRegion const & subRegion,
                                             DofManager const & dofManager,
                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                             arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( time );

  // get a reference to the degree-of-freedom numbers
  string const wellDofKey = dofManager.getKey( wellElementDofName() );

  if( isThermal() )
  {
    string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
    SingleFluidBase const & fluid = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );
    thermalSinglePhaseWellKernels::
      FaceBasedAssemblyKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( dt,
                                                 dofManager.rankOffset(),
                                                 wellDofKey,
                                                 *this,
                                                 subRegion,
                                                 fluid,
                                                 localMatrix,
                                                 localRhs );
  }
  else
  {
    singlePhaseWellKernels::
      FaceBasedAssemblyKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( dt,
                                                 dofManager.rankOffset(),
                                                 wellDofKey,
                                                 *this,
                                                 subRegion,
                                                 localMatrix,
                                                 localRhs );
  }

}


void SinglePhaseWell::assembleWellConstraintTerms( real64 const & time_n,
                                                   real64 const & GEOS_UNUSED_PARAM( dt ),
                                                   WellElementSubRegion const & subRegion,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;


  // the rank that owns the reference well element is responsible for the calculations below.

  if( !subRegion.isLocallyOwned() || !(  getWellStatus() == WellControls::Status::OPEN ))
  {
    return;
  }
  {
    forSubGroups< MinimumBHPConstraint, MaximumBHPConstraint, InjectionConstraint< VolumeRateConstraint >, ProductionConstraint< VolumeRateConstraint > >( [&]( auto & constraint )
    {
      if( constraint.getName() == getCurrentConstraint()->getName())
      {
        // found limiting constraint

        // fluid data
        constitutive::SingleFluidBase & fluidSeparator =   getSingleFluidSeparator();
        integer isThermal = fluidSeparator.isThermal();

        geos::internal::kernelLaunchSelectorThermalSwitch( isThermal, [&] ( auto ISTHERMAL )
        {
          integer constexpr IS_THERMAL = ISTHERMAL();

          singlePhaseWellConstraintKernels::ConstraintHelper< IS_THERMAL >::assembleConstraintEquation( time_n,
                                                                                                        *this,
                                                                                                        constraint,
                                                                                                        subRegion,
                                                                                                        dofManager.getKey( wellElementDofName() ),
                                                                                                        dofManager.rankOffset(),
                                                                                                        localMatrix,
                                                                                                        localRhs );
        } );
      }

    } );
  }

}

void SinglePhaseWell::assembleWellPressureRelations( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                                     real64 const & GEOS_UNUSED_PARAM( dt ),
                                                     WellElementSubRegion const & subRegion,
                                                     DofManager const & dofManager,
                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                     arrayView1d< real64 > const & localRhs )
{

  // get the degrees of freedom numbers, depth, next well elem index
  string const wellDofKey = dofManager.getKey( wellElementDofName() );
  arrayView1d< globalIndex const > const & wellElemDofNumber =
    subRegion.getReference< array1d< globalIndex > >( wellDofKey );
  arrayView1d< real64 const > const & wellElemGravCoef =
    subRegion.getField< well::gravityCoefficient >();
  arrayView1d< localIndex const > const & nextWellElemIndex =
    subRegion.getReference< array1d< localIndex > >( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString() );

  // get primary variables on well elements
  arrayView1d< real64 const > const & wellElemPressure =
    subRegion.getField< well::pressure >();

  // get well constitutive data
  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  SingleFluidBase const & fluid = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & wellElemDensity = fluid.density();
  arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const & dWellElemDensity = fluid.dDensity();

  geos::internal::kernelLaunchSelectorThermalSwitch( isThermal(), [&] ( auto ISTHERMAL )
  {
    PressureRelationKernel::launch< ISTHERMAL >( subRegion.size(),
                                                 dofManager.rankOffset(),
                                                 wellElemDofNumber,
                                                 wellElemGravCoef,
                                                 nextWellElemIndex,
                                                 wellElemPressure,
                                                 wellElemDensity,
                                                 dWellElemDensity,
                                                 localMatrix,
                                                 localRhs );
  } );

}

void SinglePhaseWell::assembleWellAccumulationTerms( real64 const & time,
                                                     real64 const & dt,
                                                     WellElementSubRegion & subRegion,
                                                     DofManager const & dofManager,
                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                     arrayView1d< real64 > const & localRhs )

{
  GEOS_UNUSED_VAR( time );
  GEOS_UNUSED_VAR( dt );
  // get a reference to the degree-of-freedom numbers
  string const wellElemDofKey = dofManager.getKey( wellElementDofName() );

  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  SingleFluidBase const & fluid = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );
  if( getWellStatus() == WellControls::Status::OPEN && !m_keepVariablesConstantDuringInitStep )
  {
    if( isThermal() )
    {
      thermalSinglePhaseWellKernels::
        ElementBasedAssemblyKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( isProducer(),
                                                   dofManager.rankOffset(),
                                                   wellElemDofKey,
                                                   subRegion,
                                                   fluid,
                                                   localMatrix,
                                                   localRhs );
    }
    else
    {
      singlePhaseWellKernels::
        ElementBasedAssemblyKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                   wellElemDofKey,
                                                   subRegion,
                                                   fluid,
                                                   localMatrix,
                                                   localRhs );
    }
    shutClosedSegments( subRegion, dofManager, localMatrix, localRhs );
  }
  else
  {
    shutEntireWell( subRegion, dofManager, localMatrix, localRhs );
  }
}

void SinglePhaseWell::resetShutInControlState()
{
  WellControls::resetShutInControlState();
  getReference< real64 >( WellControls::viewKeyStruct::currentVolRateString() ) = 0.0;
}

void SinglePhaseWell::computeWellPerforationRates( real64 const & time_n,
                                                   real64 const & GEOS_UNUSED_PARAM( dt ),
                                                   ElementRegionManager & elemManager,
                                                   WellElementSubRegion & subRegion )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( time_n );
  PerforationData * const perforationData = subRegion.getPerforationData();
  if( isWellOpen() && !m_keepVariablesConstantDuringInitStep )
  {
    string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
    SingleFluidBase const & fluid = getConstitutiveModel< SingleFluidBase >( subRegion, fluidName );

    if( isThermal() )
    {
      thermalSinglePhasePerforationFluxKernels::
        PerforationFluxKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( getFlowSolverName(),
                                                   perforationData,
                                                   subRegion,
                                                   fluid,
                                                   elemManager );
    }
    else
    {
      isothermalSinglePhasePerforationFluxKernels::
        PerforationFluxKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( getFlowSolverName(),
                                                   perforationData,
                                                   subRegion,
                                                   fluid,
                                                   elemManager );
    }
  }
  else
  {
    // Zero completion flow rate.
    arrayView1d< real64 > const perfRate = perforationData->getField< fields::well::perforationRate >();
    for( integer iperf=0; iperf<perforationData->size(); iperf++ )
    {
      perfRate[iperf] = 0.0;
    }
  }
}
real64
SinglePhaseWell::scalingForWellSystemSolution( WellElementSubRegion & subRegion,
                                               DofManager const & dofManager,
                                               arrayView1d< real64 const > const & localSolution )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( subRegion );
  GEOS_UNUSED_VAR( dofManager );
  GEOS_UNUSED_VAR( localSolution );

  return 1.0;
}

array1d< real64 >
SinglePhaseWell::calculateLocalWellResidualNorm( real64 const & time_n,
                                                 real64 const & dt,
                                                 NonlinearSolverParameters const & nonlinearSolverParameters,
                                                 WellElementSubRegion const & subRegion,
                                                 DofManager const & dofManager,
                                                 arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  integer numNorm = 1; // mass balance
  array1d< real64 > localResidualNorm;
  array1d< real64 > localResidualNormalizer;

  if( isThermal() )
  {
    numNorm = 2; // mass balance and energy balance
  }
  localResidualNorm.resize( numNorm );
  localResidualNormalizer.resize( numNorm );


  globalIndex const rankOffset = dofManager.rankOffset();
  string const wellDofKey = dofManager.getKey( wellElementDofName() );



  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  SingleFluidBase const & fluid = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );



  if( isWellOpen() )
  {
    // step 1: compute the norm in the subRegion
    if( isThermal() )
    {
      real64 subRegionResidualNorm[2]{};
      thermalSinglePhaseWellKernels::ResidualNormKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( rankOffset,
                                                   wellDofKey,
                                                   localRhs,
                                                   subRegion,
                                                   fluid,
                                                   *this,
                                                   time_n,
                                                   dt,
                                                   nonlinearSolverParameters.m_minNormalizer,
                                                   subRegionResidualNorm );
      // step 2: reduction across meshBodies/regions/subRegions

      for( integer i=0; i<numNorm; i++ )
      {
        if( subRegionResidualNorm[i] > localResidualNorm[i] )
        {
          localResidualNorm[i] = subRegionResidualNorm[i];
        }
      }
    }
    else
    {
      real64 subRegionResidualNorm[1]{};
      ResidualNormKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( rankOffset,
                                                   wellDofKey,
                                                   localRhs,
                                                   subRegion,
                                                   fluid,
                                                   *this,
                                                   time_n,
                                                   dt,
                                                   nonlinearSolverParameters.m_minNormalizer,
                                                   subRegionResidualNorm );

      // step 2: reduction across meshBodies/regions/subRegions
      if( subRegionResidualNorm[0] > localResidualNorm[0] )
      {
        localResidualNorm[0] = subRegionResidualNorm[0];
      }
    }
  }

  return localResidualNorm;
}

real64
SinglePhaseWell::calculateWellResidualNorm( real64 const & time_n,
                                            real64 const & dt,
                                            NonlinearSolverParameters const & nonlinearSolverParameters,
                                            WellElementSubRegion const & subRegion,
                                            DofManager const & dofManager,
                                            arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  integer numNorm = 1;     // mass balance
  array1d< real64 > localResidualNorm;
  array1d< real64 > localResidualNormalizer;

  if( isThermal() )
  {
    numNorm = 2;     // mass balance and energy balance
  }
  localResidualNorm.resize( numNorm );
  localResidualNormalizer.resize( numNorm );


  //globalIndex const rankOffset = dofManager.rankOffset();
  string const wellDofKey = dofManager.getKey( wellElementDofName() );

  if( isWellOpen() )
  {
    localResidualNorm = calculateLocalWellResidualNorm( time_n,
                                                        dt,
                                                        nonlinearSolverParameters,
                                                        subRegion,
                                                        dofManager,
                                                        localRhs );
  }
  real64 resNorm=localResidualNorm[0];
  if( isThermal() )
  {
    real64 globalResidualNorm[2]{};
    globalResidualNorm[0] = MpiWrapper::max( localResidualNorm[0] );
    globalResidualNorm[1] = MpiWrapper::max( localResidualNorm[1] );
    resNorm= std::sqrt( globalResidualNorm[0] * globalResidualNorm[0] + globalResidualNorm[1] * globalResidualNorm[1] );

    GEOS_LOG_LEVEL_RANK_0_NLR( logInfo::ResidualNorm, GEOS_FMT( "        ( R{} ) = ( {:4.2e} )        ( Renergy ) = ( {:4.2e} )",
                                                                coupledSolverAttributePrefix(), globalResidualNorm[0], globalResidualNorm[1] ));

    //getConvergenceStats().setResidualValue( GEOS_FMT( "R{}", coupledSolverAttributePrefix()), globalResidualNorm[0] );
    //getConvergenceStats().setResidualValue( "Renergy", globalResidualNorm[1] );
  }
  else
  {
    resNorm = MpiWrapper::max( resNorm );

    GEOS_LOG_LEVEL_RANK_0_NLR( logInfo::ResidualNorm, GEOS_FMT( "        ( R{} ) = ( {:4.2e} )",
                                                                coupledSolverAttributePrefix(), resNorm ));
    //getConvergenceStats().setResidualValue( GEOS_FMT( "R{}", coupledSolverAttributePrefix()), resNorm );
  }

  return resNorm;
}
bool SinglePhaseWell::checkWellSystemSolution( WellElementSubRegion & subRegion,
                                               DofManager const & dofManager,
                                               arrayView1d< real64 const > const & localSolution,
                                               real64 const scalingFactor,
                                               real64 & minPressure,
                                               real64 & minDensity,
                                               real64 & minTotalDensity,
                                               ElementsReporterBuffer & negPressureIds,
                                               ElementsReporterBuffer & negDensityIds,
                                               ElementsReporterBuffer & negTotalDensityIds )
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( minDensity, minTotalDensity );
  GEOS_UNUSED_VAR( negDensityIds, negTotalDensityIds );

  string const wellDofKey = dofManager.getKey( wellElementDofName() );

  globalIndex const rankOffset = dofManager.rankOffset();
  // get the degree of freedom numbers on well elements
  arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( wellDofKey );
  arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();

  // get a reference to the primary variables on well elements
  arrayView1d< real64 const > const & pressure = subRegion.getField< well::pressure >();
  auto const negPresCollector = negPressureIds.createCollector( subRegion.localToGlobalMap().toViewConst() );

  using Kernel = singlePhaseBaseKernels::SolutionCheckKernel;

  auto const results = Kernel::launch< parallelDevicePolicy<> >( localSolution,
                                                                 rankOffset,
                                                                 dofNumber,
                                                                 ghostRank,
                                                                 pressure,
                                                                 scalingFactor,
                                                                 negPresCollector );

  minPressure = std::min( minPressure, results.minNegPres );

  return (m_allowNegativePressure || 0.0 < results.minNegPres) ? 1 : 0;
}

void
SinglePhaseWell::applyWellSystemSolution( DofManager const & dofManager,
                                          arrayView1d< real64 const > const & localSolution,
                                          real64 const scalingFactor,
                                          real64 const dt,
                                          DomainPartition & domain,
                                          MeshLevel & mesh,
                                          WellElementSubRegion & subRegion )
{
  GEOS_UNUSED_VAR( dt );
  GEOS_UNUSED_VAR( subRegion );
  DofManager::CompMask pressureMask( m_numDofPerWellElement, 0, 1 );
  DofManager::CompMask connRateMask( m_numDofPerWellElement, 1, 2 );
  dofManager.addVectorToField( localSolution,
                               wellElementDofName(),
                               well::pressure::key(),
                               scalingFactor,
                               pressureMask );

  dofManager.addVectorToField( localSolution,
                               wellElementDofName(),
                               well::connectionRate::key(),
                               scalingFactor,
                               connRateMask );

  if( isThermal() )
  {
    DofManager::CompMask temperatureMask( m_numDofPerWellElement, 2, 3 );

    dofManager.addVectorToField( localSolution,
                                 wellElementDofName(),
                                 fields::well::temperature::key(),
                                 scalingFactor,
                                 temperatureMask );

  }


  FieldIdentifiers fieldsToBeSync;
  if( isThermal() )
  {
    fieldsToBeSync.addElementFields( { well::pressure::key(),
                                       well::connectionRate::key(),
                                       well::temperature::key() },
                                     getTargetRegionNames() );
  }
  else
  {
    fieldsToBeSync.addElementFields( { well::pressure::key(),
                                       well::connectionRate::key() },
                                     getTargetRegionNames() );
  }
  CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                       mesh,
                                                       domain.getNeighbors(),
                                                       true );


}

void SinglePhaseWell::applyWellBoundaryConditions( real64 const time_n,
                                                   real64 const dt,
                                                   ElementRegionManager & elemManager,
                                                   WellElementSubRegion & subRegion,
                                                   DofManager const & dofManager,
                                                   arrayView1d< real64 > const & localRhs,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix )
{
  GEOS_UNUSED_VAR( time_n );
  GEOS_UNUSED_VAR( elemManager );

  if( !isWellOpen() )
  {
    return;
  }

  string const wellDofKey = dofManager.getKey( wellElementDofName() );
  PerforationData const * const perforationData = subRegion.getPerforationData();

  if( isThermal() )
  {
    coupledReservoirAndSinglePhaseWellKernels::
      ThermalSinglePhaseWellFluxKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( isProducer(),
                                                 dt,
                                                 dofManager.rankOffset(),
                                                 wellDofKey,
                                                 subRegion,
                                                 perforationData,
                                                 localMatrix,
                                                 localRhs );
  }
  else
  {
    coupledReservoirAndSinglePhaseWellKernels::
      IsothermalSinglePhaseWellFluxKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( dt,
                                                 dofManager.rankOffset(),
                                                 wellDofKey,
                                                 subRegion,
                                                 perforationData,
                                                 localMatrix,
                                                 localRhs );
  }
}


void SinglePhaseWell::resetStateToBeginningOfStep( DomainPartition & domain,
                                                   string const & meshBodyName, ElementRegionManager const & elemManager, WellElementSubRegion & subRegion )
{


  // get a reference to the primary variables on well elements
  arrayView1d< real64 > const & wellElemPressure =
    subRegion.getField< well::pressure >();
  arrayView1d< real64 const > const & wellElemPressure_n =
    subRegion.getField< well::pressure_n >();
  wellElemPressure.setValues< parallelDevicePolicy<> >( wellElemPressure_n );

  if( isThermal() )
  {
    arrayView1d< real64 > const & wellElemTemperature =
      subRegion.getField< fields::well::temperature >();
    arrayView1d< real64 const > const & wellElemTemperature_n =
      subRegion.getField< fields::well::temperature_n >();
    wellElemTemperature.setValues< parallelDevicePolicy<> >( wellElemTemperature_n );
  }
  arrayView1d< real64 > const & connRate =
    subRegion.getField< well::connectionRate >();
  arrayView1d< real64 const > const & connRate_n =
    subRegion.getField< well::connectionRate_n >();
  connRate.setValues< parallelDevicePolicy<> >( connRate_n );
  updateSubRegionState( -1.0, domain.getMeshBody( meshBodyName ), elemManager, subRegion );

}

void SinglePhaseWell::saveState( WellElementSubRegion & subRegion )
{
  arrayView1d< real64 const > const wellElemPressure = subRegion.getField< well::pressure >();
  arrayView1d< real64 > const wellElemPressure_n = subRegion.getField< well::pressure_n >();
  wellElemPressure_n.setValues< parallelDevicePolicy<> >( wellElemPressure );

  if( isThermal() )
  {
    arrayView1d< real64 const > const wellElemTemperature = subRegion.getField< well::temperature >();
    arrayView1d< real64 > const wellElemTemperature_n = subRegion.getField< well::temperature_n >();
    wellElemTemperature_n.setValues< parallelDevicePolicy<> >( wellElemTemperature );
  }
  arrayView1d< real64 const > const connRate = subRegion.getField< well::connectionRate >();
  arrayView1d< real64 > const connRate_n = subRegion.getField< well::connectionRate_n >();
  connRate_n.setValues< parallelDevicePolicy<> >( connRate );

  SingleFluidBase const & fluid =
    getConstitutiveModel< SingleFluidBase >( subRegion, subRegion.getReference< string >( viewKeyStruct::fluidNamesString() ) );
  fluid.saveConvergedState();
}

void SinglePhaseWell::implicitStepSetup( real64 const & time_n,
                                         real64 const & dt,
                                         DomainPartition & domain,
                                         string const & meshBodyName,
                                         ElementRegionManager & elemManager,
                                         WellElementSubRegion & subRegion )
{
  GEOS_MARK_FUNCTION;
  WellControls::implicitStepSetup( time_n, dt, domain, meshBodyName, elemManager, subRegion );
  arrayView1d< real64 const > const wellElemPressure = subRegion.getField< well::pressure >();
  arrayView1d< real64 > const wellElemPressure_n = subRegion.getField< well::pressure_n >();
  wellElemPressure_n.setValues< parallelDevicePolicy<> >( wellElemPressure );

  if( isThermal() )
  {
    arrayView1d< real64 const > const wellElemTemperature = subRegion.getField< well::temperature >();
    arrayView1d< real64 > const wellElemTemperature_n = subRegion.getField< well::temperature_n >();
    wellElemTemperature_n.setValues< parallelDevicePolicy<> >( wellElemTemperature );
  }
  arrayView1d< real64 const > const connRate = subRegion.getField< well::connectionRate >();
  arrayView1d< real64 > const connRate_n = subRegion.getField< well::connectionRate_n >();
  connRate_n.setValues< parallelDevicePolicy<> >( connRate );

  SingleFluidBase const & fluid =
    getConstitutiveModel< SingleFluidBase >( subRegion, subRegion.getReference< string >( viewKeyStruct::fluidNamesString() ) );
  fluid.saveConvergedState();

  validateWellConstraints( time_n, dt, subRegion );

  updateSubRegionState( time_n, domain.getMeshBody( meshBodyName ), elemManager, subRegion );

}

void SinglePhaseWell::implicitStepComplete( real64 const & time_n,
                                            real64 const & dt,
                                            WellElementSubRegion const & subRegion )
{
  printRates( time_n, dt, subRegion );
}

void SinglePhaseWell::printRates( real64 const & time_n,
                                  real64 const & dt,
                                  WellElementSubRegion const & subRegion )
{

  GEOS_UNUSED_VAR( dt );
  // the rank that owns the reference well element is responsible for the calculations below.
  if( !subRegion.isLocallyOwned() )
  {
    return;
  }

  localIndex const iwelemRef = subRegion.getTopWellElementIndex();

  // subRegion data

  arrayView1d< real64 const > const & connRate =
    subRegion.getField< well::connectionRate >();

  // control data


  string const wellControlsName =  getName();

  // format: time,total_rate,total_vol_rate
  std::ofstream outputFile;
  if( m_writeCSV > 0 )
  {
    outputFile.open( m_ratesOutputDir + "/" + wellControlsName + ".csv", std::ios_base::app );
    outputFile << time_n;
  }

  if( !isWellOpen() )
  {
    GEOS_LOG( GEOS_FMT( "{}: well is shut", wellControlsName ) );
    if( outputFile.is_open())
    {
      // print all zeros in the rates file
      outputFile << ",0.0,0.0,0.0" << std::endl;
      outputFile.close();
    }
    return;
  }

  integer const useSurfaceCond =  useSurfaceConditions();

  real64 const & currentBHP =
    getReference< real64 >( WellControls::viewKeyStruct::currentBHPString() );
  real64 const & currentTotalVolRate =
    getReference< real64 >( WellControls::viewKeyStruct::currentVolRateString() );

  // bring everything back to host, capture the scalars by reference
  forAll< serialPolicy >( 1, [&useSurfaceCond,
                              &currentBHP,
                              connRate,
                              &currentTotalVolRate,
                              &iwelemRef,
                              &wellControlsName,
                              &outputFile] ( localIndex const )
  {
    string const conditionKey = useSurfaceCond ? "surface" : "reservoir";
    string const unitKey = useSurfaceCond ? "s" : "r";

    real64 const currentTotalRate = connRate[iwelemRef];
    GEOS_LOG( GEOS_FMT( "{}: BHP (at the specified reference elevation): {} Pa",
                        wellControlsName, currentBHP ) );
    GEOS_LOG( GEOS_FMT( "{}: Total rate: {} kg/s; total {} volumetric rate: {} {}m3/s",
                        wellControlsName, currentTotalRate, conditionKey, currentTotalVolRate, unitKey ) );
    if( outputFile.is_open())
    {
      outputFile << "," << currentBHP;
      outputFile << "," << currentTotalRate << "," << currentTotalVolRate << std::endl;
      outputFile.close();
    }
  } );
}



}// namespace geos
