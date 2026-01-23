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
#include "mesh/WellElementSubRegion.hpp"
#include "mesh/PerforationFields.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "physicsSolvers/fluidFlow/wells/LogLevelsInfo.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWellFields.hpp"
#include "physicsSolvers/fluidFlow/wells/kernels/SinglePhaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/kernels/ThermalSinglePhaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/kernels/SinglePhasePerforationFluxKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/FluidUpdateKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/SolutionCheckKernel.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseStatistics.hpp"
#include "physicsSolvers/fluidFlow/wells/kernels/SinglePhaseWellConstraintKernels.hpp"
namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace fields;
using namespace singlePhaseWellKernels;

SinglePhaseWell::SinglePhaseWell( const string & name,
                                  Group * const parent ):
  WellControls( name, parent )
{
  m_numDofPerWellElement = 2;
  m_numDofPerResElement = 1;
  m_numPhases = 1;
  m_numComponents = 1;

  this->registerWrapper( FlowSolverBase::viewKeyStruct::allowNegativePressureString(), &m_allowNegativePressure ).
    setApplyDefaultValue( 1 ). // negative pressure is allowed by default
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag indicating if negative pressure is allowed" );
}

void SinglePhaseWell::registerWellDataOnMesh( WellElementSubRegion & subRegion )
{
  WellControls::registerDataOnMesh( subRegion );

  setConstitutiveNames ( subRegion );

  //DomainPartition const & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  //ConstitutiveManager const & cm = domain.getConstitutiveManager();
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

  registerWrapper< array1d< real64 > >( viewKeyStruct::dCurrentBHPString() ).
    setSizedFromParent( 0 ).
    reference().resizeDimension< 0 >( 2 + isThermal() );       // dP, dT , dQ

  registerWrapper< real64 >( viewKeyStruct::currentVolRateString() );
  registerWrapper< array1d< real64 > >( viewKeyStruct::dCurrentVolRateString() ).
    setSizedFromParent( 0 ).
    reference().resizeDimension< 0 >( 2 + isThermal() );       // dP, dT, dQ

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
  GEOS_UNUSED_VAR( time_n ); // tjb this will be needed with validation against tables
  GEOS_UNUSED_VAR( subRegion );
  if( useSurfaceConditions() )
  {
    bool useSeg =  getReferenceReservoirRegion().empty();
    GEOS_WARNING_IF( useSeg,
                     "WellControls referenceReservoirRegion not set and well constraint fluid property calculations will use top segement pressure and temp " );
    if( useSeg )
    {
      setRegionAveragePressure( -1 );
    }
    else
    {
      // Check if region name exists in list of Reservoir's target regions
      string const regionName =  getReferenceReservoirRegion();
      SinglePhaseBase const & flowSolver = getParent().getGroup< SinglePhaseBase >( getFlowSolverName() );
      string_array const & targetRegionsNames = flowSolver.getTargetRegionNames();
      auto const pos = std::find( targetRegionsNames.begin(), targetRegionsNames.end(), regionName );
      GEOS_ERROR_IF( pos == targetRegionsNames.end(),
                     GEOS_FMT( "{}: Region {} is not a target of the reservoir solver and cannot be used for referenceReservoirRegion in WellControl {}.",
                               getDataContext(), regionName, getName() ) );

    }
  }
}

void SinglePhaseWell::initializeWellPostInitialConditionsPreSubGroups( WellElementSubRegion & subRegion )
{

  // set gravity coefficient
  setGravCoef( subRegion, getParent().getParent().getReference< R1Tensor >( PhysicsSolverManager::viewKeyStruct::gravityVectorString() ));

  // setup fluid model
  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  constitutive::SingleFluidBase & fluid = subRegion.getConstitutiveModel< constitutive::SingleFluidBase >( fluidName );
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

  // control data

  string const wellControlsName =  getName();
  real64 const & refGravCoef =  getReferenceGravityCoef();

  real64 & currentBHP =
    getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentBHPString() );
  arrayView1d< real64 > const & dCurrentBHP =
    getReference< array1d< real64 > >( SinglePhaseWell::viewKeyStruct::dCurrentBHPString() );

  geos::internal::kernelLaunchSelectorThermalSwitch( isThermal(), [&] ( auto ISTHERMAL )
  {
    integer constexpr IS_THERMAL = ISTHERMAL();
    // bring everything back to host, capture the scalars by reference
    forAll< serialPolicy >( 1, [pres,
                                dens,
                                dDens,
                                wellElemGravCoef,
                                &currentBHP,
                                &dCurrentBHP,
                                &iwelemRef,
                                &refGravCoef] ( localIndex const )
    {
      real64 const diffGravCoef = refGravCoef - wellElemGravCoef[iwelemRef];
      currentBHP = pres[iwelemRef] + dens[iwelemRef][0] * diffGravCoef;
      dCurrentBHP[DerivOffset::dP] = 1.0 + dDens[iwelemRef][0][DerivOffset::dP] *diffGravCoef;
      if constexpr ( IS_THERMAL )
      {
        dCurrentBHP[DerivOffset::dT] =  dDens[iwelemRef][0][DerivOffset::dT] * diffGravCoef;
      }
    } );
  } );

  GEOS_LOG_LEVEL_BY_RANK( logInfo::WellControl,
                          GEOS_FMT( "{}: The BHP (at the specified reference elevation) = {} Pa",
                                    wellControlsName, currentBHP ) );

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


  // control data


  // fluid data
  constitutive::SingleFluidBase & fluidSeparator =   getSingleFluidSeparator();
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & dens = fluidSeparator.density();

  real64 & currentVolRate =
    getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentVolRateString() );

  // bring everything back to host, capture the scalars by reference
  forAll< serialPolicy >( 1, [connRate,
                              dens,
                              &currentVolRate,
                              &iwelemRef] ( localIndex const )
  {
    real64 const densInv = 1.0 / dens[iwelemRef][0];
    currentVolRate = connRate[iwelemRef] * densInv;
    // tjb compute mass
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

  // the rank that owns the reference well element is responsible for the calculations below.
  if( !subRegion.isLocallyOwned() )
  {
    return;
  }

  localIndex const iwelemRef = subRegion.getTopWellElementIndex();

  // subRegion data
  arrayView1d< real64 const > const pres =
    subRegion.getField< well::pressure >();

  // control data

  string const wellControlsName =  getName();
  bool const logSurfaceCondition = isLogLevelActive< logInfo::WellControl >( getLogLevel());
  integer const useSurfaceCond =  useSurfaceConditions();

  // fluid data
  constitutive::SingleFluidBase & fluidSeparator =   getSingleFluidSeparator();
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & dens = fluidSeparator.density();

  real64 flashPressure;
  if( useSurfaceCond )
  {
    // use surface conditions
    flashPressure =  getSurfacePressure();
  }
  else
  {
    if( !getReferenceReservoirRegion().empty() )
    {
      ElementRegionBase const & region = elemManager.getRegion( getReferenceReservoirRegion());
      GEOS_ERROR_IF ( !region.hasWrapper( SinglePhaseStatistics::regionStatisticsName() ),
                      GEOS_FMT( "{}: WellControl {} referenceReservoirRegion field requires SinglePhaseStatistics to be configured for region {} ",
                                getDataContext(), getName(), getReferenceReservoirRegion() ) );

      SinglePhaseStatistics::RegionStatistics const & stats = region.getReference< SinglePhaseStatistics::RegionStatistics >(
        SinglePhaseStatistics::regionStatisticsName() );
      setRegionAveragePressure( stats.averagePressure );
      setRegionAverageTemperature( stats.averageTemperature );
      GEOS_ERROR_IF( stats.averagePressure <= 0.0,
                     GEOS_FMT( "{}: No region average quantities computed.  WellControl {} referenceReservoirRegion field requires CompositionalMultiphaseStatistics to be configured for region {} ",
                               getDataContext(), getName(), getReferenceReservoirRegion() ));

    }
    // use region conditions
    flashPressure =        getRegionAveragePressure();
    if( flashPressure < 0.0 )
    {
      // use segment conditions
      flashPressure   = pres[iwelemRef];
    }
  }

  constitutiveUpdatePassThru( fluidSeparator, [&]( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidSeparatorWrapper = castedFluid.createKernelWrapper();

    geos::internal::kernelLaunchSelectorThermalSwitch( isThermal(), [&] ( auto ISTHERMAL )
    {
      integer constexpr IS_THERMAL = ISTHERMAL();
      GEOS_UNUSED_VAR( IS_THERMAL );
      // bring everything back to host, capture the scalars by reference
      forAll< serialPolicy >( 1, [fluidSeparatorWrapper,
                                  pres,
                                  dens,
                                  logSurfaceCondition,
                                  &useSurfaceCond,
                                  &flashPressure,
                                  &iwelemRef,
                                  &wellControlsName] ( localIndex const )
      {
        //    We need to evaluate the density as follows:
        //      - Surface conditions: using the surface pressure provided by the user
        //      - Reservoir conditions: using the pressure in the top element

        if( useSurfaceCond )
        {
          // we need to compute the surface density
          fluidSeparatorWrapper.update( iwelemRef, 0, flashPressure );
          if( logSurfaceCondition )
          {

            GEOS_LOG_RANK( GEOS_FMT( "{}: surface density computed with P_surface = {} Pa",
                                     wellControlsName, flashPressure ) );
          }

#ifdef GEOS_USE_HIP
          GEOS_UNUSED_VAR( wellControlsName );
#endif

        }
        else
        {
          fluidSeparatorWrapper.update( iwelemRef, 0, flashPressure );
        }
      } );
    } );
  } );
}

real64 SinglePhaseWell::updateSubRegionState( ElementRegionManager const & elemManager, WellElementSubRegion & subRegion )
{

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


  }
  return 0.0; // change in phasevolume fraction doesnt apply
}
void SinglePhaseWell::initializeWell( DomainPartition & domain, MeshLevel & mesh, WellElementSubRegion & subRegion, real64 const & time_n )
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
            setControl( static_cast< WellControls::Control >(constraint.getControl()) );     // tjb old
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
            setControl( static_cast< WellControls::Control >(constraint.getControl()) );    // tjb old
          }
        } );
      }
    }

    PresTempInitializationKernel::SinglePhaseFlowAccessors resSinglePhaseFlowAccessors( elemManager, getFlowSolverName());
    PresTempInitializationKernel::SingleFluidAccessors resSingleFluidAccessors( elemManager, getFlowSolverName() );

    // 1) Loop over all perforations to compute an average density
    // 2) Initialize the reference pressure
    // 3) Estimate the pressures in the well elements using the average density
    PresTempInitializationKernel::
      launch( isThermal(),
              perforationData.size(),
              subRegion.size(),
              perforationData.getNumPerforationsGlobal(),
              *this,
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
    updateSubRegionState( elemManager, subRegion );

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
    constraint->setBHP ( getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentBHPString() ));
    constraint->setTotalVolumeRate ( getReference< real64 >(
                                       SinglePhaseWell::viewKeyStruct::currentVolRateString() ));
    //constraint->setMassRate( wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentMassRateString() ));
    // 7) Copy well / fluid dofs to "prop"_n variables
    saveState( subRegion );
  }
  else if( !hasNonZeroRate )
  {
    setWellState( false );
    GEOS_LOG_RANK_0( "tjb shut wells "<< subRegion.getName());
  }
  else
  {
    setWellState( true );
    // setup for restart
    if( getCurrentConstraint() == nullptr )
    {
      updateSubRegionState( elemManager, subRegion );
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
    }

  }

}

void SinglePhaseWell::shutDownWell( WellElementSubRegion & subRegion,
                                    DofManager const & dofManager,
                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                    arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  string const wellDofKey = dofManager.getKey( wellElementDofName() );



  if( isWellOpen(  ) )
  {
    return;
  }

  globalIndex const rankOffset = dofManager.rankOffset();

  arrayView1d< integer const > const ghostRank =
    subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString() );
  arrayView1d< globalIndex const > const dofNumber =
    subRegion.getReference< array1d< globalIndex > >( wellDofKey );

  arrayView1d< real64 const > const pres =
    subRegion.getField< fields::well::pressure >();
  arrayView1d< real64 const > const connRate =
    subRegion.getField< fields::well::connectionRate >();

  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    if( ghostRank[ei] >= 0 )
    {
      return;
    }

    globalIndex const dofIndex = dofNumber[ei];
    localIndex const localRow = dofIndex - rankOffset;
    real64 rhsValue;

    // 4.1. Apply pressure value to the matrix/rhs
    FieldSpecificationEqual::SpecifyFieldValue( dofIndex,
                                                rankOffset,
                                                localMatrix,
                                                rhsValue,
                                                pres[ei],           // freeze the current pressure value
                                                pres[ei] );
    localRhs[localRow] = rhsValue;

    // 4.2. Apply rate value to the matrix/rhs
    FieldSpecificationEqual::SpecifyFieldValue( dofIndex + 1,
                                                rankOffset,
                                                localMatrix,
                                                rhsValue,
                                                connRate[ei],           // freeze the current pressure value
                                                connRate[ei] );
    localRhs[localRow + 1] = rhsValue;

  } );

}
real64 SinglePhaseWell::updateWellState( ElementRegionManager const & elemManager, WellElementSubRegion & subRegion )
{
  GEOS_MARK_FUNCTION;

  updateSubRegionState( elemManager, subRegion );
  return 0.0;
}
#if 0
void SinglePhaseWell::assembleSystem( real64 const time,
                                      real64 const dt,
                                      DomainPartition & domain,
                                      DofManager const & dofManager,
                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                      arrayView1d< real64 > const & localRhs )
{
  string const wellDofKey = dofManager.getKey( wellElementDofName());


  // selects constraints one of 2 ways
  //  wellEstimator flag set to 0 => orginal logic rates are computed during update state and constraints are selected every newton
  // iteration
  //  wellEstimator flag > 0 =>   well esitmator solved for each constraint and then selects the constraint
  //                         =>   estimator solve only performed first "wellEstimator" iterations
  NonlinearSolverParameters const & nonlinearParams =  getNonlinearSolverParameters();
  selectWellConstraint( time, dt, nonlinearParams.m_numNewtonIterations, domain );


  // assemble the accumulation term in the mass balance equations
  assembleAccumulationTerms( time, dt, domain, dofManager, localMatrix, localRhs );

  // then assemble the pressure relations between well elements
  //assemblePressureRelations( time, dt, domain, dofManager, localMatrix, localRhs );
  {
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 string_array const & regionNames )
    {
      ElementRegionManager & elementRegionManager = mesh.getElemManager();
      elementRegionManager.forElementRegions< WellElementRegion >( regionNames,
                                                                   [&]( localIndex const,
                                                                        WellElementRegion & region )
      {
        WellElementSubRegion & subRegion = region.getGroup( ElementRegionBase::viewKeyStruct::elementSubRegions() )
                                             .getGroup< WellElementSubRegion >( region.getSubRegionName() );

        assembleWellConstraintTerms( time, dt, subRegion, dofManager, localMatrix.toViewConstSizes(), localRhs );
      } );
    } );
  }
  // then compute the perforation rates (later assembled by the coupled solver)
  computePerforationRates( time, dt, domain );

  // then assemble the flux terms in the mass balance equations
  // get a reference to the degree-of-freedom numbers
  // then assemble the flux terms in the mass balance equations
  assembleFluxTerms( time, dt, domain, dofManager, localMatrix, localRhs );

  // then apply a special treatment to the wells that are shut
  shutDownWell( time, domain, dofManager, localMatrix, localRhs );
}
#endif

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
    // tjb wellControls.forSubGroups< BHPConstraint, MassConstraint, VolumeRateConstraint >( [&]( auto & constraint )
    forSubGroups< BHPConstraint, InjectionConstraint< VolumeRateConstraint >, ProductionConstraint< VolumeRateConstraint > >( [&]( auto & constraint )
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
}



void SinglePhaseWell::computeWellPerforationRates( real64 const & time_n,
                                                   real64 const & GEOS_UNUSED_PARAM( dt ),
                                                   ElementRegionManager const & elemManager,
                                                   WellElementSubRegion & subRegion )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( time_n );
  // get the well data
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
    // Zero completion flow rate
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
  integer numNorm = 1;     // mass balance
  array1d< real64 > localResidualNorm;
  array1d< real64 > localResidualNormalizer;

  if( isThermal() )
  {
    numNorm = 2;     // mass balance and energy balance
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
                                               real64 const scalingFactor )
{
  GEOS_MARK_FUNCTION;

  string const wellDofKey = dofManager.getKey( wellElementDofName() );
  integer numNegativePressures = 0;
  real64 minPressure = 0.0;


  globalIndex const rankOffset = dofManager.rankOffset();
  // get the degree of freedom numbers on well elements
  arrayView1d< globalIndex const > const & dofNumber =
    subRegion.getReference< array1d< globalIndex > >( wellDofKey );
  arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();

  // get a reference to the primary variables on well elements
  arrayView1d< real64 const > const & pres =
    subRegion.getField< well::pressure >();

  auto const statistics =
    singlePhaseBaseKernels::SolutionCheckKernel::
      launch< parallelDevicePolicy<> >( localSolution, rankOffset, dofNumber, ghostRank, pres, scalingFactor );

  numNegativePressures += statistics.first;
  minPressure = std::min( minPressure, statistics.second );


  numNegativePressures = MpiWrapper::sum( numNegativePressures );

  if( numNegativePressures > 0 )
  {
    GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                           GEOS_FMT( "        {}: Number of negative pressure values: {}, minimum value: {} Pa",
                                     getName(), numNegativePressures, fmt::format( "{:.{}f}", minPressure, 3 ) ) );
  }

  return (m_allowNegativePressure || numNegativePressures == 0) ?  1 : 0;
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


void SinglePhaseWell::resetStateToBeginningOfStep( ElementRegionManager const & elemManager, WellElementSubRegion & subRegion )
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

  updateSubRegionState( elemManager, subRegion );

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
                                         ElementRegionManager & elemManager,
                                         WellElementSubRegion & subRegion )
{
  GEOS_MARK_FUNCTION;
  WellControls::implicitStepSetup( time_n, dt, elemManager, subRegion );
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

  updateSubRegionState( elemManager, subRegion );

}

void SinglePhaseWell::implicitStepComplete( real64 const & time_n,
                                            real64 const & dt,
                                            WellElementSubRegion const & subRegion )
{
//  WellSolverBase::implicitStepComplete( time_n, dt, domain );

  if( getLogLevel() > 0 )
  {
    printRates( time_n, dt, subRegion );
  }
}

void SinglePhaseWell::printRates( real64 const & time_n,
                                  real64 const & dt,
                                  WellElementSubRegion const & subRegion )
{

  GEOS_UNUSED_VAR( dt );// FIX THIS tjb
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
    getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentBHPString() );
  real64 const & currentTotalVolRate =
    getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentVolRateString() );

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

bool SinglePhaseWell::evaluateConstraints( real64 const & time_n,
                                           WellElementSubRegion & subRegion )
{

  // create list of all constraints to process
  std::vector< WellConstraintBase * > constraintList;
  if( isProducer() )
  {
    constraintList =  getProdRateConstraints();
    // Solve minimum bhp constraint first
    constraintList.insert( constraintList.begin(), getMinBHPConstraint() );
  }
  else
  {
    constraintList = getInjRateConstraints();
    // Solve maximum bhp constraint first;
    constraintList.insert( constraintList.begin(), getMaxBHPConstraint() );
  }
  // Get current constraint
  WellConstraintBase *  limitingConstraint = nullptr;
  for( auto & constraint : constraintList )
  {
    if( constraint->getName() == getCurrentConstraint()->getName())
    {
      limitingConstraint =  constraint;
      // tjb. this is likely not needed. set in update state
      constraint->setBHP ( getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentBHPString() ));
      constraint->setTotalVolumeRate ( getReference< real64 >(
                                         SinglePhaseWell::viewKeyStruct::currentVolRateString() ));

      GEOS_LOG_RANK_IF ( getLogLevel() > 4 && subRegion.isLocallyOwned(),
                         " Well " << subRegion.getName() << " Limiting Constraint " << limitingConstraint->getName() << " "  << limitingConstraint->bottomHolePressure() << " "   <<
                         limitingConstraint->totalVolumeRate() );
    }
  }

  constraintList.erase( std::find( constraintList.begin(), constraintList.end(), limitingConstraint ) );

  // Check current against other constraints
  for( auto & constraint : constraintList )
  {

    if( limitingConstraint->getName() != constraint->getName())
    {
      // limitingConstraint->getName() << std::endl;
      if( constraint->checkViolation( *limitingConstraint, time_n ) )
      {
        setControl( static_cast< WellControls::Control >(constraint->getControl()) );      // tjb old
        setCurrentConstraint( constraint );
        GEOS_LOG_RANK_IF ( getLogLevel() > 4 && subRegion.isLocallyOwned(),
                           " Well " << subRegion.getName() << " New Limiting Constraint " << constraint->getName() << " "  << constraint->getConstraintValue( time_n )  );
      }
    }
  }
  GEOS_LOG_RANK_IF ( getLogLevel() > 4 && subRegion.isLocallyOwned(),
                     " Well " << subRegion.getName() << " Limiting Constraint " << limitingConstraint->getName() << " "  << limitingConstraint->bottomHolePressure() << " " <<
                     limitingConstraint->phaseVolumeRates() << " " <<
                     limitingConstraint->totalVolumeRate() << " " << limitingConstraint->massRate());

  return true;
}


}// namespace geos
