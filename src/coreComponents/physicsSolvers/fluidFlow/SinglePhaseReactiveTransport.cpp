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
 * @file SinglePhaseReactiveTransport.cpp
 */

#include "SinglePhaseReactiveTransport.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/diffusion/DiffusionFields.hpp"
#include "constitutive/diffusion/DiffusionSelector.hpp"
#include "constitutive/fluid/reactivefluid/ReactiveFluidLayouts.hpp"
#include "constitutive/fluid/reactivefluid/ReactiveSinglePhaseFluid.hpp"
#include "constitutive/fluid/reactivefluid/ReactiveSinglePhaseFluid.cpp"
#include "constitutive/fluid/reactivefluid/ReactiveFluidSelector.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/solid/ReactiveSolid.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/SourceFluxBoundaryCondition.hpp"
#include "physicsSolvers/fluidFlow/SourceFluxStatistics.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/reactive/KernelLaunchSelectors.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"

/**
 * @namespace the geos namespace that encapsulates the majority of the code
 */
namespace geos
{

using namespace dataRepository;
using namespace constitutive;

template< typename POROUSWRAPPER_TYPE >
void updatePorosityAndPermeabilityFromPressureAndReactions( POROUSWRAPPER_TYPE porousWrapper,
                                                            ElementSubRegionBase & subRegion,
                                                            arrayView1d< real64 const > const & pressure,
                                                            arrayView2d< real64 const, compflow::USD_COMP > const & kineticReactionMolarIncrements )
{
  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < porousWrapper.numGauss(); ++q )
    {
      porousWrapper.updateStateFromPressureAndReactions( k, q,
                                                         pressure[k],
                                                         kineticReactionMolarIncrements[k] );
    }
  } );
}

template< typename POROUSWRAPPER_TYPE >
void updateSurfaceAreaFromReactions( POROUSWRAPPER_TYPE porousWrapper,
                                     ElementSubRegionBase & subRegion,
                                     arrayView2d< real64 const, compflow::USD_COMP > const & initialSurfaceArea,
                                     arrayView2d< real64, compflow::USD_COMP > const & surfaceArea )
{
  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < porousWrapper.numGauss(); ++q )
    {
      porousWrapper.updateSurfaceArea( k, q,
                                       initialSurfaceArea[k],
                                       surfaceArea[k] );
    }
  } );
}

SinglePhaseReactiveTransport::SinglePhaseReactiveTransport( const string & name,
                                                            Group * const parent ):
  SinglePhaseBase( name, parent ),
  m_numPrimarySpecies( 0 ),
  m_numKineticReactions( 0 ),
  m_hasDiffusion( 0 ),
  m_isUpdateReactivePorosity( 0 ),
  m_isUpdateSurfaceArea( 0 )
{
  // To add modeling parameters we want to add here

  this->registerWrapper( viewKeyStruct::isUpdateReactivePorosityString(), &m_isUpdateReactivePorosity ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag indicating whether use the reactive porosity or not" );

  this->registerWrapper( viewKeyStruct::isUpdateSurfaceAreaString(), &m_isUpdateSurfaceArea ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag indicating whether to update the surface area or not" );

  this->registerWrapper( viewKeyStruct::immobilePrimarySpeciesIndicesString(), &m_immobilePrimarySpeciesIndices ).
    setApplyDefaultValue( { } ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Array to store the indices of immobile species. Default is {}, which indicates no immobile species." );

  addLogLevel< logInfo::BoundaryConditions >();
}

void SinglePhaseReactiveTransport::registerDataOnMesh( Group & meshBodies )
{
  using namespace fields::flow;

  SinglePhaseBase::registerDataOnMesh( meshBodies );

  DomainPartition const & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  ConstitutiveManager const & cm = domain.getConstitutiveManager();

  // 0. Find a reactive fluid model name (at this point, models are already attached to subregions)
  forDiscretizationOnMeshTargets( meshBodies, [&]( string const &,
                                                   MeshLevel & mesh,
                                                   string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      if( m_reactiveFluidModelName.empty() )
      {
        m_reactiveFluidModelName = m_isThermal? getConstitutiveName< reactivefluid::ReactiveThermalCompressibleSinglePhaseFluid >( subRegion ):
                                   getConstitutiveName< reactivefluid::ReactiveCompressibleSinglePhaseFluid >( subRegion );
      }

      // If at least one region has a diffusion model, consider it enabled for all
      string const diffusionName = getConstitutiveName< DiffusionBase >( subRegion );
      if( !diffusionName.empty() )
      {
        m_hasDiffusion = true;
      }
    } );
  } );

  // 1. Set key dimensions of the problem
  // Check needed to avoid errors when running in schema generation mode.
  if( !m_reactiveFluidModelName.empty() )
  {
    if( m_isThermal )
    {
      reactivefluid::ReactiveThermalCompressibleSinglePhaseFluid const & reactiveFluid = cm.getConstitutiveRelation< reactivefluid::ReactiveThermalCompressibleSinglePhaseFluid >(
        m_reactiveFluidModelName );
      m_numPrimarySpecies = reactiveFluid.numPrimarySpecies();
      m_numKineticReactions = reactiveFluid.numKineticReactions();
    }
    else
    {
      reactivefluid::ReactiveCompressibleSinglePhaseFluid const & reactiveFluid = cm.getConstitutiveRelation< reactivefluid::ReactiveCompressibleSinglePhaseFluid >( m_reactiveFluidModelName );
      m_numPrimarySpecies = reactiveFluid.numPrimarySpecies();
      m_numKineticReactions = reactiveFluid.numKineticReactions();
    }
  }

  // n_c components + one pressure ( + one temperature if needed )
  m_numDofPerCell = m_isThermal ? m_numPrimarySpecies + 2 : m_numPrimarySpecies + 1;

  // 2. Register and resize all fields as necessary
  forDiscretizationOnMeshTargets( meshBodies, [&]( string const &,
                                                   MeshLevel & mesh,
                                                   string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      if( m_hasDiffusion )
      {
        subRegion.registerWrapper< string >( viewKeyStruct::diffusionNamesString() ).
          setPlotLevel( PlotLevel::NOPLOT ).
          setRestartFlags( RestartFlags::NO_WRITE ).
          setSizedFromParent( 0 ).
          setDescription( "Name of the diffusion constitutive model to use" );

        string & diffusionName = subRegion.getReference< string >( viewKeyStruct::diffusionNamesString() );
        diffusionName = getConstitutiveName< DiffusionBase >( subRegion );
        GEOS_THROW_IF( diffusionName.empty(),
                       GEOS_FMT( "Diffusion model not found on subregion {}", subRegion.getName() ),
                       InputError );
      }

      subRegion.registerField< logPrimarySpeciesConcentration >( getName() ).
        reference().resizeDimension< 1 >( m_numPrimarySpecies );

      subRegion.registerField< logPrimarySpeciesConcentration_n >( getName() ).
        reference().resizeDimension< 1 >( m_numPrimarySpecies );

      subRegion.registerField< primarySpeciesAggregateMole >( getName() ).
        reference().resizeDimension< 1 >( m_numPrimarySpecies );

      subRegion.registerField< primarySpeciesAggregateMole_n >( getName() ).
        reference().resizeDimension< 1 >( m_numPrimarySpecies );

      subRegion.registerField< bcLogPrimarySpeciesConcentration >( getName() ).
        reference().resizeDimension< 1 >( m_numPrimarySpecies );

      subRegion.registerField< kineticReactionMolarIncrements >( getName() ).
        reference().resizeDimension< 1 >( m_numKineticReactions );

      subRegion.registerField< surfaceArea >( getName() ).
        reference().resizeDimension< 1 >( m_numKineticReactions );

      subRegion.registerField< initialSurfaceArea >( getName() ).
        reference().resizeDimension< 1 >( m_numKineticReactions );
    } );
  } );
}

void SinglePhaseReactiveTransport::validateConstitutiveModels( DomainPartition & domain ) const
{
  GEOS_MARK_FUNCTION;

  SinglePhaseBase::validateConstitutiveModels( domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion )
    {
      string const & porosityModelName = getConstitutiveName< PorosityBase >( subRegion );

      PorosityBase const & porosity = getConstitutiveModel< PorosityBase >( subRegion, porosityModelName );

      GEOS_THROW_IF( m_isUpdateReactivePorosity && (porosity.getCatalogName() != "ReactivePorosity"),
                     GEOS_FMT( "SinglePhaseReactiveTransport {}: the reaction porosity update option is enabled in the solver, but the porosity model {} is not for reactive porosity",
                               getDataContext(), porosity.getDataContext() ),
                     InputError );

      if( m_isUpdateReactivePorosity )
      {
        ReactivePorosity const & reactivePorosity = getConstitutiveModel< ReactivePorosity >( subRegion, porosityModelName );

        GEOS_THROW_IF_NE_MSG( reactivePorosity.numKineticReactions(), m_numKineticReactions,
                              GEOS_FMT( "Mismatch in number of kinetic reactions, check the number of components input in porosity model {}",
                                        reactivePorosity.getDataContext() ),
                              InputError );
      }
    } );
  } );
}

void SinglePhaseReactiveTransport::setupDofs( DomainPartition const & domain,
                                              DofManager & dofManager ) const
{
  // add a field for the cell-centered degrees of freedom
  dofManager.addField( viewKeyStruct::elemDofFieldString(),
                       FieldLocation::Elem,
                       m_numDofPerCell,
                       getMeshTargets() );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  dofManager.addCoupling( viewKeyStruct::elemDofFieldString(), fluxApprox );
}

void SinglePhaseReactiveTransport::resetStateToBeginningOfStep( DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                                   auto & subRegion )
    {
      arrayView2d< real64, compflow::USD_COMP > const logPrimarySpeciesConc = subRegion.template getField< fields::flow::logPrimarySpeciesConcentration >();
      arrayView2d< real64 const, compflow::USD_COMP > const logPrimarySpeciesConc_n = subRegion.template getField< fields::flow::logPrimarySpeciesConcentration_n >();
      logPrimarySpeciesConc.setValues< parallelDevicePolicy<> >( logPrimarySpeciesConc_n );
    } );
  } );

  SinglePhaseBase::resetStateToBeginningOfStep( domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                                   auto & subRegion )
    {
      updateMixedReactionSystem( subRegion );
      updateSpeciesAmount( subRegion );
    } );
  } );
}

void SinglePhaseReactiveTransport::implicitStepSetup( real64 const & time_n,
                                                      real64 const & dt,
                                                      DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  SinglePhaseBase::implicitStepSetup( time_n, dt, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                                   auto & subRegion )
    {
      updateMixedReactionSystem( subRegion );
      updateSpeciesAmount( subRegion );
    } );
  } );
}

void SinglePhaseReactiveTransport::implicitStepComplete( real64 const & time,
                                                         real64 const & dt,
                                                         DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  SinglePhaseBase::implicitStepComplete( time, dt, domain );

  // Update molar increments of kinetic reactions for porosity update
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion )
    {
      updateKineticReactionMolarIncrements( dt, subRegion );
    } );
  } );

  if( m_hasDiffusion )
  {
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 string_array const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions( regionNames,
                                                  [&]( localIndex const,
                                                       ElementSubRegionBase & subRegion )
      {
        string const & diffusionName = subRegion.getReference< string >( viewKeyStruct::diffusionNamesString() );
        DiffusionBase const & diffusionMaterial = getConstitutiveModel< DiffusionBase >( subRegion, diffusionName );
        arrayView1d< real64 const > const temperature = subRegion.template getField< fields::flow::temperature >();
        diffusionMaterial.saveConvergedTemperatureState( temperature );
      } );
    } );
  }
}

void SinglePhaseReactiveTransport::assembleSystem( real64 const GEOS_UNUSED_PARAM( time_n ),
                                                   real64 const dt,
                                                   DomainPartition & domain,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  assembleAccumulationTermsInMassBalanceAndSpeciesAmountEqs( dt,
                                                             domain,
                                                             dofManager,
                                                             localMatrix,
                                                             localRhs );
  assembleFluxTerms( dt,
                     domain,
                     dofManager,
                     localMatrix,
                     localRhs );
}

void SinglePhaseReactiveTransport::assembleAccumulationTermsInMassBalanceAndSpeciesAmountEqs( real64 const dt,
                                                                                              DomainPartition & domain,
                                                                                              DofManager const & dofManager,
                                                                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                              arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      geos::constitutive::CoupledSolidBase const & solid =
        getConstitutiveModel< geos::constitutive::CoupledSolidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::solidNamesString() ) );

      string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

      if( m_isThermal )
      {
        reactivefluid::ReactiveThermalCompressibleSinglePhaseFluid const & fluid =
          getConstitutiveModel< reactivefluid::ReactiveThermalCompressibleSinglePhaseFluid >( subRegion, subRegion.template getReference< string >( viewKeyStruct::fluidNamesString() ) );

        thermalSinglePhaseReactiveBaseKernels::
          AccumulationKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_numPrimarySpecies,
                                                     dt,
                                                     dofManager.rankOffset(),
                                                     dofKey,
                                                     subRegion,
                                                     fluid,
                                                     solid,
                                                     localMatrix,
                                                     localRhs );
      }
      else
      {
        reactivefluid::ReactiveCompressibleSinglePhaseFluid const & fluid =
          getConstitutiveModel< reactivefluid::ReactiveCompressibleSinglePhaseFluid >( subRegion, subRegion.template getReference< string >( viewKeyStruct::fluidNamesString() ) );

        singlePhaseReactiveBaseKernels::
          AccumulationKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_numPrimarySpecies,
                                                     dt,
                                                     dofManager.rankOffset(),
                                                     dofKey,
                                                     subRegion,
                                                     fluid,
                                                     solid,
                                                     localMatrix,
                                                     localRhs );
      }
    } );
  } );
}

void SinglePhaseReactiveTransport::assembleFluxTerms( real64 const dt,
                                                      DomainPartition const & domain,
                                                      DofManager const & dofManager,
                                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                      arrayView1d< real64 > const & localRhs )
{
  array1d< integer > mobilePrimarySpeciesFlags;
  mobilePrimarySpeciesFlags.resize( m_numPrimarySpecies );

  for( integer i=0; i<mobilePrimarySpeciesFlags.size(); ++i )
  {
    mobilePrimarySpeciesFlags[i] = 1;
  }

  if( m_immobilePrimarySpeciesIndices.size() > 0 )
  {
    for( integer i = 0; i < m_immobilePrimarySpeciesIndices.size(); ++i )
    {
      localIndex const immobileSpeciesIndex = m_immobilePrimarySpeciesIndices[i];
      mobilePrimarySpeciesFlags[immobileSpeciesIndex] = 0;
    }
  }

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               string_array const & )
  {
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

    string const & dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

    fluxApprox.forAllStencils( mesh, [&] ( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

      if( m_isThermal )
      {
        thermalSinglePhaseReactiveFVMKernels::
          FluxComputeKernelFactory::createAndLaunch< parallelDevicePolicy<> >( m_numPrimarySpecies,
                                                                               m_hasDiffusion,
                                                                               mobilePrimarySpeciesFlags.toViewConst(),
                                                                               dofManager.rankOffset(),
                                                                               dofKey,
                                                                               getName(),
                                                                               mesh.getElemManager(),
                                                                               stencilWrapper,
                                                                               dt,
                                                                               localMatrix.toViewConstSizes(),
                                                                               localRhs.toView() );
      }
      else
      {
        singlePhaseReactiveFVMKernels::
          FluxComputeKernelFactory::createAndLaunch< parallelDevicePolicy<> >( m_numPrimarySpecies,
                                                                               m_hasDiffusion,
                                                                               mobilePrimarySpeciesFlags.toViewConst(),
                                                                               dofManager.rankOffset(),
                                                                               dofKey,
                                                                               getName(),
                                                                               mesh.getElemManager(),
                                                                               stencilWrapper,
                                                                               dt,
                                                                               localMatrix.toViewConstSizes(),
                                                                               localRhs.toView() );
      }
    } );
  } );
}

void SinglePhaseReactiveTransport::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  SinglePhaseBase::updateState( domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                                   auto & subRegion )
    {
      updateMixedReactionSystem( subRegion );

      updateSpeciesAmount( subRegion );
    } );
  } );
}

void SinglePhaseReactiveTransport::updateSpeciesAmount( ElementSubRegionBase & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  arrayView2d< real64, compflow::USD_COMP > const primarySpeciesAggregateMole = subRegion.getField< fields::flow::primarySpeciesAggregateMole >();
  arrayView2d< real64, compflow::USD_COMP > const primarySpeciesAggregateMole_n = subRegion.getField< fields::flow::primarySpeciesAggregateMole_n >();

  CoupledSolidBase const & porousSolid =
    getConstitutiveModel< CoupledSolidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::solidNamesString() ) );
  arrayView2d< real64 const > const porosity = porousSolid.getPorosity();
  arrayView2d< real64 const > const porosity_n = porousSolid.getPorosity_n();

  arrayView1d< real64 const > const volume = subRegion.getElementVolume();
  arrayView1d< real64 > const deltaVolume = subRegion.getField< fields::flow::deltaVolume >();

  integer const numPrimarySpecies = m_numPrimarySpecies;

  if( m_isThermal )
  {
    reactivefluid::ReactiveThermalCompressibleSinglePhaseFluid & fluid =
      getConstitutiveModel< reactivefluid::ReactiveThermalCompressibleSinglePhaseFluid >( subRegion, subRegion.getReference< string >( viewKeyStruct::fluidNamesString() ) );
    arrayView3d< real64 const, reactivefluid::USD_SPECIES > const primarySpeciesAggregateConcentration = fluid.primarySpeciesAggregateConcentration();
    arrayView3d< real64 const, reactivefluid::USD_SPECIES > const primarySpeciesAggregateConcentration_n = fluid.primarySpeciesAggregateConcentration_n();

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      for( integer is = 0; is < numPrimarySpecies; ++is )
      {
        primarySpeciesAggregateMole[ei][is] = porosity[ei][0] * ( volume[ei] + deltaVolume[ei] ) * primarySpeciesAggregateConcentration[ei][0][is];

        if( isZero( primarySpeciesAggregateMole_n[ei][is] ) )
          primarySpeciesAggregateMole_n[ei][is] = porosity_n[ei][0] * volume[ei] * primarySpeciesAggregateConcentration_n[ei][0][is];
      }
    } );
  }
  else
  {
    reactivefluid::ReactiveCompressibleSinglePhaseFluid & fluid =
      getConstitutiveModel< reactivefluid::ReactiveCompressibleSinglePhaseFluid >( subRegion, subRegion.getReference< string >( viewKeyStruct::fluidNamesString() ) );
    arrayView3d< real64 const, reactivefluid::USD_SPECIES > const primarySpeciesAggregateConcentration = fluid.primarySpeciesAggregateConcentration();
    arrayView3d< real64 const, reactivefluid::USD_SPECIES > const primarySpeciesAggregateConcentration_n = fluid.primarySpeciesAggregateConcentration_n();

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      for( integer is = 0; is < numPrimarySpecies; ++is )
      {
        primarySpeciesAggregateMole[ei][is] = porosity[ei][0] * ( volume[ei] + deltaVolume[ei] ) * primarySpeciesAggregateConcentration[ei][0][is];

        if( isZero( primarySpeciesAggregateMole_n[ei][is] ) )
          primarySpeciesAggregateMole_n[ei][is] = porosity_n[ei][0] * volume[ei] * primarySpeciesAggregateConcentration_n[ei][0][is];
      }
    } );
  }
}

void SinglePhaseReactiveTransport::updateKineticReactionMolarIncrements( real64 const dt,
                                                                         ElementSubRegionBase & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  arrayView2d< real64, compflow::USD_COMP > const kineticReactionMolarIncrements = subRegion.getField< fields::flow::kineticReactionMolarIncrements >();

  integer const numKineticReactions = m_numKineticReactions;

  if( m_isThermal )
  {
    reactivefluid::ReactiveThermalCompressibleSinglePhaseFluid & fluid =
      getConstitutiveModel< reactivefluid::ReactiveThermalCompressibleSinglePhaseFluid >( subRegion, subRegion.getReference< string >( viewKeyStruct::fluidNamesString() ) );
    arrayView3d< real64 const, reactivefluid::USD_SPECIES > const kineticReactionRates = fluid.kineticReactionRates();

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      for( integer r = 0; r < numKineticReactions; ++r )
      {
        kineticReactionMolarIncrements[ei][r] = dt* kineticReactionRates[ei][0][r];
      }
    } );
  }
  else
  {
    reactivefluid::ReactiveCompressibleSinglePhaseFluid & fluid =
      getConstitutiveModel< reactivefluid::ReactiveCompressibleSinglePhaseFluid >( subRegion, subRegion.getReference< string >( viewKeyStruct::fluidNamesString() ) );
    arrayView3d< real64 const, reactivefluid::USD_SPECIES > const kineticReactionRates = fluid.kineticReactionRates();

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      for( integer r = 0; r < numKineticReactions; ++r )
      {
        kineticReactionMolarIncrements[ei][r] = dt* kineticReactionRates[ei][0][r];
      }
    } );
  }
}

void SinglePhaseReactiveTransport::updateFluidModel( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  arrayView1d< real64 const > const pres = dataGroup.getField< fields::flow::pressure >();
  arrayView1d< real64 const > const temp = dataGroup.getField< fields::flow::temperature >();
  arrayView2d< real64 const, compflow::USD_COMP > const logPrimaryConc = dataGroup.getField< fields::flow::logPrimarySpeciesConcentration >();

  if( m_isThermal )
  {
    reactivefluid::ReactiveThermalCompressibleSinglePhaseFluid & fluid =
      getConstitutiveModel< reactivefluid::ReactiveThermalCompressibleSinglePhaseFluid >( dataGroup, dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() ) );

    constitutive::constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
    {
      typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
      singlePhaseReactiveBaseKernels::FluidUpdateKernel::launch( fluidWrapper, pres, temp, logPrimaryConc );
    } );
  }
  else
  {
    reactivefluid::ReactiveCompressibleSinglePhaseFluid & fluid =
      getConstitutiveModel< reactivefluid::ReactiveCompressibleSinglePhaseFluid >( dataGroup, dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() ) );

    constitutive::constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
    {
      typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
      singlePhaseReactiveBaseKernels::FluidUpdateKernel::launch( fluidWrapper, pres, temp, logPrimaryConc );
    } );
  }
}

void SinglePhaseReactiveTransport::updatePorosityAndPermeability( CellElementSubRegion & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  if( m_isUpdateReactivePorosity )
  {
    arrayView1d< real64 const > const & pressure = subRegion.getField< fields::flow::pressure >();
    arrayView2d< real64 const, compflow::USD_COMP > const kineticReactionMolarIncrements = subRegion.getField< fields::flow::kineticReactionMolarIncrements >();

    string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
    CoupledSolidBase & porousSolid = subRegion.template getConstitutiveModel< CoupledSolidBase >( solidName );

    constitutive::ConstitutivePassThru< ReactiveSolidBase >::execute( porousSolid, [=, &subRegion] ( auto & castedPorousSolid )
    {
      typename TYPEOFREF( castedPorousSolid ) ::KernelWrapper porousWrapper = castedPorousSolid.createKernelUpdates();
      updatePorosityAndPermeabilityFromPressureAndReactions( porousWrapper, subRegion, pressure, kineticReactionMolarIncrements );
    } );
  }
  else
  {
    FlowSolverBase::updatePorosityAndPermeability( subRegion );
  }
}

void SinglePhaseReactiveTransport::updateMixedReactionSystem( ElementSubRegionBase & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  updateSurfaceArea( subRegion );

  arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
  arrayView1d< real64 const > const temp = subRegion.getField< fields::flow::temperature >();
  arrayView2d< real64 const, compflow::USD_COMP > const logPrimaryConc = subRegion.getField< fields::flow::logPrimarySpeciesConcentration >();
  arrayView2d< real64 const, compflow::USD_COMP > const surfaceArea = subRegion.getField< fields::flow::surfaceArea >();

  if( m_isThermal )
  {
    reactivefluid::ReactiveThermalCompressibleSinglePhaseFluid & fluid =
      getConstitutiveModel< reactivefluid::ReactiveThermalCompressibleSinglePhaseFluid >( subRegion, subRegion.getReference< string >( viewKeyStruct::fluidNamesString() ) );

    constitutive::constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
    {
      singlePhaseReactiveBaseKernels::MixedSystemReactionUpdateKernel::launch( castedFluid, pres, temp, logPrimaryConc, surfaceArea );
    } );
  }
  else
  {
    reactivefluid::ReactiveCompressibleSinglePhaseFluid & fluid =
      getConstitutiveModel< reactivefluid::ReactiveCompressibleSinglePhaseFluid >( subRegion, subRegion.getReference< string >( viewKeyStruct::fluidNamesString() ) );

    constitutive::constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
    {
      singlePhaseReactiveBaseKernels::MixedSystemReactionUpdateKernel::launch( castedFluid, pres, temp, logPrimaryConc, surfaceArea );
    } );
  }
}

void SinglePhaseReactiveTransport::updateSurfaceArea( ElementSubRegionBase & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  arrayView2d< real64 const, compflow::USD_COMP > const initialSurfaceArea = subRegion.getField< fields::flow::initialSurfaceArea >();
  arrayView2d< real64, compflow::USD_COMP > const surfaceArea = subRegion.getField< fields::flow::surfaceArea >();

  if( m_isUpdateReactivePorosity && m_isUpdateSurfaceArea )
  {
    string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
    CoupledSolidBase & porousSolid = subRegion.template getConstitutiveModel< CoupledSolidBase >( solidName );

    constitutive::ConstitutivePassThru< ReactiveSolidBase >::execute( porousSolid, [=, &subRegion] ( auto & castedPorousSolid )
    {
      typename TYPEOFREF( castedPorousSolid ) ::KernelWrapper porousWrapper = castedPorousSolid.createKernelUpdates();
      updateSurfaceAreaFromReactions( porousWrapper, subRegion, initialSurfaceArea, surfaceArea );
    } );
  }
  else
  {
    integer const numKineticReactions = m_numKineticReactions;

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      for( integer ir = 0; ir < numKineticReactions; ++ir )
      {
        surfaceArea[ei][ir] = initialSurfaceArea[ei][ir];
      }
    } );
  }
}

void SinglePhaseReactiveTransport::initializeFluidState( MeshLevel & mesh, string_array const & regionNames )
{
  mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                                 auto & subRegion )
  {
    string const & fluidName = subRegion.template getReference< string >( viewKeyStruct::fluidNamesString() );
    updateFluidState( subRegion );

    // 2. save the initial density (for use in the single-phase poromechanics solver to compute the deltaBodyForce)
    if( m_isThermal )
    {
      reactivefluid::ReactiveThermalCompressibleSinglePhaseFluid const & fluid =
        getConstitutiveModel< reactivefluid::ReactiveThermalCompressibleSinglePhaseFluid >( subRegion, fluidName );

      fluid.initializeState();
    }
    else
    {
      reactivefluid::ReactiveCompressibleSinglePhaseFluid const & fluid =
        getConstitutiveModel< reactivefluid::ReactiveCompressibleSinglePhaseFluid >( subRegion, fluidName );

      fluid.initializeState();
    }

    initializeEquilibriumReaction( subRegion );

    updateMixedReactionSystem( subRegion );

    SinglePhaseBase::updateMass( subRegion );
    updateSpeciesAmount( subRegion );

    saveConvergedState( subRegion );

    // If the diffusion is supported, initialize the model
    if( m_hasDiffusion )
    {
      string const & diffusionName = subRegion.template getReference< string >( viewKeyStruct::diffusionNamesString() );
      DiffusionBase const & diffusionMaterial = getConstitutiveModel< DiffusionBase >( subRegion, diffusionName );
      arrayView1d< real64 const > const temperature = subRegion.template getField< fields::flow::temperature >();
      diffusionMaterial.initializeTemperatureState( temperature );
    }
  } );
}

void SinglePhaseReactiveTransport::initializeEquilibriumReaction( ElementSubRegionBase & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
  arrayView1d< real64 const > const temp = subRegion.getField< fields::flow::temperature >();
  arrayView2d< real64, compflow::USD_COMP > const logPrimaryConc = subRegion.getField< fields::flow::logPrimarySpeciesConcentration >();

  if( m_isThermal )
  {
    reactivefluid::ReactiveThermalCompressibleSinglePhaseFluid & fluid =
      getConstitutiveModel< reactivefluid::ReactiveThermalCompressibleSinglePhaseFluid >( subRegion, subRegion.getReference< string >( viewKeyStruct::fluidNamesString() ) );

    constitutive::constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
    {
      singlePhaseReactiveBaseKernels::EquilibriumReactionUpdateKernel::launch( castedFluid, pres, temp, logPrimaryConc );
    } );

    fluid.saveConvergedState();
  }
  else
  {
    reactivefluid::ReactiveCompressibleSinglePhaseFluid & fluid =
      getConstitutiveModel< reactivefluid::ReactiveCompressibleSinglePhaseFluid >( subRegion, subRegion.getReference< string >( viewKeyStruct::fluidNamesString() ) );

    constitutive::constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
    {
      singlePhaseReactiveBaseKernels::EquilibriumReactionUpdateKernel::launch( castedFluid, pres, temp, logPrimaryConc );
    } );

    fluid.saveConvergedState();
  }
}

void SinglePhaseReactiveTransport::initializePostInitialConditionsPreSubGroups()
{
  GEOS_MARK_FUNCTION;

  FlowSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( { fields::flow::logPrimarySpeciesConcentration::key() },
                                     regionNames );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), false );
  } );

  FlowSolverBase::initializeState( domain );
}

namespace
{
char const bcLogMessage[] =
  "SinglePhaseReactiveTransport {}: at time {}s, "
  "the <{}> boundary condition '{}' is applied to the element set '{}' in subRegion '{}'. "
  "\nThe scale of this boundary condition is {} and multiplies the value of the provided function (if any). "
  "\nThe total number of target elements (including ghost elements) is {}. "
  "\nNote that if this number is equal to zero for all subRegions, the boundary condition will not be applied on this element set.";
}

void SinglePhaseReactiveTransport::applyBoundaryConditions( real64 const time_n,
                                                            real64 const dt,
                                                            DomainPartition & domain,
                                                            DofManager const & dofManager,
                                                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                            arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  // apply Dirichlet boundary conditions.
  applyDirichletBC( time_n, dt, domain, dofManager, localMatrix.toViewConstSizes(), localRhs.toView() );

  // apply flux boundary conditions
  applySourceFluxBC( time_n, dt, domain, dofManager, localMatrix.toViewConstSizes(), localRhs.toView() );

  // apply face Dirichlet boundary conditions.
  applyFaceDirichletBC( time_n, dt, dofManager, domain, localMatrix, localRhs );
}

void SinglePhaseReactiveTransport::applySourceFluxBC( real64 const time_n,
                                                      real64 const dt,
                                                      DomainPartition & domain,
                                                      DofManager const & dofManager,
                                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                      arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  // Step 1: count individual source flux boundary conditions

  stdMap< string, localIndex > bcNameToBcId;
  localIndex bcCounter = 0;

  fsManager.forSubGroups< SourceFluxBoundaryCondition >( [&] ( SourceFluxBoundaryCondition const & bc )
  {
    // collect all the bc names to idx
    bcNameToBcId[bc.getName()] = bcCounter;
    bcCounter++;
  } );

  if( bcCounter == 0 )
  {
    return;
  }

  // Step 2: count the set size for each source flux (each source flux may have multiple target sets)

  array1d< globalIndex > bcAllSetsSize( bcNameToBcId.size() );

  computeSourceFluxSizeScalingFactor( time_n,
                                      dt,
                                      domain,
                                      bcNameToBcId,
                                      bcAllSetsSize.toView() );

  // Step 3: we are ready to impose the boundary condition, normalized by the set size

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & )
  {
    fsManager.apply< ElementSubRegionBase,
                     SourceFluxBoundaryCondition >( time_n + dt,
                                                    mesh,
                                                    SourceFluxBoundaryCondition::catalogName(),
                                                    [&]( SourceFluxBoundaryCondition const & fs,
                                                         string const & setName,
                                                         SortedArrayView< localIndex const > const & targetSet,
                                                         ElementSubRegionBase & subRegion,
                                                         string const & )
    {
      if( m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
      {
        globalIndex const numTargetElems = MpiWrapper::sum< globalIndex >( targetSet.size() );
        GEOS_LOG_LEVEL_RANK_0_ON_GROUP( logInfo::BoundaryConditions,
                                        GEOS_FMT( bcLogMessage,
                                                  getName(), time_n+dt, fs.getCatalogName(), fs.getName(),
                                                  setName, subRegion.getName(), fs.getScale(), numTargetElems ),
                                        fs );
      }

      if( targetSet.size() == 0 )
      {
        return;
      }
      if( !subRegion.hasWrapper( dofKey ) )
      {
        GEOS_LOG_LEVEL_BY_RANK_ON_GROUP( logInfo::BoundaryConditions,
                                         GEOS_FMT( "{}: trying to apply {}, but its targetSet named '{}' intersects with non-simulated region named '{}'.",
                                                   getDataContext(), SourceFluxBoundaryCondition::catalogName(), setName, subRegion.getName() ),
                                         fs );
        return;
      }

      arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

      // Step 3.1: get the values of the source boundary condition that need to be added to the rhs

      array1d< globalIndex > dofArray( targetSet.size() );
      array1d< real64 > rhsContributionArray( targetSet.size() );
      arrayView1d< real64 > rhsContributionArrayView = rhsContributionArray.toView();
      localIndex const rankOffset = dofManager.rankOffset();

      RAJA::ReduceSum< parallelDeviceReduce, real64 > massProd( 0.0 );

      // note that the dofArray will not be used after this step (simpler to use dofNumber instead)
      fs.computeRhsContribution< FieldSpecificationAdd,
                                 parallelDevicePolicy<> >( targetSet.toViewConst(),
                                                           time_n + dt,
                                                           dt,
                                                           subRegion,
                                                           dofNumber,
                                                           rankOffset,
                                                           localMatrix,
                                                           dofArray.toView(),
                                                           rhsContributionArrayView,
                                                           [] GEOS_HOST_DEVICE ( localIndex const )
      {
        return 0.0;
      } );

      // Step 3.2: we are ready to add the right-hand side contributions, taking into account our equation layout

      // get the normalizer
      real64 const sizeScalingFactor = bcAllSetsSize[bcNameToBcId.at( fs.getName())];

      if( m_isThermal )
      {
        reactivefluid::ReactiveThermalCompressibleSinglePhaseFluid const & fluid =
          getConstitutiveModel< reactivefluid::ReactiveThermalCompressibleSinglePhaseFluid >( subRegion, subRegion.template getReference< string >( viewKeyStruct::fluidNamesString() ) );

        thermalSinglePhaseReactiveBaseKernels::
          SourceFluxComputeKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_numPrimarySpecies,
                                                     rankOffset,
                                                     dofNumber,
                                                     ghostRank,
                                                     targetSet,
                                                     rhsContributionArrayView,
                                                     sizeScalingFactor,
                                                     fluid,
                                                     localMatrix,
                                                     localRhs,
                                                     massProd );
      }
      else
      {
        reactivefluid::ReactiveCompressibleSinglePhaseFluid const & fluid =
          getConstitutiveModel< reactivefluid::ReactiveCompressibleSinglePhaseFluid >( subRegion, subRegion.template getReference< string >( viewKeyStruct::fluidNamesString() ) );

        singlePhaseReactiveBaseKernels::
          SourceFluxComputeKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_numPrimarySpecies,
                                                     rankOffset,
                                                     dofNumber,
                                                     ghostRank,
                                                     targetSet,
                                                     rhsContributionArrayView,
                                                     sizeScalingFactor,
                                                     fluid,
                                                     localMatrix,
                                                     localRhs,
                                                     massProd );
      }

      SourceFluxStatsAggregator::forAllFluxStatWrappers( subRegion, fs.getName(),
                                                         [&]( SourceFluxStatsAggregator::WrappedStats & wrapper )
      {
        // set the new sub-region statistics for this timestep
        array1d< real64 > massProdArr{ 1 };
        massProdArr[0] = massProd.get();
        wrapper.gatherTimeStepStats( time_n, dt, massProdArr.toViewConst(), targetSet.size() );
      } );
    } );
  } );
}

void SinglePhaseReactiveTransport::applyDirichletBC( real64 const time_n,
                                                     real64 const dt,
                                                     DomainPartition & domain,
                                                     DofManager const & dofManager,
                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                     arrayView1d< real64 > const & localRhs ) const
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & )
  {
    // 1. Apply pressure Dirichlet BCs, store in a separate field
    applyFieldValue< ElementSubRegionBase >( time_n, dt, mesh, bcLogMessage,
                                             fields::flow::pressure::key(), fields::flow::bcPressure::key() );
    // 2. Apply primary species BC (log promary species concentration)
    applyFieldValue< ElementSubRegionBase >( time_n, dt, mesh, bcLogMessage,
                                             fields::flow::logPrimarySpeciesConcentration::key(), fields::flow::bcLogPrimarySpeciesConcentration::key() );
    // 3. Apply temperature Dirichlet BCs, store in a separate field
    if( m_isThermal )
    {
      applyFieldValue< ElementSubRegionBase >( time_n, dt, mesh, bcLogMessage,
                                               fields::flow::temperature::key(), fields::flow::bcTemperature::key() );
    }

    globalIndex const rankOffset = dofManager.rankOffset();
    string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

    // 4. Apply pressure to the system
    fsManager.apply< ElementSubRegionBase >( time_n + dt,
                                             mesh,
                                             fields::flow::pressure::key(),
                                             [&] ( FieldSpecificationBase const &,
                                                   string const &,
                                                   SortedArrayView< localIndex const > const & targetSet,
                                                   ElementSubRegionBase & subRegion,
                                                   string const & )
    {
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
      arrayView1d< globalIndex const > const dofNumber =
        subRegion.getReference< array1d< globalIndex > >( dofKey );

      arrayView1d< real64 const > const bcPres =
        subRegion.getReference< array1d< real64 > >( fields::flow::bcPressure::key() );
      arrayView1d< real64 const > const pres =
        subRegion.getReference< array1d< real64 > >( fields::flow::pressure::key() );

      forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        localIndex const ei = targetSet[a];
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
                                                    bcPres[ei],
                                                    pres[ei] );
        localRhs[localRow] = rhsValue;
      } );
    } );

    // 5. Apply log primary species concentration to the system
    fsManager.apply< ElementSubRegionBase >( time_n + dt,
                                             mesh,
                                             fields::flow::logPrimarySpeciesConcentration::key(),
                                             [&] ( FieldSpecificationBase const &,
                                                   string const &,
                                                   SortedArrayView< localIndex const > const & targetSet,
                                                   ElementSubRegionBase & subRegion,
                                                   string const & )
    {
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
      arrayView1d< globalIndex const > const dofNumber =
        subRegion.getReference< array1d< globalIndex > >( dofKey );

      arrayView2d< real64 const, compflow::USD_COMP > const bcLogPrimaryConc =
        subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( fields::flow::bcLogPrimarySpeciesConcentration::key() );
      arrayView2d< real64 const, compflow::USD_COMP > const logPrimaryConc =
        subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( fields::flow::logPrimarySpeciesConcentration::key() );

      integer const numPrimarySpecies = m_numPrimarySpecies;
      integer const isThermal = m_isThermal;

      forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        localIndex const ei = targetSet[a];
        if( ghostRank[ei] >= 0 )
        {
          return;
        }

        globalIndex const dofIndex = dofNumber[ei];
        localIndex const localRow = dofIndex - rankOffset;
        real64 rhsValue;

        integer const speciesDofBeginIndex = isThermal? 2:1;

        // 5.1. For each component, apply target global density value
        for( integer is = 0; is < numPrimarySpecies; ++is )
        {
          FieldSpecificationEqual::SpecifyFieldValue( dofIndex + is + speciesDofBeginIndex,
                                                      rankOffset,
                                                      localMatrix,
                                                      rhsValue,
                                                      bcLogPrimaryConc[ei][is],
                                                      logPrimaryConc[ei][is] );
          localRhs[localRow + is + speciesDofBeginIndex] = rhsValue;
        }
      } );
    } );

    // 6. Apply temperature to the system
    if( m_isThermal )
    {
      fsManager.apply< ElementSubRegionBase >( time_n + dt,
                                               mesh,
                                               fields::flow::temperature::key(),
                                               [&] ( FieldSpecificationBase const &,
                                                     string const &,
                                                     SortedArrayView< localIndex const > const & targetSet,
                                                     ElementSubRegionBase & subRegion,
                                                     string const & )
      {
        arrayView1d< integer const > const ghostRank =
          subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString() );
        arrayView1d< globalIndex const > const dofNumber =
          subRegion.getReference< array1d< globalIndex > >( dofKey );
        arrayView1d< real64 const > const bcTemp =
          subRegion.getReference< array1d< real64 > >( fields::flow::bcTemperature::key() );
        arrayView1d< real64 const > const temp =
          subRegion.getReference< array1d< real64 > >( fields::flow::temperature::key() );

        forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
        {
          localIndex const ei = targetSet[a];
          if( ghostRank[ei] >= 0 )
          {
            return;
          }

          globalIndex const dofIndex = dofNumber[ei];
          localIndex const localRow = dofIndex - rankOffset;
          real64 rhsValue;

          // 4.2. Apply temperature value to the matrix/rhs
          FieldSpecificationEqual::SpecifyFieldValue( dofIndex + 1,
                                                      rankOffset,
                                                      localMatrix,
                                                      rhsValue,
                                                      bcTemp[ei],
                                                      temp[ei] );
          localRhs[localRow + 1] = rhsValue;
        } );
      } );
    }
  } );
}

real64 SinglePhaseReactiveTransport::calculateResidualNorm( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                                            real64 const & GEOS_UNUSED_PARAM( dt ),
                                                            DomainPartition const & domain,
                                                            DofManager const & dofManager,
                                                            arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  integer constexpr numNorm = 3; // total mass balance, energy balance, and species mass balance
  array1d< real64 > localResidualNorm;
  array1d< real64 > localResidualNormalizer;
  localResidualNorm.resize( numNorm );
  localResidualNormalizer.resize( numNorm );

  physicsSolverBaseKernels::NormType const normType = SinglePhaseBase::getNonlinearSolverParameters().normType();

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel const & mesh,
                                                                      string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      real64 subRegionResidualNorm[numNorm]{};
      real64 subRegionResidualNormalizer[numNorm]{};

      // step 1: compute the norm in the subRegion

      if( m_isThermal )
      {
        singlePhaseReactiveBaseKernels::
          ResidualNormKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( normType,
                                                     m_numPrimarySpecies,
                                                     rankOffset,
                                                     dofKey,
                                                     localRhs,
                                                     subRegion,
                                                     m_nonlinearSolverParameters.m_minNormalizer,
                                                     subRegionResidualNorm,
                                                     subRegionResidualNormalizer );
      }
      else
      {
        real64 subRegionFlowResidualNorm[2]{};
        real64 subRegionFlowResidualNormalizer[2]{};

        singlePhaseReactiveBaseKernels::
          ResidualNormKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( normType,
                                                     m_numPrimarySpecies,
                                                     rankOffset,
                                                     dofKey,
                                                     localRhs,
                                                     subRegion,
                                                     m_nonlinearSolverParameters.m_minNormalizer,
                                                     subRegionFlowResidualNorm,
                                                     subRegionFlowResidualNormalizer );
        subRegionResidualNorm[0] = subRegionFlowResidualNorm[0];
        subRegionResidualNorm[1] = subRegionFlowResidualNorm[1];
        subRegionResidualNormalizer[0] = subRegionFlowResidualNormalizer[0];
        subRegionResidualNormalizer[1] = subRegionFlowResidualNormalizer[1];
      }

      // step 2: first reduction across meshBodies/regions/subRegions

      if( normType == physicsSolverBaseKernels::NormType::Linf )
      {
        physicsSolverBaseKernels::LinfResidualNormHelper::
          updateLocalNorm< numNorm >( subRegionResidualNorm, localResidualNorm );
      }
      else
      {
        physicsSolverBaseKernels::L2ResidualNormHelper::
          updateLocalNorm< numNorm >( subRegionResidualNorm, subRegionResidualNormalizer, localResidualNorm, localResidualNormalizer );
      }
    } );
  } );

  // step 3: second reduction across MPI ranks

  real64 residualNorm = 0.0;
  array1d< real64 > globalResidualNorm;
  globalResidualNorm.resize( numNorm );
  if( m_isThermal )
  {
    if( normType == physicsSolverBaseKernels::NormType::Linf )
    {
      physicsSolverBaseKernels::LinfResidualNormHelper::
        computeGlobalNorm( localResidualNorm, globalResidualNorm );
    }
    else
    {
      physicsSolverBaseKernels::L2ResidualNormHelper::
        computeGlobalNorm( localResidualNorm, localResidualNormalizer, globalResidualNorm );
    }
    residualNorm = sqrt( globalResidualNorm[0] * globalResidualNorm[0] + globalResidualNorm[1] * globalResidualNorm[1]  + globalResidualNorm[2] * globalResidualNorm[2] );

    GEOS_LOG_LEVEL_RANK_0_NLR( logInfo::Convergence, GEOS_FMT( "        ( RtotalMass RspeciesAmount ) = ( {:4.2e} {:4.2e} )        ( Renergy ) = ( {:4.2e} )",
                                                               globalResidualNorm[0], globalResidualNorm[2], globalResidualNorm[1] ));
  }
  else
  {
    if( normType == physicsSolverBaseKernels::NormType::Linf )
    {
      physicsSolverBaseKernels::LinfResidualNormHelper::
        computeGlobalNorm( localResidualNorm, globalResidualNorm );
    }
    else
    {
      physicsSolverBaseKernels::L2ResidualNormHelper::
        computeGlobalNorm( localResidualNorm, localResidualNormalizer, globalResidualNorm );
    }
    residualNorm = sqrt( globalResidualNorm[0] * globalResidualNorm[0] + globalResidualNorm[1] * globalResidualNorm[1] );

    GEOS_LOG_LEVEL_RANK_0_NLR( logInfo::Convergence, GEOS_FMT( "        ( RtotalMass RspeciesAmount ) = ( {:4.2e} {:4.2e} )",
                                                               globalResidualNorm[0], globalResidualNorm[1] ) );
  }
  return residualNorm;
}

void SinglePhaseReactiveTransport::applySystemSolution( DofManager const & dofManager,
                                                        arrayView1d< real64 const > const & localSolution,
                                                        real64 const scalingFactor,
                                                        real64 const dt,
                                                        DomainPartition & domain )
{
  GEOS_UNUSED_VAR( dt );

  if( m_isThermal )
  {
    DofManager::CompMask pressureMask( m_numDofPerCell, 0, 1 );
    DofManager::CompMask temperatureMask( m_numDofPerCell, 1, 2 );
    DofManager::CompMask speciesMask( m_numDofPerCell, 2, m_numPrimarySpecies+2 );

    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 fields::flow::pressure::key(),
                                 scalingFactor,
                                 pressureMask );

    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 fields::flow::temperature::key(),
                                 scalingFactor,
                                 temperatureMask );

    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 fields::flow::logPrimarySpeciesConcentration::key(),
                                 scalingFactor,
                                 speciesMask );
  }
  else
  {
    DofManager::CompMask pressureMask( m_numDofPerCell, 0, 1 );
    DofManager::CompMask speciesMask( m_numDofPerCell, 1, m_numPrimarySpecies+1 );

    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 fields::flow::pressure::key(),
                                 scalingFactor,
                                 pressureMask );

    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 fields::flow::logPrimarySpeciesConcentration::key(),
                                 scalingFactor,
                                 speciesMask );
  }

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel & mesh,
                                                                      string_array const & regionNames )
  {
    stdVector< string > fields{ fields::flow::pressure::key() };

    if( m_isThermal )
    {
      fields.emplace_back( fields::flow::temperature::key() );
    }

    fields.emplace_back( fields::flow::logPrimarySpeciesConcentration::key() );

    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( fields, regionNames );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), true );
  } );
}

void SinglePhaseReactiveTransport::saveConvergedState( ElementSubRegionBase & subRegion ) const
{
  SinglePhaseBase::saveConvergedState( subRegion );

  arrayView2d< real64 const, compflow::USD_COMP > const logPrimaryConc = subRegion.template getField< fields::flow::logPrimarySpeciesConcentration >();
  arrayView2d< real64, compflow::USD_COMP > const logPrimarySpeciesConc_n = subRegion.template getField< fields::flow::logPrimarySpeciesConcentration_n >();
  logPrimarySpeciesConc_n.setValues< parallelDevicePolicy<> >( logPrimaryConc );

  arrayView2d< real64 const, compflow::USD_COMP > const primarySpeciesAggregateMole = subRegion.template getField< fields::flow::primarySpeciesAggregateMole >();
  arrayView2d< real64, compflow::USD_COMP > const primarySpeciesAggregateMole_n = subRegion.template getField< fields::flow::primarySpeciesAggregateMole_n >();
  primarySpeciesAggregateMole_n.setValues< parallelDevicePolicy<> >( primarySpeciesAggregateMole );
}

namespace
{
char const faceBcLogMessage[] =
  "SinglePhaseReactiveTransport {}: at time {}s, "
  "the <{}> boundary condition '{}' is applied to the face set '{}' in '{}'. "
  "\nThe scale of this boundary condition is {} and multiplies the value of the provided function (if any). "
  "\nThe total number of target faces (including ghost faces) is {}."
  "\nNote that if this number is equal to zero, the boundary condition will not be applied on this face set.";
}

void SinglePhaseReactiveTransport::applyFaceDirichletBC( real64 const time_n,
                                                         real64 const dt,
                                                         DofManager const & dofManager,
                                                         DomainPartition & domain,
                                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                         arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  array1d< integer > mobilePrimarySpeciesFlags;
  mobilePrimarySpeciesFlags.resize( m_numPrimarySpecies );

  for( integer i=0; i<mobilePrimarySpeciesFlags.size(); ++i )
  {
    mobilePrimarySpeciesFlags[i] = 1;
  }

  if( m_immobilePrimarySpeciesIndices.size() > 0 )
  {
    for( integer i = 0; i < m_immobilePrimarySpeciesIndices.size(); ++i )
    {
      localIndex const immobileSpeciesIndex = m_immobilePrimarySpeciesIndices[i];
      mobilePrimarySpeciesFlags[immobileSpeciesIndex] = 0;
    }
  }

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const & dofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel & mesh,
                                                                      string_array const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    if( m_isThermal )
    {
      // Take BCs defined for "pressure" field and apply values to "facePressure"
      applyFieldValue< FaceManager >( time_n, dt, mesh, faceBcLogMessage,
                                      fields::flow::pressure::key(), fields::flow::facePressure::key() );

      // Take BCs defined for "temperature" field and apply values to "faceTemperature"
      applyFieldValue< FaceManager >( time_n, dt, mesh, faceBcLogMessage,
                                      fields::flow::temperature::key(), fields::flow::faceTemperature::key() );

      // Then launch the face Dirichlet kernel
      fsManager.apply< FaceManager >( time_n + dt,
                                      mesh,
                                      fields::flow::pressure::key(), // we have required that pressure is always present
                                      [&] ( FieldSpecificationBase const &,
                                            string const & setName,
                                            SortedArrayView< localIndex const > const &,
                                            FaceManager &,
                                            string const & )
      {
        BoundaryStencil const & stencil = fluxApprox.getStencil< BoundaryStencil >( mesh, setName );
        if( stencil.size() == 0 )
        {
          return;
        }

        // TODO: same issue as in the single-phase case
        //       currently we just use model from the first cell in this stencil
        //       since it's not clear how to create fluid kernel wrappers for arbitrary models.
        //       Can we just use cell properties for an approximate flux computation?
        //       Then we can forget about capturing the fluid model.
        localIndex const er = stencil.getElementRegionIndices()( 0, 0 );
        localIndex const esr = stencil.getElementSubRegionIndices()( 0, 0 );
        ElementSubRegionBase & subRegion = elemManager.getRegion( er ).getSubRegion( esr );
        string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
        reactivefluid::ReactiveThermalCompressibleSinglePhaseFluid & reactiveFluid = subRegion.getConstitutiveModel< reactivefluid::ReactiveThermalCompressibleSinglePhaseFluid >( fluidName );

        BoundaryStencilWrapper const stencilWrapper = stencil.createKernelWrapper();

        thermalSinglePhaseReactiveFVMKernels::
          DirichletFluxComputeKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_numPrimarySpecies,
                                                     mobilePrimarySpeciesFlags,
                                                     dofManager.rankOffset(),
                                                     dofKey,
                                                     this->getName(),
                                                     faceManager,
                                                     elemManager,
                                                     stencilWrapper,
                                                     reactiveFluid,
                                                     dt,
                                                     localMatrix,
                                                     localRhs );
      } );
    }
    else
    {
      // Take BCs defined for "pressure" field and apply values to "facePressure"
      applyFieldValue< FaceManager >( time_n, dt, mesh, faceBcLogMessage,
                                      fields::flow::pressure::key(), fields::flow::facePressure::key() );

      // Then launch the face Dirichlet kernel
      fsManager.apply< FaceManager >( time_n + dt,
                                      mesh,
                                      fields::flow::pressure::key(), // we have required that pressure is always present
                                      [&] ( FieldSpecificationBase const &,
                                            string const & setName,
                                            SortedArrayView< localIndex const > const &,
                                            FaceManager &,
                                            string const & )
      {
        BoundaryStencil const & stencil = fluxApprox.getStencil< BoundaryStencil >( mesh, setName );
        if( stencil.size() == 0 )
        {
          return;
        }

        // TODO: same issue as in the single-phase case
        //       currently we just use model from the first cell in this stencil
        //       since it's not clear how to create fluid kernel wrappers for arbitrary models.
        //       Can we just use cell properties for an approximate flux computation?
        //       Then we can forget about capturing the fluid model.
        localIndex const er = stencil.getElementRegionIndices()( 0, 0 );
        localIndex const esr = stencil.getElementSubRegionIndices()( 0, 0 );
        ElementSubRegionBase & subRegion = elemManager.getRegion( er ).getSubRegion( esr );
        string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
        reactivefluid::ReactiveCompressibleSinglePhaseFluid & reactiveFluid = subRegion.getConstitutiveModel< reactivefluid::ReactiveCompressibleSinglePhaseFluid >( fluidName );

        BoundaryStencilWrapper const stencilWrapper = stencil.createKernelWrapper();

        singlePhaseReactiveFVMKernels::
          DirichletFluxComputeKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_numPrimarySpecies,
                                                     mobilePrimarySpeciesFlags,
                                                     dofManager.rankOffset(),
                                                     dofKey,
                                                     this->getName(),
                                                     faceManager,
                                                     elemManager,
                                                     stencilWrapper,
                                                     reactiveFluid,
                                                     dt,
                                                     localMatrix,
                                                     localRhs );
      } );
    }
  } );
}

void SinglePhaseReactiveTransport::assembleEDFMFluxTerms( real64 const GEOS_UNUSED_PARAM( time_n ),
                                                          real64 const GEOS_UNUSED_PARAM( dt ),
                                                          DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                                          DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                                          CRSMatrixView< real64, globalIndex const > const & GEOS_UNUSED_PARAM( localMatrix ),
                                                          arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ),
                                                          string const & GEOS_UNUSED_PARAM( jumpDofKey ) )
{}

void SinglePhaseReactiveTransport::applyAquiferBC( real64 const GEOS_UNUSED_PARAM( time ),
                                                   real64 const GEOS_UNUSED_PARAM( dt ),
                                                   DomainPartition & GEOS_UNUSED_PARAM( domain ),
                                                   DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                                   CRSMatrixView< real64, globalIndex const > const & GEOS_UNUSED_PARAM( localMatrix ),
                                                   arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ) ) const
{}

void SinglePhaseReactiveTransport::assembleStabilizedFluxTerms( real64 const GEOS_UNUSED_PARAM( dt ),
                                                                DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                                                DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                                                CRSMatrixView< real64, globalIndex const > const & GEOS_UNUSED_PARAM( localMatrix ),
                                                                arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ) )
{}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SinglePhaseReactiveTransport, string const &, Group * const )

} /* namespace geos */
