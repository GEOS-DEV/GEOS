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
 * @file ImmiscibleMultiphaseFlow.cpp
 */

#include "ImmiscibleMultiphaseFlow.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureFields.hpp"
#include "constitutive/capillaryPressure/capillaryPressureSelector.hpp"
#include "constitutive/ConstitutivePassThru.hpp"

#if defined( __INTEL_COMPILER )
#pragma GCC optimize "O0"
#endif

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

ImmiscibleMultiphaseFlow::ImmiscibleMultiphaseFlow( const string & name,
                                                          Group * const parent )
  :
  FlowSolverBase( name, parent ),
  m_numPhases( 2 ),
  m_hasCapPressure( 0 )
{
  this->registerWrapper( viewKeyStruct::inputTemperatureString(), &m_inputTemperature ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Temperature" );
}

void ImmiscibleMultiphaseFlow::postInputInitialization()
{
  FlowSolverBase::postInputInitialization();
}

void ImmiscibleMultiphaseFlow::registerDataOnMesh( Group & meshBodies )
{
  using namespace fields::flow;

  FlowSolverBase::registerDataOnMesh( meshBodies );

  DomainPartition const & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  ConstitutiveManager const & cm = domain.getConstitutiveManager();

  // 0. Find a "reference" fluid model name (at this point, models are already attached to subregions)
  forDiscretizationOnMeshTargets( meshBodies, [&]( string const &,
                                                   MeshLevel & mesh,
                                                   arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      // If at least one region has a capillary pressure model, consider it enabled for all
      string const capPresName = getConstitutiveName< CapillaryPressureBase >( subRegion );
      if( !capPresName.empty() )
      {
        m_hasCapPressure = true;
      }
    } );
  } );

  m_numDofPerCell = m_numPhases;

  // 2. Register and resize all fields as necessary
  forDiscretizationOnMeshTargets( meshBodies, [&]( string const &,
                                                   MeshLevel & mesh,
                                                   arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      if( m_hasCapPressure )
      {
        subRegion.registerWrapper< string >( viewKeyStruct::capPressureNamesString() ).
          setPlotLevel( PlotLevel::NOPLOT ).
          setRestartFlags( RestartFlags::NO_WRITE ).
          setSizedFromParent( 0 ).
          setDescription( "Name of the capillary pressure constitutive model to use" ).
          reference();

        string & capPresName = subRegion.getReference< string >( viewKeyStruct::capPressureNamesString() );
        capPresName = getConstitutiveName< CapillaryPressureBase >( subRegion );
        GEOS_THROW_IF( capPresName.empty(),
                       GEOS_FMT( "{}: Capillary pressure model not found on subregion {}",
                                 getDataContext(), subRegion.getDataContext() ),
                       InputError );
      }

    
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );

      // The resizing of the arrays needs to happen here, before the call to initializePreSubGroups,
      // to make sure that the dimensions are properly set before the timeHistoryOutput starts its initialization.
      subRegion.registerField< phaseVolumeFraction >( getName() ).
        setDimLabels( 1, fluid.phaseNames() ).
        reference().resizeDimension< 1 >( m_numPhases );

      subRegion.registerField< phaseMobility >( getName() ).
        setDimLabels( 1, fluid.phaseNames() ).
        reference().resizeDimension< 1 >( m_numPhases );

      subRegion.registerField< dPhaseMobility >( getName() ).
        reference().resizeDimension< 1, 2 >( m_numPhases, m_numPhases ); // dP, dT, dC

      subRegion.registerField< phaseVolumeFraction_n >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );
    } );

    FaceManager & faceManager = mesh.getFaceManager();
    {
      faceManager.registerField< totalMassFlux >( getName() );
      // TODO: add conditional registration later, this is only needed when there is a face-based Dirichlet BC
      faceManager.registerField< facePressure >( getName() );
      faceManager.registerField< faceTemperature >( getName() );
    }
  } );
}

void ImmiscibleMultiphaseFlow::setConstitutiveNames( ElementSubRegionBase & subRegion ) const
{
  string & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  fluidName = getConstitutiveName< MultiFluidBase >( subRegion );
  GEOS_THROW_IF( fluidName.empty(),
                 GEOS_FMT( "{}: Fluid model not found on subregion {}",
                           getDataContext(), subRegion.getDataContext() ),
                 InputError );

  string & relPermName = subRegion.registerWrapper< string >( viewKeyStruct::relPermNamesString() ).
                           setPlotLevel( PlotLevel::NOPLOT ).
                           setRestartFlags( RestartFlags::NO_WRITE ).
                           setSizedFromParent( 0 ).
                           setDescription( "Name of the relative permeability constitutive model to use" ).
                           reference();

  relPermName = getConstitutiveName< RelativePermeabilityBase >( subRegion );

  GEOS_THROW_IF( relPermName.empty(),
                 GEOS_FMT( "{}: Relative permeability model not found on subregion {}",
                           getDataContext(), subRegion.getDataContext() ),
                 InputError );
}

void ImmiscibleMultiphaseFlow::initializePreSubGroups()
{
  FlowSolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  ConstitutiveManager const & cm = domain.getConstitutiveManager();

  // 1. Validate various models against each other (must have same phases and components)
  validateConstitutiveModels( domain );

  // 2. Set the value of temperature
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 > const temp = subRegion.getField< fields::flow::temperature >();
      temp.setValues< parallelHostPolicy >( m_inputTemperature );
    } );
  } );
}

void ImmiscibleMultiphaseFlow::validateConstitutiveModels( DomainPartition const & domain ) const
{
  GEOS_MARK_FUNCTION;

  ConstitutiveManager const & cm = domain.getConstitutiveManager();
  MultiFluidBase const & referenceFluid = cm.getConstitutiveRelation< MultiFluidBase >( m_referenceFluidModelName );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {

      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
      compareMultiphaseModels( fluid, referenceFluid );
      compareMulticomponentModels( fluid, referenceFluid );

      string const & relpermName = subRegion.getReference< string >( viewKeyStruct::relPermNamesString() );
      RelativePermeabilityBase const & relPerm = getConstitutiveModel< RelativePermeabilityBase >( subRegion, relpermName );
      compareMultiphaseModels( relPerm, referenceFluid );

    } );
  } );
}

void ImmiscibleMultiphaseFlow::updateFluidModel( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  arrayView1d< real64 const > const pres = dataGroup.getField< fields::flow::pressure >();
  arrayView1d< real64 const > const temp = dataGroup.getField< fields::flow::temperature >();

  string const & fluidName = dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( dataGroup, fluidName );
}

void ImmiscibleMultiphaseFlow::updateRelPermModel( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
    dataGroup.getField< fields::flow::phaseVolumeFraction >();

  string const & relPermName = dataGroup.getReference< string >( viewKeyStruct::relPermNamesString() );
  RelativePermeabilityBase & relPerm = getConstitutiveModel< RelativePermeabilityBase >( dataGroup, relPermName );

  // constitutive::constitutiveUpdatePassThru( relPerm, [&] ( auto & castedRelPerm )
  // {
  //   typename TYPEOFREF( castedRelPerm ) ::KernelWrapper relPermWrapper = castedRelPerm.createKernelWrapper();

  //   isothermalImmiscibleMultiphaseFlowKernels::
  //     RelativePermeabilityUpdateKernel::
  //     launch< parallelDevicePolicy<> >( dataGroup.size(),
  //                                       relPermWrapper,
  //                                       phaseVolFrac );
  // } );
}

void ImmiscibleMultiphaseFlow::updateCapPressureModel( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  if( m_hasCapPressure )
  {
    arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
      dataGroup.getField< fields::flow::phaseVolumeFraction >();

    string const & cappresName = dataGroup.getReference< string >( viewKeyStruct::capPressureNamesString() );
    CapillaryPressureBase & capPressure = getConstitutiveModel< CapillaryPressureBase >( dataGroup, cappresName );

    constitutive::constitutiveUpdatePassThru( capPressure, [&] ( auto & castedCapPres )
    {
      typename TYPEOFREF( castedCapPres ) ::KernelWrapper capPresWrapper = castedCapPres.createKernelWrapper();

      isothermalImmiscibleMultiphaseFlowKernels::
        CapillaryPressureUpdateKernel::
        launch< parallelDevicePolicy<> >( dataGroup.size(),
                                          capPresWrapper,
                                          phaseVolFrac );
    } );
  }
}

real64 ImmiscibleMultiphaseFlow::updateFluidState( ElementSubRegionBase & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  updateFluidModel( subRegion );
  updateRelPermModel( subRegion );
  updatePhaseMobility( subRegion );
  updateCapPressureModel( subRegion );

  return maxDeltaPhaseVolFrac;
}

void ImmiscibleMultiphaseFlow::initializeFluidState( MeshLevel & mesh,
                                                        DomainPartition & domain,
                                                        arrayView1d< string const > const & regionNames )
{
  GEOS_MARK_FUNCTION;

  mesh.getElemManager().forElementSubRegions( regionNames,
                                              [&]( localIndex const,
                                                   ElementSubRegionBase & subRegion )
  {
    // 2. Assume global component fractions have been prescribed.
    // Initialize constitutive state to get fluid density.
    updateFluidModel( subRegion );

    // 3. Back-calculate global component densities from fractions and total fluid density
    // in order to initialize the primary solution variables
    string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
    MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
    arrayView2d< real64 const, multifluid::USD_FLUID > const totalDens = fluid.totalDensity();

    arrayView2d< real64 const, compflow::USD_COMP > const compFrac =
      subRegion.getField< fields::flow::globalCompFraction >();
    arrayView2d< real64, compflow::USD_COMP > const compDens =
      subRegion.getField< fields::flow::globalCompDensity >();

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      for( integer ic = 0; ic < numComp; ++ic )
      {
        compDens[ei][ic] = totalDens[ei][0] * compFrac[ei][ic];
      }
    } );
  } );

  // with initial component densities defined - check if they need to be corrected to avoid zero diags etc
  chopNegativeDensities( domain );

  // for some reason CUDA does not want the host_device lambda to be defined inside the generic lambda
  // I need the exact type of the subRegion for updateSolidflowProperties to work well.
  mesh.getElemManager().forElementSubRegions< CellElementSubRegion,
                                              SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                           auto & subRegion )
  {
    // 4. Initialize/update dependent state quantities

    // 4.1 Update the constitutive models that only depend on
    //      - the primary variables
    //      - the fluid constitutive quantities (as they have already been updated)
    // We postpone the other constitutive models for now
    // In addition, to avoid multiplying permeability/porosity bay netToGross in the assembly kernel, we do it once and for all here
    arrayView1d< real64 const > const netToGross = subRegion.template getField< fields::flow::netToGross >();
    CoupledSolidBase const & porousSolid =
      getConstitutiveModel< CoupledSolidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::solidNamesString() ) );
    PermeabilityBase const & permeabilityModel =
      getConstitutiveModel< PermeabilityBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::permeabilityNamesString() ) );
    permeabilityModel.scaleHorizontalPermeability( netToGross );
    porousSolid.scaleReferencePorosity( netToGross );
    saveConvergedState( subRegion ); // necessary for a meaningful porosity update in sequential schemes
    updatePorosityAndPermeability( subRegion );
    updateCompAmount( subRegion );
    updatePhaseVolumeFraction( subRegion );

    // Now, we initialize and update each constitutive model one by one

    // 4.2 Save the computed porosity into the old porosity
    //
    // Note:
    // - This must be called after updatePorosityAndPermeability
    // - This step depends on porosity
    string const & solidName = subRegion.template getReference< string >( viewKeyStruct::solidNamesString() );
    CoupledSolidBase const & porousMaterial = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );
    porousMaterial.initializeState();

    // 4.3 Initialize/update the relative permeability model using the initial phase volume fraction
    //     This is needed to handle relative permeability hysteresis
    //     Also, initialize the fluid model
    //
    // Note:
    // - This must be called after updatePhaseVolumeFraction
    // - This step depends on phaseVolFraction

    // initialized phase volume fraction
    arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
      subRegion.template getField< fields::flow::phaseVolumeFraction >();

    string const & relpermName = subRegion.template getReference< string >( viewKeyStruct::relPermNamesString() );
    RelativePermeabilityBase & relPermMaterial =
      getConstitutiveModel< RelativePermeabilityBase >( subRegion, relpermName );
    relPermMaterial.saveConvergedPhaseVolFractionState( phaseVolFrac ); // this needs to happen before calling updateRelPermModel
    updateRelPermModel( subRegion );
    relPermMaterial.saveConvergedState(); // this needs to happen after calling updateRelPermModel

    string const & fluidName = subRegion.template getReference< string >( viewKeyStruct::fluidNamesString() );
    MultiFluidBase & fluidMaterial = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
    fluidMaterial.initializeState();

    // 4.4 Then, we initialize/update the capillary pressure model
    //
    // Note:
    // - This must be called after updatePorosityAndPermeability
    // - This step depends on porosity and permeability
    if( m_hasCapPressure )
    {
      // initialized porosity
      arrayView2d< real64 const > const porosity = porousMaterial.getPorosity();

      string const & permName = subRegion.template getReference< string >( viewKeyStruct::permeabilityNamesString() );
      PermeabilityBase const & permeabilityMaterial =
        getConstitutiveModel< PermeabilityBase >( subRegion, permName );
      // initialized permeability
      arrayView3d< real64 const > const permeability = permeabilityMaterial.permeability();

      string const & capPressureName = subRegion.template getReference< string >( viewKeyStruct::capPressureNamesString() );
      CapillaryPressureBase const & capPressureMaterial =
        getConstitutiveModel< CapillaryPressureBase >( subRegion, capPressureName );
      capPressureMaterial.initializeRockState( porosity, permeability ); // this needs to happen before calling updateCapPressureModel
      updateCapPressureModel( subRegion );
    }

    // 4.5 Update the phase mobility
    //
    // Note:
    // - This must be called after updateRelPermModel
    // - This step depends phaseRelPerm
    updatePhaseMobility( subRegion );

  } );

  // 5. Save initial pressure
  mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
    arrayView1d< real64 > const initPres = subRegion.getField< fields::flow::initialPressure >();
    arrayView1d< real64 const > const temp = subRegion.getField< fields::flow::temperature >();
    arrayView1d< real64 > const initTemp = subRegion.template getField< fields::flow::initialTemperature >();
    initPres.setValues< parallelDevicePolicy<> >( pres );
    initTemp.setValues< parallelDevicePolicy<> >( temp );
  } );
}

void ImmiscibleMultiphaseFlow::initializePostInitialConditionsPreSubGroups()
{
  GEOS_MARK_FUNCTION;

  FlowSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( { fields::flow::pressure::key(),
                                       fields::flow::globalCompDensity::key() },
                                     regionNames );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), false );

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames,
                                                                                                 [&]( localIndex const,
                                                                                                      auto & subRegion )
    {
      // set mass fraction flag on fluid models
      string const & fluidName = subRegion.template getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
      fluid.setMassFlag( m_useMass );

      saveConvergedState( subRegion ); // necessary for a meaningful porosity update in sequential schemes
      updatePorosityAndPermeability( subRegion );

      CoupledSolidBase const & porousSolid =
        getConstitutiveModel< CoupledSolidBase >( subRegion,
                                                  subRegion.template getReference< string >( viewKeyStruct::solidNamesString() ) );
      porousSolid.initializeState();
    } );

    // Initialize primary variables from applied initial conditions
    initializeFluidState( mesh, domain, regionNames );

    mesh.getElemManager().forElementRegions< SurfaceElementRegion >( regionNames,
                                                                     [&]( localIndex const,
                                                                          SurfaceElementRegion & region )
    {
      region.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
      {
        subRegion.getWrapper< real64_array >( fields::flow::hydraulicAperture::key() ).
          setApplyDefaultValue( region.getDefaultAperture() );
      } );
    } );

  } );

  // report to the user if some pore volumes are very small
  // note: this function is here because: 1) porosity has been initialized and 2) NTG has been applied
  validatePoreVolumes( domain );
}

void
ImmiscibleMultiphaseFlow::implicitStepSetup( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                                real64 const & GEOS_UNUSED_PARAM( dt ),
                                                DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion,
                                                SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                auto & subRegion )
    {
      saveConvergedState( subRegion );

      // update porosity, permeability
      updatePorosityAndPermeability( subRegion );
      // update all fluid properties
      updateFluidState( subRegion );

      // after the update, save the new saturation
      arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
        subRegion.template getField< fields::flow::phaseVolumeFraction >();

      arrayView2d< real64, compflow::USD_PHASE > const phaseVolFrac_n =
        subRegion.template getField< fields::flow::phaseVolumeFraction_n >();
      phaseVolFrac_n.setValues< parallelDevicePolicy<> >( phaseVolFrac );

    } );
  } );
}

void ImmiscibleMultiphaseFlow::assembleSystem( real64 const GEOS_UNUSED_PARAM( time_n ),
                                                  real64 const dt,
                                                  DomainPartition & domain,
                                                  DofManager const & dofManager,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  assembleAccumulationTerm( domain,
                                             dofManager,
                                             localMatrix,
                                             localRhs );

 
  assembleFluxTerms( dt,
                       domain,
                       dofManager,
                       localMatrix,
                       localRhs );
}

void ImmiscibleMultiphaseFlow::assembleAccumulationTerm( DomainPartition & domain,
                                                                             DofManager const & dofManager,
                                                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                             arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );

      MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
    } );
  } );
}

void ImmiscibleMultiphaseFlow::assembleFluxTerms( real64 const dt,
                     DomainPartition const & domain,
                     DofManager const & dofManager,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs ) const
{}

void ImmiscibleMultiphaseFlow::applyBoundaryConditions( real64 const time_n,
                                                           real64 const dt,
                                                           DomainPartition & domain,
                                                           DofManager const & dofManager,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                           arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  // apply pressure boundary conditions.
  applyDirichletBC( time_n, dt, dofManager, domain, localMatrix.toViewConstSizes(), localRhs.toView() );

  // apply flux boundary conditions
  applySourceFluxBC( time_n, dt, dofManager, domain, localMatrix.toViewConstSizes(), localRhs.toView() );
}

namespace
{
char const bcLogMessage[] =
  "ImmiscibleMultiphaseFlow {}: at time {}s, "
  "the <{}> boundary condition '{}' is applied to the element set '{}' in subRegion '{}'. "
  "\nThe scale of this boundary condition is {} and multiplies the value of the provided function (if any). "
  "\nThe total number of target elements (including ghost elements) is {}. "
  "\nNote that if this number is equal to zero for all subRegions, the boundary condition will not be applied on this element set.";
}

void ImmiscibleMultiphaseFlow::applySourceFluxBC( real64 const time,
                                                     real64 const dt,
                                                     DofManager const & dofManager,
                                                     DomainPartition & domain,
                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                     arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  // Step 1: count individual source flux boundary conditions

  std::map< string, localIndex > bcNameToBcId;
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

  computeSourceFluxSizeScalingFactor( time,
                                      dt,
                                      domain,
                                      bcNameToBcId,
                                      bcAllSetsSize.toView() );

  // Step 3: we are ready to impose the boundary condition, normalized by the set size

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & )
  {
    fsManager.apply< ElementSubRegionBase,
                     SourceFluxBoundaryCondition >( time + dt,
                                                    mesh,
                                                    SourceFluxBoundaryCondition::catalogName(),
                                                    [&]( SourceFluxBoundaryCondition const & fs,
                                                         string const & setName,
                                                         SortedArrayView< localIndex const > const & targetSet,
                                                         ElementSubRegionBase & subRegion,
                                                         string const & )
    {
      if( fs.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
      {
        globalIndex const numTargetElems = MpiWrapper::sum< globalIndex >( targetSet.size() );
        GEOS_LOG_RANK_0( GEOS_FMT( bcLogMessage,
                                   getName(), time+dt, fs.getCatalogName(), fs.getName(),
                                   setName, subRegion.getName(), fs.getScale(), numTargetElems ) );
      }

      if( targetSet.size() == 0 )
      {
        return;
      }
      if( !subRegion.hasWrapper( dofKey ) )
      {
        if( fs.getLogLevel() >= 1 )
        {
          GEOS_LOG_RANK( GEOS_FMT( "{}: trying to apply SourceFlux, but its targetSet named '{}' intersects with non-simulated region named '{}'.",
                                   getDataContext(), setName, subRegion.getName() ) );
        }
        return;
      }

      arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

      // Step 3.1: get the values of the source boundary condition that need to be added to the rhs
      // We don't use FieldSpecificationBase::applyConditionToSystem here because we want to account for the row permutation used in the
      // compositional solvers

      array1d< globalIndex > dofArray( targetSet.size() );
      array1d< real64 > rhsContributionArray( targetSet.size() );
      arrayView1d< real64 > rhsContributionArrayView = rhsContributionArray.toView();
      localIndex const rankOffset = dofManager.rankOffset();

      RAJA::ReduceSum< parallelDeviceReduce, real64 > massProd( 0.0 );

      // note that the dofArray will not be used after this step (simpler to use dofNumber instead)
      fs.computeRhsContribution< FieldSpecificationAdd,
                                 parallelDevicePolicy<> >( targetSet.toViewConst(),
                                                           time + dt,
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

      integer const fluidComponentId = fs.getComponent();
      integer const numFluidComponents = m_numComponents;
      integer const useTotalMassEquation = m_useTotalMassEquation;
      forAll< parallelDevicePolicy<> >( targetSet.size(), [sizeScalingFactor,
                                                           targetSet,
                                                           rankOffset,
                                                           ghostRank,
                                                           fluidComponentId,
                                                           numFluidComponents,
                                                           useTotalMassEquation,
                                                           dofNumber,
                                                           rhsContributionArrayView,
                                                           localRhs,
                                                           massProd] GEOS_HOST_DEVICE ( localIndex const a )
      {
        // we need to filter out ghosts here, because targetSet may contain them
        localIndex const ei = targetSet[a];
        if( ghostRank[ei] >= 0 )
        {
          return;
        }

        real64 const rhsValue = rhsContributionArrayView[a] / sizeScalingFactor; // scale the contribution by the sizeScalingFactor here!
        massProd += rhsValue;
        if( useTotalMassEquation > 0 )
        {
          // for all "fluid components", we add the value to the total mass balance equation
          globalIndex const totalMassBalanceRow = dofNumber[ei] - rankOffset;
          localRhs[totalMassBalanceRow] += rhsValue;
          if( fluidComponentId < numFluidComponents - 1 )
          {
            globalIndex const compMassBalanceRow = totalMassBalanceRow + fluidComponentId + 1; // component mass bal equations are shifted
            localRhs[compMassBalanceRow] += rhsValue;
          }
        }
        else
        {
          globalIndex const compMassBalanceRow = dofNumber[ei] - rankOffset + fluidComponentId;
          localRhs[compMassBalanceRow] += rhsValue;
        }
      } );

      SourceFluxStatsAggregator::forAllFluxStatWrappers( subRegion, fs.getName(),
                                                         [&]( SourceFluxStatsAggregator::WrappedStats & wrapper )
      {
        // set the new sub-region statistics for this timestep
        array1d< real64 > massProdArr{ m_numComponents };
        massProdArr[fluidComponentId] = massProd.get();
        wrapper.gatherTimeStepStats( time, dt, massProdArr.toViewConst(), targetSet.size() );
      } );
    } );
  } );
}

bool ImmiscibleMultiphaseFlow::validateDirichletBC( DomainPartition & domain,
                                                       real64 const time ) const
{
  constexpr integer MAX_NC = MultiFluidBase::MAX_NUM_COMPONENTS;
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  bool bcConsistent = true;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & )
  {
    // map: regionName -> subRegionName -> setName -> numComps to check pressure/comp are present consistent
    map< string, map< string, map< string, ComponentMask< MAX_NC > > > > bcPresCompStatusMap;
    // map: regionName -> subRegionName -> setName check to that temperature is present/consistent
    map< string, map< string, set< string > > > bcTempStatusMap;

    // 1. Check pressure Dirichlet BCs
    fsManager.apply< ElementSubRegionBase >( time,
                                             mesh,
                                             fields::flow::pressure::key(),
                                             [&]( FieldSpecificationBase const &,
                                                  string const & setName,
                                                  SortedArrayView< localIndex const > const &,
                                                  ElementSubRegionBase & subRegion,
                                                  string const & )
    {
      // Check whether pressure has already been applied to this set
      string const & subRegionName = subRegion.getName();
      string const & regionName = subRegion.getParent().getParent().getName();

      auto & subRegionSetMap = bcPresCompStatusMap[regionName][subRegionName];
      if( subRegionSetMap.count( setName ) > 0 )
      {
        bcConsistent = false;
        GEOS_WARNING( GEOS_FMT( "Conflicting pressure boundary conditions on set {}/{}/{}", regionName, subRegionName, setName ) );
      }
      subRegionSetMap[setName].setNumComp( m_numComponents );
    } );
  } );
  return bcConsistent;
}

void ImmiscibleMultiphaseFlow::applyDirichletBC( real64 const time_n,
                                                    real64 const dt,
                                                    DofManager const & dofManager,
                                                    DomainPartition & domain,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  // Only validate BC at the beginning of Newton loop
  if( m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
  {
    bool const bcConsistent = validateDirichletBC( domain, time_n + dt );
    GEOS_ERROR_IF( !bcConsistent, GEOS_FMT( "ImmiscibleMultiphaseFlow {}: inconsistent boundary conditions", getDataContext() ) );
  }

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & )
  {

    // 1. Apply pressure Dirichlet BCs, store in a separate field
    applyFieldValue< ElementSubRegionBase >( time_n, dt, mesh, bcLogMessage,
                                             fields::flow::pressure::key(), fields::flow::bcPressure::key() );
    // 2. Apply composition BC (global component fraction) and store them for constitutive call
    applyFieldValue< ElementSubRegionBase >( time_n, dt, mesh, bcLogMessage,
                                             fields::flow::globalCompFraction::key(), fields::flow::globalCompFraction::key() );
    // 3. Apply temperature Dirichlet BCs, store in a separate field
    if( m_isThermal )
    {
      applyFieldValue< ElementSubRegionBase >( time_n, dt, mesh, bcLogMessage,
                                               fields::flow::temperature::key(), fields::flow::bcTemperature::key() );
    }

    globalIndex const rankOffset = dofManager.rankOffset();
    string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

    // 4. Call constitutive update, back-calculate target global component densities and apply to the system
    fsManager.apply< ElementSubRegionBase >( time_n + dt,
                                             mesh,
                                             fields::flow::pressure::key(),
                                             [&] ( FieldSpecificationBase const &,
                                                   string const &,
                                                   SortedArrayView< localIndex const > const & targetSet,
                                                   ElementSubRegionBase & subRegion,
                                                   string const & )
    {
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );

      // in the isothermal case, we use the reservoir temperature to enforce the boundary condition
      // in the thermal case, the validation function guarantees that temperature has been provided
      string const temperatureKey = m_isThermal ? fields::flow::bcTemperature::key() : fields::flow::temperature::key();

      arrayView1d< real64 const > const bcPres =
        subRegion.getReference< array1d< real64 > >( fields::flow::bcPressure::key() );
      arrayView1d< real64 const > const bcTemp =
        subRegion.getReference< array1d< real64 > >( temperatureKey );
      arrayView2d< real64 const, compflow::USD_COMP > const compFrac =
        subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( fields::flow::globalCompFraction::key() );

      constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
      {
        using FluidType = TYPEOFREF( castedFluid );
        using ExecPolicy = typename FluidType::exec_policy;
        typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

        thermalImmiscibleMultiphaseFlowKernels::
          FluidUpdateKernel::
          launch< ExecPolicy >( targetSet,
                                fluidWrapper,
                                bcPres,
                                bcTemp,
                                compFrac );
      } );

      arrayView1d< integer const > const ghostRank =
        subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString() );
      arrayView1d< globalIndex const > const dofNumber =
        subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< real64 const > const pres =
        subRegion.getReference< array1d< real64 > >( fields::flow::pressure::key() );
      arrayView2d< real64 const, compflow::USD_COMP > const compDens =
        subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( fields::flow::globalCompDensity::key() );
      arrayView2d< real64 const, multifluid::USD_FLUID > const totalDens = fluid.totalDensity();

      integer const numComp = m_numComponents;
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

        // 4.2. For each component, apply target global density value
        for( integer ic = 0; ic < numComp; ++ic )
        {
          FieldSpecificationEqual::SpecifyFieldValue( dofIndex + ic + 1,
                                                      rankOffset,
                                                      localMatrix,
                                                      rhsValue,
                                                      totalDens[ei][0] * compFrac[ei][ic],
                                                      compDens[ei][ic] );
          localRhs[localRow + ic + 1] = rhsValue;
        }
      } );
    } );

    // 5. Apply temperature to the system
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

        integer const numComp = m_numComponents;
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
          FieldSpecificationEqual::SpecifyFieldValue( dofIndex + numComp + 1,
                                                      rankOffset,
                                                      localMatrix,
                                                      rhsValue,
                                                      bcTemp[ei],
                                                      temp[ei] );
          localRhs[localRow + numComp + 1] = rhsValue;
        } );
      } );
    }
  } );
}



void ImmiscibleMultiphaseFlow::resetStateToBeginningOfStep( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion,
                                                SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                auto & subRegion )
    {
      arrayView1d< real64 > const & pres =
        subRegion.template getField< fields::flow::pressure >();
      arrayView1d< real64 const > const & pres_n =
        subRegion.template getField< fields::flow::pressure_n >();
      pres.setValues< parallelDevicePolicy<> >( pres_n );

      arrayView2d< real64, compflow::USD_COMP > const & compDens =
        subRegion.template getField< fields::flow::globalCompDensity >();
      arrayView2d< real64 const, compflow::USD_COMP > const & compDens_n =
        subRegion.template getField< fields::flow::globalCompDensity_n >();
      compDens.setValues< parallelDevicePolicy<> >( compDens_n );

      if( m_isThermal )
      {
        arrayView1d< real64 > const & temp =
          subRegion.template getField< fields::flow::temperature >();
        arrayView1d< real64 const > const & temp_n =
          subRegion.template getField< fields::flow::temperature_n >();
        temp.setValues< parallelDevicePolicy<> >( temp_n );
      }

      // update porosity, permeability
      updatePorosityAndPermeability( subRegion );
      // update all fluid properties
      updateFluidState( subRegion );
    } );
  } );
}

void ImmiscibleMultiphaseFlow::implicitStepComplete( real64 const & time,
                                                        real64 const & dt,
                                                        DomainPartition & domain )
{
  // Step 1: save the converged aquifer state
  // note: we have to save the aquifer state **before** updating the pressure,
  // otherwise the aquifer flux is saved with the wrong pressure time level
  saveAquiferConvergedState( time, dt, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      // update deltaPressure
      arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
      arrayView1d< real64 const > const initPres = subRegion.getField< fields::flow::initialPressure >();
      arrayView1d< real64 > const deltaPres = subRegion.getField< fields::flow::deltaPressure >();
      isothermalImmiscibleMultiphaseFlowKernels::StatisticsKernel::
        saveDeltaPressure< parallelDevicePolicy<> >( subRegion.size(), pres, initPres, deltaPres );

      // Step 2: save the converged fluid state
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluidMaterial = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
      fluidMaterial.saveConvergedState();

      // Step 3: save the converged solid state
      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
      CoupledSolidBase const & porousMaterial = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );
      if( m_keepFlowVariablesConstantDuringInitStep )
      {
        porousMaterial.ignoreConvergedState(); // newPorosity <- porosity_n
      }
      else
      {
        porousMaterial.saveConvergedState(); // porosity_n <- porosity
      }

      // Step 4: save converged state for the relperm model to handle hysteresis
      arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
        subRegion.getField< fields::flow::phaseVolumeFraction >();
      string const & relPermName = subRegion.getReference< string >( viewKeyStruct::relPermNamesString() );
      RelativePermeabilityBase const & relPermMaterial =
        getConstitutiveModel< RelativePermeabilityBase >( subRegion, relPermName );
      relPermMaterial.saveConvergedPhaseVolFractionState( phaseVolFrac );

      // Step 5: if capillary pressure is supported, send the converged porosity and permeability to the capillary pressure model
      // note: this is needed when the capillary pressure depends on porosity and permeability (Leverett J-function for instance)
      if( m_hasCapPressure )
      {
        arrayView2d< real64 const > const porosity = porousMaterial.getPorosity();

        string const & permName = subRegion.getReference< string >( viewKeyStruct::permeabilityNamesString() );
        PermeabilityBase const & permeabilityMaterial =
          getConstitutiveModel< PermeabilityBase >( subRegion, permName );
        arrayView3d< real64 const > const permeability = permeabilityMaterial.permeability();

        string const & capPressName = subRegion.getReference< string >( viewKeyStruct::capPressureNamesString() );
        CapillaryPressureBase const & capPressureMaterial =
          getConstitutiveModel< CapillaryPressureBase >( subRegion, capPressName );
        capPressureMaterial.saveConvergedRockState( porosity, permeability );
      }
    } );
  } );
}

void ImmiscibleMultiphaseFlow::saveConvergedState( ElementSubRegionBase & subRegion ) const
{
  FlowSolverBase::saveConvergedState( subRegion );

  arrayView2d< real64 const, compflow::USD_COMP > const & compDens =
    subRegion.template getField< fields::flow::globalCompDensity >();
  arrayView2d< real64, compflow::USD_COMP > const & compDens_n =
    subRegion.template getField< fields::flow::globalCompDensity_n >();
  compDens_n.setValues< parallelDevicePolicy<> >( compDens );

  arrayView2d< real64 const, compflow::USD_COMP > const & compAmount =
    subRegion.template getField< fields::flow::compAmount >();
  arrayView2d< real64, compflow::USD_COMP > const & compAmount_n =
    subRegion.template getField< fields::flow::compAmount_n >();
  compAmount_n.setValues< parallelDevicePolicy<> >( compAmount );

  if( m_isFixedStressPoromechanicsUpdate )
  {
    arrayView2d< real64, compflow::USD_COMP > const & compDens_k =
      subRegion.template getField< fields::flow::globalCompDensity_k >();
    compDens_k.setValues< parallelDevicePolicy<> >( compDens );
  }
}

void ImmiscibleMultiphaseFlow::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  real64 maxDeltaPhaseVolFrac = 0.0;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion,
                                                SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                             auto & subRegion )
    {
      // update porosity, permeability, and solid internal energy
      updatePorosityAndPermeability( subRegion );
      // update all fluid properties
      real64 const deltaPhaseVolFrac = updateFluidState( subRegion );
      maxDeltaPhaseVolFrac = LvArray::math::max( maxDeltaPhaseVolFrac, deltaPhaseVolFrac );
      // for thermal, update solid internal energy
      if( m_isThermal )
      {
        updateSolidInternalEnergyModel( subRegion );
        updateEnergy( subRegion );
      }
    } );
  } );

  maxDeltaPhaseVolFrac = MpiWrapper::max( maxDeltaPhaseVolFrac );

  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "        {}: Max phase volume fraction change = {}", getName(), fmt::format( "{:.{}f}", maxDeltaPhaseVolFrac, 4 ) ) );
}

real64 ImmiscibleMultiphaseFlow::setNextDt( const geos::real64 & currentDt, geos::DomainPartition & domain )
{
  return SolverBase::setNextDt( currentDt, domain );
}

REGISTER_CATALOG_ENTRY( SolverBase, ImmiscibleMultiphaseFlow, string const &, Group * const )

} // namespace geos
