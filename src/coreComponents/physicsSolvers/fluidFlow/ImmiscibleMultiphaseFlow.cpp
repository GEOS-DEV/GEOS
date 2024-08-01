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

#include "FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/ImmiscibleMultiphaseFlowFields.hpp"
#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseFVMKernels.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureFields.hpp"
#include "constitutive/capillaryPressure/capillaryPressureSelector.hpp"
#include "constitutive/relativePermeability/RelativePermeabilitySelector.hpp"

#include "fieldSpecification/SourceFluxBoundaryCondition.hpp"


#include "constitutive/ConstitutivePassThru.hpp"

#if defined( __INTEL_COMPILER )
#pragma GCC optimize "O0"
#endif

namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace fields::immiscibleMultiphaseFlow;

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
  FlowSolverBase::registerDataOnMesh( meshBodies );

  // DomainPartition const & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  // ConstitutiveManager const & cm = domain.getConstitutiveManager();

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

      // The resizing of the arrays needs to happen here, before the call to initializePreSubGroups,
      // to make sure that the dimensions are properly set before the timeHistoryOutput starts its initialization.
      subRegion.registerField< phaseVolumeFraction >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );

      subRegion.registerField< phaseVolumeFraction_n >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );

      subRegion.registerField< phaseMobility >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );

      subRegion.registerField< dPhaseMobility >( getName() ).
        reference().resizeDimension< 1, 2 >( m_numPhases, m_numPhases ); // dP, dS

    } );

    // FaceManager & faceManager = mesh.getFaceManager();
    // {
    //   faceManager.registerField< totalMassFlux >( getName() );
    // }
  } );
}

void ImmiscibleMultiphaseFlow::setConstitutiveNames( ElementSubRegionBase & subRegion ) const
{
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
  // ConstitutiveManager const & cm = domain.getConstitutiveManager();

  // // 1. Validate various models against each other (must have same phases and components)
  // validateConstitutiveModels( domain );

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

void ImmiscibleMultiphaseFlow::updateFluidModel( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  // arrayView1d< real64 const > const pres = dataGroup.getField< fields::flow::pressure >();
  // arrayView1d< real64 const > const temp = dataGroup.getField< fields::flow::temperature >();
  GEOS_UNUSED_VAR( dataGroup );
}

void ImmiscibleMultiphaseFlow::updateRelPermModel( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( dataGroup );

  arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac =
    dataGroup.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();

  string const & relPermName = dataGroup.getReference< string >( viewKeyStruct::relPermNamesString() );
  RelativePermeabilityBase & relPerm = getConstitutiveModel< RelativePermeabilityBase >( dataGroup, relPermName );

  constitutive::constitutiveUpdatePassThru( relPerm, [&] ( auto & castedRelPerm )
  {
    typename TYPEOFREF( castedRelPerm ) ::KernelWrapper relPermWrapper = castedRelPerm.createKernelWrapper();

    isothermalCompositionalMultiphaseBaseKernels::
      RelativePermeabilityUpdateKernel::
      launch< parallelDevicePolicy<> >( dataGroup.size(),
                                        relPermWrapper,
                                        phaseVolFrac );
  } );
}

void ImmiscibleMultiphaseFlow::updateCapPressureModel( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  if( m_hasCapPressure )
  {
    arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac =
      dataGroup.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();

    string const & cappresName = dataGroup.getReference< string >( viewKeyStruct::capPressureNamesString() );
    CapillaryPressureBase & capPressure = getConstitutiveModel< CapillaryPressureBase >( dataGroup, cappresName );

    constitutive::constitutiveUpdatePassThru( capPressure, [&] ( auto & castedCapPres )
    {
      typename TYPEOFREF( castedCapPres ) ::KernelWrapper capPresWrapper = castedCapPres.createKernelWrapper();

      isothermalCompositionalMultiphaseBaseKernels::
        CapillaryPressureUpdateKernel::
        launch< parallelDevicePolicy<> >( dataGroup.size(),
                                          capPresWrapper,
                                          phaseVolFrac );
    } );
  }
}

void ImmiscibleMultiphaseFlow::updateFluidState( ElementSubRegionBase & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  updateFluidModel( subRegion );
  updateRelPermModel( subRegion );
  updatePhaseMobility( subRegion );
  updateCapPressureModel( subRegion );
}

void ImmiscibleMultiphaseFlow::updatePhaseMobility( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;
  
  GEOS_UNUSED_VAR( dataGroup );

  /// Matteo: looks like you will to create a new update function for the mobility. I have left the code as an example.
  // // note that the phase mobility computed here also includes phase density
  // string const & fluidName = dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() );
  // MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( dataGroup, fluidName );

  // string const & relpermName = dataGroup.getReference< string >( viewKeyStruct::relPermNamesString() );
  // RelativePermeabilityBase const & relperm = getConstitutiveModel< RelativePermeabilityBase >( dataGroup, relpermName );

  // isothermalCompositionalMultiphaseFVMKernels::
  //   PhaseMobilityKernelFactory::
  //   createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
  //                                              m_numPhases,
  //                                              dataGroup,
  //                                              fluid,
  //                                              relperm );
}

void ImmiscibleMultiphaseFlow::initializeFluidState( MeshLevel & mesh,
                                                     DomainPartition & GEOS_UNUSED_PARAM( domain ),
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

  } );

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
    arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac =
      subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();

    string const & relpermName = subRegion.template getReference< string >( viewKeyStruct::relPermNamesString() );
    RelativePermeabilityBase & relPermMaterial =
      getConstitutiveModel< RelativePermeabilityBase >( subRegion, relpermName );
    relPermMaterial.saveConvergedPhaseVolFractionState( phaseVolFrac ); // this needs to happen before calling updateRelPermModel
    updateRelPermModel( subRegion );
    relPermMaterial.saveConvergedState(); // this needs to happen after calling updateRelPermModel

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
      PermeabilityBase const & permeabilityMaterial = getConstitutiveModel< PermeabilityBase >( subRegion, permName );
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
                                       fields::immiscibleMultiphaseFlow::phaseVolumeFraction::key() },
                                     regionNames );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), false );

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames,
                                                                                                 [&]( localIndex const,
                                                                                                      auto & subRegion )
    {
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
      arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac =
        subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();

      arrayView2d< real64, immiscibleFlow::USD_PHASE > const phaseVolFrac_n =
        subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction_n >();
      phaseVolFrac_n.setValues< parallelDevicePolicy<> >( phaseVolFrac );

      arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const & phaseMass =
        subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseMass >();

      arrayView2d< real64, immiscibleFlow::USD_PHASE > const & phaseMass_n =
        subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseMass_n >();
      phaseMass_n.setValues< parallelDevicePolicy<> >( phaseMass );

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


  GEOS_UNUSED_VAR( domain, dofManager, localMatrix, localRhs );

  // forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
  //                                                              MeshLevel const & mesh,
  //                                                              arrayView1d< string const > const & regionNames )
  // {
  //   mesh.getElemManager().forElementSubRegions( regionNames,
  //                                               [&]( localIndex const,
  //                                                    ElementSubRegionBase const & subRegion )
  //   {
  //     string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  //     string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
  //   } );
  // } );
}

void ImmiscibleMultiphaseFlow::assembleFluxTerms( real64 const dt,
                                                  DomainPartition const & domain,
                                                  DofManager const & dofManager,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs ) const
{
  GEOS_UNUSED_VAR( dt, domain, dofManager, localMatrix, localRhs );
}

// Ryan: Looks like this will need to be overwritten as well...
// I have left the CompositionalMultiphaseFVM implementation for reference
void ImmiscibleMultiphaseFlow::setupDofs( DomainPartition const & domain,
                                            DofManager & dofManager ) const
{
  GEOS_UNUSED_VAR(domain, dofManager);
  //// add a field for the cell-centered degrees of freedom
  //dofManager.addField( viewKeyStruct::elemDofFieldString(),
  //                     FieldLocation::Elem,
  //                     m_numDofPerCell,
  //                     getMeshTargets() );

  //// this call with instruct GEOS to reorder the dof numbers
  //dofManager.setLocalReorderingType( viewKeyStruct::elemDofFieldString(),
  //                                   DofManager::LocalReorderingType::ReverseCutHillMcKee );

  //// for the volume balance equation, disable global coupling
  //// this equation is purely local (not coupled to neighbors or other physics)
  //dofManager.disableGlobalCouplingForEquation( viewKeyStruct::elemDofFieldString(),
  //                                             m_numComponents );


  //NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  //FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  //FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
  //dofManager.addCoupling( viewKeyStruct::elemDofFieldString(), fluxApprox );
}

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
}

void ImmiscibleMultiphaseFlow::applyDirichletBC( real64 const time,
                                                 real64 const dt,
                                                 DofManager const & dofManager,
                                                 DomainPartition & domain,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs ) const
{
  GEOS_UNUSED_VAR( time, dt, dofManager, domain, localMatrix, localRhs );
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

      // after the update, save the new saturation
      arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac_n =
        subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction_n >();

      arrayView2d< real64, immiscibleFlow::USD_PHASE > const phaseVolFrac =
        subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();
      phaseVolFrac.setValues< parallelDevicePolicy<> >( phaseVolFrac );

      arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const & phaseMass_n =
        subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseMass_n >();

      arrayView2d< real64, immiscibleFlow::USD_PHASE > const & phaseMass =
        subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseMass >();
      phaseMass.setValues< parallelDevicePolicy<> >( phaseMass_n );

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
      arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac =
        subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();
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

  arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac =
    subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();
  arrayView2d< real64, immiscibleFlow::USD_PHASE > const phaseVolFrac_n =
    subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction_n >();
  phaseVolFrac_n.setValues< parallelDevicePolicy<> >( phaseVolFrac );



}

void ImmiscibleMultiphaseFlow::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

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
      updateFluidState( subRegion );
    } );
  } );
}

real64 ImmiscibleMultiphaseFlow::setNextDt( const geos::real64 & currentDt, geos::DomainPartition & domain )
{
  return SolverBase::setNextDt( currentDt, domain );
}

REGISTER_CATALOG_ENTRY( SolverBase, ImmiscibleMultiphaseFlow, string const &, Group * const )

} // namespace geos
