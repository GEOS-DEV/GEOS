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
 * @file ProppantTransport.cpp
 */

#include "ProppantTransport.hpp"

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/fluid/SingleFluidExtrinsicData.hpp"
#include "constitutive/fluid/slurryFluidSelector.hpp"
#include "constitutive/fluid/SlurryFluidExtrinsicData.hpp"
#include "constitutive/fluid/particleFluidSelector.hpp"
#include "constitutive/fluid/ParticleFluidExtrinsicData.hpp"
#include "constitutive/permeability/PermeabilityExtrinsicData.hpp"
#include "constitutive/permeability/ProppantPermeability.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/proppantTransport/ProppantTransportExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/proppantTransport/ProppantTransportKernels.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace proppantTransportKernels;

ProppantTransport::ProppantTransport( const string & name,
                                      Group * const parent ):
  FlowSolverBase( name, parent )
{
  registerWrapper( viewKeyStruct::bridgingFactorString(), &m_bridgingFactor ).setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Bridging factor used for bridging/screen-out calculation" );

  registerWrapper( viewKeyStruct::maxProppantConcentrationString(), &m_maxProppantConcentration ).setApplyDefaultValue( 0.6 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum proppant concentration" );

  registerWrapper( viewKeyStruct::proppantDiameterString(), &m_proppantDiameter ).setApplyDefaultValue( 0.4e-3 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Proppant diameter" );

  registerWrapper( viewKeyStruct::proppantDensityString(), &m_proppantDensity ).setApplyDefaultValue( 2500.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Proppant density" );

  registerWrapper( viewKeyStruct::criticalShieldsNumberString(), &m_criticalShieldsNumber ).setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Critical Shields number" );

  registerWrapper( viewKeyStruct::frictionCoefficientString(), &m_frictionCoefficient ).setApplyDefaultValue( 0.03 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Friction coefficient" );

  registerWrapper( viewKeyStruct::updateProppantPackingString(), &m_updateProppantPacking ).setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag that enables/disables proppant-packing update" );

}

void ProppantTransport::postProcessInput()
{
  FlowSolverBase::postProcessInput();
}

void ProppantTransport::registerDataOnMesh( Group & meshBodies )
{
  using namespace extrinsicMeshData::proppant;

  FlowSolverBase::registerDataOnMesh( meshBodies );

  forMeshTargets( meshBodies, [&]( string const &,
                                   MeshLevel & mesh,
                                   arrayView1d< string const > const & regionNames )
  {

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames,
                                                                        [&]( localIndex const,
                                                                             CellElementSubRegion & subRegion )
    {
      subRegion.registerExtrinsicData< proppantConcentration >( getName() );
      subRegion.registerExtrinsicData< deltaProppantConcentration >( getName() );
      subRegion.registerExtrinsicData< componentConcentration >( getName() );
      subRegion.registerExtrinsicData< deltaComponentConcentration >( getName() );
      subRegion.registerExtrinsicData< bcComponentConcentration >( getName() );
      subRegion.registerExtrinsicData< cellBasedFlux >( getName() ).
        reference().resizeDimension< 1 >( 3 );
      subRegion.registerExtrinsicData< isProppantBoundary >( getName() );

      setConstitutiveNames( subRegion );
    } );

    mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          FaceElementSubRegion & subRegion )
    {
      subRegion.registerExtrinsicData< proppantConcentration >( getName() );
      subRegion.registerExtrinsicData< deltaProppantConcentration >( getName() );
      subRegion.registerExtrinsicData< componentConcentration >( getName() );
      subRegion.registerExtrinsicData< deltaComponentConcentration >( getName() );
      subRegion.registerExtrinsicData< bcComponentConcentration >( getName() );
      subRegion.registerExtrinsicData< oldComponentDensity >( getName() );
      subRegion.registerExtrinsicData< cellBasedFlux >( getName() ).
        reference().resizeDimension< 1 >( 3 );

      subRegion.registerExtrinsicData< isProppantBoundary >( getName() );
      subRegion.registerExtrinsicData< isProppantMobile >( getName() );
      subRegion.registerExtrinsicData< proppantPackVolumeFraction >( getName() );
      subRegion.registerExtrinsicData< proppantExcessPackVolume >( getName() );
      subRegion.registerExtrinsicData< proppantLiftFlux >( getName() );

      setConstitutiveNames( subRegion );

    } );
  } );
}


void ProppantTransport::setConstitutiveNames( ElementSubRegionBase & subRegion ) const
{
  string & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  fluidName = getConstitutiveName< SlurryFluidBase >( subRegion );
  GEOSX_THROW_IF( fluidName.empty(),
                  GEOSX_FMT( "Fluid model not found on subregion {}", subRegion.getName() ),
                  InputError );

  subRegion.registerWrapper< string >( viewKeyStruct::proppantNamesString() );
  string & proppantName = subRegion.getReference< string >( viewKeyStruct::proppantNamesString() );
  proppantName = getConstitutiveName< ParticleFluidBase >( subRegion );
  GEOSX_THROW_IF( proppantName.empty(),
                  GEOSX_FMT( "Proppant model not found on subregion {}", subRegion.getName() ),
                  InputError );

}


void ProppantTransport::initializePreSubGroups()
{
  FlowSolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  ConstitutiveManager & cm = domain.getConstitutiveManager();

  // Validate proppant models in regions
  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                                   auto & subRegion )

    {
      if( m_numDofPerCell < 1 )
      {
        SlurryFluidBase const & fluid0 = cm.getConstitutiveRelation< SlurryFluidBase >( subRegion.template getReference< string >( viewKeyStruct::fluidNamesString() ) );

        m_numComponents = fluid0.numFluidComponents();

        m_numDofPerCell = m_numComponents + 1;
      }
    } );

    if( m_numComponents > 0 )
    {
      mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                            CellElementSubRegion & subRegion )

      {
        subRegion.getExtrinsicData< extrinsicMeshData::proppant::componentConcentration >().resizeDimension< 1 >( m_numComponents );
        subRegion.getExtrinsicData< extrinsicMeshData::proppant::deltaComponentConcentration >().resizeDimension< 1 >( m_numComponents );
      } );
    }
  } );
}

void ProppantTransport::resizeFractureFields( MeshLevel & mesh, arrayView1d< string const > const & regionNames )
{
  if( m_numComponents > 0 )
  {
    mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          FaceElementSubRegion & subRegion )
    {
      subRegion.getExtrinsicData< extrinsicMeshData::proppant::componentConcentration >().resizeDimension< 1 >( m_numComponents );
      subRegion.getExtrinsicData< extrinsicMeshData::proppant::deltaComponentConcentration >().resizeDimension< 1 >( m_numComponents );
      subRegion.getExtrinsicData< extrinsicMeshData::proppant::oldComponentDensity >().resizeDimension< 1 >( m_numComponents );
      subRegion.getExtrinsicData< extrinsicMeshData::proppant::bcComponentConcentration >().resizeDimension< 1 >( m_numComponents );
    } );
  }
}

void ProppantTransport::updateFluidModel( ObjectManagerBase & dataGroup )
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const pres  = dataGroup.getExtrinsicData< extrinsicMeshData::flow::pressure >();
  arrayView1d< real64 const > const dPres = dataGroup.getExtrinsicData< extrinsicMeshData::flow::deltaPressure >();

  arrayView2d< real64 const > const componentConc  = dataGroup.getExtrinsicData< extrinsicMeshData::proppant::componentConcentration >();
  arrayView2d< real64 const > const dComponentConc = dataGroup.getExtrinsicData< extrinsicMeshData::proppant::deltaComponentConcentration >();

  SlurryFluidBase & fluid = getConstitutiveModel< SlurryFluidBase >( dataGroup, dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() ) );

  constitutive::constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    FluidUpdateKernel::launch( fluidWrapper,
                               pres,
                               dPres,
                               componentConc,
                               dComponentConc );
  } );
}

void ProppantTransport::updateComponentDensity( ObjectManagerBase & dataGroup )
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const pres  = dataGroup.getExtrinsicData< extrinsicMeshData::flow::pressure >();
  arrayView1d< real64 const > const dPres = dataGroup.getExtrinsicData< extrinsicMeshData::flow::deltaPressure >();

  arrayView2d< real64 const > const componentConc  = dataGroup.getExtrinsicData< extrinsicMeshData::proppant::componentConcentration >();
  arrayView2d< real64 const > const dComponentConc = dataGroup.getExtrinsicData< extrinsicMeshData::proppant::deltaComponentConcentration >();

  SlurryFluidBase & fluid = getConstitutiveModel< SlurryFluidBase >( dataGroup, dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() ) );

  constitutive::constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    ComponentDensityUpdateKernel::launch( fluidWrapper,
                                          pres,
                                          dPres,
                                          componentConc,
                                          dComponentConc );
  } );
}


void ProppantTransport::updateProppantModel( ObjectManagerBase & dataGroup )
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const proppantConc  = dataGroup.getExtrinsicData< extrinsicMeshData::proppant::proppantConcentration >();
  arrayView1d< real64 const > const dProppantConc = dataGroup.getExtrinsicData< extrinsicMeshData::proppant::deltaProppantConcentration >();

  SlurryFluidBase const & fluid = getConstitutiveModel< SlurryFluidBase >( dataGroup, dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() ) );

  arrayView2d< real64 const > const fluidDens            = fluid.fluidDensity();
  arrayView2d< real64 const > const dFluidDens_dPres     = fluid.dFluidDensity_dPressure();
  arrayView3d< real64 const > const dFluidDens_dCompConc = fluid.dFluidDensity_dComponentConcentration();
  arrayView2d< real64 const > const fluidVisc            = fluid.fluidViscosity();
  arrayView2d< real64 const > const dFluidVisc_dPres     = fluid.dFluidViscosity_dPressure();
  arrayView3d< real64 const > const dFluidVisc_dCompConc = fluid.dFluidViscosity_dComponentConcentration();

  ParticleFluidBase & proppant = getConstitutiveModel< ParticleFluidBase >( dataGroup, dataGroup.getReference< string >( viewKeyStruct::proppantNamesString() ) );

  constitutiveUpdatePassThru( proppant, [&]( auto & castedProppant )
  {
    typename TYPEOFREF( castedProppant ) ::KernelWrapper proppantWrapper = castedProppant.createKernelWrapper();
    ProppantUpdateKernel::launch( proppantWrapper,
                                  proppantConc,
                                  dProppantConc,
                                  fluidDens,
                                  dFluidDens_dPres,
                                  dFluidDens_dCompConc,
                                  fluidVisc,
                                  dFluidVisc_dPres,
                                  dFluidVisc_dCompConc );
  } );
}

void ProppantTransport::updateProppantMobility( ObjectManagerBase & dataGroup )
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const conc = dataGroup.getExtrinsicData< extrinsicMeshData::proppant::proppantConcentration >();
  arrayView1d< real64 const > const aperture = dataGroup.getReference< array1d< real64 > >( FaceElementSubRegion::viewKeyStruct::elementApertureString() );
  arrayView1d< integer > const isProppantMobile = dataGroup.getExtrinsicData< extrinsicMeshData::proppant::isProppantMobile >();

  real64 const minAperture = m_minAperture;
  real64 const maxProppantConcentration = m_maxProppantConcentration;

  forAll< parallelDevicePolicy<> >( dataGroup.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
  {
    isProppantMobile[a] = aperture[a] > minAperture && conc[a] < maxProppantConcentration;
  } );

}

void ProppantTransport::updateState( ObjectManagerBase & dataGroup )
{
  GEOSX_MARK_FUNCTION;

  updateFluidModel( dataGroup );
  updateProppantModel( dataGroup );
}

void ProppantTransport::initializePostInitialConditionsPreSubGroups()
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( extrinsicMeshData::proppant::proppantConcentration::key() );
  fieldNames["elems"].emplace_back( extrinsicMeshData::proppant::componentConcentration::key() );

  integer const numComponents = m_numComponents;

  // We have to redo the below loop after fractures are generated
  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )
  {

    CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), true );

    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion )
    {
      // We have to redo the below loop after fractures are generated
      updateState( subRegion );

      SlurryFluidBase const & fluid =
        getConstitutiveModel< SlurryFluidBase >( subRegion, subRegion.getReference< string >( viewKeyStruct::fluidNamesString() ) );
      arrayView3d< real64 const > const componentDens = fluid.componentDensity();
      arrayView2d< real64 > const componentDensOld = subRegion.getExtrinsicData< extrinsicMeshData::proppant::oldComponentDensity >();
      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
        for( localIndex c = 0; c < numComponents; ++c )
        {
          componentDensOld[ei][c] = componentDens[ei][0][c];
        }
      } );
    } );
  } );
  m_minAperture = m_bridgingFactor * m_proppantDiameter;
}

void ProppantTransport::preStepUpdate( real64 const & time,
                                       real64 const & GEOSX_UNUSED_PARAM( dt ),
                                       DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )
  {

    FlowSolverBase::precomputeData( mesh, regionNames );

    if( time <= 0 )
    {
      mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                    ElementSubRegionBase & subRegion )
      {
        updateProppantMobility( subRegion );
      } );
    }

    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion )
    {
      SlurryFluidBase const & fluid = getConstitutiveModel< SlurryFluidBase >( subRegion, subRegion.getReference< string >( viewKeyStruct::fluidNamesString() ) );

      arrayView3d< real64 const > const componentDens = fluid.componentDensity();
      arrayView2d< real64 > const componentDensOld = subRegion.getExtrinsicData< extrinsicMeshData::proppant::oldComponentDensity >();

      arrayView1d< real64 > const excessPackVolume = subRegion.getExtrinsicData< extrinsicMeshData::proppant::proppantExcessPackVolume >();
      arrayView2d< real64 > const cellBasedFlux = subRegion.getExtrinsicData< extrinsicMeshData::proppant::cellBasedFlux >();

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
        for( localIndex c = 0; c < m_numComponents; ++c )
        {
          componentDensOld[ei][c] = componentDens[ei][0][c];
        }
        excessPackVolume[ei] = 0.0;
        LvArray::tensorOps::fill< 3 >( cellBasedFlux[ei], 0.0 );
      } );
    } );

  } );

  updateCellBasedFlux( time, domain );
}

void ProppantTransport::postStepUpdate( real64 const & time_n,
                                        real64 const & dt_return,
                                        DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )
  {

    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion )
    {
      updateProppantMobility( subRegion );
    } );

    real64 const maxProppantConcentration = m_maxProppantConcentration;
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 > const & packVolFrac = subRegion.getExtrinsicData< extrinsicMeshData::proppant::proppantPackVolumeFraction >();
      arrayView1d< real64 > const & proppantConc = subRegion.getExtrinsicData< extrinsicMeshData::proppant::proppantConcentration >();

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
        if( proppantConc[ei] >= maxProppantConcentration || packVolFrac[ei] >= 1.0 )
        {
          packVolFrac[ei] = 1.0;
          proppantConc[ei] = maxProppantConcentration;
        }
      } );
    } );
  } );
  if( m_updateProppantPacking == 1 )
  {
    updateProppantPackVolume( time_n, dt_return, domain );
  }
}

void ProppantTransport::implicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                           real64 const & GEOSX_UNUSED_PARAM( dt ),
                                           DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

void ProppantTransport::implicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                              real64 const & GEOSX_UNUSED_PARAM( dt ),
                                              DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  integer const & numComponents = m_numComponents;

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion )
    {

      arrayView1d< real64 > const proppantConc =
        subRegion.getExtrinsicData< extrinsicMeshData::proppant::proppantConcentration >();
      arrayView1d< real64 const > const dProppantConc =
        subRegion.getExtrinsicData< extrinsicMeshData::proppant::deltaProppantConcentration >();

      arrayView2d< real64 > const componentConc =
        subRegion.getExtrinsicData< extrinsicMeshData::proppant::componentConcentration >();
      arrayView2d< real64 const > const dComponentConc =
        subRegion.getExtrinsicData< extrinsicMeshData::proppant::deltaComponentConcentration >();

      arrayView1d< real64 > const proppantLiftFlux =
        subRegion.getExtrinsicData< extrinsicMeshData::proppant::proppantLiftFlux >();

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
        proppantConc[ei] += dProppantConc[ei];
        proppantLiftFlux[ei] = 0.0;

        for( localIndex c = 0; c < numComponents; ++c )
        {
          componentConc[ei][c] += dComponentConc[ei][c];
        }
      } );
    } );
  } );
}

void ProppantTransport::setupDofs( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                   DofManager & dofManager ) const
{
  dofManager.addField( extrinsicMeshData::proppant::proppantConcentration::key(),
                       DofManager::Location::Elem,
                       m_numDofPerCell,
                       m_meshTargets );

  dofManager.addCoupling( extrinsicMeshData::proppant::proppantConcentration::key(),
                          extrinsicMeshData::proppant::proppantConcentration::key(),
                          DofManager::Connector::Face );
}


void ProppantTransport::assembleSystem( real64 const time,
                                        real64 const dt,
                                        DomainPartition & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  assembleAccumulationTerms( dt,
                             domain,
                             dofManager,
                             localMatrix,
                             localRhs );

  assembleFluxTerms( time,
                     dt,
                     domain,
                     dofManager,
                     localMatrix,
                     localRhs );
}

void ProppantTransport::assembleAccumulationTerms( real64 const dt,
                                                   DomainPartition const & domain,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  string const dofKey = dofManager.getKey( extrinsicMeshData::proppant::proppantConcentration::key() );

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel const & mesh,
                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase const & subRegion )
    {
      arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const & volume = subRegion.getElementVolume();

      arrayView2d< real64 const > const componentDensOld =
        subRegion.getExtrinsicData< extrinsicMeshData::proppant::oldComponentDensity >();
      arrayView1d< real64 const > const proppantConc =
        subRegion.getExtrinsicData< extrinsicMeshData::proppant::proppantConcentration >();
      arrayView1d< real64 const > const dProppantConc =
        subRegion.getExtrinsicData< extrinsicMeshData::proppant::deltaProppantConcentration >();
      arrayView1d< real64 const > const proppantPackVolFrac =
        subRegion.getExtrinsicData< extrinsicMeshData::proppant::proppantPackVolumeFraction >();
      arrayView1d< real64 const > const proppantLiftFlux =
        subRegion.getExtrinsicData< extrinsicMeshData::proppant::proppantLiftFlux >();

      SlurryFluidBase const & fluid =
        getConstitutiveModel< SlurryFluidBase >( subRegion, subRegion.getReference< string >( viewKeyStruct::fluidNamesString() ) );

      arrayView3d< real64 const > const componentDens = fluid.componentDensity();
      arrayView3d< real64 const > const dCompDens_dPres = fluid.dComponentDensity_dPressure();
      arrayView4d< real64 const > const dCompDens_dCompConc = fluid.dComponentDensity_dComponentConcentration();

      AccumulationKernel::launch( subRegion.size(),
                                  m_numComponents,
                                  m_numDofPerCell,
                                  dofManager.rankOffset(),
                                  dofNumber,
                                  elemGhostRank,
                                  proppantConc,
                                  dProppantConc,
                                  componentDensOld,
                                  componentDens,
                                  dCompDens_dPres,
                                  dCompDens_dCompConc,
                                  volume,
                                  proppantPackVolFrac,
                                  proppantLiftFlux,
                                  dt,
                                  m_maxProppantConcentration,
                                  localMatrix,
                                  localRhs );
    } );
  } );
}


void ProppantTransport::assembleFluxTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                           real64 const dt,
                                           DomainPartition const & domain,
                                           DofManager const & dofManager,
                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                           arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  R1Tensor downVector = gravityVector();
  LvArray::tensorOps::normalize< 3 >( downVector );

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel const & mesh,
                                               arrayView1d< string const > const & )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();

    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );


    string const dofKey = dofManager.getKey( extrinsicMeshData::proppant::proppantConcentration::key() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
      elemManager.constructViewAccessor< array1d< globalIndex >, arrayView1d< globalIndex const > >( dofKey );

    typename FluxKernel::FlowAccessors flowAccessors( elemManager, getName() );
    typename FluxKernel::ParticleFluidAccessors particleFluidAccessors( elemManager, getName() );
    typename FluxKernel::SlurryFluidAccessors slurryFluidAccessors( elemManager, getName() );
    typename FluxKernel::PermeabilityAccessors permAccessors( elemManager, getName() );

    fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( auto const & stencil )
    {

      SurfaceElementStencilWrapper stencilWrapper = stencil.createStencilWrapper();

      FluxKernel::launch( stencilWrapper,
                          m_numDofPerCell,
                          dt,
                          dofManager.rankOffset(),
                          m_updateProppantPacking,
                          downVector,
                          dofNumberAccessor.toNestedViewConst(),
                          flowAccessors.get< extrinsicMeshData::ghostRank >(),
                          flowAccessors.get< extrinsicMeshData::flow::pressure >(),
                          flowAccessors.get< extrinsicMeshData::flow::deltaPressure >(),
                          flowAccessors.get< extrinsicMeshData::proppant::proppantConcentration >(),
                          flowAccessors.get< extrinsicMeshData::proppant::deltaProppantConcentration >(),
                          slurryFluidAccessors.get< extrinsicMeshData::slurryfluid::componentDensity >(),
                          slurryFluidAccessors.get< extrinsicMeshData::slurryfluid::dComponentDensity_dPressure >(),
                          slurryFluidAccessors.get< extrinsicMeshData::slurryfluid::dComponentDensity_dComponentConcentration >(),
                          flowAccessors.get< extrinsicMeshData::flow::gravityCoefficient >(),
                          slurryFluidAccessors.get< extrinsicMeshData::singlefluid::density >(),
                          slurryFluidAccessors.get< extrinsicMeshData::singlefluid::dDensity_dPressure >(),
                          slurryFluidAccessors.get< extrinsicMeshData::slurryfluid::dDensity_dProppantConcentration >(),
                          slurryFluidAccessors.get< extrinsicMeshData::slurryfluid::dDensity_dComponentConcentration >(),
                          slurryFluidAccessors.get< extrinsicMeshData::singlefluid::viscosity >(),
                          slurryFluidAccessors.get< extrinsicMeshData::singlefluid::dViscosity_dPressure >(),
                          slurryFluidAccessors.get< extrinsicMeshData::slurryfluid::dViscosity_dProppantConcentration >(),
                          slurryFluidAccessors.get< extrinsicMeshData::slurryfluid::dViscosity_dComponentConcentration >(),
                          slurryFluidAccessors.get< extrinsicMeshData::slurryfluid::fluidDensity >(),
                          slurryFluidAccessors.get< extrinsicMeshData::slurryfluid::dFluidDensity_dPressure >(),
                          slurryFluidAccessors.get< extrinsicMeshData::slurryfluid::dFluidDensity_dComponentConcentration >(),
                          particleFluidAccessors.get< extrinsicMeshData::particlefluid::settlingFactor >(),
                          particleFluidAccessors.get< extrinsicMeshData::particlefluid::dSettlingFactor_dPressure >(),
                          particleFluidAccessors.get< extrinsicMeshData::particlefluid::dSettlingFactor_dProppantConcentration >(),
                          particleFluidAccessors.get< extrinsicMeshData::particlefluid::dSettlingFactor_dComponentConcentration >(),
                          particleFluidAccessors.get< extrinsicMeshData::particlefluid::collisionFactor >(),
                          particleFluidAccessors.get< extrinsicMeshData::particlefluid::dCollisionFactor_dProppantConcentration >(),
                          flowAccessors.get< extrinsicMeshData::proppant::isProppantMobile >(),
                          permAccessors.permeability(),
                          permAccessors.permeabilityMultiplier(),
                          flowAccessors.get< extrinsicMeshData::elementAperture >(),
                          localMatrix,
                          localRhs );
    } );
  } );
}

void ProppantTransport::applyBoundaryConditions( real64 const time_n,
                                                 real64 const dt,
                                                 DomainPartition & domain,
                                                 DofManager const & dofManager,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  string const dofKey = dofManager.getKey( extrinsicMeshData::proppant::proppantConcentration::key() );
  globalIndex const rankOffset = dofManager.rankOffset();

  //  Apply Dirichlet BC for proppant concentration

  fsManager.apply( time_n + dt,
                   domain,
                   "ElementRegions",
                   extrinsicMeshData::proppant::proppantConcentration::key(),
                   [&]( FieldSpecificationBase const & fs,
                        string const &,
                        SortedArrayView< localIndex const > const & lset,
                        Group & subRegion,
                        string const & )
  {
    arrayView1d< globalIndex const > const
    dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

    arrayView1d< real64 const > const
    proppantConc = subRegion.getReference< array1d< real64 > >( extrinsicMeshData::proppant::proppantConcentration::key() );

    arrayView1d< real64 const > const
    dProppantConc = subRegion.getReference< array1d< real64 > >( extrinsicMeshData::proppant::deltaProppantConcentration::key() );

    fs.applyBoundaryConditionToSystem< FieldSpecificationEqual,
                                       parallelDevicePolicy<> >( lset,
                                                                 time_n + dt,
                                                                 subRegion,
                                                                 dofNumber,
                                                                 rankOffset,
                                                                 localMatrix,
                                                                 localRhs,
                                                                 [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      return proppantConc[a] + dProppantConc[a];
    } );
  } );

  //  Apply Dirichlet BC for component concentration
  if( m_numComponents > 0 )
  {
    map< string, map< string, array1d< bool > > > bcStatusMap; // map to check consistent application of BC

    fsManager.apply( time_n + dt,
                     domain,
                     "ElementRegions",
                     extrinsicMeshData::proppant::proppantConcentration::key(),
                     [&]( FieldSpecificationBase const &,
                          string const & setName,
                          SortedArrayView< localIndex const > const &,
                          Group & subRegion,
                          string const & )
    {

      string const & subRegionName = subRegion.getName();
      GEOSX_ERROR_IF( bcStatusMap[subRegionName].count( setName ) > 0, "Conflicting proppant boundary conditions on set " << setName );
      bcStatusMap[subRegionName][setName].resize( m_numComponents );
      bcStatusMap[subRegionName][setName].setValues< serialPolicy >( false );

    } );

    fsManager.apply( time_n + dt,
                     domain,
                     "ElementRegions",
                     extrinsicMeshData::proppant::componentConcentration::key(),
                     [&] ( FieldSpecificationBase const & fs,
                           string const & setName,
                           SortedArrayView< localIndex const > const & targetSet,
                           Group & subRegion,
                           string const & )
    {

      string const & subRegionName = subRegion.getName();
      localIndex const comp = fs.getComponent();

      GEOSX_ERROR_IF( bcStatusMap[subRegionName].count( setName ) == 0, "Proppant boundary condition not prescribed on set '" << setName << "'" );
      GEOSX_ERROR_IF( bcStatusMap[subRegionName][setName][comp], "Conflicting composition[" << comp << "] boundary conditions on set '" << setName << "'" );
      bcStatusMap[subRegionName][setName][comp] = true;

      fs.applyFieldValue< FieldSpecificationEqual >( targetSet,
                                                     time_n + dt,
                                                     subRegion,
                                                     extrinsicMeshData::proppant::bcComponentConcentration::key() );

    } );

    bool bcConsistent = true;
    for( auto const & bcStatusEntryOuter : bcStatusMap )
    {
      for( auto const & bcStatusEntryInner : bcStatusEntryOuter.second )
      {
        for( localIndex ic = 0; ic < m_numComponents; ++ic )
        {
          bcConsistent &= bcStatusEntryInner.second[ic];
          GEOSX_WARNING_IF( !bcConsistent, "Composition boundary condition not applied to component " << ic
                                                                                                      << " on region '" << bcStatusEntryOuter.first << "',"
                                                                                                      << " set '" << bcStatusEntryInner.first << "'" );
        }
      }
    }

    GEOSX_ERROR_IF( !bcConsistent, "Inconsistent composition boundary conditions" );

    fsManager.apply( time_n + dt,
                     domain,
                     "ElementRegions",
                     extrinsicMeshData::proppant::proppantConcentration::key(),
                     [&] ( FieldSpecificationBase const &,
                           string const &,
                           SortedArrayView< localIndex const > const & targetSet,
                           Group & subRegion,
                           string const & )
    {
      arrayView1d< integer const > const ghostRank =
        subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString() );
      arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

      arrayView2d< real64 const > const compConc =
        subRegion.getReference< array2d< real64 > >( extrinsicMeshData::proppant::componentConcentration::key() );
      arrayView2d< real64 const > const deltaCompConc =
        subRegion.getReference< array2d< real64 > >( extrinsicMeshData::proppant::deltaComponentConcentration::key() );
      arrayView2d< real64 const > const bcCompConc =
        subRegion.getReference< array2d< real64 > >( extrinsicMeshData::proppant::bcComponentConcentration::key() );

      forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
      {
        localIndex const ei = targetSet[a];
        if( ghostRank[ei] >= 0 )
          return;

        globalIndex const dofIndex = dofNumber[ei];
        localIndex const localRow = dofIndex - rankOffset;
        real64 rhsValue;

        for( localIndex ic = 0; ic < m_numComponents; ++ic )
        {
          FieldSpecificationEqual::SpecifyFieldValue( dofIndex + ic + 1,
                                                      rankOffset,
                                                      localMatrix,
                                                      rhsValue,
                                                      bcCompConc[ei][ic],
                                                      compConc[ei][ic] + deltaCompConc[ei][ic] );
          localRhs[localRow + ic + 1] = rhsValue;
        }
      } );
    } );
  }
}

real64
ProppantTransport::calculateResidualNorm( DomainPartition const & domain,
                                          DofManager const & dofManager,
                                          arrayView1d< real64 const > const & localRhs )
{
  localIndex const NDOF = m_numDofPerCell;

  localIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( extrinsicMeshData::proppant::proppantConcentration::key() );

  // compute the norm of local residual scaled by cell pore volume
  real64 localResidualNorm = 0.0;

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel const & mesh,
                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const volume = subRegion.getElementVolume();

      RAJA::ReduceSum< parallelDeviceReduce, real64 > localSum( 0.0 );

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
        if( elemGhostRank[ei] < 0 )
        {
          localIndex const lid = dofNumber[ei] - rankOffset;
          for( localIndex idof = 0; idof < NDOF; ++idof )
          {
            real64 const val = localRhs[lid] / volume[ei];
            localSum += val * val;
          }
        }
      } );

      localResidualNorm += localSum.get();
    } );
  } );

  // compute global residual norm
  real64 const globalResidualNorm = MpiWrapper::sum( localResidualNorm, MPI_COMM_GEOSX );

  return sqrt( globalResidualNorm );
}

void ProppantTransport::applySystemSolution( DofManager const & dofManager,
                                             arrayView1d< real64 const > const & localSolution,
                                             real64 const scalingFactor,
                                             DomainPartition & domain )
{
  dofManager.addVectorToField( localSolution,
                               extrinsicMeshData::proppant::proppantConcentration::key(),
                               extrinsicMeshData::proppant::deltaProppantConcentration::key(),
                               scalingFactor,
                               { m_numDofPerCell, 0, 1 } );


  if( m_numDofPerCell > 1 )
  {
    dofManager.addVectorToField( localSolution,
                                 extrinsicMeshData::proppant::proppantConcentration::key(),
                                 extrinsicMeshData::proppant::deltaComponentConcentration::key(),
                                 scalingFactor,
                                 { m_numDofPerCell, 1, m_numDofPerCell } );
  }

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( extrinsicMeshData::proppant::deltaProppantConcentration::key() );
  fieldNames["elems"].emplace_back( extrinsicMeshData::proppant::deltaComponentConcentration::key() );

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )
  {


    CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), true );

    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      updateComponentDensity( subRegion );
    } );

  } );

}

void ProppantTransport::solveSystem( DofManager const & dofManager,
                                     ParallelMatrix & matrix,
                                     ParallelVector & rhs,
                                     ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::solveSystem( dofManager, matrix, rhs, solution );
}

void ProppantTransport::resetStateToBeginningOfStep( DomainPartition & domain )
{
  integer const numComponents = m_numComponents;

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 > const & dProppantConc =
        subRegion.getExtrinsicData< extrinsicMeshData::proppant::deltaProppantConcentration >();
      arrayView2d< real64 > const & dComponentConc =
        subRegion.getExtrinsicData< extrinsicMeshData::proppant::deltaComponentConcentration >();
      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
        dProppantConc[ei] = 0.0;
        for( localIndex c = 0; c < numComponents; ++c )
        {
          dComponentConc[ei][c] = 0.0;
        }
      } );

      updateState( subRegion );
    } );
  } );
}



void ProppantTransport::updateCellBasedFlux( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                             DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  R1Tensor downVector = gravityVector();
  LvArray::tensorOps::normalize< 3 >( downVector );

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager & elemManager = mesh.getElemManager();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > const & cellBasedFluxAccessor =
    elemManager.constructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( extrinsicMeshData::proppant::cellBasedFlux::key() );

  typename FluxKernel::CellBasedFluxFlowAccessors flowAccessors( elemManager, getName() );
  typename FluxKernel::CellBasedFluxSlurryFluidAccessors slurryFluidAccessors( elemManager, getName() );
  typename FluxKernel::PermeabilityAccessors permAccessors( elemManager, getName() );

  fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( auto const & stencil )
  {
    SurfaceElementStencilWrapper stencilWrapper = stencil.createStencilWrapper();

    FluxKernel::launchCellBasedFluxCalculation( stencilWrapper,
                                                downVector,
                                                flowAccessors.pressure(),
                                                flowAccessors.gravityCoefficient(),
                                                slurryFluidAccessors.density(),
                                                slurryFluidAccessors.viscosity(),
                                                permAccessors.permeability(),
                                                permAccessors.permeabilityMultiplier(),
                                                flowAccessors.elementAperture(),
                                                cellBasedFluxAccessor.toNestedView() );
  } );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( extrinsicMeshData::proppant::cellBasedFlux::key() );

  CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), true );
}

void ProppantTransport::updateProppantPackVolume( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                  real64 const dt,
                                                  DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  R1Tensor downVector = gravityVector();
  LvArray::tensorOps::normalize< 3 >( downVector );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );


  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    // For data modified through an accessor, we must create the view accessor
    // every time in order to ensure the data gets properly touched on device
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const proppantConc =
      elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( extrinsicMeshData::proppant::proppantConcentration::key() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const proppantPackVolFrac =
      elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( extrinsicMeshData::proppant::proppantPackVolumeFraction::key() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const proppantExcessPackVolume =
      elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( extrinsicMeshData::proppant::proppantExcessPackVolume::key() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const proppantLiftFlux =
      elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( extrinsicMeshData::proppant::proppantLiftFlux::key() );

    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const
    aperture = elemManager.constructArrayViewAccessor< real64, 1 >( FaceElementSubRegion::viewKeyStruct::elementApertureString() );

    typename ProppantPackVolumeKernel::FlowAccessors flowAccessors( elemManager, getName() );
    typename ProppantPackVolumeKernel::SlurryFluidAccessors slurryFluidAccessors( elemManager, getName() );
    typename ProppantPackVolumeKernel::ParticleFluidAccessors particleFluidAccessors( elemManager, getName() );

    fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( auto const & stencil )
    {
      ProppantPackVolumeKernel::launchProppantPackVolumeCalculation( stencil,
                                                                     dt,
                                                                     m_proppantDensity,
                                                                     m_proppantDiameter,
                                                                     m_maxProppantConcentration,
                                                                     downVector,
                                                                     m_criticalShieldsNumber,
                                                                     m_frictionCoefficient,
                                                                     particleFluidAccessors.get< extrinsicMeshData::particlefluid::settlingFactor >(),
                                                                     slurryFluidAccessors.get< extrinsicMeshData::singlefluid::density >(),
                                                                     slurryFluidAccessors.get< extrinsicMeshData::slurryfluid::fluidDensity >(),
                                                                     slurryFluidAccessors.get< extrinsicMeshData::slurryfluid::fluidViscosity >(),
                                                                     flowAccessors.get< extrinsicMeshData::proppant::isProppantMobile >(),
                                                                     flowAccessors.get< extrinsicMeshData::proppant::isProppantBoundary >(),
                                                                     flowAccessors.get< extrinsicMeshData::elementAperture >(),
                                                                     flowAccessors.get< extrinsicMeshData::elementVolume >(),
                                                                     flowAccessors.get< extrinsicMeshData::ghostRank >(),
                                                                     flowAccessors.get< extrinsicMeshData::proppant::cellBasedFlux >(),
                                                                     proppantConc.toNestedView(),
                                                                     proppantPackVolFrac.toNestedView(),
                                                                     proppantExcessPackVolume.toNestedView(),
                                                                     proppantLiftFlux.toNestedView() );
    } );

    {
      std::map< string, string_array > fieldNames;
      fieldNames["elems"].emplace_back( extrinsicMeshData::proppant::proppantConcentration::key() );
      fieldNames["elems"].emplace_back( extrinsicMeshData::proppant::proppantPackVolumeFraction::key() );
      fieldNames["elems"].emplace_back( extrinsicMeshData::proppant::proppantExcessPackVolume::key() );
      fieldNames["elems"].emplace_back( extrinsicMeshData::proppant::proppantLiftFlux::key() );


      CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), true );
    }

    elemManager.forElementSubRegions( regionNames,
                                      [&]( localIndex const,
                                           ElementSubRegionBase & subRegion )
    {
      updateProppantMobility( subRegion );
    } );


    fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( auto const & stencil )
    {
      ProppantPackVolumeKernel::launchProppantPackVolumeUpdate( stencil,
                                                                downVector,
                                                                m_maxProppantConcentration,
                                                                flowAccessors.get< extrinsicMeshData::proppant::isProppantMobile >(),
                                                                proppantExcessPackVolume.toNestedViewConst(),
                                                                proppantConc.toNestedView(),
                                                                proppantPackVolFrac.toNestedView() );
    } );

    {
      std::map< string, string_array > fieldNames;
      fieldNames["elems"].emplace_back( extrinsicMeshData::proppant::proppantConcentration::key() );
      fieldNames["elems"].emplace_back( extrinsicMeshData::proppant::proppantPackVolumeFraction::key() );

      CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), true );
    }

    elemManager.forElementSubRegions( regionNames,
                                      [&]( localIndex const,
                                           ElementSubRegionBase & subRegion )
    {
      updateProppantMobility( subRegion );
    } );

  } );

}


REGISTER_CATALOG_ENTRY( SolverBase, ProppantTransport, string const &, Group * const )
} /* namespace geosx */
