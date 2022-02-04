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
 * @file SinglePhaseFVM.cpp
 */

#include "SinglePhaseFVM.hpp"

#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/fluid/singleFluidSelector.hpp"
#include "constitutive/permeability/PermeabilityExtrinsicData.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseFVMKernels.hpp"
#include "physicsSolvers/multiphysics/SinglePhasePoromechanicsFluxKernels.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace SinglePhaseBaseKernels;
using namespace SinglePhaseFVMKernels;
using namespace SinglePhasePoromechanicsFluxKernels;


template< typename BASE >
SinglePhaseFVM< BASE >::SinglePhaseFVM( const string & name,
                                        Group * const parent ):
  BASE( name, parent )
{
  m_numDofPerCell = 1;
}

template< typename BASE >
void SinglePhaseFVM< BASE >::initializePreSubGroups()
{
  BASE::initializePreSubGroups();

  DomainPartition & domain = this->template getGroupByPath< DomainPartition >( "/Problem/domain" );
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();

  if( !fvManager.hasGroup< FluxApproximationBase >( m_discretizationName ) )
  {
    GEOSX_ERROR( "A discretization deriving from FluxApproximationBase must be selected with SinglePhaseFVM" );
  }
}

template< typename BASE >
void SinglePhaseFVM< BASE >::setupDofs( DomainPartition const & domain,
                                        DofManager & dofManager ) const
{
  dofManager.addField( extrinsicMeshData::flow::pressure::key(),
                       DofManager::Location::Elem,
                       1,
                       targetRegionNames() );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  dofManager.addCoupling( extrinsicMeshData::flow::pressure::key(), fluxApprox );
}

template< typename BASE >
void SinglePhaseFVM< BASE >::setupSystem( DomainPartition & domain,
                                          DofManager & dofManager,
                                          CRSMatrix< real64, globalIndex > & localMatrix,
                                          ParallelVector & rhs,
                                          ParallelVector & solution,
                                          bool const setSparsity )
{
  GEOSX_MARK_FUNCTION;
  BASE::setupSystem( domain,
                     dofManager,
                     localMatrix,
                     rhs,
                     solution,
                     setSparsity );

}

template< typename BASE >
real64 SinglePhaseFVM< BASE >::calculateResidualNorm( DomainPartition const & domain,
                                                      DofManager const & dofManager,
                                                      arrayView1d< real64 const > const & localRhs )
{
  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  string const dofKey = dofManager.getKey( extrinsicMeshData::flow::pressure::key() );
  globalIndex const rankOffset = dofManager.rankOffset();

  // compute the norm of local residual scaled by cell pore volume
  real64 localResidualNorm[3] = { 0.0, 0.0, 0.0 };
  forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                  auto const & subRegion )
  {
    arrayView1d< globalIndex const > const & dofNumber = subRegion.template getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
    arrayView1d< real64 const > const & volume         = subRegion.getElementVolume();
    arrayView1d< real64 const > const & densOld        = subRegion.template getExtrinsicData< extrinsicMeshData::flow::densityOld >();

    CoupledSolidBase const & solidModel = subRegion.template getConstitutiveModel< CoupledSolidBase >( m_solidModelNames[targetIndex] );

    arrayView2d< real64 const > const & porosityOld = solidModel.getOldPorosity();

    ResidualNormKernel::launch< parallelDevicePolicy<>, parallelDeviceReduce >( localRhs,
                                                                                rankOffset,
                                                                                dofNumber,
                                                                                elemGhostRank,
                                                                                volume,
                                                                                densOld,
                                                                                porosityOld,
                                                                                localResidualNorm );
  } );

  // compute global residual norm
  real64 globalResidualNorm[3] = {0, 0, 0};
  MpiWrapper::allReduce( localResidualNorm,
                         globalResidualNorm,
                         3,
                         MPI_SUM,
                         MPI_COMM_GEOSX );


  real64 const residual = sqrt( globalResidualNorm[0] ) / ( ( globalResidualNorm[1] + m_fluxEstimate ) / (globalResidualNorm[2]+1) );
  return residual;
}


template< typename BASE >
void SinglePhaseFVM< BASE >::applySystemSolution( DofManager const & dofManager,
                                                  arrayView1d< real64 const > const & localSolution,
                                                  real64 const scalingFactor,
                                                  DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  dofManager.addVectorToField( localSolution,
                               extrinsicMeshData::flow::pressure::key(),
                               extrinsicMeshData::flow::deltaPressure::key(),
                               scalingFactor );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( extrinsicMeshData::flow::deltaPressure::key() ) );

  CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), true );
}

template<>
void SinglePhaseFVM< SinglePhaseBase >::assembleFluxTerms( real64 const GEOSX_UNUSED_PARAM ( time_n ),
                                                           real64 const dt,
                                                           DomainPartition const & domain,
                                                           DofManager const & dofManager,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                           arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
  ElementRegionManager const & elemManager = mesh.getElemManager();

  string const & dofKey = dofManager.getKey( extrinsicMeshData::flow::pressure::key() );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  elemDofNumber = elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
  elemDofNumber.setName( this->getName() + "/accessors/" + dofKey );

  fluxApprox.forAllStencils( mesh, [&]( auto & stencil )
  {
    typename TYPEOFREF( stencil ) ::StencilWrapper stencilWrapper = stencil.createStencilWrapper();


    typename FluxKernel::SinglePhaseFlowAccessors flowAccessors( elemManager, getName() );
    typename FluxKernel::SinglePhaseFluidAccessors fluidAccessors( elemManager, getName(), targetRegionNames(), fluidModelNames() );
    typename FluxKernel::PermeabilityAccessors permAccessors( elemManager, getName(), targetRegionNames(), permeabilityModelNames() );


    FluxKernel::launch( stencilWrapper,
                        dt,
                        dofManager.rankOffset(),
                        elemDofNumber.toNestedViewConst(),
                        flowAccessors.get< extrinsicMeshData::ghostRank >(),
                        flowAccessors.get< extrinsicMeshData::flow::pressure >(),
                        flowAccessors.get< extrinsicMeshData::flow::deltaPressure >(),
                        flowAccessors.get< extrinsicMeshData::flow::gravityCoefficient >(),
                        fluidAccessors.get< extrinsicMeshData::singlefluid::density >(),
                        fluidAccessors.get< extrinsicMeshData::singlefluid::dDensity_dPressure >(),
                        flowAccessors.get< extrinsicMeshData::flow::mobility >(),
                        flowAccessors.get< extrinsicMeshData::flow::dMobility_dPressure >(),
                        permAccessors.get< extrinsicMeshData::permeability::permeability >(),
                        permAccessors.get< extrinsicMeshData::permeability::dPerm_dPressure >(),
                        localMatrix,
                        localRhs );

  } );
}


template<>
void SinglePhaseFVM< SinglePhaseProppantBase >::assembleFluxTerms( real64 const GEOSX_UNUSED_PARAM ( time_n ),
                                                                   real64 const dt,
                                                                   DomainPartition const & domain,
                                                                   DofManager const & dofManager,
                                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                   arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
  ElementRegionManager const & elemManager = mesh.getElemManager();

  string const & dofKey = dofManager.getKey( extrinsicMeshData::flow::pressure::key() );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  elemDofNumber = elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
  elemDofNumber.setName( this->getName() + "/accessors/" + dofKey );

  fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( auto & stencil )
  {
    typename TYPEOFREF( stencil ) ::StencilWrapper stencilWrapper = stencil.createStencilWrapper();

    typename FluxKernel::SinglePhaseFlowAccessors flowAccessors( elemManager, getName() );
    typename FluxKernel::SlurryFluidAccessors fluidAccessors( elemManager, getName(), targetRegionNames(), fluidModelNames() );
    typename FluxKernel::ProppantPermeabilityAccessors permAccessors( elemManager, getName(), targetRegionNames(), permeabilityModelNames() );


    FaceElementFluxKernel::launch( stencilWrapper,
                                   dt,
                                   dofManager.rankOffset(),
                                   elemDofNumber.toNestedViewConst(),
                                   flowAccessors.get< extrinsicMeshData::ghostRank >(),
                                   flowAccessors.get< extrinsicMeshData::flow::pressure >(),
                                   flowAccessors.get< extrinsicMeshData::flow::deltaPressure >(),
                                   flowAccessors.get< extrinsicMeshData::flow::gravityCoefficient >(),
                                   fluidAccessors.get< extrinsicMeshData::slurryfluid::density >(),
                                   fluidAccessors.get< extrinsicMeshData::slurryfluid::dDensity_dPressure >(),
                                   flowAccessors.get< extrinsicMeshData::flow::mobility >(),
                                   flowAccessors.get< extrinsicMeshData::flow::dMobility_dPressure >(),
                                   permAccessors.get< extrinsicMeshData::permeability::permeability >(),
                                   permAccessors.get< extrinsicMeshData::permeability::dPerm_dPressure >(),
                                   permAccessors.get< extrinsicMeshData::permeability::dPerm_dDispJump >(),
                                   permAccessors.get< extrinsicMeshData::permeability::permeabilityMultiplier >(),
                                   this->gravityVector(),
                                   localMatrix,
                                   localRhs );
  } );
}


template< typename BASE >
void SinglePhaseFVM< BASE >::assembleFluxTermsExplicit( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                        real64 const dt ,
                                                        DomainPartition & domain )
{
	  GEOSX_MARK_FUNCTION;
	  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

	  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
	  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
	  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
	  ElementRegionManager const & elemManager = mesh.getElemManager();

	  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > >
	  	  m_referencePressure = elemManager.constructArrayViewAccessor< real64, 1 >( SinglePhaseBase::viewKeyStruct::referencePressureString() );
	  m_referencePressure.setName( this->getName() + "/accessors/" + SinglePhaseBase::viewKeyStruct::referencePressureString() );

	  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
	  ElementRegionManager & elemManager0 = meshLevel.getElemManager();
	  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > >
	  	  m_fluidMass = elemManager0.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( SinglePhaseBase::viewKeyStruct::fluidMassString() );
	  //m_fluidMass = elemManager0.ConstructArrayViewAccessor1< real64, 1 >( viewKeyStruct::fluidMassString );
	  m_fluidMass.setName( this->getName() + "/accessors/" + SinglePhaseBase::viewKeyStruct::fluidMassString() );

	  fluxApprox.forAllStencils( mesh, [&]( auto & stencil )
	  {
	    typename TYPEOFREF( stencil ) ::StencilWrapper stencilWrapper = stencil.createStencilWrapper();


	    typename FluxKernel::SinglePhaseFlowAccessors flowAccessors( elemManager, this->getName() );
	    typename FluxKernel::SinglePhaseFluidAccessors fluidAccessors( elemManager, this->getName(), this->targetRegionNames(), this->fluidModelNames() );
	    typename FluxKernel::PermeabilityAccessors permAccessors( elemManager, this->getName(), this->targetRegionNames(), this->permeabilityModelNames() );

	    FluxKernel::launch( stencilWrapper,
	                        dt,
	                        flowAccessors.get< extrinsicMeshData::flow::pressure >(),
	                        flowAccessors.get< extrinsicMeshData::flow::gravityCoefficient >(),
	                        fluidAccessors.get< extrinsicMeshData::singlefluid::density >(),
	                        flowAccessors.get< extrinsicMeshData::flow::mobility >(),
	                        permAccessors.get< extrinsicMeshData::permeability::permeability >(),
							permAccessors.get< extrinsicMeshData::permeability::dPerm_dPressure >(),
							m_referencePressure.toNestedViewConst(),
							&m_fluidMass );

	  } );

}

template< typename BASE >
void SinglePhaseFVM< BASE >::assemblePoroelasticFluxTerms( real64 const GEOSX_UNUSED_PARAM ( time_n ),
                                                           real64 const dt,
                                                           DomainPartition const & domain,
                                                           DofManager const & dofManager,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                           arrayView1d< real64 > const & localRhs,
                                                           string const & jumpDofKey )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
  ElementRegionManager const & elemManager = mesh.getElemManager();

  string const & pressureDofKey = dofManager.getKey( extrinsicMeshData::flow::pressure::key() );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  pressureDofNumber = mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( pressureDofKey );
  pressureDofNumber.setName( this->getName() + "/accessors/" + pressureDofKey );

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  jumpDofNumber = mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( jumpDofKey );
  jumpDofNumber.setName( this->getName() + "/accessors/" + jumpDofKey );

  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const > > dPerm_dDispJump =
    mesh.getElemManager().constructMaterialArrayViewAccessor< real64, 4 >( extrinsicMeshData::permeability::dPerm_dDispJump::key(),
                                                                           targetRegionNames(),
                                                                           m_permeabilityModelNames,
                                                                           true );

  fluxApprox.forStencils< CellElementStencilTPFA, SurfaceElementStencil, EmbeddedSurfaceToCellStencil >( mesh, [&]( auto & stencil )
  {
    typename TYPEOFREF( stencil ) ::StencilWrapper stencilWrapper = stencil.createStencilWrapper();

    typename FluxKernel::SinglePhaseFlowAccessors flowAccessors( elemManager, this->getName() );
    typename FluxKernel::SinglePhaseFluidAccessors fluidAccessors( elemManager, this->getName(), this->targetRegionNames(), this->fluidModelNames() );
    typename FluxKernel::PermeabilityAccessors permAccessors( elemManager, this->getName(), this->targetRegionNames(), this->permeabilityModelNames() );

    EmbeddedSurfaceFluxKernel::launch( stencilWrapper,
                                       dt,
                                       dofManager.rankOffset(),
                                       pressureDofNumber.toNestedViewConst(),
                                       jumpDofNumber.toNestedViewConst(),
                                       flowAccessors.get< extrinsicMeshData::ghostRank >(),
                                       flowAccessors.get< extrinsicMeshData::flow::pressure >(),
                                       flowAccessors.get< extrinsicMeshData::flow::deltaPressure >(),
                                       flowAccessors.get< extrinsicMeshData::flow::gravityCoefficient >(),
                                       fluidAccessors.get< extrinsicMeshData::singlefluid::density >(),
                                       fluidAccessors.get< extrinsicMeshData::singlefluid::dDensity_dPressure >(),
                                       flowAccessors.get< extrinsicMeshData::flow::mobility >(),
                                       flowAccessors.get< extrinsicMeshData::flow::dMobility_dPressure >(),
                                       permAccessors.get< extrinsicMeshData::permeability::permeability >(),
                                       permAccessors.get< extrinsicMeshData::permeability::dPerm_dPressure >(),
                                       dPerm_dDispJump.toNestedViewConst(),
                                       localMatrix,
                                       localRhs );
  } );
}

template< typename BASE >
void SinglePhaseFVM< BASE >::assembleHydrofracFluxTerms( real64 const GEOSX_UNUSED_PARAM ( time_n ),
                                                         real64 const dt,
                                                         DomainPartition const & domain,
                                                         DofManager const & dofManager,
                                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                         arrayView1d< real64 > const & localRhs,
                                                         CRSMatrixView< real64, localIndex const > const & dR_dAper )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
  ElementRegionManager const & elemManager = mesh.getElemManager();

  string const & dofKey = dofManager.getKey( extrinsicMeshData::flow::pressure::key() );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  elemDofNumber = elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
  elemDofNumber.setName( this->getName() + "/accessors/" + dofKey );

  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const > > dPerm_dDispJump =
    mesh.getElemManager().constructMaterialArrayViewAccessor< real64, 4 >( extrinsicMeshData::permeability::dPerm_dDispJump::key(),
                                                                           targetRegionNames(),
                                                                           m_permeabilityModelNames,
                                                                           true );

  fluxApprox.forStencils< CellElementStencilTPFA, SurfaceElementStencil, FaceElementToCellStencil >( mesh, [&]( auto & stencil )
  {
    typename TYPEOFREF( stencil ) ::StencilWrapper stencilWrapper = stencil.createStencilWrapper();

    typename FluxKernel::SinglePhaseFlowAccessors flowAccessors( elemManager, this->getName() );
    typename FluxKernel::SinglePhaseFluidAccessors fluidAccessors( elemManager, this->getName(), this->targetRegionNames(), this->fluidModelNames() );
    typename FluxKernel::PermeabilityAccessors permAccessors( elemManager, this->getName(), this->targetRegionNames(), this->permeabilityModelNames() );

    FaceElementFluxKernel::launch( stencilWrapper,
                                   dt,
                                   dofManager.rankOffset(),
                                   elemDofNumber.toNestedViewConst(),
                                   flowAccessors.get< extrinsicMeshData::ghostRank >(),
                                   flowAccessors.get< extrinsicMeshData::flow::pressure >(),
                                   flowAccessors.get< extrinsicMeshData::flow::deltaPressure >(),
                                   flowAccessors.get< extrinsicMeshData::flow::gravityCoefficient >(),
                                   fluidAccessors.get< extrinsicMeshData::singlefluid::density >(),
                                   fluidAccessors.get< extrinsicMeshData::singlefluid::dDensity_dPressure >(),
                                   flowAccessors.get< extrinsicMeshData::flow::mobility >(),
                                   flowAccessors.get< extrinsicMeshData::flow::dMobility_dPressure >(),
                                   permAccessors.get< extrinsicMeshData::permeability::permeability >(),
                                   permAccessors.get< extrinsicMeshData::permeability::dPerm_dPressure >(),
                                   dPerm_dDispJump.toNestedViewConst(),
                                   localMatrix,
                                   localRhs,
                                   dR_dAper );
  } );
}

template< typename BASE >
void
SinglePhaseFVM< BASE >::applyBoundaryConditions( real64 const time_n,
                                                 real64 const dt,
                                                 DomainPartition & domain,
                                                 DofManager const & dofManager,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  BASE::applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );
  applyFaceDirichletBC( time_n, dt, dofManager, domain, localMatrix, localRhs );
}

namespace internal
{
string const faceBcLogMessage = string( "SinglePhaseFVM {}: at time {}s, " )
                                + string( "the <{}> boundary condition '{}' is applied to the face set '{}' in '{}'. " )
                                + string( "\nThe total number of target faces (including ghost faces) is {}. " )
                                + string( "\nNote that if this number is equal to zero, the boundary condition will not be applied on this face set." );
}


template< typename BASE >
void SinglePhaseFVM< BASE >::applyFaceDirichletBC( real64 const time_n,
                                                   real64 const dt,
                                                   DofManager const & dofManager,
                                                   DomainPartition & domain,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  FaceManager & faceManager = mesh.getFaceManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  ConstitutiveManager & constitutiveManager = domain.getConstitutiveManager();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  // make a list of region indices to be included
  map< localIndex, localIndex > regionFluidMap;
  forTargetRegionsComplete( mesh, [&]( localIndex const targetIndex, localIndex const er, ElementRegionBase & )
  {
    localIndex const modelIndex = constitutiveManager.getSubGroups().getIndex( m_fluidModelNames[targetIndex] );
    regionFluidMap.emplace( er, modelIndex );
  } );

  arrayView1d< real64 const > const presFace =
    faceManager.getExtrinsicData< extrinsicMeshData::flow::facePressure >();

  arrayView1d< real64 const > const gravCoefFace =
    faceManager.getExtrinsicData< extrinsicMeshData::flow::gravityCoefficient >();

  string const & dofKey = dofManager.getKey( extrinsicMeshData::flow::pressure::key() );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  elemDofNumber = elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
  elemDofNumber.setName( this->getName() + "/accessors/" + dofKey );

  // Take BCs defined for "pressure" field and apply values to "facePressure"
  fsManager.apply( time_n + dt,
                   domain,
                   "faceManager",
                   extrinsicMeshData::flow::pressure::key(),
                   [&] ( FieldSpecificationBase const & fs,
                         string const & setName,
                         SortedArrayView< localIndex const > const & targetSet,
                         Group & targetGroup,
                         string const & )
  {
    BoundaryStencil const & stencil = fluxApprox.getStencil< BoundaryStencil >( mesh, setName );
    if( fs.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
    {
      globalIndex const numTargetFaces = MpiWrapper::sum< globalIndex >( stencil.size() );
      GEOSX_LOG_RANK_0( GEOSX_FMT( geosx::internal::faceBcLogMessage,
                                   this->getName(), time_n+dt, AquiferBoundaryCondition::catalogName(),
                                   fs.getName(), setName, targetGroup.getName(), numTargetFaces ) );
    }

    if( stencil.size() == 0 )
    {
      return;
    }

    // first, evaluate BC to get primary field values (pressure)
    fs.applyFieldValue< FieldSpecificationEqual, parallelDevicePolicy<> >( targetSet,
                                                                           time_n + dt,
                                                                           targetGroup,
                                                                           extrinsicMeshData::flow::facePressure::key() );

    // Now run the actual kernel
    BoundaryStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
    BoundaryStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
    BoundaryStencil::IndexContainerViewConstType const & sefi = stencil.getElementIndices();
    BoundaryStencil::WeightContainerViewConstType const & trans = stencil.getWeights();

    // TODO: currently we just use model from the first cell in this stencil
    //       since it's not clear how to create fluid kernel wrappers for arbitrary models.
    //       Can we just use cell properties for an approximate flux computation?
    //       Then we can forget about capturing the fluid model.
    SingleFluidBase & fluidBase = constitutiveManager.getConstitutiveRelation< SingleFluidBase >( regionFluidMap[seri( 0, 0 )] );

    constitutiveUpdatePassThru( fluidBase, [&]( auto & fluid )
    {
      // create the fluid compute wrapper suitable for capturing in a kernel lambda
      typename TYPEOFREF( fluid ) ::KernelWrapper fluidWrapper = fluid.createKernelWrapper();

      typename FluxKernel::SinglePhaseFlowAccessors flowAccessors( elemManager, this->getName() );
      typename FluxKernel::SinglePhaseFluidAccessors fluidAccessors( elemManager, this->getName(), this->targetRegionNames(), this->fluidModelNames() );

      FaceDirichletBCKernel::launch( seri, sesri, sefi, trans,
                                     flowAccessors.get< extrinsicMeshData::ghostRank >(),
                                     elemDofNumber.toNestedViewConst(),
                                     dofManager.rankOffset(),
                                     flowAccessors.get< extrinsicMeshData::flow::pressure >(),
                                     flowAccessors.get< extrinsicMeshData::flow::deltaPressure >(),
                                     flowAccessors.get< extrinsicMeshData::flow::gravityCoefficient >(),
                                     fluidAccessors.get< extrinsicMeshData::singlefluid::density >(),
                                     fluidAccessors.get< extrinsicMeshData::singlefluid::dDensity_dPressure >(),
                                     flowAccessors.get< extrinsicMeshData::flow::mobility >(),
                                     flowAccessors.get< extrinsicMeshData::flow::dMobility_dPressure >(),
                                     presFace,
                                     gravCoefFace,
                                     fluidWrapper,
                                     dt,
                                     localMatrix,
                                     localRhs );
    } );
  } );
}

template<>
void SinglePhaseFVM< SinglePhaseProppantBase >::applyAquiferBC( real64 const GEOSX_UNUSED_PARAM( time ),
                                                                real64 const GEOSX_UNUSED_PARAM( dt ),
                                                                DomainPartition & GEOSX_UNUSED_PARAM( domain ),
                                                                DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                                                                CRSMatrixView< real64, globalIndex const > const & GEOSX_UNUSED_PARAM( localMatrix ),
                                                                arrayView1d< real64 > const & GEOSX_UNUSED_PARAM( localRhs ) ) const
{
  // Aquifer does not make sense for proppant flow in fractures
}

template<>
void SinglePhaseFVM< SinglePhaseBase >::applyAquiferBC( real64 const time,
                                                        real64 const dt,
                                                        DomainPartition & domain,
                                                        DofManager const & dofManager,
                                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                        arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager const & elemManager = mesh.getElemManager();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const & elemDofKey = dofManager.getKey( extrinsicMeshData::flow::pressure::key() );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > elemDofNumber =
    elemManager.constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
  elemDofNumber.setName( this->getName() + "/accessors/" + elemDofKey );

  typename FluxKernel::SinglePhaseFlowAccessors flowAccessors( elemManager, this->getName() );
  typename FluxKernel::SinglePhaseFluidAccessors fluidAccessors( elemManager, this->getName(), this->targetRegionNames(), this->fluidModelNames() );

  fsManager.apply< AquiferBoundaryCondition >( time + dt,
                                               domain,
                                               "faceManager",
                                               AquiferBoundaryCondition::catalogName(),
                                               [&] ( AquiferBoundaryCondition const & bc,
                                                     string const & setName,
                                                     SortedArrayView< localIndex const > const &,
                                                     Group & targetGroup,
                                                     string const & )
  {
    BoundaryStencil const & stencil = fluxApprox.getStencil< BoundaryStencil >( mesh, setName );
    if( bc.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
    {
      globalIndex const numTargetFaces = MpiWrapper::sum< globalIndex >( stencil.size() );
      GEOSX_LOG_RANK_0( GEOSX_FMT( geosx::internal::faceBcLogMessage,
                                   this->getName(), time+dt, AquiferBoundaryCondition::catalogName(),
                                   bc.getName(), setName, targetGroup.getName(), numTargetFaces ) );
    }

    if( stencil.size() == 0 )
    {
      return;
    }

    AquiferBoundaryCondition::KernelWrapper aquiferBCWrapper = bc.createKernelWrapper();
    real64 const & aquiferDens = bc.getWaterPhaseDensity();

    SinglePhaseFVMKernels::AquiferBCKernel::launch( stencil,
                                                    dofManager.rankOffset(),
                                                    elemDofNumber.toNestedViewConst(),
                                                    flowAccessors.get< extrinsicMeshData::ghostRank >(),
                                                    aquiferBCWrapper,
                                                    aquiferDens,
                                                    flowAccessors.get< extrinsicMeshData::flow::pressure >(),
                                                    flowAccessors.get< extrinsicMeshData::flow::deltaPressure >(),
                                                    flowAccessors.get< extrinsicMeshData::flow::gravityCoefficient >(),
                                                    fluidAccessors.get< extrinsicMeshData::singlefluid::density >(),
                                                    fluidAccessors.get< extrinsicMeshData::singlefluid::dDensity_dPressure >(),
                                                    time,
                                                    dt,
                                                    localMatrix.toViewConstSizes(),
                                                    localRhs.toView() );

  } );
}


template< typename BASE >
void
SinglePhaseFVM< BASE >::calculateAndApplyMassFlux( real64 const & time_n ,
                                                   real64 const & dt ,
                                                   DomainPartition & domain )
{
	FunctionManager const & functionManager = FunctionManager::getInstance();

	assembleFluxTermsExplicit( time_n, dt, domain );

  // apply mass flux boundary condition
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  fsManager.apply( time_n + dt,
 	 	 	 	   domain,
                   "ElementRegions",
                   FieldSpecificationBase::viewKeyStruct::fluxBoundaryConditionString(),
                   [&]( FieldSpecificationBase const & fs,
                        string const &,
                        SortedArrayView< localIndex const > const & lset,
						Group & subRegion,
                        string const & ) -> void
  {
    arrayView1d< integer const > const
    ghostRank = subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString() );

    SortedArray< localIndex > localSet;
    for( localIndex const a : lset )
    {
      if( ghostRank[a] < 0 )
      {
        localSet.insert( a );
      }
    }
    string const & functionName = fs.getFunctionName();
    arrayView1d< real64 > const mass = subRegion.getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::fluidMassString() );

	  if( functionName.empty() )
	  {
		  real64 value = fs.getScale();
		  forAll< parallelDevicePolicy<> >( localSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
		    {
		      localIndex const kf = localSet[ i ];
		      mass( kf ) += value * dt;
		    } );
	  }
	  else
	  {
		  FunctionBase const & function = functionManager.getGroup< FunctionBase >( functionName );
		  if( function.isFunctionOfTime() == 2 )
		  {
			  real64 value = fs.getScale() * function.evaluate( &time_n );
			  forAll< parallelDevicePolicy<> >( localSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
			    {
			      localIndex const kf = localSet[ i ];
			      mass( kf ) += value * dt;
			    } );
		  }
	  }


  } );

  // apply pressure boundary condition in the explicit solver
  //FieldSpecificationManager & fsManager = FieldSpecificationManager::get();
  fsManager.apply( time_n + dt,
                   domain,
                   "ElementRegions",
				   extrinsicMeshData::flow::pressure::key(),
                    [&]( FieldSpecificationBase const & fs,
                         string const &,
                         SortedArrayView<localIndex const> const & lset,
						 Group & subRegion,
                         string const & ) -> void
  {
    fs.applyFieldValue<FieldSpecificationEqual>( lset,
                                                  time_n + dt,
                                                  subRegion,
												  extrinsicMeshData::flow::pressure::key() );
    arrayView1d< real64 const > const vol = subRegion.getReference< array1d< real64 > >(  CellElementSubRegion::viewKeyStruct::elementVolumeString() );
    arrayView1d< real64 > const mass = subRegion.getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::fluidMassString() );

    //string a = m_fluidModelNames[0];
    Group & region = subRegion.getParent().getParent();
    localIndex targetIndex = SolverBase::targetRegionIndex( region.getName() );

    CompressibleSinglePhaseFluid & fluid = SolverBase::getConstitutiveModel< CompressibleSinglePhaseFluid>( subRegion, m_fluidModelNames[ targetIndex ] );
    real64 referenceDensity = fluid.referenceDensity();
    real64 referencePressure = fluid.referencePressure();
    real64 compressibility = fluid.compressibility();

    CoupledSolidBase const & solidModel = SolverBase::getConstitutiveModel< CoupledSolidBase > ( subRegion, m_solidModelNames[targetIndex] );
    arrayView2d< real64 const > const & porosity = solidModel.getPorosity();
    for( localIndex const a : lset )
    {
    	mass[a] = referenceDensity / (1.0 - (fs.getScale() -referencePressure) * compressibility) * vol[a] * porosity[a][0];
    }
  });
  // synchronize element fields
  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( SinglePhaseBase::viewKeyStruct::fluidMassString() ) );
  CommunicationTools::getInstance().synchronizeFields( fieldNames,
                                                       domain.getMeshBody( 0 ).getMeshLevel( 0 ),
                                                       domain.getNeighbors(),
                                                       true );

}

namespace
{
typedef SinglePhaseFVM< SinglePhaseBase > NoProppant;
typedef SinglePhaseFVM< SinglePhaseProppantBase > Proppant;
REGISTER_CATALOG_ENTRY( SolverBase, NoProppant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( SolverBase, Proppant, string const &, Group * const )
}
} /* namespace geosx */
