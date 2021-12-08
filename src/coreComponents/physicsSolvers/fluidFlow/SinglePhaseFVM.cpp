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
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
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
  dofManager.addField( BASE::viewKeyStruct::pressureString(),
                       DofManager::Location::Elem,
                       1,
                       BASE::m_meshTargets );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  dofManager.addCoupling( BASE::viewKeyStruct::pressureString(), fluxApprox );
}

template< typename BASE >
void SinglePhaseFVM< BASE >::setupSystem( DomainPartition & domain,
                                          DofManager & dofManager,
                                          CRSMatrix< real64, globalIndex > & localMatrix,
                                          array1d< real64 > & localRhs,
                                          array1d< real64 > & localSolution,
                                          bool const setSparsity )
{
  GEOSX_MARK_FUNCTION;
  BASE::setupSystem( domain,
                     dofManager,
                     localMatrix,
                     localRhs,
                     localSolution,
                     setSparsity );

}

template< typename BASE >
real64 SinglePhaseFVM< BASE >::calculateResidualNorm( DomainPartition const & domain,
                                                      DofManager const & dofManager,
                                                      arrayView1d< real64 const > const & localRhs )
{
  real64 residual = 0;
  integer numMeshTargets = 0;
  string const dofKey = dofManager.getKey( BASE::viewKeyStruct::pressureString() );
  globalIndex const rankOffset = dofManager.rankOffset();
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    real64 localResidualNorm[3] = { 0.0, 0.0, 0.0 };
    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&] ( localIndex const,
                                                                    auto const & subRegion )
    {
      arrayView1d< globalIndex const > const & dofNumber = subRegion.template getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const & volume         = subRegion.getElementVolume();
      arrayView1d< real64 const > const & densOld        = subRegion.template getReference< array1d< real64 > >( BASE::viewKeyStruct::densityOldString() );

      string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
      CoupledSolidBase const & solidModel = getConstitutiveModel< CoupledSolidBase >( subRegion, solidMaterialName );

      arrayView2d< real64 const > const & porosityOld = solidModel.getOldPorosity();

      ResidualNormKernel::launch< parallelDevicePolicy<>, parallelDeviceReduce >( localRhs,
                                                                                  rankOffset,
                                                                                  dofNumber,
                                                                                  elemGhostRank,
                                                                                  volume,
                                                                                  densOld,
                                                                                  porosityOld,
                                                                                  localResidualNorm );

    } ) :

      // compute global residual norm
      real64 globalResidualNorm[3] = {0, 0, 0};
    MpiWrapper::allReduce( localResidualNorm,
                           globalResidualNorm,
                           3,
                           MPI_SUM,
                           MPI_COMM_GEOSX );

    residual += sqrt( globalResidualNorm[0] ) / ( ( globalResidualNorm[1] + m_fluxEstimate ) / (globalResidualNorm[2]+1) );
    numMeshTargets++
  } );

  return residual / numMeshTargets;
}


template< typename BASE >
void SinglePhaseFVM< BASE >::applySystemSolution( DofManager const & dofManager,
                                                  arrayView1d< real64 const > const & localSolution,
                                                  real64 const scalingFactor,
                                                  DomainPartition & domain )
{
  dofManager.addVectorToField( localSolution,
                               BASE::viewKeyStruct::pressureString(),
                               BASE::viewKeyStruct::deltaPressureString(),
                               scalingFactor );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & )
  {
    std::map< string, string_array > fieldNames;
    fieldNames["elems"].emplace_back( string( BASE::viewKeyStruct::deltaPressureString() ) );

    CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), true );
  } );
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

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const & dofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::pressureString() );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & )
  {
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
    elemDofNumber = mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    elemDofNumber.setName( this->getName() + "/accessors/" + dofKey );

    fluxApprox.forAllStencils( mesh, [&]( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::StencilWrapper stencilWrapper = stencil.createStencilWrapper();

      FluxKernel::launch( stencilWrapper,
                          dt,
                          dofManager.rankOffset(),
                          elemDofNumber.toNestedViewConst(),
                          m_elemGhostRank.toNestedViewConst(),
                          m_pressure.toNestedViewConst(),
                          m_deltaPressure.toNestedViewConst(),
                          m_gravCoef.toNestedViewConst(),
                          m_density.toNestedViewConst(),
                          m_dDens_dPres.toNestedViewConst(),
                          m_mobility.toNestedViewConst(),
                          m_dMobility_dPres.toNestedViewConst(),
                          m_permeability.toNestedViewConst(),
                          m_dPerm_dPressure.toNestedViewConst(),
                          localMatrix,
                          localRhs );

    } );

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

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const & dofKey = dofManager.getKey( SinglePhaseProppantBase::viewKeyStruct::pressureString() );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & )
  {
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
    elemDofNumber = mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    elemDofNumber.setName( this->getName() + "/accessors/" + dofKey );

    ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > > dPerm_dAper =
      mesh.getElemManager().constructMaterialArrayViewAccessor< PermeabilityBase, real64, 3 >( PermeabilityBase::viewKeyStruct::dPerm_dApertureString() );

    fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::StencilWrapper stencilWrapper = stencil.createStencilWrapper();

      FaceElementFluxKernel::launch( stencilWrapper,
                                     dt,
                                     dofManager.rankOffset(),
                                     elemDofNumber.toNestedViewConst(),
                                     m_elemGhostRank.toNestedViewConst(),
                                     m_pressure.toNestedViewConst(),
                                     m_deltaPressure.toNestedViewConst(),
                                     m_gravCoef.toNestedViewConst(),
                                     m_density.toNestedViewConst(),
                                     m_dDens_dPres.toNestedViewConst(),
                                     m_mobility.toNestedViewConst(),
                                     m_dMobility_dPres.toNestedViewConst(),
                                     m_permeability.toNestedViewConst(),
                                     m_dPerm_dPressure.toNestedViewConst(),
                                     dPerm_dAper.toNestedViewConst(),
                                     SinglePhaseProppantBase::m_permeabilityMultiplier.toNestedViewConst(),
                                     this->gravityVector(),
                                     localMatrix,
                                     localRhs );
    } );
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

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const & pressureDofKey = dofManager.getKey( BASE::viewKeyStruct::pressureString() );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & )
  {
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
    pressureDofNumber = mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( pressureDofKey );
    pressureDofNumber.setName( this->getName() + "/accessors/" + pressureDofKey );

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
    jumpDofNumber = mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( jumpDofKey );
    jumpDofNumber.setName( this->getName() + "/accessors/" + jumpDofKey );

    ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > > dPerm_dAper =
      mesh.getElemManager().constructMaterialArrayViewAccessor< PermeabilityBase, real64, 3 >( PermeabilityBase::viewKeyStruct::dPerm_dApertureString() );


    fluxApprox.forStencils< CellElementStencilTPFA, SurfaceElementStencil, EmbeddedSurfaceToCellStencil >( mesh, [&]( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::StencilWrapper stencilWrapper = stencil.createStencilWrapper();

      EmbeddedSurfaceFluxKernel::launch( stencilWrapper,
                                         dt,
                                         dofManager.rankOffset(),
                                         pressureDofNumber.toNestedViewConst(),
                                         jumpDofNumber.toNestedViewConst(),
                                         m_elemGhostRank.toNestedViewConst(),
                                         m_pressure.toNestedViewConst(),
                                         m_deltaPressure.toNestedViewConst(),
                                         m_gravCoef.toNestedViewConst(),
                                         m_density.toNestedViewConst(),
                                         m_dDens_dPres.toNestedViewConst(),
                                         m_mobility.toNestedViewConst(),
                                         m_dMobility_dPres.toNestedViewConst(),
                                         m_permeability.toNestedViewConst(),
                                         m_dPerm_dPressure.toNestedViewConst(),
                                         dPerm_dAper.toNestedViewConst(),
                                         localMatrix,
                                         localRhs );
    } );

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

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const & dofKey = dofManager.getKey( BASE::viewKeyStruct::pressureString() );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & )
  {
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
    elemDofNumber = mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    elemDofNumber.setName( this->getName() + "/accessors/" + dofKey );

    ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > > dPerm_dAper =
      mesh.getElemManager().constructMaterialArrayViewAccessor< PermeabilityBase, real64, 3 >( PermeabilityBase::viewKeyStruct::dPerm_dApertureString() );

    fluxApprox.forStencils< CellElementStencilTPFA, SurfaceElementStencil, FaceElementToCellStencil >( mesh, [&]( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::StencilWrapper stencilWrapper = stencil.createStencilWrapper();

      FaceElementFluxKernel::launch( stencilWrapper,
                                     dt,
                                     dofManager.rankOffset(),
                                     elemDofNumber.toNestedViewConst(),
                                     m_elemGhostRank.toNestedViewConst(),
                                     m_pressure.toNestedViewConst(),
                                     m_deltaPressure.toNestedViewConst(),
                                     m_gravCoef.toNestedViewConst(),
                                     m_density.toNestedViewConst(),
                                     m_dDens_dPres.toNestedViewConst(),
                                     m_mobility.toNestedViewConst(),
                                     m_dMobility_dPres.toNestedViewConst(),
                                     m_permeability.toNestedViewConst(),
                                     m_dPerm_dPressure.toNestedViewConst(),
                                     dPerm_dAper.toNestedViewConst(),
                                     localMatrix,
                                     localRhs,
                                     dR_dAper );
    } );
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
  ConstitutiveManager & constitutiveManager = domain.getConstitutiveManager();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const & dofKey = dofManager.getKey( BASE::viewKeyStruct::pressureString() );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();

    arrayView1d< real64 const > const presFace =
      faceManager.getReference< array1d< real64 > >( BASE::viewKeyStruct::facePressureString() );

    arrayView1d< real64 const > const gravCoefFace =
      faceManager.getReference< array1d< real64 > >( BASE::viewKeyStruct::gravityCoefString() );

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
    elemDofNumber = mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    elemDofNumber.setName( this->getName() + "/accessors/" + dofKey );

    // Take BCs defined for "pressure" field and apply values to "facePressure"
    fsManager.apply( time_n + dt,
                     domain,
                     "faceManager",
                     BASE::viewKeyStruct::pressureString(),
                     [&] ( FieldSpecificationBase const & fs,
                           string const & setName,
                           SortedArrayView< localIndex const > const & targetSet,
                           Group & targetGroup,
                           string const & )
    {
      BoundaryStencil const & stencil = fluxApprox.getStencil< BoundaryStencil >( mesh, setName );
      if( stencil.size() == 0 )
      {
        return;
      }

      // first, evaluate BC to get primary field values (pressure)
      fs.applyFieldValue< FieldSpecificationEqual, parallelDevicePolicy<> >( targetSet,
                                                                             time_n + dt,
                                                                             targetGroup,
                                                                             BASE::viewKeyStruct::facePressureString() );

      // Now run the actual kernel
      BoundaryStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
      BoundaryStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
      BoundaryStencil::IndexContainerViewConstType const & sefi = stencil.getElementIndices();
      BoundaryStencil::WeightContainerViewConstType const & trans = stencil.getWeights();

      // TODO: currently we just use model from the first cell in this stencil
      //       since it's not clear how to create fluid kernel wrappers for arbitrary models.
      //       Can we just use cell properties for an approximate flux computation?
      //       Then we can forget about capturing the fluid model.
      ElementSubRegionBase & subRegion = mesh.getElemManager().getRegion( seri( 0, 0 ) ).getSubRegion( sesri( 0, 0 ) );

      string & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );

      SingleFluidBase & fluidBase = subRegion.template getConstitutiveModel< SingleFluidBase >( fluidName );

      constitutiveUpdatePassThru( fluidBase, [&]( auto & fluid )
      {
        // create the fluid compute wrapper suitable for capturing in a kernel lambda
        typename TYPEOFREF( fluid ) ::KernelWrapper fluidWrapper = fluid.createKernelWrapper();

        FaceDirichletBCKernel::launch( seri, sesri, sefi, trans,
                                       m_elemGhostRank.toNestedViewConst(),
                                       elemDofNumber.toNestedViewConst(),
                                       dofManager.rankOffset(),
                                       m_pressure.toNestedViewConst(),
                                       m_deltaPressure.toNestedViewConst(),
                                       m_gravCoef.toNestedViewConst(),
                                       m_density.toNestedViewConst(),
                                       m_dDens_dPres.toNestedViewConst(),
                                       m_mobility.toNestedViewConst(),
                                       m_dMobility_dPres.toNestedViewConst(),
                                       presFace,
                                       gravCoefFace,
                                       fluidWrapper,
                                       dt,
                                       localMatrix,
                                       localRhs );
      } );
    } );
  } );


}

template< typename BASE >
void SinglePhaseFVM< BASE >::applyAquiferBC( real64 const time,
                                             real64 const dt,
                                             DomainPartition & domain,
                                             DofManager const & dofManager,
                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                             arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const & elemDofKey = dofManager.getKey( BASE::viewKeyStruct::pressureString() );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & )
  {
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > elemDofNumber =
      mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
    elemDofNumber.setName( this->getName() + "/accessors/" + elemDofKey );

    fsManager.apply< AquiferBoundaryCondition >( time + dt,
                                                 domain,
                                                 "faceManager",
                                                 AquiferBoundaryCondition::catalogName(),
                                                 [&] ( AquiferBoundaryCondition const & bc,
                                                       string const & setName,
                                                       SortedArrayView< localIndex const > const &,
                                                       Group &,
                                                       string const & )
    {
      BoundaryStencil const & stencil = fluxApprox.getStencil< BoundaryStencil >( mesh, setName );
      if( stencil.size() == 0 )
      {
        return;
      }

      AquiferBoundaryCondition::KernelWrapper aquiferBCWrapper = bc.createKernelWrapper();
      real64 const & aquiferDens = bc.getWaterPhaseDensity();

      SinglePhaseFVMKernels::AquiferBCKernel::launch( stencil,
                                                      dofManager.rankOffset(),
                                                      elemDofNumber.toNestedViewConst(),
                                                      m_elemGhostRank.toNestedViewConst(),
                                                      aquiferBCWrapper,
                                                      aquiferDens,
                                                      m_pressure.toNestedViewConst(),
                                                      m_deltaPressure.toNestedViewConst(),
                                                      m_gravCoef.toNestedViewConst(),
                                                      m_density.toNestedViewConst(),
                                                      m_dDens_dPres.toNestedViewConst(),
                                                      time,
                                                      dt,
                                                      localMatrix.toViewConstSizes(),
                                                      localRhs.toView() );

    } );
  } );


}

namespace
{
typedef SinglePhaseFVM< SinglePhaseBase > NoProppant;
typedef SinglePhaseFVM< SinglePhaseProppantBase > Proppant;
REGISTER_CATALOG_ENTRY( SolverBase, NoProppant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( SolverBase, Proppant, string const &, Group * const )
}
} /* namespace geosx */
