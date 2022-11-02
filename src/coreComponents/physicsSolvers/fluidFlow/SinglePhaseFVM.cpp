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

#include "common/TimingMacros.hpp"
#include "constitutive/fluid/singleFluidSelector.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/ThermalSinglePhaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseFVMKernels.hpp"
#include "physicsSolvers/fluidFlow/ThermalSinglePhaseFVMKernels.hpp"
#include "physicsSolvers/multiphysics/SinglePhasePoromechanicsFluxKernels.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace singlePhaseBaseKernels;
using namespace singlePhaseFVMKernels;
using namespace singlePhasePoromechanicsFluxKernels;

template< typename BASE >
SinglePhaseFVM< BASE >::SinglePhaseFVM( const string & name,
                                        Group * const parent ):
  BASE( name, parent )
{}

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
  dofManager.addField( BASE::viewKeyStruct::elemDofFieldString(),
                       FieldLocation::Elem,
                       m_numDofPerCell,
                       BASE::getMeshTargets() );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  dofManager.addCoupling( BASE::viewKeyStruct::elemDofFieldString(), fluxApprox );
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
real64 SinglePhaseFVM< BASE >::calculateResidualNorm( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                                      real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                      DomainPartition const & domain,
                                                      DofManager const & dofManager,
                                                      arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  integer const numNorm = m_isThermal ? 2 : 1;
  array1d< real64 > localResidualNorm;
  array1d< real64 > localResidualNormalizer;
  localResidualNorm.resize( numNorm );
  localResidualNormalizer.resize( numNorm );

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( BASE::viewKeyStruct::elemDofFieldString() );

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel const & mesh,
                                                                      arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      real64 subRegionResidualNorm[2]{};
      real64 subRegionResidualNormalizer[2]{};

      string const & fluidName = subRegion.template getReference< string >( BASE::viewKeyStruct::fluidNamesString() );
      SingleFluidBase const & fluid = SolverBase::getConstitutiveModel< SingleFluidBase >( subRegion, fluidName );

      string const & solidName = subRegion.template getReference< string >( BASE::viewKeyStruct::solidNamesString() );
      CoupledSolidBase const & solid = SolverBase::getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

      // step 1: compute the norm in the subRegion

      if( m_isThermal )
      {
        string const & solidInternalEnergyName = subRegion.template getReference< string >( BASE::viewKeyStruct::solidInternalEnergyNamesString() );
        SolidInternalEnergy const & solidInternalEnergy = SolverBase::getConstitutiveModel< SolidInternalEnergy >( subRegion, solidInternalEnergyName );

        thermalSinglePhaseBaseKernels::
          ResidualNormKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( BASE::m_normType,
                                                     rankOffset,
                                                     dofKey,
                                                     localRhs,
                                                     subRegion,
                                                     fluid,
                                                     solid,
                                                     solidInternalEnergy,
                                                     subRegionResidualNorm,
                                                     subRegionResidualNormalizer );
      }
      else
      {
        real64 subRegionFlowResidualNorm[1]{};
        real64 subRegionFlowResidualNormalizer[1]{};
        singlePhaseBaseKernels::
          ResidualNormKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( BASE::m_normType,
                                                     rankOffset,
                                                     dofKey,
                                                     localRhs,
                                                     subRegion,
                                                     fluid,
                                                     solid,
                                                     subRegionFlowResidualNorm,
                                                     subRegionFlowResidualNormalizer );
        subRegionResidualNorm[0] = subRegionFlowResidualNorm[0];
        subRegionResidualNormalizer[0] = subRegionFlowResidualNormalizer[0];
      }

      // step 2: first reduction across meshBodies/regions/subRegions

      if( BASE::m_normType == solverBaseKernels::NormType::Linf )
      {
        solverBaseKernels::LinfResidualNormHelper::updateLocalNorm< 2 >( subRegionResidualNorm, localResidualNorm );
      }
      else
      {
        solverBaseKernels::L2ResidualNormHelper::
          updateLocalNorm< 2 >( subRegionResidualNorm, subRegionResidualNormalizer, localResidualNorm, localResidualNormalizer );
      }
    } );
  } );

  // step 3: second reduction across MPI ranks

  real64 residualNorm = 0.0;
  if( m_isThermal )
  {

    array1d< real64 > globalResidualNorm;
    if( BASE::m_normType == solverBaseKernels::NormType::Linf )
    {
      solverBaseKernels::LinfResidualNormHelper::computeGlobalNorm( localResidualNorm, globalResidualNorm );
    }
    else
    {
      solverBaseKernels::L2ResidualNormHelper::computeGlobalNorm( localResidualNorm, localResidualNormalizer, globalResidualNorm );
    }
    residualNorm = sqrt( globalResidualNorm[0] * globalResidualNorm[0] + globalResidualNorm[1] * globalResidualNorm[1] );

    if( getLogLevel() >= 1 && logger::internal::rank == 0 )
    {
      std::cout << GEOSX_FMT( "    ( R{} ) = ( {:4.2e} ) ; ( Renergy ) = ( {:4.2e} ) ; ",
                              FlowSolverBase::coupledSolverAttributePrefix(), globalResidualNorm[0], globalResidualNorm[1] );
    }
  }
  else
  {

    if( BASE::m_normType == solverBaseKernels::NormType::Linf )
    {
      solverBaseKernels::LinfResidualNormHelper::computeGlobalNorm( localResidualNorm[0], residualNorm );
    }
    else
    {
      solverBaseKernels::L2ResidualNormHelper::computeGlobalNorm( localResidualNorm[0], localResidualNormalizer[0], residualNorm );
    }

    if( getLogLevel() >= 1 && logger::internal::rank == 0 )
    {
      std::cout << GEOSX_FMT( "    ( R{} ) = ( {:4.2e} ) ; ", FlowSolverBase::coupledSolverAttributePrefix(), residualNorm );
    }
  }
  return residualNorm;
}


template< typename BASE >
void SinglePhaseFVM< BASE >::applySystemSolution( DofManager const & dofManager,
                                                  arrayView1d< real64 const > const & localSolution,
                                                  real64 const scalingFactor,
                                                  DomainPartition & domain )
{
  if( m_isThermal )
  {
    DofManager::CompMask pressureMask( m_numDofPerCell, 0, 1 );
    DofManager::CompMask temperatureMask( m_numDofPerCell, 1, 2 );

    dofManager.addVectorToField( localSolution,
                                 BASE::viewKeyStruct::elemDofFieldString(),
                                 fields::flow::pressure::key(),
                                 scalingFactor,
                                 pressureMask );

    dofManager.addVectorToField( localSolution,
                                 BASE::viewKeyStruct::elemDofFieldString(),
                                 fields::flow::temperature::key(),
                                 scalingFactor,
                                 temperatureMask );
  }
  else
  {
    dofManager.addVectorToField( localSolution,
                                 BASE::viewKeyStruct::elemDofFieldString(),
                                 fields::flow::pressure::key(),
                                 scalingFactor );
  }

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel & mesh,
                                                                      arrayView1d< string const > const & regionNames )
  {
    std::vector< string > fields{ fields::flow::pressure::key() };

    if( m_isThermal )
    {
      fields.emplace_back( fields::flow::temperature::key() );
    }

    FieldIdentifiers fieldsToBeSync;

    fieldsToBeSync.addElementFields( fields, regionNames );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), true );
  } );
}

template< >
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

  string const & dofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & )
  {
    fluxApprox.forAllStencils( mesh, [&]( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();


      if( m_isThermal )
      {
        thermalSinglePhaseFVMKernels::
          FaceBasedAssemblyKernelFactory::createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
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
        singlePhaseFVMKernels::
          FaceBasedAssemblyKernelFactory::createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
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

  string const & dofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
    elemDofNumber = elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    elemDofNumber.setName( this->getName() + "/accessors/" + dofKey );

    fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

      typename FluxKernel::SinglePhaseFlowAccessors flowAccessors( elemManager, getName() );
      typename FluxKernel::SlurryFluidAccessors fluidAccessors( elemManager, getName() );
      typename FluxKernel::ProppantPermeabilityAccessors permAccessors( elemManager, getName() );

      FaceElementFluxKernel::launch( stencilWrapper,
                                     dt,
                                     dofManager.rankOffset(),
                                     elemDofNumber.toNestedViewConst(),
                                     flowAccessors.get< fields::ghostRank >(),
                                     flowAccessors.get< fields::flow::pressure >(),
                                     flowAccessors.get< fields::flow::gravityCoefficient >(),
                                     fluidAccessors.get< fields::singlefluid::density >(),
                                     fluidAccessors.get< fields::singlefluid::dDensity_dPressure >(),
                                     flowAccessors.get< fields::flow::mobility >(),
                                     flowAccessors.get< fields::flow::dMobility_dPressure >(),
                                     permAccessors.get< fields::permeability::permeability >(),
                                     permAccessors.get< fields::permeability::dPerm_dPressure >(),
                                     permAccessors.get< fields::permeability::dPerm_dDispJump >(),
                                     permAccessors.get< fields::permeability::permeabilityMultiplier >(),
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

  string const & pressureDofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel const & mesh,
                                                                      arrayView1d< string const > const & )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
    pressureDofNumber = elemManager.constructArrayViewAccessor< globalIndex, 1 >( pressureDofKey );
    pressureDofNumber.setName( this->getName() + "/accessors/" + pressureDofKey );

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
    jumpDofNumber = elemManager.constructArrayViewAccessor< globalIndex, 1 >( jumpDofKey );
    jumpDofNumber.setName( this->getName() + "/accessors/" + jumpDofKey );

    ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const > > dPerm_dDispJump =
      elemManager.constructMaterialArrayViewAccessor< PermeabilityBase, real64, 4 >( fields::permeability::dPerm_dDispJump::key() );

    fluxApprox.forStencils< CellElementStencilTPFA, SurfaceElementStencil, EmbeddedSurfaceToCellStencil >( mesh, [&]( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

      typename FluxKernel::SinglePhaseFlowAccessors flowAccessors( elemManager, this->getName() );
      typename FluxKernel::SinglePhaseFluidAccessors fluidAccessors( elemManager, this->getName() );
      typename FluxKernel::PermeabilityAccessors permAccessors( elemManager, this->getName() );

      EmbeddedSurfaceFluxKernel::launch( stencilWrapper,
                                         dt,
                                         dofManager.rankOffset(),
                                         pressureDofNumber.toNestedViewConst(),
                                         jumpDofNumber.toNestedViewConst(),
                                         flowAccessors.get< fields::ghostRank >(),
                                         flowAccessors.get< fields::flow::pressure >(),
                                         flowAccessors.get< fields::flow::gravityCoefficient >(),
                                         fluidAccessors.get< fields::singlefluid::density >(),
                                         fluidAccessors.get< fields::singlefluid::dDensity_dPressure >(),
                                         flowAccessors.get< fields::flow::mobility >(),
                                         flowAccessors.get< fields::flow::dMobility_dPressure >(),
                                         permAccessors.get< fields::permeability::permeability >(),
                                         permAccessors.get< fields::permeability::dPerm_dPressure >(),
                                         dPerm_dDispJump.toNestedViewConst(),
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

  string const & dofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );


  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel const & mesh,
                                                                      arrayView1d< string const > const & )
  {
    std::cout<<mesh.getName()<<std::endl;

    ElementRegionManager const & elemManager = mesh.getElemManager();
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
    elemDofNumber = elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    elemDofNumber.setName( this->getName() + "/accessors/" + dofKey );

    ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const > > dPerm_dDispJump =
      elemManager.constructMaterialArrayViewAccessor< PermeabilityBase, real64, 4 >( fields::permeability::dPerm_dDispJump::key() );

    fluxApprox.forStencils< CellElementStencilTPFA, SurfaceElementStencil, FaceElementToCellStencil >( mesh, [&]( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

      typename FluxKernel::SinglePhaseFlowAccessors flowAccessors( elemManager, this->getName() );
      typename FluxKernel::SinglePhaseFluidAccessors fluidAccessors( elemManager, this->getName() );
      typename FluxKernel::PermeabilityAccessors permAccessors( elemManager, this->getName() );

      FaceElementFluxKernel::launch( stencilWrapper,
                                     dt,
                                     dofManager.rankOffset(),
                                     elemDofNumber.toNestedViewConst(),
                                     flowAccessors.get< fields::ghostRank >(),
                                     flowAccessors.get< fields::flow::pressure >(),
                                     flowAccessors.get< fields::flow::gravityCoefficient >(),
                                     fluidAccessors.get< fields::singlefluid::density >(),
                                     fluidAccessors.get< fields::singlefluid::dDensity_dPressure >(),
                                     flowAccessors.get< fields::flow::mobility >(),
                                     flowAccessors.get< fields::flow::dMobility_dPressure >(),
                                     permAccessors.get< fields::permeability::permeability >(),
                                     permAccessors.get< fields::permeability::dPerm_dPressure >(),
                                     dPerm_dDispJump.toNestedViewConst(),
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

namespace
{
char const faceBcLogMessage[] =
  "SinglePhaseFVM {}: at time {}s, "
  "the <{}> boundary condition '{}' is applied to the face set '{}' in '{}'. "
  "\nThe total number of target faces (including ghost faces) is {}. "
  "\nNote that if this number is equal to zero, the boundary condition will not be applied on this face set.";
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

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const & dofKey = dofManager.getKey( BASE::viewKeyStruct::elemDofFieldString() );

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel & mesh,
                                                                      arrayView1d< string const > const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    arrayView1d< real64 const > const presFace =
      faceManager.getField< fields::flow::facePressure >();

    arrayView1d< real64 const > const gravCoefFace =
      faceManager.getField< fields::flow::gravityCoefficient >();

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
    elemDofNumber = elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    elemDofNumber.setName( this->getName() + "/accessors/" + dofKey );

    // Take BCs defined for "pressure" field and apply values to "facePressure"
    fsManager.apply< FaceManager >( time_n + dt,
                                    mesh,
                                    fields::flow::pressure::key(),
                                    [&] ( FieldSpecificationBase const & fs,
                                          string const & setName,
                                          SortedArrayView< localIndex const > const & targetSet,
                                          FaceManager & targetGroup,
                                          string const & )
    {
      BoundaryStencil const & stencil = fluxApprox.getStencil< BoundaryStencil >( mesh, setName );

      if( fs.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
      {
        globalIndex const numTargetFaces = MpiWrapper::sum< globalIndex >( stencil.size() );
        GEOSX_LOG_RANK_0( GEOSX_FMT( faceBcLogMessage,
                                     this->getName(), time_n+dt, FieldSpecificationBase::catalogName(),
                                     fs.getName(), setName, targetGroup.getName(), numTargetFaces ) );
      }

      if( stencil.size() == 0 )
      {
        return;
      }

      // first, evaluate BC to get primary field values (pressure)
      fs.applyFieldValue< FieldSpecificationEqual,
                          parallelDevicePolicy<> >( targetSet,
                                                    time_n + dt,
                                                    targetGroup,
                                                    fields::flow::facePressure::key() );

      // TODO: currently we just use model from the first cell in this stencil
      //       since it's not clear how to create fluid kernel wrappers for arbitrary models.
      //       Can we just use cell properties for an approximate flux computation?
      //       Then we can forget about capturing the fluid model.
      localIndex const er = stencil.getElementRegionIndices()( 0, 0 );
      localIndex const esr = stencil.getElementSubRegionIndices()( 0, 0 );
      ElementSubRegionBase & subRegion = mesh.getElemManager().getRegion( er ).getSubRegion( esr );
      string const & fluidName = subRegion.getReference< string >( BASE::viewKeyStruct::fluidNamesString() );
      SingleFluidBase & fluidBase = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );

      constitutiveUpdatePassThru( fluidBase, [&]( auto & fluid )
      {
        typename TYPEOFREF( fluid ) ::KernelWrapper fluidWrapper = fluid.createKernelWrapper();

        typename FluxKernel::SinglePhaseFlowAccessors flowAccessors( elemManager, this->getName() );
        typename FluxKernel::SinglePhaseFluidAccessors fluidAccessors( elemManager, this->getName() );
        typename FluxKernel::PermeabilityAccessors permAccessors( elemManager, this->getName() );

        FaceDirichletBCKernel::launch( stencil.createKernelWrapper(),
                                       flowAccessors.get< fields::ghostRank >(),
                                       elemDofNumber.toNestedViewConst(),
                                       dofManager.rankOffset(),
                                       permAccessors.get< fields::permeability::permeability >(),
                                       permAccessors.get< fields::permeability::dPerm_dPressure >(),
                                       flowAccessors.get< fields::flow::pressure >(),
                                       flowAccessors.get< fields::flow::gravityCoefficient >(),
                                       fluidAccessors.get< fields::singlefluid::density >(),
                                       fluidAccessors.get< fields::singlefluid::dDensity_dPressure >(),
                                       flowAccessors.get< fields::flow::mobility >(),
                                       flowAccessors.get< fields::flow::dMobility_dPressure >(),
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

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const & elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > elemDofNumber =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
    elemDofNumber.setName( this->getName() + "/accessors/" + elemDofKey );

    typename FluxKernel::SinglePhaseFlowAccessors flowAccessors( elemManager, this->getName() );
    typename FluxKernel::SinglePhaseFluidAccessors fluidAccessors( elemManager, this->getName() );

    fsManager.apply< FaceManager,
                     AquiferBoundaryCondition >( time + dt,
                                                 mesh,
                                                 AquiferBoundaryCondition::catalogName(),
                                                 [&] ( AquiferBoundaryCondition const & bc,
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

      AquiferBoundaryCondition::KernelWrapper aquiferBCWrapper = bc.createKernelWrapper();
      real64 const & aquiferDens = bc.getWaterPhaseDensity();

      singlePhaseFVMKernels::AquiferBCKernel::launch( stencil,
                                                      dofManager.rankOffset(),
                                                      elemDofNumber.toNestedViewConst(),
                                                      flowAccessors.get< fields::ghostRank >(),
                                                      aquiferBCWrapper,
                                                      aquiferDens,
                                                      flowAccessors.get< fields::flow::pressure >(),
                                                      flowAccessors.get< fields::flow::pressure_n >(),
                                                      flowAccessors.get< fields::flow::gravityCoefficient >(),
                                                      fluidAccessors.get< fields::singlefluid::density >(),
                                                      fluidAccessors.get< fields::singlefluid::dDensity_dPressure >(),
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
