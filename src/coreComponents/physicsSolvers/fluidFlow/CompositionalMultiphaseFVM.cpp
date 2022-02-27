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
 * @file CompositionalMultiphaseFVM.cpp
 */

#include "CompositionalMultiphaseFVM.hpp"

#include "common/DataTypes.hpp"
#include "common/MpiWrapper.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "dataRepository/Group.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/ThermalCompositionalMultiphaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseFVMKernels.hpp"
#include "physicsSolvers/fluidFlow/ThermalCompositionalMultiphaseFVMKernels.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

CompositionalMultiphaseFVM::CompositionalMultiphaseFVM( const string & name,
                                                        Group * const parent )
  :
  CompositionalMultiphaseBase( name, parent )
{
  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseFVM;
}

void CompositionalMultiphaseFVM::initializePreSubGroups()
{
  CompositionalMultiphaseBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  if( !fvManager.hasGroup< FluxApproximationBase >( m_discretizationName ) )
  {
    GEOSX_ERROR( "A discretization deriving from FluxApproximationBase must be selected with CompositionalMultiphaseFlow" );
  }

}

void CompositionalMultiphaseFVM::setupDofs( DomainPartition const & domain,
                                            DofManager & dofManager ) const
{
  dofManager.addField( viewKeyStruct::elemDofFieldString(),
                       DofManager::Location::Elem,
                       m_numDofPerCell,
                       m_meshTargets );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
  dofManager.addCoupling( viewKeyStruct::elemDofFieldString(), fluxApprox );
}


void CompositionalMultiphaseFVM::assembleFluxTerms( real64 const dt,
                                                    DomainPartition const & domain,
                                                    DofManager const & dofManager,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel const & mesh,
                                               arrayView1d< string const > const & )
  {
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

    string const & elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > elemDofNumber =
      mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
    elemDofNumber.setName( getName() + "/accessors/" + elemDofKey );

    fluxApprox.forAllStencils( mesh, [&] ( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

      if( m_isThermal )
      {
        thermalCompositionalMultiphaseFVMKernels::
          FaceBasedAssemblyKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                     m_numPhases,
                                                     dofManager.rankOffset(),
                                                     elemDofKey,
                                                     m_capPressureFlag,
                                                     getName(),
                                                     mesh.getElemManager(),
                                                     stencilWrapper,
                                                     dt,
                                                     localMatrix.toViewConstSizes(),
                                                     localRhs.toView() );
      }
      else
      {
        isothermalCompositionalMultiphaseFVMKernels::
          FaceBasedAssemblyKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                     m_numPhases,
                                                     dofManager.rankOffset(),
                                                     elemDofKey,
                                                     m_capPressureFlag,
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

void CompositionalMultiphaseFVM::implicitStepComplete( real64 const & time,
                                                       real64 const & dt,
                                                       DomainPartition & domain )
{
  CompositionalMultiphaseBase::implicitStepComplete( time, dt, domain );

  if( m_computeCFLNumbers )
  {
    computeCFLNumbers( dt, domain );
  }
}

void CompositionalMultiphaseFVM::computeCFLNumbers( real64 const & dt,
                                                    DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  // Step 1: reset the arrays involved in the computation of CFL numbers
  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView2d< real64, compflow::USD_PHASE > const & phaseOutflux =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::phaseOutflux >();
      arrayView2d< real64, compflow::USD_COMP > const & compOutflux =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::componentOutflux >();
      phaseOutflux.zero();
      compOutflux.zero();
    } );

    // Step 2: compute the total volumetric outflux for each reservoir cell by looping over faces
    NumericalMethodsManager & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

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
                                                   arrayView2d< real64, compflow::USD_PHASE > >( extrinsicMeshData::flow::phaseOutflux::key() );

    ElementRegionManager::ElementViewAccessor< arrayView2d< real64, compflow::USD_COMP > > const compOutfluxAccessor =
      mesh.getElemManager().constructViewAccessor< array2d< real64, compflow::LAYOUT_COMP >,
                                                   arrayView2d< real64, compflow::USD_COMP > >( extrinsicMeshData::flow::componentOutflux::key() );


    fluxApprox.forAllStencils( mesh, [&] ( auto & stencil )
    {

      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

      // While this kernel is waiting for a factory class, pass all the accessors here
      isothermalCompositionalMultiphaseBaseKernels::KernelLaunchSelector1
      < isothermalCompositionalMultiphaseFVMKernels::CFLFluxKernel >( m_numComponents,
                                                                      m_numPhases,
                                                                      dt,
                                                                      stencilWrapper,
                                                                      compFlowAccessors.get( extrinsicMeshData::flow::pressure{} ),
                                                                      compFlowAccessors.get( extrinsicMeshData::flow::gravityCoefficient{} ),
                                                                      compFlowAccessors.get( extrinsicMeshData::flow::phaseVolumeFraction{} ),
                                                                      permeabilityAccessors.get( extrinsicMeshData::permeability::permeability{} ),
                                                                      permeabilityAccessors.get( extrinsicMeshData::permeability::dPerm_dPressure{} ),
                                                                      relPermAccessors.get( extrinsicMeshData::relperm::phaseRelPerm{} ),
                                                                      multiFluidAccessors.get( extrinsicMeshData::multifluid::phaseViscosity{} ),
                                                                      multiFluidAccessors.get( extrinsicMeshData::multifluid::phaseDensity{} ),
                                                                      multiFluidAccessors.get( extrinsicMeshData::multifluid::phaseMassDensity{} ),
                                                                      multiFluidAccessors.get( extrinsicMeshData::multifluid::phaseCompFraction{} ),
                                                                      phaseOutfluxAccessor.toNestedView(),
                                                                      compOutfluxAccessor.toNestedView() );
    } );
  } );

  // Step 3: finalize the (cell-based) computation of the CFL numbers
  real64 localMaxPhaseCFLNumber = 0.0;
  real64 localMaxCompCFLNumber = 0.0;

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView2d< real64 const, compflow::USD_PHASE > const & phaseOutflux =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::phaseOutflux >();
      arrayView2d< real64 const, compflow::USD_COMP > const & compOutflux =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::componentOutflux >();

      arrayView1d< real64 > const & phaseCFLNumber = subRegion.getExtrinsicData< extrinsicMeshData::flow::phaseCFLNumber >();
      arrayView1d< real64 > const & compCFLNumber = subRegion.getExtrinsicData< extrinsicMeshData::flow::componentCFLNumber >();

      arrayView1d< real64 const > const & volume = subRegion.getElementVolume();

      arrayView2d< real64 const, compflow::USD_COMP > const & compDens =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompDensity >();
      arrayView2d< real64 const, compflow::USD_COMP > const compFrac =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompFraction >();
      arrayView2d< real64, compflow::USD_PHASE > const phaseVolFrac =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::phaseVolumeFraction >();

      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
      arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseVisc = fluid.phaseViscosity();

      string const & relpermName = subRegion.getReference< string >( viewKeyStruct::relPermNamesString() );
      RelativePermeabilityBase const & relperm = getConstitutiveModel< RelativePermeabilityBase >( subRegion, relpermName );
      arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseRelPerm = relperm.phaseRelPerm();
      arrayView4d< real64 const, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac = relperm.dPhaseRelPerm_dPhaseVolFraction();

      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
      CoupledSolidBase const & solidModel = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

      real64 subRegionMaxPhaseCFLNumber = 0.0;
      real64 subRegionMaxCompCFLNumber = 0.0;

      arrayView2d< real64 const > const & porosity    = solidModel.getPorosity();

      isothermalCompositionalMultiphaseBaseKernels::KernelLaunchSelector2
      < isothermalCompositionalMultiphaseFVMKernels::CFLKernel >( m_numComponents, m_numPhases,
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

  GEOSX_LOG_LEVEL_RANK_0( 1, getName() << ": Max phase CFL number: " << globalMaxPhaseCFLNumber );
  GEOSX_LOG_LEVEL_RANK_0( 1, getName() << ": Max component CFL number: " << globalMaxCompCFLNumber );
}

real64 CompositionalMultiphaseFVM::calculateResidualNorm( DomainPartition const & domain,
                                                          DofManager const & dofManager,
                                                          arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  real64 localFlowResidualNorm = 0.0;
  real64 localEnergyResidualNorm = 0.0;

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel const & mesh,
                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {

      arrayView1d< globalIndex const > dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
      arrayView1d< real64 const > const & totalDensOld = subRegion.getExtrinsicData< extrinsicMeshData::flow::totalDensityOld >();

      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
      CoupledSolidBase const & solidModel = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

      arrayView1d< real64 const > const & referencePorosity = solidModel.getReferencePorosity();

      real64 subRegionFlowResidualNorm = 0.0;
      real64 subRegionEnergyResidualNorm = 0.0;

      if( m_isThermal )
      {
        thermalCompositionalMultiphaseBaseKernels::
          ResidualNormKernel::
          launch< parallelDevicePolicy<>,
                  parallelDeviceReduce >( localRhs,
                                          rankOffset,
                                          numFluidComponents(),
                                          dofNumber,
                                          elemGhostRank,
                                          referencePorosity,
                                          volume,
                                          totalDensOld,
                                          subRegionFlowResidualNorm,
                                          subRegionEnergyResidualNorm );
      }
      else
      {
        isothermalCompositionalMultiphaseBaseKernels::
          ResidualNormKernel::
          launch< parallelDevicePolicy<>,
                  parallelDeviceReduce >( localRhs,
                                          rankOffset,
                                          numFluidComponents(),
                                          dofNumber,
                                          elemGhostRank,
                                          referencePorosity,
                                          volume,
                                          totalDensOld,
                                          subRegionFlowResidualNorm );
      }
      localFlowResidualNorm   += subRegionFlowResidualNorm;
      localEnergyResidualNorm += subRegionEnergyResidualNorm;
    } );
  } );

  // compute global residual norms
  real64 residual = 0.0;
  if( m_isThermal )
  {
    real64 const flowResidual = std::sqrt( MpiWrapper::sum( localFlowResidualNorm ) );
    real64 const energyResidual = std::sqrt( MpiWrapper::sum( localEnergyResidualNorm ) );
    residual = std::sqrt( flowResidual*flowResidual + energyResidual*energyResidual );
    if( getLogLevel() >= 1 && logger::internal::rank == 0 )
    {
      std::cout << GEOSX_FMT( "    ( Rfluid ) = ( {:4.2e} ) ; ( Renergy ) = ( {:4.2e} ) ; ", flowResidual, energyResidual );
    }
  }
  else
  {
    residual = std::sqrt( MpiWrapper::sum( localFlowResidualNorm ) );
    if( getLogLevel() >= 1 && logger::internal::rank == 0 )
    {
      std::cout << GEOSX_FMT( "    ( Rfluid ) = ( {:4.2e} ) ; ", residual );
    }
  }
  return residual;
}

real64 CompositionalMultiphaseFVM::scalingForSystemSolution( DomainPartition const & domain,
                                                             DofManager const & dofManager,
                                                             arrayView1d< real64 const > const & localSolution )
{
  GEOSX_MARK_FUNCTION;

  // check if we want to rescale the Newton update
  if( m_maxCompFracChange >= 1.0 )
  {
    // no rescaling wanted, we just return 1.0;
    return 1.0;
  }

  real64 constexpr eps = isothermalCompositionalMultiphaseBaseKernels::minDensForDivision;
  real64 const maxCompFracChange = m_maxCompFracChange;

  localIndex const NC = m_numComponents;

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  real64 scalingFactor = 1.0;

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel const & mesh,
                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();


      arrayView2d< real64 const, compflow::USD_COMP > const & compDens =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompDensity >();
      arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaGlobalCompDensity >();

      RAJA::ReduceMin< parallelDeviceReduce, real64 > minVal( 1.0 );

      forAll< parallelDevicePolicy<> >( dofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
        if( elemGhostRank[ei] < 0 )
        {
          real64 prevTotalDens = 0;
          for( localIndex ic = 0; ic < NC; ++ic )
          {
            prevTotalDens += compDens[ei][ic] + dCompDens[ei][ic];
          }

          // compute the change in component densities and component fractions
          for( localIndex ic = 0; ic < NC; ++ic )
          {
            localIndex const lid = dofNumber[ei] + ic + 1 - rankOffset;

            // compute scaling factor based on relative change in component densities
            real64 const absCompDensChange = fabs( localSolution[lid] );
            real64 const maxAbsCompDensChange = maxCompFracChange * prevTotalDens;

            // This actually checks the change in component fraction, using a lagged total density
            // Indeed we can rewrite the following check as:
            //    | prevCompDens / prevTotalDens - newCompDens / prevTotalDens | > maxCompFracChange
            // Note that the total density in the second term is lagged (i.e, we use prevTotalDens)
            // because I found it more robust than using directly newTotalDens (which can vary also
            // wildly when the compDens change is large)
            if( absCompDensChange > maxAbsCompDensChange && absCompDensChange > eps )
            {
              minVal.min( maxAbsCompDensChange / absCompDensChange );
            }
          }
        }
      } );

      if( minVal.get() < scalingFactor )
      {
        scalingFactor = minVal.get();
      }
    } );
  } );

  return LvArray::math::max( MpiWrapper::min( scalingFactor, MPI_COMM_GEOSX ), m_minScalingFactor );
}

bool CompositionalMultiphaseFVM::checkSystemSolution( DomainPartition const & domain,
                                                      DofManager const & dofManager,
                                                      arrayView1d< real64 const > const & localSolution,
                                                      real64 const scalingFactor )
{
  GEOSX_MARK_FUNCTION;

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  localIndex localCheck = 1;

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel const & mesh,
                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

      arrayView1d< real64 const > const & pres = subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::pressure::key() );
      arrayView1d< real64 const > const & dPres = subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::deltaPressure::key() );
      arrayView2d< real64 const, compflow::USD_COMP > const & compDens =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompDensity >();
      arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaGlobalCompDensity >();

      localIndex const subRegionSolutionCheck =
        isothermalCompositionalMultiphaseBaseKernels::
          SolutionCheckKernel::launch< parallelDevicePolicy<>,
                                       parallelDeviceReduce >( localSolution,
                                                               dofManager.rankOffset(),
                                                               numFluidComponents(),
                                                               dofNumber,
                                                               elemGhostRank,
                                                               pres,
                                                               dPres,
                                                               compDens,
                                                               dCompDens,
                                                               m_allowCompDensChopping,
                                                               scalingFactor );

      localCheck = std::min( localCheck, subRegionSolutionCheck );
    } );
  } );

  return MpiWrapper::min( localCheck, MPI_COMM_GEOSX );
}

void CompositionalMultiphaseFVM::applySystemSolution( DofManager const & dofManager,
                                                      arrayView1d< real64 const > const & localSolution,
                                                      real64 const scalingFactor,
                                                      DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  DofManager::CompMask pressureMask( m_numDofPerCell, 0, 1 );

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::elemDofFieldString(),
                               extrinsicMeshData::flow::deltaPressure::key(),
                               scalingFactor,
                               pressureMask );

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::elemDofFieldString(),
                               extrinsicMeshData::flow::deltaGlobalCompDensity::key(),
                               scalingFactor,
                               ~pressureMask );

  // if component density chopping is allowed, some component densities may be negative after the update
  // these negative component densities are set to zero in this function
  if( m_allowCompDensChopping )
  {
    chopNegativeDensities( domain );
  }

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & )
  {
    std::map< string, string_array > fieldNames;
    fieldNames["elems"].emplace_back( extrinsicMeshData::flow::deltaPressure::key() );
    fieldNames["elems"].emplace_back( extrinsicMeshData::flow::deltaGlobalCompDensity::key() );

    if( m_isThermal )
    {
      fieldNames["elems"].emplace_back( extrinsicMeshData::flow::deltaTemperature::key() );
    }

    CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), true );
  } );
}

void CompositionalMultiphaseFVM::updatePhaseMobility( ObjectManagerBase & dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  // note that the phase mobility computed here also includes phase density
  string const & fluidName = dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( dataGroup, fluidName );

  string const & relpermName = dataGroup.getReference< string >( viewKeyStruct::relPermNamesString() );
  RelativePermeabilityBase const & relperm = getConstitutiveModel< RelativePermeabilityBase >( dataGroup, relpermName );

  if( m_isThermal )
  {
    thermalCompositionalMultiphaseFVMKernels::
      PhaseMobilityKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                 m_numPhases,
                                                 dataGroup,
                                                 fluid,
                                                 relperm );
  }
  else
  {
    isothermalCompositionalMultiphaseFVMKernels::
      PhaseMobilityKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                 m_numPhases,
                                                 dataGroup,
                                                 fluid,
                                                 relperm );
  }
}

void CompositionalMultiphaseFVM::applyAquiferBC( real64 const time,
                                                 real64 const dt,
                                                 DofManager const & dofManager,
                                                 DomainPartition & domain,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & )
  {
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

    string const & elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > elemDofNumber =
      mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
    elemDofNumber.setName( getName() + "/accessors/" + elemDofKey );

    isothermalCompositionalMultiphaseFVMKernels::
      AquiferBCKernel::CompFlowAccessors compFlowAccessors( mesh.getElemManager(), getName() );
    isothermalCompositionalMultiphaseFVMKernels::
      AquiferBCKernel::MultiFluidAccessors multiFluidAccessors( mesh.getElemManager(), getName() );

    fsManager.apply< AquiferBoundaryCondition >( time + dt,
                                                 domain,
                                                 "faceManager",
                                                 AquiferBoundaryCondition::catalogName(),
                                                 [&] ( AquiferBoundaryCondition const & bc,
                                                       string const & setName,
                                                       SortedArrayView< localIndex const > const &,
                                                       Group & subRegion,
                                                       string const & )
    {
      BoundaryStencil const & stencil = fluxApprox.getStencil< BoundaryStencil >( mesh, setName );
      if( bc.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
      {
        globalIndex const numTargetFaces = MpiWrapper::sum< globalIndex >( stencil.size() );
        string const logMessage = string( "CompositionalMultiphaseFVM {}: at time {}s, " )
                                  + string( "the <{}> boundary condition '{}' is applied to the face set '{}' in '{}'. " )
                                  + string( "\nThe total number of target faces (including ghost faces) is {}. " )
                                  + string( "\nNote that if this number is equal to zero, the boundary condition will not be applied on this face set." );
        GEOSX_LOG_RANK_0( GEOSX_FMT( logMessage,
                                     getName(), time+dt, AquiferBoundaryCondition::catalogName(),
                                     bc.getName(), setName, subRegion.getName(), numTargetFaces ) );
      }

      if( stencil.size() == 0 )
      {
        return;
      }

      AquiferBoundaryCondition::KernelWrapper aquiferBCWrapper = bc.createKernelWrapper();
      bool const allowAllPhasesIntoAquifer = bc.allowAllPhasesIntoAquifer();
      localIndex const waterPhaseIndex = bc.getWaterPhaseIndex();
      real64 const & aquiferWaterPhaseDens = bc.getWaterPhaseDensity();
      arrayView1d< real64 const > const & aquiferWaterPhaseCompFrac = bc.getWaterPhaseComponentFraction();

      // While this kernel is waiting for a factory class, pass all the accessors here
      isothermalCompositionalMultiphaseBaseKernels::KernelLaunchSelector1
      < isothermalCompositionalMultiphaseFVMKernels::AquiferBCKernel >( m_numComponents,
                                                                        m_numPhases,
                                                                        waterPhaseIndex,
                                                                        allowAllPhasesIntoAquifer,
                                                                        stencil,
                                                                        dofManager.rankOffset(),
                                                                        elemDofNumber.toNestedViewConst(),
                                                                        aquiferBCWrapper,
                                                                        aquiferWaterPhaseDens,
                                                                        aquiferWaterPhaseCompFrac,
                                                                        compFlowAccessors.get( extrinsicMeshData::ghostRank{} ),
                                                                        compFlowAccessors.get( extrinsicMeshData::flow::pressure{} ),
                                                                        compFlowAccessors.get( extrinsicMeshData::flow::deltaPressure{} ),
                                                                        compFlowAccessors.get( extrinsicMeshData::flow::gravityCoefficient{} ),
                                                                        compFlowAccessors.get( extrinsicMeshData::flow::phaseVolumeFraction{} ),
                                                                        compFlowAccessors.get( extrinsicMeshData::flow::dPhaseVolumeFraction_dPressure{} ),
                                                                        compFlowAccessors.get( extrinsicMeshData::flow::dPhaseVolumeFraction_dGlobalCompDensity{} ),
                                                                        compFlowAccessors.get( extrinsicMeshData::flow::dGlobalCompFraction_dGlobalCompDensity{} ),
                                                                        multiFluidAccessors.get( extrinsicMeshData::multifluid::phaseDensity{} ),
                                                                        multiFluidAccessors.get( extrinsicMeshData::multifluid::dPhaseDensity{} ),
                                                                        multiFluidAccessors.get( extrinsicMeshData::multifluid::phaseCompFraction{} ),
                                                                        multiFluidAccessors.get( extrinsicMeshData::multifluid::dPhaseCompFraction{} ),
                                                                        time,
                                                                        dt,
                                                                        localMatrix.toViewConstSizes(),
                                                                        localRhs.toView() );
    } );
  } );

}

//START_SPHINX_INCLUDE_01
REGISTER_CATALOG_ENTRY( SolverBase, CompositionalMultiphaseFVM, string const &, Group * const )
//END_SPHINX_INCLUDE_01
}// namespace geosx
