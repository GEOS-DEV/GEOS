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
{}

void CompositionalMultiphaseFVM::initializePreSubGroups()
{
  CompositionalMultiphaseBase::initializePreSubGroups();

  m_linearSolverParameters.get().mgr.strategy = m_isThermal
    ? LinearSolverParameters::MGR::StrategyType::thermalCompositionalMultiphaseFVM
    : LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseFVM;

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
                       FieldLocation::Elem,
                       m_numDofPerCell,
                       getMeshTargets() );

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

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & )
  {
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

    string const & elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

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
                                                     m_hasCapPressure,
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
                                                     m_hasCapPressure,
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

real64 CompositionalMultiphaseFVM::calculateResidualNorm( DomainPartition const & domain,
                                                          DofManager const & dofManager,
                                                          arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  real64 localFlowResidualNorm = 0.0;
  real64 localEnergyResidualNorm = 0.0;

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {

      arrayView1d< globalIndex const > dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const volume = subRegion.getElementVolume();

      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
      arrayView2d< real64 const, multifluid::USD_FLUID > const totalDens_n = fluid.totalDensity_n();

      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
      CoupledSolidBase const & solidModel = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );
      arrayView1d< real64 const > const referencePorosity = solidModel.getReferencePorosity();

      real64 subRegionFlowResidualNorm = 0.0;
      real64 subRegionEnergyResidualNorm = 0.0;

      if( m_isThermal )
      {
        arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac_n =
          subRegion.getExtrinsicData< extrinsicMeshData::flow::phaseVolumeFraction_n >();
        arrayView3d< real64 const, multifluid::USD_PHASE > const phaseDens_n = fluid.phaseDensity_n();
        arrayView3d< real64 const, multifluid::USD_PHASE > const phaseInternalEnergy_n = fluid.phaseInternalEnergy_n();

        string const & solidInternalEnergyName = subRegion.getReference< string >( viewKeyStruct::solidInternalEnergyNamesString() );
        SolidInternalEnergy const & solidInternalEnergy = getConstitutiveModel< SolidInternalEnergy >( subRegion, solidInternalEnergyName );
        arrayView2d< real64 const > const solidInternalEnergy_n = solidInternalEnergy.getInternalEnergy_n();

        thermalCompositionalMultiphaseBaseKernels::
          ResidualNormKernel::
          launch< parallelDevicePolicy<> >( localRhs,
                                            rankOffset,
                                            numFluidPhases(),
                                            numFluidComponents(),
                                            dofNumber,
                                            elemGhostRank,
                                            referencePorosity,
                                            volume,
                                            solidInternalEnergy_n,
                                            phaseVolFrac_n,
                                            totalDens_n,
                                            phaseDens_n,
                                            phaseInternalEnergy_n,
                                            subRegionFlowResidualNorm,
                                            subRegionEnergyResidualNorm );
      }
      else
      {
        isothermalCompositionalMultiphaseBaseKernels::
          ResidualNormKernel::
          launch< parallelDevicePolicy<> >( localRhs,
                                            rankOffset,
                                            numFluidComponents(),
                                            dofNumber,
                                            elemGhostRank,
                                            referencePorosity,
                                            volume,
                                            totalDens_n,
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
      std::cout << GEOSX_FMT( "    ( R{} ) = ( {:4.2e} ) ; ( Renergy ) = ( {:4.2e} ) ; ", coupledSolverAttributePrefix(), flowResidual, energyResidual );
    }
  }
  else
  {
    residual = std::sqrt( MpiWrapper::sum( localFlowResidualNorm ) );
    if( getLogLevel() >= 1 && logger::internal::rank == 0 )
    {
      std::cout << GEOSX_FMT( "    ( R{} ) = ( {:4.2e} ) ; ", coupledSolverAttributePrefix(), residual );
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

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
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

      RAJA::ReduceMin< parallelDeviceReduce, real64 > minVal( 1.0 );

      forAll< parallelDevicePolicy<> >( dofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
        if( elemGhostRank[ei] < 0 )
        {
          real64 prevTotalDens = 0;
          for( localIndex ic = 0; ic < NC; ++ic )
          {
            prevTotalDens += compDens[ei][ic];
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

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
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
      arrayView2d< real64 const, compflow::USD_COMP > const & compDens =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompDensity >();

      localIndex const subRegionSolutionCheck =
        isothermalCompositionalMultiphaseBaseKernels::
          SolutionCheckKernel::launch< parallelDevicePolicy<> >( localSolution,
                                                                 dofManager.rankOffset(),
                                                                 numFluidComponents(),
                                                                 dofNumber,
                                                                 elemGhostRank,
                                                                 pres,
                                                                 compDens,
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
  DofManager::CompMask componentMask( m_numDofPerCell, 1, m_numComponents+1 );

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::elemDofFieldString(),
                               extrinsicMeshData::flow::pressure::key(),
                               scalingFactor,
                               pressureMask );

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::elemDofFieldString(),
                               extrinsicMeshData::flow::globalCompDensity::key(),
                               scalingFactor,
                               componentMask );

  if( m_isThermal )
  {
    DofManager::CompMask temperatureMask( m_numDofPerCell, m_numComponents+1, m_numComponents+2 );
    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 extrinsicMeshData::flow::temperature::key(),
                                 scalingFactor,
                                 temperatureMask );
  }

  // if component density chopping is allowed, some component densities may be negative after the update
  // these negative component densities are set to zero in this function
  if( m_allowCompDensChopping )
  {
    chopNegativeDensities( domain );
  }

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    std::vector< string > fields{ extrinsicMeshData::flow::pressure::key(), extrinsicMeshData::flow::globalCompDensity::key() };
    if( m_isThermal )
    {
      fields.emplace_back( extrinsicMeshData::flow::temperature::key() );
    }
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( fields, regionNames );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), true );
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

void CompositionalMultiphaseFVM::applyBoundaryConditions( real64 time_n,
                                                          real64 dt,
                                                          DomainPartition & domain,
                                                          DofManager const & dofManager,
                                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                          arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  CompositionalMultiphaseBase::applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );
  applyFaceDirichletBC( time_n, dt, dofManager, domain, localMatrix, localRhs );
}

bool CompositionalMultiphaseFVM::validateFaceDirichletBC( DomainPartition & domain,
                                                          real64 const time ) const
{
  constexpr integer MAX_NC = MultiFluidBase::MAX_NUM_COMPONENTS;
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  bool bcConsistent = true;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & )
  {

    // maps to check consistent application of BC
    // maps: setName (-> numComps)
    map< string, ComponentMask< MAX_NC > > bcPresCompStatusMap; // check that pressure/comp are present/consistent
    set< string > bcTempStatusMap; // check that temperature is present/consistent

    // 1. Check pressure Dirichlet BCs
    fsManager.apply< FaceManager >( time,
                                    mesh,
                                    extrinsicMeshData::flow::pressure::key(),
                                    [&]( FieldSpecificationBase const &,
                                         string const & setName,
                                         SortedArrayView< localIndex const > const &,
                                         FaceManager &,
                                         string const & )
    {
      // Check whether pressure has already been applied to this set
      if( bcPresCompStatusMap.count( setName ) > 0 )
      {
        bcConsistent = false;
        GEOSX_WARNING( GEOSX_FMT( "Conflicting pressure boundary conditions on set {}", setName ) );
      }
      bcPresCompStatusMap[setName].setNumComp( m_numComponents );
    } );

    // 2. Check temperature Dirichlet BCs (we always require a temperature for face-based BCs)
    fsManager.apply< FaceManager >( time,
                                    mesh,
                                    extrinsicMeshData::flow::temperature::key(),
                                    [&]( FieldSpecificationBase const &,
                                         string const & setName,
                                         SortedArrayView< localIndex const > const &,
                                         FaceManager &,
                                         string const & )
    {
      // 2.1 Check whether temperature has already been applied to this set
      if( bcTempStatusMap.count( setName ) > 0 )
      {
        bcConsistent = false;
        GEOSX_WARNING( GEOSX_FMT( "Conflicting temperature boundary conditions on set {}", setName ) );
      }
      bcTempStatusMap.insert( setName );

      // 2.2 Check that there is pressure bc applied to this set
      if( bcPresCompStatusMap.count( setName ) == 0 )
      {
        bcConsistent = false;
        GEOSX_WARNING( GEOSX_FMT( "Pressure boundary condition not prescribed on set {}", setName ) );
      }
    } );

    // 3. Check composition BC (global component fraction)
    fsManager.apply< FaceManager >( time,
                                    mesh,
                                    extrinsicMeshData::flow::globalCompFraction::key(),
                                    [&] ( FieldSpecificationBase const & fs,
                                          string const & setName,
                                          SortedArrayView< localIndex const > const &,
                                          FaceManager &,
                                          string const & )
    {
      // 3.1 Check pressure, temperature, and record composition bc application
      integer const comp = fs.getComponent();

      if( bcPresCompStatusMap.count( setName ) == 0 )
      {
        bcConsistent = false;
        GEOSX_WARNING( GEOSX_FMT( "Pressure boundary condition not prescribed on set {}", setName ) );
      }
      if( bcTempStatusMap.count( setName ) == 0 )
      {
        bcConsistent = false;
        GEOSX_WARNING( GEOSX_FMT( "Temperature boundary condition not prescribed on set {}. \n"
                                  "Note that for face boundary conditions, you must provide a temperature", setName ) );
      }
      if( comp < 0 || comp >= m_numComponents )
      {
        bcConsistent = false;
        GEOSX_WARNING( GEOSX_FMT( "Invalid component index [{}] in composition boundary condition {}", comp, fs.getName() ) );
        return; // can't check next part with invalid component id
      }

      ComponentMask< MAX_NC > & compMask = bcPresCompStatusMap[setName];
      if( compMask[comp] )
      {
        bcConsistent = false;
        GEOSX_WARNING( GEOSX_FMT( "Conflicting composition[{}] boundary conditions on set {}", comp, setName ) );
      }
      compMask.set( comp );
    } );

    // 3.2 Check consistency between composition BC applied to sets
    for( auto const & setEntry : bcPresCompStatusMap )
    {
      ComponentMask< MAX_NC > const & compMask = setEntry.second;
      for( integer ic = 0; ic < m_numComponents; ++ic )
      {
        if( !compMask[ic] )
        {
          bcConsistent = false;
          GEOSX_WARNING( GEOSX_FMT( "Boundary condition not applied to composition[{}] on set {}", ic, setEntry.first ) );
        }
      }
    }
  } );

  return bcConsistent;
}

namespace
{
char const faceBcLogMessage[] =
  "CompositionalMultiphaseFVM {}: at time {}s, "
  "the <{}> boundary condition '{}' is applied to the face set '{}' in '{}'. "
  "\nThe scale of this boundary condition is {} and multiplies the value of the provided function (if any). "
  "\nThe total number of target faces (including ghost faces) is {}."
  "\nNote that if this number is equal to zero, the boundary condition will not be applied on this face set.";
}

void CompositionalMultiphaseFVM::applyFaceDirichletBC( real64 const time_n,
                                                       real64 const dt,
                                                       DofManager const & dofManager,
                                                       DomainPartition & domain,
                                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                       arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  // Only validate BC at the beginning of Newton loop
  if( m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
  {
    bool const bcConsistent = validateFaceDirichletBC( domain, time_n + dt );
    GEOSX_ERROR_IF( !bcConsistent, GEOSX_FMT( "CompositionalMultiphaseBase {}: inconsistent boundary conditions", getName() ) );
  }

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    // Take BCs defined for "pressure" field and apply values to "facePressure"
    applyFieldValue< FaceManager >( time_n, dt, mesh, faceBcLogMessage,
                                    extrinsicMeshData::flow::pressure::key(), extrinsicMeshData::flow::facePressure::key() );
    // Take BCs defined for "globalCompFraction" field and apply values to "faceGlobalCompFraction"
    applyFieldValue< FaceManager >( time_n, dt, mesh, faceBcLogMessage,
                                    extrinsicMeshData::flow::globalCompFraction::key(), extrinsicMeshData::flow::faceGlobalCompFraction::key() );
    // Take BCs defined for "temperature" field and apply values to "faceTemperature"
    applyFieldValue< FaceManager >( time_n, dt, mesh, faceBcLogMessage,
                                    extrinsicMeshData::flow::temperature::key(), extrinsicMeshData::flow::faceTemperature::key() );

    // Then launch the face Dirichlet kernel
    fsManager.apply< FaceManager >( time_n + dt,
                                    mesh,
                                    extrinsicMeshData::flow::pressure::key(), // we have required that pressure is always present
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
      MultiFluidBase & multiFluidBase = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

      BoundaryStencilWrapper const stencilWrapper = stencil.createKernelWrapper();

      string const & elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

      if( m_isThermal )
      {
        thermalCompositionalMultiphaseFVMKernels::
          DirichletFaceBasedAssemblyKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                     m_numPhases,
                                                     dofManager.rankOffset(),
                                                     elemDofKey,
                                                     getName(),
                                                     faceManager,
                                                     elemManager,
                                                     stencilWrapper,
                                                     multiFluidBase,
                                                     dt,
                                                     localMatrix,
                                                     localRhs );
      }
      else
      {
        isothermalCompositionalMultiphaseFVMKernels::
          DirichletFaceBasedAssemblyKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                     m_numPhases,
                                                     dofManager.rankOffset(),
                                                     elemDofKey,
                                                     getName(),
                                                     faceManager,
                                                     elemManager,
                                                     stencilWrapper,
                                                     multiFluidBase,
                                                     dt,
                                                     localMatrix,
                                                     localRhs );
      }

    } );
  } );
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

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
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

    fsManager.apply< FaceManager,
                     AquiferBoundaryCondition >( time + dt,
                                                 mesh,
                                                 AquiferBoundaryCondition::catalogName(),
                                                 [&] ( AquiferBoundaryCondition const & bc,
                                                       string const & setName,
                                                       SortedArrayView< localIndex const > const &,
                                                       FaceManager & faceManager,
                                                       string const & )
    {
      BoundaryStencil const & stencil = fluxApprox.getStencil< BoundaryStencil >( mesh, setName );
      if( bc.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
      {
        globalIndex const numTargetFaces = MpiWrapper::sum< globalIndex >( stencil.size() );
        GEOSX_LOG_RANK_0( GEOSX_FMT( faceBcLogMessage,
                                     getName(), time+dt, AquiferBoundaryCondition::catalogName(),
                                     bc.getName(), setName, faceManager.getName(), bc.getScale(), numTargetFaces ) );
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
                                                                        compFlowAccessors.get( extrinsicMeshData::flow::pressure_n{} ),
                                                                        compFlowAccessors.get( extrinsicMeshData::flow::gravityCoefficient{} ),
                                                                        compFlowAccessors.get( extrinsicMeshData::flow::phaseVolumeFraction{} ),
                                                                        compFlowAccessors.get( extrinsicMeshData::flow::dPhaseVolumeFraction{} ),
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
