/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CompositionalMultiphaseHybridFVM.cpp
 */

#include "CompositionalMultiphaseHybridFVM.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseHybridFVMKernels.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseHybridFVMKernels.hpp"


/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace CompositionalMultiphaseBaseKernels;
using namespace CompositionalMultiphaseHybridFVMKernels;

CompositionalMultiphaseHybridFVM::CompositionalMultiphaseHybridFVM( const std::string & name,
                                                                    Group * const parent ):
  CompositionalMultiphaseBase( name, parent ),
  m_lengthTolerance( 0 )
{

  this->registerWrapper( viewKeyStruct::maxRelativePresChangeString, &m_maxRelativePresChange )->
    setSizedFromParent( 0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setApplyDefaultValue( 1.0 )->
    setDescription( "Maximum (relative) change in (face) pressure between two Newton iterations" );

  m_linearSolverParameters.get().mgr.strategy = "CompositionalMultiphaseHybridFVM";

}

void CompositionalMultiphaseHybridFVM::RegisterDataOnMesh( Group * const MeshBodies )
{
  GEOSX_MARK_FUNCTION;

  // 1) Register the elem-centered data
  CompositionalMultiphaseBase::RegisterDataOnMesh( MeshBodies );

  // 2) Register the face data
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel & meshLevel = *Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    FaceManager & faceManager = *meshLevel.getFaceManager();

    // primary variables: face pressure changes
    faceManager.registerWrapper< array1d< real64 > >( viewKeyStruct::deltaFacePressureString )->
      setPlotLevel( PlotLevel::LEVEL_0 )->
      setRegisteringObjects( this->getName())->
      setDescription( "An array that holds the accumulated phase pressure updates at the faces." );

    // auxiliary data for the buoyancy coefficient
    faceManager.registerWrapper< array1d< real64 > >( viewKeyStruct::mimGravityCoefString );
  }
}

void CompositionalMultiphaseHybridFVM::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
{
  GEOSX_MARK_FUNCTION;

  DomainPartition & domain = *rootGroup->GetGroup< DomainPartition >( keys::domain );
  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager const & elemManager = *mesh.getElemManager();
  FaceManager const & faceManager = *mesh.getFaceManager();
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  m_lengthTolerance = domain.getMeshBody( 0 )->getGlobalLengthScale() * 1e-8;
  string const & coeffName = fluxApprox.template getReference< string >( FluxApproximationBase::viewKeyStruct::coeffNameString );
  m_transMultName = coeffName + FluxApproximationBase::viewKeyStruct::transMultiplierString;

  CompositionalMultiphaseBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  // in the flux kernel, we need to make sure that we act only on the target regions
  // for that, we need the following region filter
  for( string const & regionName : targetRegionNames() )
  {
    m_regionFilter.insert( elemManager.GetRegions().getIndex( regionName ) );
  }

  // check that multipliers are stricly larger than 0, which would work with SinglePhaseFVM, but not with SinglePhaseHybridFVM.
  // To deal with a 0 multiplier, we would just have to skip the corresponding face in the FluxKernel
  arrayView1d< real64 const > const & transMultiplier =
    faceManager.getReference< array1d< real64 > >( m_transMultName );

  RAJA::ReduceMin< parallelDeviceReduce, real64 > minVal( 1.0 );
  forAll< parallelDevicePolicy<> >( faceManager.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iface )
  {
    minVal.min( transMultiplier[iface] );
  } );

  GEOSX_ERROR_IF_LE_MSG( minVal.get(), 0.0,
                         "The transmissibility multipliers used in SinglePhaseHybridFVM must strictly larger than 0.0" );

}

void CompositionalMultiphaseHybridFVM::PrecomputeData( MeshLevel & mesh )
{
  FlowSolverBase::PrecomputeData( mesh );

  NodeManager const & nodeManager = *mesh.getNodeManager();
  FaceManager & faceManager = *mesh.getFaceManager();

  array1d< RAJA::ReduceSum< parallelDeviceReduce, real64 > > mimFaceGravCoefNumerator;
  array1d< RAJA::ReduceSum< parallelDeviceReduce, real64 > > mimFaceGravCoefDenominator;
  mimFaceGravCoefNumerator.resize( faceManager.size() );
  mimFaceGravCoefDenominator.resize( faceManager.size() );

  // node data

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  // face data

  arrayView1d< real64 const > const & transMultiplier =
    faceManager.getReference< array1d< real64 > >( m_transMultName );

  arrayView1d< real64 > const mimFaceGravCoef =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::mimGravityCoefString );

  ArrayOfArraysView< localIndex const > const & faceToNodes = faceManager.nodeList().toViewConst();

  real64 const lengthTolerance = m_lengthTolerance;

  forTargetSubRegionsComplete< CellElementSubRegion >( mesh,
                                                       [&]( localIndex const,
                                                            localIndex const,
                                                            localIndex const,
                                                            ElementRegionBase const &,
                                                            auto const & subRegion )
  {
    arrayView2d< real64 const > const & elemCenter =
      subRegion.template getReference< array2d< real64 > >( CellBlock::viewKeyStruct::elementCenterString );
    arrayView1d< R1Tensor const > const & elemPerm =
      subRegion.template getReference< array1d< R1Tensor > >( CompositionalMultiphaseBase::viewKeyStruct::permeabilityString );
    arrayView1d< real64 const > const elemGravCoef =
      subRegion.template getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );
    arrayView2d< localIndex const > const & elemToFaces = subRegion.faceList();

    SinglePhaseHybridFVMKernels::KernelLaunchSelector< PrecomputeKernel >( subRegion.numFacesPerElement(),
                                                                           subRegion.size(),
                                                                           faceManager.size(),
                                                                           nodePosition,
                                                                           faceToNodes,
                                                                           elemCenter,
                                                                           elemPerm,
                                                                           elemGravCoef,
                                                                           elemToFaces,
                                                                           transMultiplier,
                                                                           lengthTolerance,
                                                                           mimFaceGravCoefNumerator.toView(),
                                                                           mimFaceGravCoefDenominator.toView(),
                                                                           mimFaceGravCoef );
  } );

}

void CompositionalMultiphaseHybridFVM::ImplicitStepSetup( real64 const & time_n,
                                                          real64 const & dt,
                                                          DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  // setup the elem-centered fields
  CompositionalMultiphaseBase::ImplicitStepSetup( time_n, dt, domain );

  // setup the face fields
  MeshLevel & meshLevel = *domain.getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  FaceManager & faceManager = *meshLevel.getFaceManager();

  // get the accumulated pressure updates
  arrayView1d< real64 > const & dFacePres =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::deltaFacePressureString );

  // zero out the face pressures
  dFacePres.setValues< parallelDevicePolicy<> >( 0.0 );
}


void CompositionalMultiphaseHybridFVM::ImplicitStepComplete( real64 const & time_n,
                                                             real64 const & dt,
                                                             DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  // increment the elem-centered fields
  CompositionalMultiphaseBase::ImplicitStepComplete( time_n, dt, domain );

  // increment the face fields
  MeshLevel & meshLevel = *domain.getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  FaceManager & faceManager = *meshLevel.getFaceManager();

  // get the face-based pressure
  arrayView1d< real64 > const facePres =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::facePressureString );
  arrayView1d< real64 > const dFacePres =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::deltaFacePressureString );

  forAll< parallelDevicePolicy<> >( faceManager.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iface )
  {
    facePres[iface] += dFacePres[iface];
    dFacePres[iface] = 0.0;
  } );
}


void CompositionalMultiphaseHybridFVM::SetupDofs( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                                  DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;

  // setup the connectivity of elem fields
  // we need Connectivity::Face because of the two-point upwinding
  // in AssembleOneSidedMassFluxes
  dofManager.addField( viewKeyStruct::elemDofFieldString,
                       DofManager::Location::Elem,
                       m_numDofPerCell,
                       targetRegionNames() );

  dofManager.addCoupling( viewKeyStruct::elemDofFieldString,
                          viewKeyStruct::elemDofFieldString,
                          DofManager::Connector::Face );

  // setup the connectivity of face fields
  dofManager.addField( viewKeyStruct::faceDofFieldString,
                       DofManager::Location::Face,
                       1,
                       targetRegionNames() );

  dofManager.addCoupling( viewKeyStruct::faceDofFieldString,
                          viewKeyStruct::faceDofFieldString,
                          DofManager::Connector::Elem );

  // setup coupling between pressure and face pressure
  dofManager.addCoupling( viewKeyStruct::faceDofFieldString,
                          viewKeyStruct::elemDofFieldString,
                          DofManager::Connector::Elem,
                          true );

}


void CompositionalMultiphaseHybridFVM::AssembleFluxTerms( real64 const dt,
                                                          DomainPartition const & domain,
                                                          DofManager const & dofManager,
                                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                          arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager const & nodeManager = *mesh.getNodeManager();
  FaceManager const & faceManager = *mesh.getFaceManager();

  // node data (for transmissibility computation)

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  // face data

  // get the face-based DOF numbers for the assembly
  string const faceDofKey = dofManager.getKey( viewKeyStruct::faceDofFieldString );
  arrayView1d< globalIndex const > const & faceDofNumber =
    faceManager.getReference< array1d< globalIndex > >( faceDofKey );
  arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();

  // get the element dof numbers for the assembly
  string const & elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > elemDofNumber =
    mesh.getElemManager()->ConstructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
  elemDofNumber.setName( getName() + "/accessors/" + elemDofKey );

  // get the face-centered pressures
  arrayView1d< real64 const > const & facePres =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::facePressureString );
  arrayView1d< real64 const > const & dFacePres =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::deltaFacePressureString );

  // get the face-centered depth
  arrayView1d< real64 const > const & faceGravCoef =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );
  arrayView1d< real64 const > const & mimFaceGravCoef =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::mimGravityCoefString );

  // get the face-centered transMultiplier
  // TODO: implement some kind of HybridFVMApprox that inherits from FluxApproximationBase
  arrayView1d< real64 const > const & transMultiplier =
    faceManager.getReference< array1d< real64 > >( m_transMultName );

  // get the face-to-nodes connectivity for the transmissibility calculation
  ArrayOfArraysView< localIndex const > const & faceToNodes = faceManager.nodeList().toViewConst();

  arrayView2d< localIndex const > const & elemRegionList    = faceManager.elementRegionList().toViewConst();
  arrayView2d< localIndex const > const & elemSubRegionList = faceManager.elementSubRegionList().toViewConst();
  arrayView2d< localIndex const > const & elemList          = faceManager.elementList().toViewConst();

  // tolerance for transmissibility calculation
  real64 const lengthTolerance = m_lengthTolerance;

  forTargetSubRegionsComplete< CellElementSubRegion >( mesh,
                                                       [&]( localIndex const,
                                                            localIndex const er,
                                                            localIndex const esr,
                                                            ElementRegionBase const &,
                                                            auto const & subRegion )
  {
    KernelLaunchSelector< FluxKernel >( subRegion.numFacesPerElement(),
                                        m_numComponents, m_numPhases,
                                        er, esr, subRegion,
                                        m_regionFilter.toViewConst(),
                                        nodePosition,
                                        elemRegionList,
                                        elemSubRegionList,
                                        elemList,
                                        faceToNodes,
                                        faceDofNumber,
                                        faceGhostRank,
                                        facePres,
                                        dFacePres,
                                        faceGravCoef,
                                        mimFaceGravCoef,
                                        transMultiplier,
                                        m_phaseDens.toNestedViewConst(),
                                        m_dPhaseDens_dPres.toNestedViewConst(),
                                        m_dPhaseDens_dComp.toNestedViewConst(),
                                        m_phaseMob.toNestedViewConst(),
                                        m_dPhaseMob_dPres.toNestedViewConst(),
                                        m_dPhaseMob_dCompDens.toNestedViewConst(),
                                        m_dCompFrac_dCompDens.toNestedViewConst(),
                                        m_phaseCompFrac.toNestedViewConst(),
                                        m_dPhaseCompFrac_dPres.toNestedViewConst(),
                                        m_dPhaseCompFrac_dComp.toNestedViewConst(),
                                        elemDofNumber.toNestedViewConst(),
                                        dofManager.rankOffset(),
                                        lengthTolerance,
                                        dt,
                                        localMatrix,
                                        localRhs );

  } );
}

real64 CompositionalMultiphaseHybridFVM::ScalingForSystemSolution( DomainPartition const & domain,
                                                                   DofManager const & dofManager,
                                                                   arrayView1d< real64 const > const & localSolution )
{
  GEOSX_MARK_FUNCTION;

  // check if we want to rescale the Newton update
  if( m_maxCompFracChange >= 1.0 && m_maxRelativePresChange >= 1.0 )
  {
    // no rescaling wanted, we just return 1.0;
    return 1.0;
  }

  real64 constexpr eps = CompositionalMultiphaseBaseKernels::minDensForDivision;
  real64 const maxCompFracChange = m_maxCompFracChange;
  real64 const maxRelativePresChange = m_maxRelativePresChange;

  localIndex const NC = m_numComponents;

  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  FaceManager const & faceManager = *mesh.getFaceManager();

  string const faceDofKey = dofManager.getKey( viewKeyStruct::faceDofFieldString );
  arrayView1d< globalIndex const > const & faceDofNumber =
    faceManager.getReference< array1d< globalIndex > >( faceDofKey );
  arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();

  arrayView1d< real64 const > const & facePressure =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::facePressureString );
  arrayView1d< real64 const > const & dFacePressure =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::deltaFacePressureString );

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString );
  real64 scalingFactor = 1.0;

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase const & subRegion )
  {
    arrayView1d< globalIndex const > const & dofNumber =
      subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

    arrayView1d< real64 const > const & pressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView1d< real64 const > const & dPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

    arrayView2d< real64 const > const & compDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString );
    arrayView2d< real64 const > const & dCompDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString );

    RAJA::ReduceMin< parallelDeviceReduce, real64 > minElemVal( 1.0 );

    forAll< parallelDevicePolicy<> >( dofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( elemGhostRank[ei] < 0 )
      {

        // compute the change in pressure
        real64 const pres = pressure[ei] + dPressure[ei];
        real64 const absPresChange = fabs( localSolution[dofNumber[ei] - rankOffset] );
        if( pres > eps )
        {
          real64 const relativePresChange = fabs( absPresChange ) / pres;
          if( relativePresChange > maxRelativePresChange )
          {
            minElemVal.min( maxRelativePresChange / relativePresChange );
          }
        }

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
            minElemVal.min( maxAbsCompDensChange / absCompDensChange );
          }
        }
      }
    } );

    if( minElemVal.get() < scalingFactor )
    {
      scalingFactor = minElemVal.get();
    }
  } );

  RAJA::ReduceMin< parallelDeviceReduce, real64 > minFaceVal( 1.0 );
  forAll< parallelDevicePolicy<> >( faceManager.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iface )
  {
    if( faceGhostRank[iface] < 0 && faceDofNumber[iface] >= 0 )
    {
      real64 const facePres = facePressure[iface] + dFacePressure[iface];
      real64 const absPresChange = fabs( localSolution[faceDofNumber[iface] - rankOffset] );
      if( facePres > eps )
      {
        real64 const relativePresChange = fabs( absPresChange ) / facePres;
        if( relativePresChange > maxRelativePresChange )
        {
          minFaceVal.min( maxRelativePresChange / relativePresChange );
        }
      }
    }
  } );

  if( minFaceVal.get() < scalingFactor )
  {
    scalingFactor = minFaceVal.get();
  }

  return LvArray::math::max( MpiWrapper::Min( scalingFactor, MPI_COMM_GEOSX ), m_minScalingFactor );
}


bool CompositionalMultiphaseHybridFVM::CheckSystemSolution( DomainPartition const & domain,
                                                            DofManager const & dofManager,
                                                            arrayView1d< real64 const > const & localSolution,
                                                            real64 const scalingFactor )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  FaceManager const & faceManager = *mesh.getFaceManager();

  // 1. check cell-centered variables for each region

  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString );
  localIndex localCheck = 1;

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase const & subRegion )
  {
    arrayView1d< globalIndex const > const & elemDofNumber = subRegion.getReference< array1d< globalIndex > >( elemDofKey );
    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

    arrayView1d< real64 const > const & elemPres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView1d< real64 const > const & dElemPres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

    arrayView2d< real64 const > const & compDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString );
    arrayView2d< real64 const > const & dCompDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString );

    localIndex const subRegionSolutionCheck =
      CompositionalMultiphaseBaseKernels::SolutionCheckKernel::Launch< parallelDevicePolicy<>,
                                                                       parallelDeviceReduce >( localSolution,
                                                                                               dofManager.rankOffset(),
                                                                                               numFluidComponents(),
                                                                                               elemDofNumber,
                                                                                               elemGhostRank,
                                                                                               elemPres,
                                                                                               dElemPres,
                                                                                               compDens,
                                                                                               dCompDens,
                                                                                               m_allowCompDensChopping,
                                                                                               scalingFactor );
    localCheck = std::min( localCheck, subRegionSolutionCheck );
  } );

  // 2. check face-centered variables in the domain

  string const faceDofKey = dofManager.getKey( viewKeyStruct::faceDofFieldString );
  arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();
  arrayView1d< globalIndex const > const & faceDofNumber =
    faceManager.getReference< array1d< globalIndex > >( faceDofKey );

  arrayView1d< real64 const > const & facePres =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::facePressureString );
  arrayView1d< real64 const > const & dFacePres =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::deltaFacePressureString );

  localIndex const faceSolutionCheck =
    CompositionalMultiphaseHybridFVMKernels::SolutionCheckKernel::Launch< parallelDevicePolicy<>,
                                                                          parallelDeviceReduce >( localSolution,
                                                                                                  dofManager.rankOffset(),
                                                                                                  faceDofNumber,
                                                                                                  faceGhostRank,
                                                                                                  facePres,
                                                                                                  dFacePres,
                                                                                                  scalingFactor );
  localCheck = std::min( localCheck, faceSolutionCheck );

  return MpiWrapper::Min( localCheck, MPI_COMM_GEOSX );
}

void CompositionalMultiphaseHybridFVM::ApplyBoundaryConditions( real64 const time_n,
                                                                real64 const dt,
                                                                DomainPartition & domain,
                                                                DofManager const & dofManager,
                                                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  CompositionalMultiphaseBase::ApplyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );
}

real64 CompositionalMultiphaseHybridFVM::CalculateResidualNorm( DomainPartition const & domain,
                                                                DofManager const & dofManager,
                                                                arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  FaceManager const & faceManager = *mesh.getFaceManager();

  // here we compute the cell-centered residual norm in the derived class
  // to avoid duplicating a synchronization point

  // get a view into local residual vector

  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString );
  string const faceDofKey = dofManager.getKey( viewKeyStruct::faceDofFieldString );

  globalIndex const rankOffset = dofManager.rankOffset();

  // local residual
  real64 localResidualNorm = 0;

  // 1. Compute the residual for the mass conservation equations

  forTargetSubRegionsComplete( mesh,
                               [&]( localIndex const,
                                    localIndex const,
                                    localIndex const,
                                    ElementRegionBase const &,
                                    ElementSubRegionBase const & subRegion )
  {
    arrayView1d< globalIndex const > const & elemDofNumber = subRegion.getReference< array1d< globalIndex > >( elemDofKey );
    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
    arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
    arrayView1d< real64 const > const & refPoro = subRegion.getReference< array1d< real64 > >( viewKeyStruct::referencePorosityString );
    arrayView1d< real64 const > const & totalDensOld = subRegion.getReference< array1d< real64 > >( viewKeyStruct::totalDensityOldString );

    CompositionalMultiphaseBaseKernels::ResidualNormKernel::Launch< parallelDevicePolicy<>,
                                                                    parallelDeviceReduce >( localRhs,
                                                                                            rankOffset,
                                                                                            numFluidComponents(),
                                                                                            elemDofNumber,
                                                                                            elemGhostRank,
                                                                                            refPoro,
                                                                                            volume,
                                                                                            totalDensOld,
                                                                                            localResidualNorm );
  } );

  arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();
  arrayView1d< globalIndex const > const & faceDofNumber =
    faceManager.getReference< array1d< globalIndex > >( faceDofKey );

  arrayView2d< localIndex const > const & elemRegionList    = faceManager.elementRegionList();
  arrayView2d< localIndex const > const & elemSubRegionList = faceManager.elementSubRegionList();
  arrayView2d< localIndex const > const & elemList          = faceManager.elementList();

  // 2. Compute the residual for the face-based constraints
  CompositionalMultiphaseHybridFVMKernels::ResidualNormKernel::Launch< parallelDevicePolicy<>,
                                                                       parallelDeviceReduce >( localRhs,
                                                                                               rankOffset,
                                                                                               numFluidPhases(),
                                                                                               faceDofNumber.toNestedViewConst(),
                                                                                               faceGhostRank.toNestedViewConst(),
                                                                                               m_regionFilter.toViewConst(),
                                                                                               elemRegionList.toNestedViewConst(),
                                                                                               elemSubRegionList.toNestedViewConst(),
                                                                                               elemList.toNestedViewConst(),
                                                                                               m_volume.toNestedViewConst(),
                                                                                               m_phaseMobOld.toNestedViewConst(),
                                                                                               localResidualNorm );

  // 3. Combine the two norms

  // compute global residual norm
  real64 const residual = std::sqrt( MpiWrapper::Sum( localResidualNorm ) );

  if( getLogLevel() >= 1 && logger::internal::rank==0 )
  {
    char output[200] = {0};
    sprintf( output, "    ( Rfluid ) = (%4.2e) ; ", residual );
    std::cout<<output;
  }

  return residual;
}

void CompositionalMultiphaseHybridFVM::ApplySystemSolution( DofManager const & dofManager,
                                                            arrayView1d< real64 const > const & localSolution,
                                                            real64 const scalingFactor,
                                                            DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  // 1. apply the elem-based update

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::elemDofFieldString,
                               viewKeyStruct::deltaPressureString,
                               scalingFactor,
                               0, 1 );

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::elemDofFieldString,
                               viewKeyStruct::deltaGlobalCompDensityString,
                               scalingFactor,
                               1, m_numDofPerCell );

  // if component density chopping is allowed, some component densities may be negative after the update
  // these negative component densities are set to zero in this function
  if( m_allowCompDensChopping )
  {
    ChopNegativeDensities( domain );
  }

  // 2. apply the face-based update

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::faceDofFieldString,
                               viewKeyStruct::deltaFacePressureString,
                               scalingFactor );

  // 3. synchronize

  std::map< string, string_array > fieldNames;
  fieldNames["face"].emplace_back( string( viewKeyStruct::deltaFacePressureString ) );
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaPressureString ) );
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaGlobalCompDensityString ) );
  CommunicationTools::SynchronizeFields( fieldNames,
                                         &mesh,
                                         domain.getNeighbors() );

  // 4. update secondary variables

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                  ElementSubRegionBase & subRegion )
  {
    UpdateState( subRegion, targetIndex );
  } );
}


void CompositionalMultiphaseHybridFVM::ResetStateToBeginningOfStep( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  // 1. Reset the cell-centered fields
  CompositionalMultiphaseBase::ResetStateToBeginningOfStep( domain );

  // 2. Reset the face-based fields
  MeshLevel & mesh = *domain.getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  FaceManager & faceManager = *mesh.getFaceManager();

  // get the accumulated face pressure updates
  arrayView1d< real64 > const & dFacePres =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::deltaFacePressureString );

  // zero out the face pressures
  dFacePres.setValues< parallelDevicePolicy<> >( 0.0 );
}


void CompositionalMultiphaseHybridFVM::UpdatePhaseMobility( Group & dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  // note that the phase mobility computed here does NOT include phase density

  // outputs

  arrayView2d< real64 > const phaseMob =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::phaseMobilityString );

  arrayView2d< real64 > const dPhaseMob_dPres =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::dPhaseMobility_dPressureString );

  arrayView3d< real64 > const dPhaseMob_dComp =
    dataGroup.getReference< array3d< real64 > >( viewKeyStruct::dPhaseMobility_dGlobalCompDensityString );

  // inputs

  arrayView2d< real64 const > const dPhaseVolFrac_dPres =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString );

  arrayView3d< real64 const > const dPhaseVolFrac_dComp =
    dataGroup.getReference< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

  arrayView3d< real64 const > const dCompFrac_dCompDens =
    dataGroup.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  MultiFluidBase const & fluid = GetConstitutiveModel< MultiFluidBase >( dataGroup, m_fluidModelNames[targetIndex] );

  arrayView3d< real64 const > const & phaseVisc = fluid.phaseViscosity();
  arrayView3d< real64 const > const & dPhaseVisc_dPres = fluid.dPhaseViscosity_dPressure();
  arrayView4d< real64 const > const & dPhaseVisc_dComp = fluid.dPhaseViscosity_dGlobalCompFraction();

  RelativePermeabilityBase const & relperm = GetConstitutiveModel< RelativePermeabilityBase >( dataGroup, m_relPermModelNames[targetIndex] );

  arrayView3d< real64 const > const & phaseRelPerm = relperm.phaseRelPerm();
  arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac = relperm.dPhaseRelPerm_dPhaseVolFraction();

  KernelLaunchSelector2< PhaseMobilityKernel >( m_numComponents, m_numPhases,
                                                dataGroup.size(),
                                                dCompFrac_dCompDens,
                                                phaseVisc,
                                                dPhaseVisc_dPres,
                                                dPhaseVisc_dComp,
                                                phaseRelPerm,
                                                dPhaseRelPerm_dPhaseVolFrac,
                                                dPhaseVolFrac_dPres,
                                                dPhaseVolFrac_dComp,
                                                phaseMob,
                                                dPhaseMob_dPres,
                                                dPhaseMob_dComp );
}


REGISTER_CATALOG_ENTRY( SolverBase, CompositionalMultiphaseHybridFVM, std::string const &, Group * const )
} /* namespace geosx */
