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
 * @file CompositionalMultiphaseHybridFVM.cpp
 */

#include "CompositionalMultiphaseHybridFVM.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "finiteVolume/HybridMimeticDiscretization.hpp"
#include "finiteVolume/MimeticInnerProductDispatch.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseHybridFVMKernels.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseHybridFVMKernels.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace isothermalCompositionalMultiphaseBaseKernels;
using namespace compositionalMultiphaseHybridFVMKernels;
using namespace mimeticInnerProduct;

CompositionalMultiphaseHybridFVM::CompositionalMultiphaseHybridFVM( const std::string & name,
                                                                    Group * const parent ):
  CompositionalMultiphaseBase( name, parent ),
  m_lengthTolerance( 0 )
{

  this->registerWrapper( viewKeyStruct::maxRelativePresChangeString(), &m_maxRelativePresChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Maximum (relative) change in (face) pressure between two Newton iterations" );

  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseHybridFVM;

}

void CompositionalMultiphaseHybridFVM::registerDataOnMesh( Group & meshBodies )
{
  GEOSX_MARK_FUNCTION;

  // 1) Register the elem-centered data
  CompositionalMultiphaseBase::registerDataOnMesh( meshBodies );

  // 2) Register the face data
  meshBodies.forSubGroups< MeshBody >( [&]( MeshBody & meshBody )
  {
    MeshLevel & meshLevel = meshBody.getMeshLevel( 0 );

    FaceManager & faceManager = meshLevel.getFaceManager();

    // primary variables: face pressure changes

    faceManager.registerExtrinsicData< extrinsicMeshData::flow::facePressure_n >( getName() );

    // auxiliary data for the buoyancy coefficient
    faceManager.registerExtrinsicData< extrinsicMeshData::flow::mimGravityCoefficient >( getName() );
  } );
}

void CompositionalMultiphaseHybridFVM::initializePreSubGroups()
{
  CompositionalMultiphaseBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();

  GEOSX_THROW_IF( !fvManager.hasGroup< HybridMimeticDiscretization >( m_discretizationName ),
                  catalogName() << " " << getName() <<
                  ": the HybridMimeticDiscretization must be selected with CompositionalMultiphaseHybridFVM",
                  InputError );

  GEOSX_THROW_IF( m_hasCapPressure,
                  catalogName() << " " << getName() <<
                  ": capillary pressure is not yet supported by CompositionalMultiphaseHybridFVM",
                  InputError );
}

void CompositionalMultiphaseHybridFVM::initializePostInitialConditionsPreSubGroups()
{
  GEOSX_MARK_FUNCTION;

  if( m_computeCFLNumbers )
  {
    GEOSX_LOG_RANK_0( catalogName() << " " << getName()
                                    << ": the computation of CFL numbers in not supported by CompositionalMultiphaseHybridFVM yet" );
  }

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  HybridMimeticDiscretization const & hmDiscretization = fvManager.getHybridMimeticDiscretization( m_discretizationName );
  MimeticInnerProductBase const & mimeticInnerProductBase =
    hmDiscretization.getReference< MimeticInnerProductBase >( HybridMimeticDiscretization::viewKeyStruct::innerProductString() );
  if( dynamicCast< QuasiRTInnerProduct const * >( &mimeticInnerProductBase )  ||
      dynamicCast< QuasiTPFAInnerProduct const * >( &mimeticInnerProductBase )  ||
      dynamicCast< SimpleInnerProduct const * >( &mimeticInnerProductBase ) )
  {
    GEOSX_ERROR( "The QuasiRT, QuasiTPFA, and Simple inner products are only available in SinglePhaseHybridFVM" );
  }

  m_transMultName = viewKeyStruct::transMultiplierString();

  m_lengthTolerance = domain.getMeshBody( 0 ).getGlobalLengthScale() * 1e-8;

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    CompositionalMultiphaseBase::initializePostInitialConditionsPreSubGroups();

    // in the flux kernel, we need to make sure that we act only on the target regions
    // for that, we need the following region filter
    for( string const & regionName : regionNames )
    {
      m_regionFilter.insert( elemManager.getRegions().getIndex( regionName ) );
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

    GEOSX_THROW_IF( minVal.get() <= 0.0,
                    catalogName() << " " << getName()
                                  << ": the transmissibility multipliers used in SinglePhaseHybridFVM must strictly larger than 0.0",
                    std::runtime_error );

    FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
    fsManager.forSubGroups< AquiferBoundaryCondition >( [&] ( AquiferBoundaryCondition const & bc )
    {
      GEOSX_LOG_RANK_0( catalogName() << " " << getName() <<
                        "An aquifer boundary condition named " << bc.getName() << " was requested in the XML file. \n"
                                                                                  "This type of boundary condition is not yet supported by CompositionalMultiphaseHybridFVM and will be ignored" );
    } );
  } );

}

void CompositionalMultiphaseHybridFVM::precomputeData( MeshLevel & mesh, arrayView1d< string const > const & regionNames )
{
  FlowSolverBase::precomputeData( mesh, regionNames );

  NodeManager const & nodeManager = mesh.getNodeManager();
  FaceManager & faceManager = mesh.getFaceManager();

  array1d< RAJA::ReduceSum< serialReduce, real64 > > mimFaceGravCoefNumerator;
  array1d< RAJA::ReduceSum< serialReduce, real64 > > mimFaceGravCoefDenominator;
  mimFaceGravCoefNumerator.resize( faceManager.size() );
  mimFaceGravCoefDenominator.resize( faceManager.size() );

  // node data

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  // face data

  arrayView1d< real64 const > const & transMultiplier =
    faceManager.getReference< array1d< real64 > >( m_transMultName );

  arrayView1d< real64 > const mimFaceGravCoef =
    faceManager.getExtrinsicData< extrinsicMeshData::flow::mimGravityCoefficient >();

  ArrayOfArraysView< localIndex const > const & faceToNodes = faceManager.nodeList().toViewConst();

  real64 const lengthTolerance = m_lengthTolerance;

  mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                        CellElementSubRegion & subRegion )
  {
    arrayView2d< real64 const > const & elemCenter =
      subRegion.template getReference< array2d< real64 > >( CellElementSubRegion::viewKeyStruct::elementCenterString() );
    string & permModelName = subRegion.getReference< string >( viewKeyStruct::permeabilityNamesString() );
    arrayView3d< real64 const > const & elemPerm =
      getConstitutiveModel< PermeabilityBase >( subRegion, permModelName ).permeability();
    arrayView1d< real64 const > const elemGravCoef =
      subRegion.template getReference< array1d< real64 > >( extrinsicMeshData::flow::gravityCoefficient::key() );
    arrayView1d< real64 const > const & elemVolume = subRegion.getElementVolume();
    arrayView2d< localIndex const > const & elemToFaces = subRegion.faceList();

    // here we precompute some quantities (mimFaceFracCoef) used in the FluxKernel to assemble the one-sided gravity term in the transport
    // scheme
    // This one-sided gravity term is currently always treated with TPFA, as in MRST.
    // In the future, I will change that (here and in the FluxKernel) to have a consistent inner product for the gravity term as well
    singlePhaseHybridFVMKernels::KernelLaunchSelector< mimeticInnerProduct::TPFAInnerProduct,
                                                       PrecomputeKernel >( subRegion.numFacesPerElement(),
                                                                           subRegion.size(),
                                                                           faceManager.size(),
                                                                           nodePosition,
                                                                           faceToNodes,
                                                                           elemCenter,
                                                                           elemVolume,
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

void CompositionalMultiphaseHybridFVM::implicitStepSetup( real64 const & time_n,
                                                          real64 const & dt,
                                                          DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  // setup the elem-centered fields
  CompositionalMultiphaseBase::implicitStepSetup( time_n, dt, domain );

  // setup the face fields
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();

    arrayView1d< real64 > const & facePres_n =
      faceManager.getExtrinsicData< extrinsicMeshData::flow::facePressure_n >();
    arrayView1d< real64 const > const & facePres =
      faceManager.getExtrinsicData< extrinsicMeshData::flow::facePressure >();
    facePres_n.setValues< parallelDevicePolicy<> >( facePres );
  } );

}


void CompositionalMultiphaseHybridFVM::setupDofs( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                                  DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;

  // setup the connectivity of elem fields
  // we need Connectivity::Face because of the two-point upwinding
  // in AssembleOneSidedMassFluxes
  dofManager.addField( viewKeyStruct::elemDofFieldString(),
                       DofManager::Location::Elem,
                       m_numDofPerCell,
                       m_meshTargets );

  dofManager.addCoupling( viewKeyStruct::elemDofFieldString(),
                          viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Face );

  // setup the connectivity of face fields
  dofManager.addField( viewKeyStruct::faceDofFieldString(),
                       DofManager::Location::Face,
                       1,
                       m_meshTargets );

  dofManager.addCoupling( viewKeyStruct::faceDofFieldString(),
                          viewKeyStruct::faceDofFieldString(),
                          DofManager::Connector::Elem );

  // setup coupling between pressure and face pressure
  dofManager.addCoupling( viewKeyStruct::faceDofFieldString(),
                          viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Elem );

}


void CompositionalMultiphaseHybridFVM::assembleFluxTerms( real64 const dt,
                                                          DomainPartition const & domain,
                                                          DofManager const & dofManager,
                                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                          arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  HybridMimeticDiscretization const & hmDiscretization = fvManager.getHybridMimeticDiscretization( m_discretizationName );
  MimeticInnerProductBase const & mimeticInnerProductBase =
    hmDiscretization.getReference< MimeticInnerProductBase >( HybridMimeticDiscretization::viewKeyStruct::innerProductString() );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel const & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    // node data (for transmissibility computation)

    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

    // face data

    // get the face-based DOF numbers for the assembly
    string const faceDofKey = dofManager.getKey( viewKeyStruct::faceDofFieldString() );
    arrayView1d< globalIndex const > const & faceDofNumber =
      faceManager.getReference< array1d< globalIndex > >( faceDofKey );
    arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();

    // get the element dof numbers for the assembly
    string const & elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > elemDofNumber =
      mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
    elemDofNumber.setName( getName() + "/accessors/" + elemDofKey );

    // get the face-centered pressures
    arrayView1d< real64 const > const & facePres =
      faceManager.getExtrinsicData< extrinsicMeshData::flow::facePressure >();

    // get the face-centered depth
    arrayView1d< real64 const > const & faceGravCoef =
      faceManager.getExtrinsicData< extrinsicMeshData::flow::gravityCoefficient >();
    arrayView1d< real64 const > const & mimFaceGravCoef =
      faceManager.getExtrinsicData< extrinsicMeshData::flow::mimGravityCoefficient >();

    // get the face-centered transMultiplier
    arrayView1d< real64 const > const & transMultiplier =
      faceManager.getReference< array1d< real64 > >( m_transMultName );

    // get the face-to-nodes connectivity for the transmissibility calculation
    ArrayOfArraysView< localIndex const > const & faceToNodes = faceManager.nodeList().toViewConst();

    arrayView2d< localIndex const > const & elemRegionList    = faceManager.elementRegionList().toViewConst();
    arrayView2d< localIndex const > const & elemSubRegionList = faceManager.elementSubRegionList().toViewConst();
    arrayView2d< localIndex const > const & elemList          = faceManager.elementList().toViewConst();


    // tolerance for transmissibility calculation
    real64 const lengthTolerance = m_lengthTolerance;

    FluxKernel::CompFlowAccessors compFlowAccessors( mesh.getElemManager(), getName() );
    FluxKernel::MultiFluidAccessors multiFluidAccessors( mesh.getElemManager(), getName() );

    mesh.getElemManager().forElementSubRegionsComplete< CellElementSubRegion >( regionNames,
                                                                                [&]( localIndex const,
                                                                                     localIndex const er,
                                                                                     localIndex const esr,
                                                                                     ElementRegionBase const &,
                                                                                     CellElementSubRegion const & subRegion )
    {
      PermeabilityBase const & permeabilityModel =
        getConstitutiveModel< PermeabilityBase >( subRegion, subRegion.getReference< string >( viewKeyStruct::permeabilityNamesString() ) );

      mimeticInnerProductReducedDispatch( mimeticInnerProductBase,
                                          [&] ( auto const mimeticInnerProduct )
      {
        using IP_TYPE = TYPEOFREF( mimeticInnerProduct );
        KernelLaunchSelector< FluxKernel,
                              IP_TYPE >( subRegion.numFacesPerElement(),
                                         m_numComponents, m_numPhases,
                                         er, esr, subRegion,
                                         permeabilityModel,
                                         m_regionFilter.toViewConst(),
                                         nodePosition,
                                         elemRegionList,
                                         elemSubRegionList,
                                         elemList,
                                         faceToNodes,
                                         faceDofNumber,
                                         faceGhostRank,
                                         facePres,
                                         faceGravCoef,
                                         mimFaceGravCoef,
                                         transMultiplier,
                                         compFlowAccessors.get( extrinsicMeshData::flow::phaseMobility{} ),
                                         compFlowAccessors.get( extrinsicMeshData::flow::dPhaseMobility_dPressure{} ),
                                         compFlowAccessors.get( extrinsicMeshData::flow::dPhaseMobility_dGlobalCompDensity{} ),
                                         compFlowAccessors.get( extrinsicMeshData::flow::dGlobalCompFraction_dGlobalCompDensity{} ),
                                         multiFluidAccessors.get( extrinsicMeshData::multifluid::phaseDensity{} ),
                                         multiFluidAccessors.get( extrinsicMeshData::multifluid::dPhaseDensity{} ),
                                         multiFluidAccessors.get( extrinsicMeshData::multifluid::phaseMassDensity{} ),
                                         multiFluidAccessors.get( extrinsicMeshData::multifluid::dPhaseMassDensity{} ),
                                         multiFluidAccessors.get( extrinsicMeshData::multifluid::phaseCompFraction{} ),
                                         multiFluidAccessors.get( extrinsicMeshData::multifluid::dPhaseCompFraction{} ),
                                         elemDofNumber.toNestedViewConst(),
                                         dofManager.rankOffset(),
                                         lengthTolerance,
                                         dt,
                                         localMatrix,
                                         localRhs );

      } );
    } );

  } );
}

void CompositionalMultiphaseHybridFVM::assembleStabilizedFluxTerms( real64 const dt,
                                                                    DomainPartition const & domain,
                                                                    DofManager const & dofManager,
                                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                    arrayView1d< real64 > const & localRhs ) const
{
  // stab not implemented
  assembleFluxTerms( dt, domain, dofManager, localMatrix, localRhs );
}

real64 CompositionalMultiphaseHybridFVM::scalingForSystemSolution( DomainPartition const & domain,
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

  real64 constexpr eps = isothermalCompositionalMultiphaseBaseKernels::minDensForDivision;
  real64 const maxCompFracChange = m_maxCompFracChange;
  real64 const maxRelativePresChange = m_maxRelativePresChange;

  localIndex const NC = m_numComponents;

  real64 solutionScaling;

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel const & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    FaceManager const & faceManager = mesh.getFaceManager();

    string const faceDofKey = dofManager.getKey( viewKeyStruct::faceDofFieldString() );
    arrayView1d< globalIndex const > const & faceDofNumber =
      faceManager.getReference< array1d< globalIndex > >( faceDofKey );
    arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();
    arrayView1d< real64 const > const & facePressure =
      faceManager.getExtrinsicData< extrinsicMeshData::flow::facePressure >();
    globalIndex const rankOffset = dofManager.rankOffset();
    string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
    real64 scalingFactor = 1.0;

    mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames, [&]( localIndex const,
                                                                                          ElementSubRegionBase const & subRegion )
    {
      arrayView1d< globalIndex const > const & dofNumber =
        subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

      arrayView1d< real64 const > const & pressure =
        subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::pressure::key() );
      arrayView2d< real64 const, compflow::USD_COMP > const & compDens =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompDensity >();

      RAJA::ReduceMin< parallelDeviceReduce, real64 > minElemVal( 1.0 );

      forAll< parallelDevicePolicy<> >( dofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
        if( elemGhostRank[ei] < 0 )
        {

          // compute the change in pressure
          real64 const pres = pressure[ei];
          real64 const absPresChange = LvArray::math::abs( localSolution[dofNumber[ei] - rankOffset] );
          if( pres > eps )
          {
            real64 const relativePresChange = LvArray::math::abs( absPresChange ) / pres;
            if( relativePresChange > maxRelativePresChange )
            {
              minElemVal.min( maxRelativePresChange / relativePresChange );
            }
          }

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
            real64 const absCompDensChange = LvArray::math::abs( localSolution[lid] );
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
        real64 const facePres = facePressure[iface];
        real64 const absPresChange = LvArray::math::abs( localSolution[faceDofNumber[iface] - rankOffset] );
        if( facePres > eps )
        {
          real64 const relativePresChange = LvArray::math::abs( absPresChange ) / facePres;
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

    solutionScaling = LvArray::math::max( MpiWrapper::min( scalingFactor, MPI_COMM_GEOSX ), m_minScalingFactor );
  } );

  return solutionScaling;
}


bool CompositionalMultiphaseHybridFVM::checkSystemSolution( DomainPartition const & domain,
                                                            DofManager const & dofManager,
                                                            arrayView1d< real64 const > const & localSolution,
                                                            real64 const scalingFactor )
{
  GEOSX_MARK_FUNCTION;

  // 1. check cell-centered variables for each region
  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  localIndex localCheck = 1;

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel const & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    FaceManager const & faceManager = mesh.getFaceManager();

    // 1. check cell-centered variables for each region
    mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames, [&]( localIndex const,
                                                                                          ElementSubRegionBase const & subRegion )
    {
      arrayView1d< globalIndex const > const & elemDofNumber = subRegion.getReference< array1d< globalIndex > >( elemDofKey );
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

      arrayView1d< real64 const > const & elemPres =
        subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::pressure::key() );
      arrayView2d< real64 const, compflow::USD_COMP > const & compDens =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompDensity >();

      localIndex const subRegionSolutionCheck =
        isothermalCompositionalMultiphaseBaseKernels::
          SolutionCheckKernel::launch< parallelDevicePolicy<> >( localSolution,
                                                                 dofManager.rankOffset(),
                                                                 numFluidComponents(),
                                                                 elemDofNumber,
                                                                 elemGhostRank,
                                                                 elemPres,
                                                                 compDens,
                                                                 m_allowCompDensChopping,
                                                                 scalingFactor );
      localCheck = std::min( localCheck, subRegionSolutionCheck );
    } );

    // 2. check face-centered variables in the domain

    string const faceDofKey = dofManager.getKey( viewKeyStruct::faceDofFieldString() );
    arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();
    arrayView1d< globalIndex const > const & faceDofNumber =
      faceManager.getReference< array1d< globalIndex > >( faceDofKey );

    arrayView1d< real64 const > const & facePres =
      faceManager.getExtrinsicData< extrinsicMeshData::flow::facePressure >();

    localIndex const faceSolutionCheck =
      compositionalMultiphaseHybridFVMKernels::
        SolutionCheckKernel::launch< parallelDevicePolicy<> >( localSolution,
                                                               dofManager.rankOffset(),
                                                               faceDofNumber,
                                                               faceGhostRank,
                                                               facePres,
                                                               scalingFactor );
    localCheck = std::min( localCheck, faceSolutionCheck );

  } );

  return MpiWrapper::min( localCheck, MPI_COMM_GEOSX );
}

void CompositionalMultiphaseHybridFVM::applyBoundaryConditions( real64 const time_n,
                                                                real64 const dt,
                                                                DomainPartition & domain,
                                                                DofManager const & dofManager,
                                                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  CompositionalMultiphaseBase::applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );

  // TODO: implement face boundary conditions here
}

void CompositionalMultiphaseHybridFVM::applyAquiferBC( real64 const time,
                                                       real64 const dt,
                                                       DofManager const & dofManager,
                                                       DomainPartition & domain,
                                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                       arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  GEOSX_UNUSED_VAR( time, dt, dofManager, domain, localMatrix, localRhs );

}

void CompositionalMultiphaseHybridFVM::saveAquiferConvergedState( real64 const & time,
                                                                  real64 const & dt,
                                                                  DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  GEOSX_UNUSED_VAR( time, dt, domain );
}


real64 CompositionalMultiphaseHybridFVM::calculateResidualNorm( DomainPartition const & domain,
                                                                DofManager const & dofManager,
                                                                arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  // local residual
  real64 localResidualNorm = 0;
  globalIndex const rankOffset = dofManager.rankOffset();

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel const & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    FaceManager const & faceManager = mesh.getFaceManager();

    // here we compute the cell-centered residual norm in the derived class
    // to avoid duplicating a synchronization point

    // get a view into local residual vector

    string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
    string const faceDofKey = dofManager.getKey( viewKeyStruct::faceDofFieldString() );


    StencilAccessors< extrinsicMeshData::elementVolume,
                      extrinsicMeshData::flow::phaseMobility_n >
    compFlowAccessors( mesh.getElemManager(), getName() );


    // 1. Compute the residual for the mass conservation equations

    mesh.getElemManager().forElementSubRegionsComplete< ElementSubRegionBase >( regionNames,
                                                                                [&]( localIndex const,
                                                                                     localIndex const,
                                                                                     localIndex const,
                                                                                     ElementRegionBase const &,
                                                                                     ElementSubRegionBase const & subRegion )
    {
      arrayView1d< globalIndex const > const & elemDofNumber = subRegion.getReference< array1d< globalIndex > >( elemDofKey );
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const & volume = subRegion.getElementVolume();

      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
      arrayView2d< real64 const, multifluid::USD_FLUID > const & totalDens_n = fluid.totalDensity_n();

      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
      CoupledSolidBase const & solid = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );
      arrayView1d< real64 const > const & referencePorosity = solid.getReferencePorosity();

      real64 subRegionResidualNorm = 0.0;
      isothermalCompositionalMultiphaseBaseKernels::
        ResidualNormKernel::launch< parallelDevicePolicy<> >( localRhs,
                                                              rankOffset,
                                                              numFluidComponents(),
                                                              elemDofNumber,
                                                              elemGhostRank,
                                                              referencePorosity,
                                                              volume,
                                                              totalDens_n,
                                                              subRegionResidualNorm );
      localResidualNorm += subRegionResidualNorm;
    } );

    arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();
    arrayView1d< globalIndex const > const & faceDofNumber =
      faceManager.getReference< array1d< globalIndex > >( faceDofKey );

    arrayView2d< localIndex const > const & elemRegionList    = faceManager.elementRegionList();
    arrayView2d< localIndex const > const & elemSubRegionList = faceManager.elementSubRegionList();
    arrayView2d< localIndex const > const & elemList          = faceManager.elementList();

    // 2. Compute the residual for the face-based constraints
    real64 faceResidualNorm = 0.0;
    compositionalMultiphaseHybridFVMKernels::
      ResidualNormKernel::launch< parallelDevicePolicy<> >( localRhs,
                                                            rankOffset,
                                                            numFluidPhases(),
                                                            faceDofNumber.toNestedViewConst(),
                                                            faceGhostRank.toNestedViewConst(),
                                                            m_regionFilter.toViewConst(),
                                                            elemRegionList.toNestedViewConst(),
                                                            elemSubRegionList.toNestedViewConst(),
                                                            elemList.toNestedViewConst(),
                                                            compFlowAccessors.get( extrinsicMeshData::elementVolume{} ),
                                                            compFlowAccessors.get( extrinsicMeshData::flow::phaseMobility_n{} ),
                                                            faceResidualNorm );
    localResidualNorm += faceResidualNorm;
  } );

  // 3. Combine the two norms

  // compute global residual norm
  real64 const residual = std::sqrt( MpiWrapper::sum( localResidualNorm ) );

  if( getLogLevel() >= 1 && logger::internal::rank == 0 )
  {
    std::cout << GEOSX_FMT( "    ( Rfluid ) = ( {:4.2e} ) ;", residual );
  }

  return residual;
}

void CompositionalMultiphaseHybridFVM::applySystemSolution( DofManager const & dofManager,
                                                            arrayView1d< real64 const > const & localSolution,
                                                            real64 const scalingFactor,
                                                            DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  // 1. apply the elem-based update
  DofManager::CompMask pressureMask( m_numDofPerCell, 0, 1 );

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::elemDofFieldString(),
                               extrinsicMeshData::flow::pressure::key(),
                               scalingFactor,
                               pressureMask );

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::elemDofFieldString(),
                               extrinsicMeshData::flow::globalCompDensity::key(),
                               scalingFactor,
                               ~pressureMask );

  // if component density chopping is allowed, some component densities may be negative after the update
  // these negative component densities are set to zero in this function
  if( m_allowCompDensChopping )
  {
    chopNegativeDensities( domain );
  }

  // 2. apply the face-based update

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::faceDofFieldString(),
                               extrinsicMeshData::flow::facePressure::key(),
                               scalingFactor );

  // 3. synchronize
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & )
  {
    std::map< string, string_array > fieldNames;
    fieldNames["face"].emplace_back( extrinsicMeshData::flow::facePressure::key() );
    fieldNames["elems"].emplace_back( extrinsicMeshData::flow::pressure::key() );
    fieldNames["elems"].emplace_back( extrinsicMeshData::flow::globalCompDensity::key() );
    CommunicationTools::getInstance().synchronizeFields( fieldNames,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );
}


void CompositionalMultiphaseHybridFVM::resetStateToBeginningOfStep( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  // 1. Reset the cell-centered fields
  CompositionalMultiphaseBase::resetStateToBeginningOfStep( domain );

  // 2. Reset the face-based fields
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();

    arrayView1d< real64 const > const & facePres_n =
      faceManager.getExtrinsicData< extrinsicMeshData::flow::facePressure_n >();
    arrayView1d< real64 > const & facePres =
      faceManager.getExtrinsicData< extrinsicMeshData::flow::facePressure >();
    facePres.setValues< parallelDevicePolicy<> >( facePres_n );
  } );
}

void CompositionalMultiphaseHybridFVM::updatePhaseMobility( ObjectManagerBase & dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  MultiFluidBase const & fluid =
    getConstitutiveModel< MultiFluidBase >( dataGroup,
                                            dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() ) );
  RelativePermeabilityBase const & relperm =
    getConstitutiveModel< RelativePermeabilityBase >( dataGroup,
                                                      dataGroup.getReference< string >( viewKeyStruct::relPermNamesString() ) );

  compositionalMultiphaseHybridFVMKernels::
    PhaseMobilityKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                               m_numPhases,
                                               dataGroup,
                                               fluid,
                                               relperm );
}

REGISTER_CATALOG_ENTRY( SolverBase, CompositionalMultiphaseHybridFVM, std::string const &, Group * const )
} /* namespace geosx */
