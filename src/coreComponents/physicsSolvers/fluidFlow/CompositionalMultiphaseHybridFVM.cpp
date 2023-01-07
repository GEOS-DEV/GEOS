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
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseHybridFVMKernels.hpp"
#include "physicsSolvers/fluidFlow/HybridFVMHelperKernels.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseHybridFVMKernels.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace isothermalCompositionalMultiphaseBaseKernels;
using namespace isothermalCompositionalMultiphaseHybridFVMKernels;
using namespace mimeticInnerProduct;

CompositionalMultiphaseHybridFVM::CompositionalMultiphaseHybridFVM( const std::string & name,
                                                                    Group * const parent ):
  CompositionalMultiphaseBase( name, parent ),
  m_lengthTolerance( 0 )
{
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
    MeshLevel & meshLevel = meshBody.getBaseDiscretization();

    FaceManager & faceManager = meshLevel.getFaceManager();

    // primary variables: face pressure changes

    faceManager.registerField< fields::flow::facePressure_n >( getName() );

    // auxiliary data for the buoyancy coefficient
    faceManager.registerField< fields::flow::mimGravityCoefficient >( getName() );
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

  m_lengthTolerance = domain.getMeshBody( 0 ).getGlobalLengthScale() * 1e-8;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
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
    arrayView1d< real64 const > const & transMultiplier = faceManager.getField< fields::flow::transMultiplier >();

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
    faceManager.getField< fields::flow::transMultiplier >();

  arrayView1d< real64 > const mimFaceGravCoef =
    faceManager.getField< fields::flow::mimGravityCoefficient >();

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
      subRegion.template getReference< array1d< real64 > >( fields::flow::gravityCoefficient::key() );
    arrayView1d< real64 const > const & elemVolume = subRegion.getElementVolume();
    arrayView2d< localIndex const > const & elemToFaces = subRegion.faceList();

    // here we precompute some quantities (mimFaceFracCoef) used in the FluxKernel to assemble the one-sided gravity term in the transport
    // scheme
    // This one-sided gravity term is currently always treated with TPFA, as in MRST.
    // In the future, I will change that (here and in the FluxKernel) to have a consistent inner product for the gravity term as well
    hybridFVMKernels::
      simpleKernelLaunchSelector< PrecomputeKernel,
                                  mimeticInnerProduct::TPFAInnerProduct >( subRegion.numFacesPerElement(),
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
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();

    arrayView1d< real64 > const & facePres_n =
      faceManager.getField< fields::flow::facePressure_n >();
    arrayView1d< real64 const > const & facePres =
      faceManager.getField< fields::flow::facePressure >();
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
                       FieldLocation::Elem,
                       m_numDofPerCell,
                       getMeshTargets() );

  dofManager.addCoupling( viewKeyStruct::elemDofFieldString(),
                          viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Face );

  // setup the connectivity of face fields
  dofManager.addField( viewKeyStruct::faceDofFieldString(),
                       FieldLocation::Face,
                       1,
                       getMeshTargets() );

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

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    // get the face-based DOF numbers for the assembly
    string const faceDofKey = dofManager.getKey( viewKeyStruct::faceDofFieldString() );
    // get the element dof numbers for the assembly
    string const & elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

    mesh.getElemManager().forElementSubRegionsComplete< CellElementSubRegion >( regionNames,
                                                                                [&]( localIndex const,
                                                                                     localIndex const er,
                                                                                     localIndex const esr,
                                                                                     ElementRegionBase const &,
                                                                                     CellElementSubRegion const & subRegion )
    {
      PermeabilityBase const & permeability =
        getConstitutiveModel< PermeabilityBase >( subRegion, subRegion.getReference< string >( viewKeyStruct::permeabilityNamesString() ) );

      isothermalCompositionalMultiphaseHybridFVMKernels::
        ElementBasedAssemblyKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                   m_numPhases,
                                                   dofManager.rankOffset(),
                                                   er,
                                                   esr,
                                                   m_lengthTolerance,
                                                   elemDofKey,
                                                   faceDofKey,
                                                   getName(),
                                                   nodeManager,
                                                   faceManager,
                                                   mesh.getElemManager(),
                                                   subRegion,
                                                   mimeticInnerProductBase,
                                                   permeability,
                                                   m_regionFilter.toViewConst(),
                                                   dt,
                                                   localMatrix,
                                                   localRhs );


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

  bool const skipCompFracDamping = m_maxCompFracChange >= 1.0;
  bool const skipPresDamping = m_maxRelativePresChange >= 1.0;

  // check if we want to rescale the Newton update
  if( skipCompFracDamping && skipPresDamping )
  {
    // no rescaling wanted, we just return 1.0;
    return 1.0;
  }

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  real64 scalingFactor = 1.0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames, [&]( localIndex const,
                                                                                          ElementSubRegionBase const & subRegion )
    {
      real64 const subRegionScalingFactor =
        isothermalCompositionalMultiphaseBaseKernels::
          ScalingForSystemSolutionKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_maxRelativePresChange,
                                                     m_maxCompFracChange,
                                                     dofManager.rankOffset(),
                                                     m_numComponents,
                                                     dofKey,
                                                     subRegion,
                                                     localSolution );

      scalingFactor = std::min( scalingFactor, subRegionScalingFactor );
    } );

    FaceManager const & faceManager = mesh.getFaceManager();

    string const faceDofKey = dofManager.getKey( viewKeyStruct::faceDofFieldString() );
    real64 const maxRelativePresChange = m_maxRelativePresChange;

    arrayView1d< globalIndex const > const & faceDofNumber =
      faceManager.getReference< array1d< globalIndex > >( faceDofKey );
    arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();
    arrayView1d< real64 const > const & facePressure =
      faceManager.getField< fields::flow::facePressure >();
    globalIndex const rankOffset = dofManager.rankOffset();

    RAJA::ReduceMin< parallelDeviceReduce, real64 > minFaceVal( 1.0 );
    forAll< parallelDevicePolicy<> >( faceManager.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iface )
    {
      if( faceGhostRank[iface] < 0 && faceDofNumber[iface] >= 0 )
      {
        real64 const facePres = facePressure[iface];
        real64 const absPresChange = LvArray::math::abs( localSolution[faceDofNumber[iface] - rankOffset] );
        if( facePres > isothermalCompositionalMultiphaseBaseKernels::minDensForDivision )
        {
          real64 const relativePresChange = absPresChange / facePres;
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
  } );
  return LvArray::math::max( MpiWrapper::min( scalingFactor ), m_minScalingFactor );
}


bool CompositionalMultiphaseHybridFVM::checkSystemSolution( DomainPartition const & domain,
                                                            DofManager const & dofManager,
                                                            arrayView1d< real64 const > const & localSolution,
                                                            real64 const scalingFactor )
{
  GEOSX_MARK_FUNCTION;

  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  integer localCheck = 1;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    // 1. check cell-centered variables for each region
    mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames, [&]( localIndex const,
                                                                                          ElementSubRegionBase const & subRegion )
    {
      // check that pressure and component densities are non-negative
      integer const subRegionSolutionCheck =
        isothermalCompositionalMultiphaseBaseKernels::
          SolutionCheckKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_allowCompDensChopping,
                                                     scalingFactor,
                                                     dofManager.rankOffset(),
                                                     m_numComponents,
                                                     elemDofKey,
                                                     subRegion,
                                                     localSolution );

      localCheck = std::min( localCheck, subRegionSolutionCheck );
    } );

    // 2. check face-centered variables in the domain

    FaceManager const & faceManager = mesh.getFaceManager();
    string const faceDofKey = dofManager.getKey( viewKeyStruct::faceDofFieldString() );
    arrayView1d< integer const > const faceGhostRank = faceManager.ghostRank();
    arrayView1d< globalIndex const > const faceDofNumber =
      faceManager.getReference< array1d< globalIndex > >( faceDofKey );
    arrayView1d< real64 const > const facePres =
      faceManager.getField< fields::flow::facePressure >();

    // check that face pressure is non-negative
    integer const faceSolutionCheck =
      isothermalCompositionalMultiphaseHybridFVMKernels::
        SolutionCheckKernel::launch< parallelDevicePolicy<> >( localSolution,
                                                               dofManager.rankOffset(),
                                                               faceDofNumber,
                                                               faceGhostRank,
                                                               facePres,
                                                               scalingFactor );
    localCheck = std::min( localCheck, faceSolutionCheck );

  } );

  return MpiWrapper::min( localCheck );
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

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    FaceManager const & faceManager = mesh.getFaceManager();

    // here we compute the cell-centered residual norm in the derived class
    // to avoid duplicating a synchronization point

    // get a view into local residual vector

    string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
    string const faceDofKey = dofManager.getKey( viewKeyStruct::faceDofFieldString() );


    StencilAccessors< fields::elementVolume,
                      fields::flow::phaseMobility_n >
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
    isothermalCompositionalMultiphaseHybridFVMKernels::
      ResidualNormKernel::launch< parallelDevicePolicy<> >( localRhs,
                                                            rankOffset,
                                                            numFluidPhases(),
                                                            faceDofNumber.toNestedViewConst(),
                                                            faceGhostRank.toNestedViewConst(),
                                                            m_regionFilter.toViewConst(),
                                                            elemRegionList.toNestedViewConst(),
                                                            elemSubRegionList.toNestedViewConst(),
                                                            elemList.toNestedViewConst(),
                                                            compFlowAccessors.get( fields::elementVolume{} ),
                                                            compFlowAccessors.get( fields::flow::phaseMobility_n{} ),
                                                            faceResidualNorm );
    localResidualNorm += faceResidualNorm;
  } );

  // 3. Combine the two norms

  // compute global residual norm
  real64 const residual = std::sqrt( MpiWrapper::sum( localResidualNorm ) );
  if( getLogLevel() >= 1 && logger::internal::rank == 0 )
  {
    GEOSX_FMT( "    ( R{} ) = ( {:4.2e} ) ; ", coupledSolverAttributePrefix(), residual );
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
                               fields::flow::pressure::key(),
                               scalingFactor,
                               pressureMask );

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::elemDofFieldString(),
                               fields::flow::globalCompDensity::key(),
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
                               fields::flow::facePressure::key(),
                               scalingFactor );

  // 3. synchronize
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    FieldIdentifiers fieldsToBeSync;

    {
      fieldsToBeSync.addElementFields( { fields::flow::pressure::key(),
                                         fields::flow::globalCompDensity::key() },
                                       regionNames );

      fieldsToBeSync.addFields( FieldLocation::Face, { fields::flow::facePressure::key() } );
    };

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
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
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();

    arrayView1d< real64 const > const & facePres_n =
      faceManager.getField< fields::flow::facePressure_n >();
    arrayView1d< real64 > const & facePres =
      faceManager.getField< fields::flow::facePressure >();
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

  isothermalCompositionalMultiphaseHybridFVMKernels::
    PhaseMobilityKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                               m_numPhases,
                                               dataGroup,
                                               fluid,
                                               relperm );
}

REGISTER_CATALOG_ENTRY( SolverBase, CompositionalMultiphaseHybridFVM, std::string const &, Group * const )
} /* namespace geosx */
