/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhaseMixedMFD.cpp
 */

#include "SinglePhaseMixedMFD.hpp"

#include "common/logger/Logger.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "fieldSpecification/FieldSpecificationImpl.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mixedMimetic/MixedMimeticDiscretization.hpp"
#include "mixedMimetic/MixedMimeticDiscretizationManager.hpp"
#include "mixedMimetic/MixedMimeticFields.hpp"
#include "mixedMimetic/adaptivity/GlobalAdaptationKernels.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/SinglePhaseMixedMFDKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/ResidualNormKernel.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace fields;
using namespace mimeticInnerProduct;

SinglePhaseMixedMFD::SinglePhaseMixedMFD( const string & name,
                                          Group * const parent ):
  SinglePhaseBase( name, parent ),
  m_areaRelTol( 1e-8 )
{
  // one cell-centered dof per cell
  m_numDofPerCell = 1;
  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhaseMixedMFD;
}

void SinglePhaseMixedMFD::registerDataOnMesh( Group & meshBodies )
{
  SinglePhaseBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    string_array const & regionNames )
  {
    // 1) Register the cell-centered adaptation data
    ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      subRegion.registerField< mixedMimetic::stencilFlag >( getName() );
      subRegion.registerField< mixedMimetic::consistencyIndicator >( getName() );
    } );

    // 2) Register the face data
    FaceManager & faceManager = mesh.getFaceManager();
    {
      // primary variables: face mass fluxes
      faceManager.registerField< mixedMimetic::faceMassFlux >( getName() );
      faceManager.registerField< mixedMimetic::faceMassFlux_n >( getName() );

      // boundary condition data
      faceManager.registerField< flow::bcPressure >( getName() );
      faceManager.registerField< flow::isBoundaryFace >( getName() );

      // Global Adaptation face residual
      faceManager.registerField< mixedMimetic::faceResidual >( getName() );
    }
  } );
}

void SinglePhaseMixedMFD::initializePreSubGroups()
{
  SinglePhaseBase::initializePreSubGroups();

  GEOS_THROW_IF( m_isThermal,
                 "The thermal option is not supported by SinglePhaseMixedMFD",
                 InputError, getDataContext() );

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  MixedMimeticDiscretizationManager const & mmManager = numericalMethodManager.getMixedMimeticDiscretizationManager();

  GEOS_THROW_IF( !mmManager.hasGroup< MixedMimeticDiscretization >( m_discretizationName ),
                 "A MixedMimeticDiscretization must be selected with SinglePhaseMixedMFD",
                 InputError, getDataContext() );
}

void SinglePhaseMixedMFD::initializePostInitialConditionsPreSubGroups()
{
  GEOS_MARK_FUNCTION;

  SinglePhaseBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    GEOS_UNUSED_VAR( regionNames );

    ElementRegionManager const & elemManager = mesh.getElemManager();

    // in the kernels, we need to make sure that we act only on the target regions
    for( string const & regionName : regionNames )
    {
      m_regionFilter.insert( elemManager.getRegions().getIndex( regionName ) );
    }

    // flag the faces on which a pressure boundary condition is imposed:
    // all the remaining domain-boundary faces are treated as no-flow faces
    FaceManager & faceManager = mesh.getFaceManager();
    arrayView1d< integer > const isPresBcFace = faceManager.getField< flow::isBoundaryFace >();

    fsManager.apply< FaceManager >( 0.0,
                                    mesh,
                                    flow::bcPressure::key(),
                                    [&] ( FieldSpecification const &,
                                          string const &,
                                          SortedArrayView< localIndex const > const & targetSet,
                                          FaceManager &,
                                          string const & )
    {
      forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        isPresBcFace[targetSet[a]] = 1;
      } );
    } );

    fsManager.forSubGroups< AquiferBoundaryCondition >( [&] ( AquiferBoundaryCondition const & bc )
    {
      GEOS_UNUSED_VAR( bc );
      GEOS_WARNING( "The aquifer boundary condition was requested in the XML file. \n"
                    "This type of boundary condition is not yet supported by SinglePhaseMixedMFD and will be ignored",
                    getDataContext(), bc.getDataContext() );
    } );
  } );

  // run the residual-based Global Adaptation pipeline (or activate the selected inner product everywhere)
  computeGlobalAdaptationIndicators( domain );
}

void SinglePhaseMixedMFD::computeGlobalAdaptationIndicators( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  MixedMimeticDiscretizationManager const & mmManager = numericalMethodManager.getMixedMimeticDiscretizationManager();
  MixedMimeticDiscretization const & discretization = mmManager.getMixedMimeticDiscretization( m_discretizationName );

  if( !discretization.isAdaptive() )
  {
    // no adaptation: activate the selected inner product in every cell
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                  MeshLevel & mesh,
                                                                  string_array const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                                          [&]( localIndex const,
                                                                               ElementSubRegionBase & subRegion )
      {
        subRegion.getField< mixedMimetic::stencilFlag >().template setValues< parallelDevicePolicy<> >( 1 );
      } );
    } );
    return;
  }

  real64 const lengthTolerance = domain.getMeshBody( 0 ).getGlobalLengthScale() * m_areaRelTol;
  real64 const tolerance = discretization.getResidualTolerance();
  R1Tensor const gradientInput = discretization.getNominalGradient();
  real64 const gradient[3] = { gradientInput[0], gradientInput[1], gradientInput[2] };

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    arrayView1d< real64 > const faceResidual = faceManager.getField< mixedMimetic::faceResidual >();
    faceResidual.zero();

    // step 1: projection of the admissible flow field induced by the nominal gradient
    array1d< real64 > projFaceFluxArray( faceManager.size() );
    arrayView1d< real64 > const projFaceFlux = projFaceFluxArray.toView();

    ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > const elemCenterAccessor =
      elemManager.constructViewAccessor< array2d< real64 >, arrayView2d< real64 const > >( ElementSubRegionBase::viewKeyStruct::elementCenterString() );

    using PermeabilityAccessors = StencilMaterialAccessors< PermeabilityBase, fields::permeability::permeability >;
    PermeabilityAccessors const permAccessors( elemManager, getName() );

    mixedMimeticKernels::FaceFluxProjectionKernel::
      launch< parallelDevicePolicy<> >( faceManager.size(),
                                        faceManager.elementRegionList(),
                                        faceManager.elementSubRegionList(),
                                        faceManager.elementList(),
                                        m_regionFilter.toViewConst(),
                                        faceManager.faceCenter(),
                                        faceManager.faceNormal(),
                                        faceManager.faceArea(),
                                        elemCenterAccessor.toNestedViewConst(),
                                        permAccessors.get( fields::permeability::permeability {} ),
                                        gradient,
                                        projFaceFlux );

    // steps 2-3: localized normalized residuals, assembled on the global face orientation
    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   CellElementSubRegion const & subRegion )
    {
      string const & permName = subRegion.getReference< string >( viewKeyStruct::permeabilityNamesString() );
      PermeabilityBase const & permeability = getConstitutiveModel< PermeabilityBase >( subRegion, permName );

      mixedMimeticKernels::internal::kernelLaunchSelectorFaceSwitch( subRegion.numFacesPerElement(), [&] ( auto NUM_FACES )
      {
        mixedMimeticKernels::LocalResidualKernel< NUM_FACES >::
          template launch< parallelDevicePolicy<> >( subRegion.size(),
                                                     nodeManager.referencePosition(),
                                                     faceManager.nodeList().toViewConst(),
                                                     subRegion.faceList().toViewConst(),
                                                     subRegion.getElementCenter(),
                                                     subRegion.getElementVolume(),
                                                     permeability.permeability(),
                                                     faceManager.faceCenter(),
                                                     faceManager.faceNormal(),
                                                     projFaceFlux.toViewConst(),
                                                     gradient,
                                                     lengthTolerance,
                                                     faceResidual );
      } );
    } );

    // step 4: thresholding
    localIndex numMfdCells = 0;
    localIndex numCells = 0;
    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   CellElementSubRegion & subRegion )
    {
      arrayView1d< real64 > const consistencyIndicator = subRegion.getField< mixedMimetic::consistencyIndicator >();
      arrayView1d< integer > const stencilFlag = subRegion.getField< mixedMimetic::stencilFlag >();

      mixedMimeticKernels::internal::kernelLaunchSelectorFaceSwitch( subRegion.numFacesPerElement(), [&] ( auto NUM_FACES )
      {
        numMfdCells += mixedMimeticKernels::MarkingKernel< NUM_FACES >::
                         template launch< parallelDevicePolicy<> >( subRegion.size(),
                                                                    subRegion.faceList().toViewConst(),
                                                                    subRegion.ghostRank(),
                                                                    faceResidual.toViewConst(),
                                                                    tolerance,
                                                                    consistencyIndicator,
                                                                    stencilFlag );
      } );
      numCells += subRegion.size() - subRegion.getNumberOfGhosts();
    } );

    // make the marking consistent on ghost cells
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( { mixedMimetic::stencilFlag::key(), mixedMimetic::consistencyIndicator::key() }, regionNames );
    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), false );

    globalIndex const globalNumMfdCells = MpiWrapper::sum< globalIndex >( numMfdCells );
    globalIndex const globalNumCells = MpiWrapper::sum< globalIndex >( numCells );
    GEOS_LOG_RANK_0( GEOS_FMT( "{}: Global Adaptation marked {} / {} cells as MFD-compatible (tolerance = {})",
                               getName(), globalNumMfdCells, globalNumCells, tolerance ) );
  } );
}

void SinglePhaseMixedMFD::implicitStepSetup( real64 const & time_n,
                                             real64 const & dt,
                                             DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  // setup the cell-centered fields
  SinglePhaseBase::implicitStepSetup( time_n, dt, domain );

  // setup the face fields
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();

    arrayView1d< real64 const > const & faceFlux =
      faceManager.getField< mixedMimetic::faceMassFlux >();
    arrayView1d< real64 > const & faceFlux_n =
      faceManager.getField< mixedMimetic::faceMassFlux_n >();
    faceFlux_n.setValues< parallelDevicePolicy<> >( faceFlux );
  } );

  // evaluate the boundary face pressure values used in the constitutive rows
  applyFacePressureBCValues( time_n + dt, domain );
}

void SinglePhaseMixedMFD::implicitStepComplete( real64 const & time,
                                                real64 const & dt,
                                                DomainPartition & domain )
{
  SinglePhaseBase::implicitStepComplete( time, dt, domain );
}

namespace
{
char const faceBcLogMessage[] =
  "SinglePhaseMixedMFD {}: at time {}s, "
  "the <{}> boundary condition '{}' is applied to the face set '{}' in '{}'. "
  "\nThe total number of target faces (including ghost faces) is {}. "
  "\nNote that if this number is equal to zero, the boundary condition will not be applied on this face set.";
}

void SinglePhaseMixedMFD::applyFacePressureBCValues( real64 const time,
                                                     DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel & mesh,
                                                                      string_array const & )
  {
    fsManager.apply< FaceManager >( time,
                                    mesh,
                                    flow::bcPressure::key(),
                                    [&] ( FieldSpecification const & fs,
                                          string const & setName,
                                          SortedArrayView< localIndex const > const & targetSet,
                                          FaceManager & targetGroup,
                                          string const & )
    {
      // provide some logging at the first nonlinear iteration
      if( m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
      {
        globalIndex const numTargetFaces = MpiWrapper::sum< globalIndex >( targetSet.size() );
        GEOS_LOG_LEVEL_RANK_0_ON_GROUP( logInfo::BoundaryConditions,
                                        GEOS_FMT( faceBcLogMessage,
                                                  this->getName(), time, fs.getCatalogName(), fs.getName(),
                                                  setName, targetGroup.getName(), numTargetFaces ),
                                        fs );
      }

      // populate the boundary face pressure values: they enter the residual directly during assembly
      FieldSpecificationImpl::applyFieldValue< FieldSpecificationEqual,
                                               parallelDevicePolicy<> >( fs,
                                                                         targetSet,
                                                                         time,
                                                                         targetGroup,
                                                                         flow::bcPressure::key() );
    } );
  } );
}

void SinglePhaseMixedMFD::setupSystem( DomainPartition & domain,
                                       DofManager & dofManager,
                                       CRSMatrix< real64, globalIndex > & localMatrix,
                                       ParallelVector & rhs,
                                       ParallelVector & solution,
                                       bool const setSparsity )
{
  SinglePhaseBase::setupSystem( domain, dofManager, localMatrix, rhs, solution, setSparsity );

  // with the dof numbering finalized, build the per-dof labels driving the
  // stencilFlag-guided three-level MGR reduction
  computeMgrPointMarkers( domain, dofManager );
}

void SinglePhaseMixedMFD::computeMgrPointMarkers( DomainPartition const & domain,
                                                  DofManager const & dofManager )
{
  GEOS_MARK_FUNCTION;

  string const faceDofKey = dofManager.getKey( mixedMimetic::faceMassFlux::key() );
  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  globalIndex const rankOffset = dofManager.rankOffset();

  // when the selected inner product is itself TPFA, the effective operator is diagonal in
  // every cell regardless of the stencil activation flag, and all face-flux rows belong to
  // the exactly-eliminable class
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  MixedMimeticDiscretizationManager const & mmManager = numericalMethodManager.getMixedMimeticDiscretizationManager();
  MixedMimeticDiscretization const & discretization = mmManager.getMixedMimeticDiscretization( m_discretizationName );
  MimeticInnerProductBase const & mimeticInnerProductBase =
    discretization.getReference< MimeticInnerProductBase >( MixedMimeticDiscretization::viewKeyStruct::innerProductString() );
  bool const effectiveTpfa = dynamic_cast< TPFAInnerProduct const * >( &mimeticInnerProductBase ) != nullptr;

  array1d< integer > & pointMarkers = m_linearSolverParameters.get().mgr.customPointMarkers;
  pointMarkers.resize( dofManager.numLocalDofs() );
  arrayView1d< integer > const markers = pointMarkers.toView();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               string_array const & regionNames )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > const stencilFlagAccessor =
      elemManager.constructArrayViewAccessor< integer, 1 >( mixedMimetic::stencilFlag::key() );
    ElementRegionManager::ElementViewConst< arrayView1d< integer const > > const stencilFlag =
      stencilFlagAccessor.toNestedViewConst();

    arrayView1d< globalIndex const > const faceDofNumber =
      faceManager.getReference< array1d< globalIndex > >( faceDofKey );
    arrayView1d< integer const > const faceGhostRank = faceManager.ghostRank();
    arrayView2d< localIndex const > const elemRegionList = faceManager.elementRegionList();
    arrayView2d< localIndex const > const elemSubRegionList = faceManager.elementSubRegionList();
    arrayView2d< localIndex const > const elemList = faceManager.elementList();
    SortedArrayView< localIndex const > const regionFilter = m_regionFilter.toViewConst();

    // face-flux dofs: label 0 if every adjacent target cell is TPFA-compatible (the
    // assembled flux row is exactly diagonal), label 1 otherwise
    forAll< parallelHostPolicy >( faceManager.size(), [=]( localIndex const kf )
    {
      if( faceGhostRank[kf] >= 0 || faceDofNumber[kf] < 0 )
      {
        return;
      }
      integer label = 0;
      if( !effectiveTpfa )
      {
        for( integer k = 0; k < elemRegionList.size( 1 ); ++k )
        {
          localIndex const er  = elemRegionList[kf][k];
          localIndex const esr = elemSubRegionList[kf][k];
          localIndex const ei  = elemList[kf][k];
          if( er >= 0 && esr >= 0 && ei >= 0 && regionFilter.contains( er ) )
          {
            label = LvArray::math::max( label, stencilFlag[er][esr][ei] );
          }
        }
      }
      markers[faceDofNumber[kf] - rankOffset] = label;
    } );

    // cell-pressure dofs: label 2
    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase const & subRegion )
    {
      arrayView1d< globalIndex const > const elemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( elemDofKey );
      arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();

      forAll< parallelHostPolicy >( subRegion.size(), [=]( localIndex const ei )
      {
        if( elemGhostRank[ei] < 0 )
        {
          markers[elemDofNumber[ei] - rankOffset] = 2;
        }
      } );
    } );
  } );
}

void SinglePhaseMixedMFD::setupDofs( DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                     DofManager & dofManager ) const
{
  // face mass-flux unknowns: the dense NF x NF local mass matrices couple
  // the fluxes of the faces of each cell
  dofManager.addField( mixedMimetic::faceMassFlux::key(),
                       FieldLocation::Face,
                       1,
                       getMeshTargets() );

  dofManager.addCoupling( mixedMimetic::faceMassFlux::key(),
                          mixedMimetic::faceMassFlux::key(),
                          DofManager::Connector::Elem );

  // cell pressure unknowns: the mass conservation equation couples the cell
  // to its own faces only (no cell-to-cell coupling in the mixed form)
  dofManager.addField( viewKeyStruct::elemDofFieldString(),
                       FieldLocation::Elem,
                       1,
                       getMeshTargets() );

  dofManager.addCoupling( viewKeyStruct::elemDofFieldString(),
                          viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::None );

  // coupling between the face fluxes and the cell pressures
  dofManager.addCoupling( mixedMimetic::faceMassFlux::key(),
                          viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Elem );
}

void SinglePhaseMixedMFD::assembleFluxTerms( real64 const dt,
                                             DomainPartition const & domain,
                                             DofManager const & dofManager,
                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                             arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  MixedMimeticDiscretizationManager const & mmManager = numericalMethodManager.getMixedMimeticDiscretizationManager();
  MixedMimeticDiscretization const & discretization = mmManager.getMixedMimeticDiscretization( m_discretizationName );
  MimeticInnerProductBase const & mimeticInnerProductBase =
    discretization.getReference< MimeticInnerProductBase >( MixedMimeticDiscretization::viewKeyStruct::innerProductString() );

  string const faceDofKey = dofManager.getKey( mixedMimetic::faceMassFlux::key() );
  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  // tolerance for the mass matrix computations
  real64 const lengthTolerance = domain.getMeshBody( 0 ).getGlobalLengthScale() * m_areaRelTol;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               string_array const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames,
                                                                        [&]( localIndex const,
                                                                             CellElementSubRegion const & subRegion )
    {
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      SingleFluidBase const & fluid = getConstitutiveModel< SingleFluidBase >( subRegion, fluidName );

      string const & permName = subRegion.getReference< string >( viewKeyStruct::permeabilityNamesString() );
      PermeabilityBase const & permeability = getConstitutiveModel< PermeabilityBase >( subRegion, permName );

      singlePhaseMixedMFDKernels::
        ElementBasedAssemblyKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                   lengthTolerance,
                                                   elemDofKey,
                                                   faceDofKey,
                                                   nodeManager,
                                                   faceManager,
                                                   mesh.getElemManager(),
                                                   subRegion,
                                                   mimeticInnerProductBase,
                                                   fluid,
                                                   permeability,
                                                   dt,
                                                   localMatrix,
                                                   localRhs );
    } );
  } );
}

void SinglePhaseMixedMFD::assembleStabilizedFluxTerms( real64 const dt,
                                                       DomainPartition const & domain,
                                                       DofManager const & dofManager,
                                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                       arrayView1d< real64 > const & localRhs )
{
  // pressure stabilization not implemented
  GEOS_UNUSED_VAR( dt, domain, dofManager, localMatrix, localRhs );
  GEOS_ERROR( "Stabilized flux not available for this flow solver" );
}

void SinglePhaseMixedMFD::assembleEDFMFluxTerms( real64 const GEOS_UNUSED_PARAM( time_n ),
                                                 real64 const dt,
                                                 DomainPartition const & domain,
                                                 DofManager const & dofManager,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs,
                                                 string const & jumpDofKey )
{
  GEOS_UNUSED_VAR( jumpDofKey );

  assembleFluxTerms( dt,
                     domain,
                     dofManager,
                     localMatrix,
                     localRhs );
}

void SinglePhaseMixedMFD::applyBoundaryConditions( real64 const time_n,
                                                   real64 const dt,
                                                   DomainPartition & domain,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  // the face pressure boundary values are applied inside the constitutive rows during assembly;
  // only the cell-centered boundary conditions (Dirichlet cells, source fluxes) remain to be applied here
  SinglePhaseBase::applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );
}

void SinglePhaseMixedMFD::applyAquiferBC( real64 const time,
                                          real64 const dt,
                                          DomainPartition & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( time, dt, dofManager, domain, localMatrix, localRhs );
}

void SinglePhaseMixedMFD::saveAquiferConvergedState( real64 const & time,
                                                     real64 const & dt,
                                                     DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( time, dt, domain );
}

real64 SinglePhaseMixedMFD::calculateResidualNorm( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                                   real64 const & GEOS_UNUSED_PARAM( dt ),
                                                   DomainPartition const & domain,
                                                   DofManager const & dofManager,
                                                   arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  real64 localResidualNorm = 0.0;
  real64 localResidualNormalizer = 0.0;

  physicsSolverBaseKernels::NormType const normType = getNonlinearSolverParameters().normType();

  globalIndex const rankOffset = dofManager.rankOffset();
  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  string const faceDofKey = dofManager.getKey( mixedMimetic::faceMassFlux::key() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                string_array const & regionNames )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    // step 1: compute the residual for the element-based mass conservation equations

    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase const & subRegion )
    {
      real64 subRegionResidualNorm[1]{};
      real64 subRegionResidualNormalizer[1]{};

      singlePhaseBaseKernels::
        ResidualNormKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( normType,
                                                   rankOffset,
                                                   elemDofKey,
                                                   localRhs,
                                                   subRegion,
                                                   m_nonlinearSolverParameters.m_minNormalizer,
                                                   subRegionResidualNorm,
                                                   subRegionResidualNormalizer );

      if( normType == physicsSolverBaseKernels::NormType::Linf )
      {
        if( subRegionResidualNorm[0] > localResidualNorm )
        {
          localResidualNorm = subRegionResidualNorm[0];
        }
      }
      else
      {
        localResidualNorm += subRegionResidualNorm[0];
        localResidualNormalizer += subRegionResidualNormalizer[0];
      }
    } );

    // step 2: compute the residual for the face-based constitutive rows

    real64 faceResidualNorm[1]{};
    real64 faceResidualNormalizer[1]{};

    singlePhaseMixedMFDKernels::
      ResidualNormKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( normType,
                                                 rankOffset,
                                                 faceDofKey,
                                                 localRhs,
                                                 m_regionFilter.toViewConst(),
                                                 getName(),
                                                 elemManager,
                                                 faceManager,
                                                 m_nonlinearSolverParameters.m_minNormalizer,
                                                 faceResidualNorm,
                                                 faceResidualNormalizer );

    if( normType == physicsSolverBaseKernels::NormType::Linf )
    {
      if( faceResidualNorm[0] > localResidualNorm )
      {
        localResidualNorm = faceResidualNorm[0];
      }
    }
    else
    {
      localResidualNorm += faceResidualNorm[0];
      localResidualNormalizer += faceResidualNormalizer[0];
    }
  } );

  // step 3: second reduction across MPI ranks

  real64 residualNorm = 0.0;
  if( normType == physicsSolverBaseKernels::NormType::Linf )
  {
    physicsSolverBaseKernels::LinfResidualNormHelper::computeGlobalNorm( localResidualNorm, residualNorm );
  }
  else
  {
    physicsSolverBaseKernels::L2ResidualNormHelper::computeGlobalNorm( localResidualNorm, localResidualNormalizer, residualNorm );
  }

  GEOS_LOG_LEVEL_RANK_0_NLR( logInfo::ResidualNorm,
                             GEOS_FMT( "        ( R{} ) = ( {:4.2e} )", coupledSolverAttributePrefix(), residualNorm ));
  getConvergenceStats().setResidualValue( GEOS_FMT( "R{}", coupledSolverAttributePrefix()), residualNorm );

  return residualNorm;
}

void SinglePhaseMixedMFD::applySystemSolution( DofManager const & dofManager,
                                               arrayView1d< real64 const > const & localSolution,
                                               real64 const scalingFactor,
                                               real64 const dt,
                                               DomainPartition & domain )
{
  GEOS_UNUSED_VAR( dt );

  // 1. apply the cell-centered update

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::elemDofFieldString(),
                               flow::pressure::key(),
                               scalingFactor );

  // 2. apply the face-based update

  dofManager.addVectorToField( localSolution,
                               mixedMimetic::faceMassFlux::key(),
                               mixedMimetic::faceMassFlux::key(),
                               scalingFactor );

  // 3. synchronize
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    FieldIdentifiers fieldsToBeSync;

    fieldsToBeSync.addElementFields( { flow::pressure::key() }, regionNames );
    fieldsToBeSync.addFields( FieldLocation::Face, { mixedMimetic::faceMassFlux::key() } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), true );
  } );
}

void SinglePhaseMixedMFD::resetStateToBeginningOfStep( DomainPartition & domain )
{
  // Reset the cell-centered fields
  SinglePhaseBase::resetStateToBeginningOfStep( domain );

  // Reset the face-based fields
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();

    arrayView1d< real64 > const & faceFlux =
      faceManager.getField< mixedMimetic::faceMassFlux >();
    arrayView1d< real64 const > const & faceFlux_n =
      faceManager.getField< mixedMimetic::faceMassFlux_n >();
    faceFlux.setValues< parallelDevicePolicy<> >( faceFlux_n );
  } );
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SinglePhaseMixedMFD, string const &, Group * const )
} /* namespace geos */
