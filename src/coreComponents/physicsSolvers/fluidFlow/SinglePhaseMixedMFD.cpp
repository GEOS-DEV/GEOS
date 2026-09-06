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

#include <algorithm>
#include <unordered_map>
#include <array>

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
      subRegion.registerField< mixedMimetic::mfdFlag >( getName() );
      subRegion.registerField< mixedMimetic::prescribedMfdFlag >( getName() );
      subRegion.registerField< mixedMimetic::consistencyIndicator >( getName() );
      subRegion.registerField< mixedMimetic::degeneracyIndicator >( getName() );
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

      // face classification driving the TPFA-face condensation and the MGR labels
      faceManager.registerField< mixedMimetic::faceStencilLabel >( getName() );
      faceManager.registerField< mixedMimetic::faceOrientationCell >( getName() );
      faceManager.registerField< mixedMimetic::faceDofScale >( getName() );
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

  computeFaceOrientation( domain );

  // run the residual-based Global Adaptation pipeline (or activate the selected inner product everywhere)
  computeGlobalAdaptationIndicators( domain );
}

void SinglePhaseMixedMFD::computeFaceOrientation( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();
    arrayView2d< localIndex const > const elemRegionList = faceManager.elementRegionList();
    arrayView2d< localIndex const > const elemSubRegionList = faceManager.elementSubRegionList();
    arrayView2d< localIndex const > const elemList = faceManager.elementList();
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const elemLocalToGlobal =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( ObjectManagerBase::viewKeyStruct::localToGlobalMapString() );
    ElementRegionManager::ElementViewConst< arrayView1d< globalIndex const > > const l2g = elemLocalToGlobal.toNestedViewConst();
    arrayView1d< globalIndex > const orientationCell = faceManager.getField< mixedMimetic::faceOrientationCell >();

    // the owner of a face sees both of its cells; ghost faces take the owner's value through the sync
    forAll< parallelHostPolicy >( faceManager.size(), [=]( localIndex const kf )
    {
      globalIndex gMin = -1;
      for( localIndex k = 0; k < 2; ++k )
      {
        localIndex const er = elemRegionList[kf][k];
        localIndex const esr = elemSubRegionList[kf][k];
        localIndex const ei = elemList[kf][k];
        if( er >= 0 && esr >= 0 && ei >= 0 )
        {
          globalIndex const g = l2g[er][esr][ei];
          gMin = ( gMin < 0 || g < gMin ) ? g : gMin;
        }
      }
      orientationCell[kf] = gMin;
    } );

    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addFields( FieldLocation::Face, { mixedMimetic::faceOrientationCell::key() } );
    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), false );
  } );
}

void SinglePhaseMixedMFD::computeGlobalAdaptationIndicators( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  MixedMimeticDiscretizationManager const & mmManager = numericalMethodManager.getMixedMimeticDiscretizationManager();
  MixedMimeticDiscretization const & discretization = mmManager.getMixedMimeticDiscretization( m_discretizationName );

  // every cell starts with the consistent (MFD) product; the layers below can only refine that choice
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                                        [&]( localIndex const,
                                                                             ElementSubRegionBase & subRegion )
    {
      subRegion.getField< mixedMimetic::mfdFlag >().template setValues< parallelDevicePolicy<> >( 1 );
    } );
  } );

  if( !discretization.isAdaptiveConsistency() )
  {
    applyLayersAndLabel( domain );
    return;
  }

  real64 const lengthTolerance = domain.getMeshBody( 0 ).getGlobalLengthScale() * m_areaRelTol;
  real64 const tolerance = discretization.getConsistencyTolerance();
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
                                        lengthTolerance,
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
      arrayView1d< integer > const mfdFlag = subRegion.getField< mixedMimetic::mfdFlag >();

      mixedMimeticKernels::internal::kernelLaunchSelectorFaceSwitch( subRegion.numFacesPerElement(), [&] ( auto NUM_FACES )
      {
        numMfdCells += mixedMimeticKernels::MarkingKernel< NUM_FACES >::
                       template launch< parallelDevicePolicy<> >( subRegion.size(),
                                                                  subRegion.faceList().toViewConst(),
                                                                  subRegion.ghostRank(),
                                                                  faceResidual.toViewConst(),
                                                                  tolerance,
                                                                  consistencyIndicator,
                                                                  mfdFlag );
      } );
      numCells += subRegion.size() - subRegion.getNumberOfGhosts();
    } );

    // make the marking consistent on ghost cells
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( { mixedMimetic::mfdFlag::key(), mixedMimetic::consistencyIndicator::key() }, regionNames );
    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), false );

    globalIndex const globalNumMfdCells = MpiWrapper::sum< globalIndex >( numMfdCells );
    globalIndex const globalNumCells = MpiWrapper::sum< globalIndex >( numCells );
    GEOS_LOG_RANK_0( GEOS_FMT( "mixedMFD Flow: consistency layer (tolerance = {}) set eta = 1 on {} / {} cells",
                               tolerance, globalNumMfdCells, globalNumCells ) );
  } );

  applyLayersAndLabel( domain );
}

void SinglePhaseMixedMFD::applyLayersAndLabel( DomainPartition & domain )
{
  // prescription, then degeneracy (admissibility, which alone may override the prescription)
  std::pair< localIndex, localIndex > const prescribed = applyPrescribedFlag( domain );
  std::pair< localIndex, localIndex > const degenerate = applyDegeneracyLayer( domain );
  globalIndex const numPrescribed0 = MpiWrapper::sum< globalIndex >( prescribed.first );
  globalIndex const numPrescribed1 = MpiWrapper::sum< globalIndex >( prescribed.second );
  globalIndex const numRejected = MpiWrapper::sum< globalIndex >( degenerate.second );
  if( numPrescribed0 + numPrescribed1 > 0 )
  {
    GEOS_LOG_RANK_0( GEOS_FMT( "mixedMFD Flow: prescribed eta = 0 on {} cells and eta = 1 on {} cells ({} rejected by the degeneracy layer)",
                               numPrescribed0, numPrescribed1, numRejected ) );
  }
  computeFaceStencilLabels( domain );
}

std::pair< localIndex, localIndex > SinglePhaseMixedMFD::applyPrescribedFlag( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  localIndex numPrescribed0 = 0;
  localIndex numPrescribed1 = 0;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    FieldIdentifiers prescribedToBeSync;
    prescribedToBeSync.addElementFields( { mixedMimetic::prescribedMfdFlag::key() }, regionNames );
    CommunicationTools::getInstance().synchronizeFields( prescribedToBeSync, mesh, domain.getNeighbors(), false );

    mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                                        [&]( localIndex const,
                                                                             ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 const > const prescribed = subRegion.getField< mixedMimetic::prescribedMfdFlag >();
      arrayView1d< integer > const mfdFlag = subRegion.getField< mixedMimetic::mfdFlag >();
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        if( prescribed[ei] < 0.0 )
        {
          continue;
        }
        integer const eta = prescribed[ei] > 0.5 ? 1 : 0;
        mfdFlag[ei] = eta;
        if( ghostRank[ei] < 0 )
        {
          numPrescribed0 += ( eta == 0 );
          numPrescribed1 += ( eta == 1 );
        }
      }
    } );
  } );
  return { numPrescribed0, numPrescribed1 };
}

std::pair< localIndex, localIndex > SinglePhaseMixedMFD::applyDegeneracyLayer( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  MixedMimeticDiscretizationManager const & mmManager = numericalMethodManager.getMixedMimeticDiscretizationManager();
  MixedMimeticDiscretization const & discretization = mmManager.getMixedMimeticDiscretization( m_discretizationName );
  real64 const tolerance = discretization.getDegeneracyTolerance();

  localIndex numDegenerate = 0;
  localIndex numRejected = 0;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    ElementRegionManager & elemManager = mesh.getElemManager();
    ArrayOfArraysView< localIndex const > const nodeToRegion = nodeManager.elementRegionList();
    ArrayOfArraysView< localIndex const > const nodeToSubRegion = nodeManager.elementSubRegionList();
    ArrayOfArraysView< localIndex const > const nodeToElem = nodeManager.elementList();
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const elemVolume =
      elemManager.constructArrayViewAccessor< real64, 1 >( ElementSubRegionBase::viewKeyStruct::elementVolumeString() );

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   CellElementSubRegion & subRegion )
    {
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemToNodes = subRegion.nodeList().toViewConst();
      arrayView1d< real64 const > const volume = subRegion.getElementVolume();
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
      arrayView1d< real64 > const indicator = subRegion.getField< mixedMimetic::degeneracyIndicator >();
      arrayView1d< integer > const mfdFlag = subRegion.getField< mixedMimetic::mfdFlag >();
      arrayView1d< real64 const > const prescribed = subRegion.getField< mixedMimetic::prescribedMfdFlag >();
      localIndex const numNodes = subRegion.numNodesPerElement();

      // the node star of a cell: every cell sharing a vertex with it, counted once
      stdVector< std::array< localIndex, 3 > > star;
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        star.clear();
        for( localIndex a = 0; a < numNodes; ++a )
        {
          localIndex const n = elemToNodes( ei, a );
          for( localIndex k = 0; k < nodeToElem.sizeOfArray( n ); ++k )
          {
            if( nodeToRegion( n, k ) >= 0 && nodeToSubRegion( n, k ) >= 0 && nodeToElem( n, k ) >= 0 )
            {
              star.push_back( { nodeToRegion( n, k ), nodeToSubRegion( n, k ), nodeToElem( n, k ) } );
            }
          }
        }
        std::sort( star.begin(), star.end() );
        star.erase( std::unique( star.begin(), star.end() ), star.end() );
        real64 sum = 0.0;
        for( auto const & c : star )
        {
          sum += elemVolume[c[0]][c[1]][c[2]];
        }
        // share of the cell in the volume of its node star, in percent
        real64 const percent = sum > 0.0 ? 100.0 * volume[ei] / sum : 0.0;
        indicator[ei] = percent;
        if( percent < tolerance && mfdFlag[ei] == 1 )
        {
          mfdFlag[ei] = 0;
          numDegenerate += ( ghostRank[ei] < 0 ) ? 1 : 0;
          numRejected += ( ghostRank[ei] < 0 && prescribed[ei] > 0.5 ) ? 1 : 0;
        }
      }
    } );

    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( { mixedMimetic::mfdFlag::key(), mixedMimetic::degeneracyIndicator::key() }, regionNames );
    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), false );
  } );
  return { numDegenerate, numRejected };
}

void SinglePhaseMixedMFD::computeFaceStencilLabels( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  MixedMimeticDiscretizationManager const & mmManager = numericalMethodManager.getMixedMimeticDiscretizationManager();
  MixedMimeticDiscretization const & discretization = mmManager.getMixedMimeticDiscretization( m_discretizationName );

  // with a TPFA inner product the effective operator is diagonal in every cell,
  // regardless of the stencil activation flags
  bool const effectiveTpfa = discretization.isTpfaInnerProduct();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > const mfdFlagAccessor =
      elemManager.constructArrayViewAccessor< integer, 1 >( mixedMimetic::mfdFlag::key() );

    mixedMimeticKernels::FaceLabelKernel::
      launch< parallelDevicePolicy<> >( faceManager.size(),
                                        faceManager.elementRegionList(),
                                        faceManager.elementSubRegionList(),
                                        faceManager.elementList(),
                                        m_regionFilter.toViewConst(),
                                        mfdFlagAccessor.toNestedViewConst(),
                                        effectiveTpfa,
                                        faceManager.getField< mixedMimetic::faceStencilLabel >() );
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
                                        GEOS_FMT_RUNTIME( faceBcLogMessage,
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

  // with the dof numbering finalized, build the per-dof labels of the MGR reduction
  computeMgrPointMarkers( domain, dofManager );
}

void SinglePhaseMixedMFD::computeMgrPointMarkers( DomainPartition const & domain,
                                                  DofManager const & dofManager )
{
  GEOS_MARK_FUNCTION;

  string const faceDofKey = dofManager.getKey( mixedMimetic::faceMassFlux::key() );
  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  globalIndex const rankOffset = dofManager.rankOffset();

  array1d< integer > & pointMarkers = m_linearSolverParameters.get().mgr.customPointMarkers;
  pointMarkers.resize( dofManager.numLocalDofs() );
  arrayView1d< integer > const markers = pointMarkers.toView();

  // flux dofs kept in the saddle point (label 1) and condensed into the pressure system (label 0)
  localIndex numLiveFaces = 0;
  localIndex numCondensedFaces = 0;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               string_array const & )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    arrayView1d< globalIndex const > const faceDofNumber =
      faceManager.getReference< array1d< globalIndex > >( faceDofKey );
    arrayView1d< integer const > const faceGhostRank = faceManager.ghostRank();
    arrayView1d< integer const > const faceStencilLabel = faceManager.getField< mixedMimetic::faceStencilLabel >();

    RAJA::ReduceSum< parallelHostReduce, localIndex > numLive( 0 );
    RAJA::ReduceSum< parallelHostReduce, localIndex > numCondensed( 0 );
    forAll< parallelHostPolicy >( faceManager.size(), [=]( localIndex const kf )
    {
      bool const owned = faceGhostRank[kf] < 0 && faceDofNumber[kf] >= 0;
      numLive += ( owned && faceStencilLabel[kf] == 1 ) ? 1 : 0;
      numCondensed += ( owned && faceStencilLabel[kf] == 0 ) ? 1 : 0;
    } );
    numLiveFaces += numLive.get();
    numCondensedFaces += numCondensed.get();
  } );
  globalIndex const globalNumLiveFaces = MpiWrapper::sum< globalIndex >( numLiveFaces );
  globalIndex const globalNumCondensedFaces = MpiWrapper::sum< globalIndex >( numCondensedFaces );
  GEOS_LOG_RANK_0( GEOS_FMT( "mixedMFD Flow: flux dofs {} non-condensed (saddle point), {} condensed (two-point closure)",
                             globalNumLiveFaces, globalNumCondensedFaces ) );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               string_array const & regionNames )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    arrayView1d< globalIndex const > const faceDofNumber =
      faceManager.getReference< array1d< globalIndex > >( faceDofKey );
    arrayView1d< integer const > const faceGhostRank = faceManager.ghostRank();
    arrayView1d< integer const > const faceStencilLabel = faceManager.getField< mixedMimetic::faceStencilLabel >();

    // face-flux dofs: the face classification is the MGR label (0 = exactly-diagonal row)
    forAll< parallelHostPolicy >( faceManager.size(), [=]( localIndex const kf )
    {
      if( faceGhostRank[kf] >= 0 || faceDofNumber[kf] < 0 )
      {
        return;
      }
      markers[faceDofNumber[kf] - rankOffset] = faceStencilLabel[kf];
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

  // cell pressure unknowns: the TPFA-face condensation writes two-point stencil
  // entries directly into the mass conservation rows, so cell-to-cell coupling
  // through the faces is required
  dofManager.addField( viewKeyStruct::elemDofFieldString(),
                       FieldLocation::Elem,
                       1,
                       getMeshTargets() );

  dofManager.addCoupling( viewKeyStruct::elemDofFieldString(),
                          viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Face );

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
                                                   subRegion,
                                                   mimeticInnerProductBase,
                                                   fluid,
                                                   permeability,
                                                   dt,
                                                   localMatrix,
                                                   localRhs );
    } );

    // condensed (label-0) faces: two-point flux contributions to the mass conservation
    // rows and one-way closure rows, assembled in a single face-based sweep
    singlePhaseMixedMFDKernels::
      TpfaCondensedFluxKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                 lengthTolerance,
                                                 elemDofKey,
                                                 faceDofKey,
                                                 getName(),
                                                 nodeManager,
                                                 faceManager,
                                                 mesh.getElemManager(),
                                                 m_regionFilter.toViewConst(),
                                                 dt,
                                                 localMatrix,
                                                 localRhs );
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

  // residual-norm weights w_i = ||A_i S||_inf = max_j |A_ij| s_j, S = diag(s_j) the
  // characteristic scales of the unknowns: each term |A_ij| s_j has the units of r_i,
  // and w_i is invariant under row rescaling and under the TPFA-face condensation
  localIndex const numRows = localMatrix.numRows();
  m_residualWeight.resize( numRows );
  m_dofScale.resize( numRows );
  arrayView1d< real64 > const residualWeight = m_residualWeight.toView();
  arrayView1d< real64 > const dofScale = m_dofScale.toView();
  residualWeight.zero();
  dofScale.zero();

  globalIndex const rankOffset = dofManager.rankOffset();
  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  string const faceDofKey = dofManager.getKey( mixedMimetic::faceMassFlux::key() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    FaceManager & faceManager = mesh.getFaceManager();

    // s_p = |p_n| for the cell-pressure dofs
    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase const & subRegion )
    {
      arrayView1d< globalIndex const > const elemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( elemDofKey );
      arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const presN = subRegion.getField< fields::flow::pressure_n >();

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        if( elemGhostRank[ei] < 0 )
        {
          dofScale[elemDofNumber[ei] - rankOffset] = LvArray::math::abs( presN[ei] );
        }
      } );
    } );

    // s_m = p_scale / |M_ff| for the face-flux dofs, p_scale the mean adjacent |p_n|
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const presNAccessor =
      elemManager.constructArrayViewAccessor< real64, 1 >( fields::flow::pressure_n::key() );
    ElementRegionManager::ElementViewConst< arrayView1d< real64 const > > const presN = presNAccessor.toNestedViewConst();

    arrayView1d< globalIndex const > const faceDofNumber =
      faceManager.getReference< array1d< globalIndex > >( faceDofKey );
    arrayView1d< integer const > const faceGhostRank = faceManager.ghostRank();
    arrayView2d< localIndex const > const elemRegionList = faceManager.elementRegionList();
    arrayView2d< localIndex const > const elemSubRegionList = faceManager.elementSubRegionList();
    arrayView2d< localIndex const > const elemList = faceManager.elementList();
    SortedArrayView< localIndex const > const regionFilter = m_regionFilter.toViewConst();
    arrayView1d< real64 > const faceScale = faceManager.getField< mixedMimetic::faceDofScale >();

    forAll< parallelDevicePolicy<> >( faceManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const kf )
    {
      if( faceGhostRank[kf] >= 0 || faceDofNumber[kf] < 0 )
      {
        return;
      }
      real64 presSum = 0.0;
      integer elemCounter = 0;
      for( integer k = 0; k < elemRegionList.size( 1 ); ++k )
      {
        localIndex const er  = elemRegionList[kf][k];
        localIndex const esr = elemSubRegionList[kf][k];
        localIndex const ei  = elemList[kf][k];
        if( er >= 0 && esr >= 0 && ei >= 0 && regionFilter.contains( er ) )
        {
          presSum += LvArray::math::abs( presN[er][esr][ei] );
          elemCounter++;
        }
      }
      real64 const pScale = elemCounter > 0 ? presSum / elemCounter : 0.0;
      localIndex const localRow = faceDofNumber[kf] - rankOffset;
      real64 diag = 0.0;
      arraySlice1d< globalIndex const > const columns = localMatrix.getColumns( localRow );
      arraySlice1d< real64 const > const entries = localMatrix.getEntries( localRow );
      for( localIndex k = 0; k < localMatrix.numNonZeros( localRow ); ++k )
      {
        if( columns[k] == faceDofNumber[kf] )
        {
          diag = LvArray::math::abs( entries[k] );
        }
      }
      dofScale[localRow] = diag > 0.0 ? pScale / diag : pScale;
      faceScale[kf] = dofScale[localRow];
    } );
  } );

  // the scales of the ghost dofs: a cell may own none of its faces, so the face scales of its
  // owners are synchronized; the pressure scale of a ghost cell is its synchronized |p_n|
  std::unordered_map< globalIndex, real64 > ghostScale;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addFields( FieldLocation::Face, { mixedMimetic::faceDofScale::key() } );
    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), false );

    FaceManager const & faceManager = mesh.getFaceManager();
    arrayView1d< globalIndex const > const faceDofNumber = faceManager.getReference< array1d< globalIndex > >( faceDofKey );
    arrayView1d< integer const > const faceGhostRank = faceManager.ghostRank();
    arrayView1d< real64 const > const faceScale = faceManager.getField< mixedMimetic::faceDofScale >();
    for( localIndex kf = 0; kf < faceManager.size(); ++kf )
    {
      if( faceGhostRank[kf] >= 0 && faceDofNumber[kf] >= 0 )
      {
        ghostScale[faceDofNumber[kf]] = faceScale[kf];
      }
    }
    mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                                        [&]( localIndex const,
                                                                             ElementSubRegionBase const & subRegion )
    {
      arrayView1d< globalIndex const > const elemDofNumber = subRegion.getReference< array1d< globalIndex > >( elemDofKey );
      arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const presN = subRegion.getField< fields::flow::pressure_n >();
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        if( elemGhostRank[ei] >= 0 && elemDofNumber[ei] >= 0 )
        {
          ghostScale[elemDofNumber[ei]] = LvArray::math::abs( presN[ei] );
        }
      }
    } );
  } );

  // w_i = max_j |A_ij| s_j over every column of row i
  forAll< serialPolicy >( numRows, [&]( localIndex const i )
  {
    real64 w = 0.0;
    arraySlice1d< globalIndex const > const columns = localMatrix.getColumns( i );
    arraySlice1d< real64 const > const entries = localMatrix.getEntries( i );
    for( localIndex k = 0; k < localMatrix.numNonZeros( i ); ++k )
    {
      globalIndex const localCol = columns[k] - rankOffset;
      real64 scale = 0.0;
      if( localCol >= 0 && localCol < numRows )
      {
        scale = dofScale[localCol];
      }
      else
      {
        auto const it = ghostScale.find( columns[k] );
        scale = it != ghostScale.end() ? it->second : 0.0;
      }
      w = LvArray::math::max( w, LvArray::math::abs( entries[k] ) * scale );
    }
    residualWeight[i] = w;
  } );
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
                                                   DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                                   DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                                   arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  // row-equilibrated residual norm ||D^{-1} r||, D = diag(max(eps, w_i)), w_i = ||A_i S||_inf
  physicsSolverBaseKernels::NormType const normType = getNonlinearSolverParameters().normType();
  real64 const minNormalizer = m_nonlinearSolverParameters.m_minNormalizer;

  arrayView1d< real64 const > const residualWeight = m_residualWeight.toViewConst();
  localIndex const numRows = localRhs.size();

  real64 residualNorm = 0.0;
  if( normType == physicsSolverBaseKernels::NormType::Linf )
  {
    RAJA::ReduceMax< parallelDeviceReduce, real64 > maxVal( 0.0 );
    forAll< parallelDevicePolicy<> >( numRows, [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      maxVal.max( LvArray::math::abs( localRhs[i] ) /
                  LvArray::math::max( minNormalizer, residualWeight[i] ) );
    } );
    real64 const localNorm = maxVal.get();
    physicsSolverBaseKernels::LinfResidualNormHelper::computeGlobalNorm( localNorm, residualNorm );
  }
  else
  {
    RAJA::ReduceSum< parallelDeviceReduce, real64 > sumVal( 0.0 );
    forAll< parallelDevicePolicy<> >( numRows, [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      real64 const r = localRhs[i] / LvArray::math::max( minNormalizer, residualWeight[i] );
      sumVal += r * r;
    } );
    real64 const localSum = sumVal.get();
    real64 const localCount = static_cast< real64 >( numRows );
    physicsSolverBaseKernels::L2ResidualNormHelper::computeGlobalNorm( localSum, localCount, residualNorm );
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
