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
 * @file SinglePhasePoromechanicsConformingFracturesALM.cpp
 */

#include "SinglePhasePoromechanicsConformingFracturesALM.hpp"

#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/solid/PorousSolid.hpp"
#include "constitutive/contact/HydraulicApertureBase.hpp"
#include "constitutive/contact/HydraulicApertureRelationSelector.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "linearAlgebra/utilities/SparsityPatternUtilities.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanicsFractures.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanicsConformingFracturesALMKernels.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"
#include "physicsSolvers/solidMechanics/contact/ContactFields.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;

template< typename FLOW_SOLVER >
SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::SinglePhasePoromechanicsConformingFracturesALM( const string & name,
                                                                                                               Group * const parent )
  : Base( name, parent )
{}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::setupCoupling( DomainPartition const & domain,
                                                                                   DofManager & dofManager ) const
{
  GEOS_MARK_FUNCTION;

  // 1. Poromechanical coupling in the bulk (from base class)
  Base::setupCoupling( domain, dofManager );

  // 2. Pressure - bubble displacement coupling in the fracture
  dofManager.addCoupling( m_pressureKey,
                          contact::totalBubbleDisplacement::key(),
                          DofManager::Connector::Elem );
}


template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::setSparsityPattern( DomainPartition & domain,
                                                                                        DofManager & dofManager,
                                                                                        CRSMatrix< real64, globalIndex > & localMatrix,
                                                                                        SparsityPattern< globalIndex > & pattern )
{
  GEOS_MARK_FUNCTION;

  // Initialize ALM contact solver internal data structures
  // These must be called before assembling the contact-dependent pattern.
  this->solidMechanicsSolver()->createFaceTypeList( domain );
  this->solidMechanicsSolver()->updateStickSlipList( domain );
  this->solidMechanicsSolver()->createBubbleCellList( domain );

  // Start from both subsolver patterns. The flow pattern may contain well and
  // flux couplings, while the mechanics pattern contains the nodal-bubble couplings.
  SparsityPattern< globalIndex > flowPattern;
  this->flowSolver()->setSparsityPattern( domain, dofManager, localMatrix, flowPattern );

  SparsityPattern< globalIndex > mechanicsPattern;
  this->solidMechanicsSolver()->setSparsityPattern( domain, dofManager, localMatrix, mechanicsPattern );
  GEOS_ERROR_IF_NE( flowPattern.numRows(), mechanicsPattern.numRows() );
  GEOS_ERROR_IF_NE( flowPattern.numColumns(), mechanicsPattern.numColumns() );

  // Count the union of the two sorted row sets. This avoids double-counting
  // shared diagonal entries without assuming either pattern contains the
  // other.
  array1d< localIndex > rowLengths( flowPattern.numRows() );
  rowLengths.zero();
  for( localIndex localRow = 0; localRow < flowPattern.numRows(); ++localRow )
  {
    arraySlice1d< globalIndex const > const flowColumns = flowPattern.getColumns( localRow );
    arraySlice1d< globalIndex const > const mechanicsColumns = mechanicsPattern.getColumns( localRow );
    localIndex flowColumn = 0;
    localIndex mechanicsColumn = 0;
    while( flowColumn < flowColumns.size() || mechanicsColumn < mechanicsColumns.size() )
    {
      if( mechanicsColumn == mechanicsColumns.size() ||
          ( flowColumn < flowColumns.size() && flowColumns[flowColumn] < mechanicsColumns[mechanicsColumn] ) )
      {
        ++flowColumn;
      }
      else if( flowColumn == flowColumns.size() || mechanicsColumns[mechanicsColumn] < flowColumns[flowColumn] )
      {
        ++mechanicsColumn;
      }
      else
      {
        ++flowColumn;
        ++mechanicsColumn;
      }
      ++rowLengths[localRow];
    }
  }

  // Add the number of nonzeros induced by coupling
  addTransmissibilityCouplingNNZ( domain, dofManager, rowLengths.toView() );
  addPressureForceCouplingNNZ( domain, dofManager, rowLengths.toView() );
  addMatrixPressureBubbleCouplingNNZ( domain, dofManager, rowLengths.toView() );

  // Allocate the coupled pattern in one pass. Growing an already populated
  // pattern row by row would shift every subsequent row of the contiguous
  // column buffer on each call, which is quadratic in the number of rows.
  pattern.resizeFromRowCapacities< parallelHostPolicy >( flowPattern.numRows(),
                                                         flowPattern.numColumns(),
                                                         rowLengths.data() );

  // Copy both subsolver patterns in. insertNonZeros discards entries that are
  // already present, so what remains is the union counted above.
  appendSparsityPattern( pattern, flowPattern );
  appendSparsityPattern( pattern, mechanicsPattern );

  // Add the nonzeros from coupling
  addTransmissibilityCouplingPattern( domain, dofManager, pattern.toView() );
  addPressureForceCouplingPattern( domain, dofManager, pattern.toView() );
  addMatrixPressureBubbleCouplingPattern( domain, dofManager, pattern.toView() );

  // Set up the derivative flux residual matrix
  setUpDflux_dApertureMatrix( domain );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::assembleSystem( real64 const time_n,
                                                                                    real64 const dt,
                                                                                    DomainPartition & domain,
                                                                                    DofManager const & dofManager,
                                                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                    arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  // Synchronize fracture state
  this->solidMechanicsSolver()->synchronizeFractureState( domain );

  // setSparsityPattern owns the contact lookup tables and this matrix. Rather
  // than lazily rebuilding a second copy of that prologue here, require that it
  // has run: setupSystem( ..., setSparsity = false ) is not supported.
  GEOS_ERROR_IF( !m_derivativeFluxResidual_dAperture,
                 GEOS_FMT( "{}: setupSystem must be called with sparsity construction enabled before assembling.",
                           this->getName() ) );

  // Move without touching: zero() below memsets the entries in this space, and
  // touching here would mark the immutable sparsity structure dirty on device.
  m_derivativeFluxResidual_dAperture->move( parallelDeviceMemorySpace, false );
  m_derivativeFluxResidual_dAperture->zero();

  // Assemble element-based contributions (mechanics + flow accumulation)
  assembleElementBasedContributions( time_n, dt, domain, dofManager, localMatrix, localRhs );

  // Assemble flux terms and get dFluidResidual/dAperture
  this->flowSolver()->assembleHydrofracFluxTerms( time_n,
                                                  dt,
                                                  domain,
                                                  dofManager,
                                                  localMatrix,
                                                  localRhs,
                                                  getDerivativeFluxResidual_dNormalJump(),
                                                  &m_derivativeFluxResidual_dApertureOffsets );

  // The flux kernel populates the derivative matrix in device memory. The
  // coupling assembly below reads it on the host, so bring it over without
  // touching it: touching would force a re-upload on the next assembly.
  m_derivativeFluxResidual_dAperture->move( hostMemorySpace, false );

  // Assemble coupling terms (must be after flux assembly)
  assembleCouplingTerms( time_n, dt, domain, dofManager, localMatrix, localRhs );

  if constexpr ( hasWells )
  {
    this->flowSolver()->wellSolver()->assembleSystem( time_n, dt, domain, dofManager, localMatrix, localRhs );
    this->flowSolver()->assembleCouplingTerms( time_n, dt, domain, dofManager, localMatrix, localRhs );
  }

}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::assembleElementBasedContributions( real64 const time_n,
                                                                                                       real64 const dt,
                                                                                                       DomainPartition & domain,
                                                                                                       DofManager const & dofManager,
                                                                                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                                       arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  // Assemble poromechanics terms (from base class)
  Base::assembleElementBasedTerms( time_n, dt, domain, dofManager, localMatrix, localRhs );

  // Flow accumulation for fractures
  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel & mesh,
                                                                      string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          FaceElementSubRegion const & subRegion )
    {
      this->flowSolver()->accumulationAssemblyLaunch( dofManager, subRegion, localMatrix, localRhs );
    } );
  } );

  // Assemble contact terms (ALM) - note: assembleContact requires time and dt
  this->solidMechanicsSolver()->assembleContact( time_n, dt, domain, dofManager, localMatrix, localRhs );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::assembleCouplingTerms( real64 const time_n,
                                                                                           real64 const dt,
                                                                                           DomainPartition const & domain,
                                                                                           DofManager const & dofManager,
                                                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                           arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( time_n, dt );

  // These steps must occur after the fluxes are assembled because that's when DerivativeFluxResidual_dAperture is filled.
  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                      MeshLevel const & mesh,
                                                                      string_array const & regionNames )
  {
    // Assemble Force Residual w.r.t. pressure (Aup) - fracture pressure contribution
    assembleForceResidualDerivativeWrtPressure( meshName, mesh, regionNames, dofManager, localMatrix, localRhs );

    // Assemble Fluid mass residual w.r.t. displacement (Apu)
    assembleFluidMassResidualDerivativeWrtDisplacement( meshName, mesh, regionNames, dofManager, localMatrix, localRhs );
  } );

  // Assemble matrix cell pressure contribution on bubble DOFs (Abp_matrix)
  // This must be outside the lambda because it uses regionBasedKernelApplication
  assembleMatrixPressureBubbleContribution( dt, const_cast< DomainPartition & >( domain ), dofManager, localMatrix, localRhs );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  // Call base poromechanics update
  Base::updateState( domain );

  // Need to call solid mechanics update separately to compute face displacement jump
  this->solidMechanicsSolver()->updateState( domain );

  // Remove the contribution of the hydraulic aperture from the stencil weights
  this->flowSolver()->prepareStencilWeights( domain );

  // Update hydraulic aperture and fracture permeability
  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel & mesh,
                                                                      string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   FaceElementSubRegion & subRegion )
    {
      arrayView2d< real64 const > const dispJump = subRegion.getField< contact::dispJump >();
      arrayView1d< real64 const > const area = subRegion.getElementArea();
      arrayView1d< real64 const > const volume = subRegion.getElementVolume();
      arrayView2d< real64 const > const fractureTraction = subRegion.getField< contact::traction >();
      arrayView1d< real64 const > const pressure = subRegion.getField< flow::pressure >();
      arrayView1d< real64 const > const oldHydraulicAperture = subRegion.getField< flow::aperture0 >();

      arrayView1d< real64 > const aperture = subRegion.getElementAperture();
      arrayView1d< real64 > const hydraulicAperture = subRegion.getField< flow::hydraulicAperture >();
      arrayView1d< real64 > const deltaVolume = subRegion.getField< flow::deltaVolume >();
      arrayView1d< integer > const & fractureState = subRegion.getField< contact::fractureState >();

      string const & hydraulicApertureRelationName = subRegion.template getReference< string >( Base::viewKeyStruct::hydraulicApertureRelationNameString() );
      HydraulicApertureBase const & hydraulicApertureModel = this->template getConstitutiveModel< HydraulicApertureBase >( subRegion, hydraulicApertureRelationName );

      string const porousSolidName = subRegion.getReference< string >( FlowSolverBase::viewKeyStruct::solidNamesString() );
      CoupledSolidBase & porousSolid = subRegion.template getConstitutiveModel< CoupledSolidBase >( porousSolidName );

      constitutiveUpdatePassThru( hydraulicApertureModel, [&] ( auto & castedHydraulicAperture )
      {
        using HydraulicApertureType = TYPEOFREF( castedHydraulicAperture );
        typename HydraulicApertureType::KernelWrapper hydraulicApertureWrapper = castedHydraulicAperture.createKernelWrapper();

        ConstitutivePassThru< CompressibleSolidBase >::execute( porousSolid, [=, &subRegion] ( auto & castedPorousSolid )
        {
          typename TYPEOFREF( castedPorousSolid ) ::KernelWrapper porousMaterialWrapper = castedPorousSolid.createKernelUpdates();

          poromechanicsFracturesKernels::StateUpdateKernel::
            launch< parallelDevicePolicy<> >( subRegion.size(),
                                              porousMaterialWrapper,
                                              hydraulicApertureWrapper,
                                              dispJump,
                                              pressure,
                                              area,
                                              volume,
                                              deltaVolume,
                                              aperture,
                                              oldHydraulicAperture,
                                              hydraulicAperture,
                                              fractureTraction,
                                              fractureState );
        } );
      } );
    } );
  } );

  // Update the stencil weights using the updated hydraulic aperture
  this->flowSolver()->updateStencilWeights( domain );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::
setUpDflux_dApertureMatrix( DomainPartition & domain )
{
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( this->flowSolver()->getDiscretizationName() );
  string const & fractureRegionName = this->solidMechanicsSolver()->getUniqueFractureRegionName();

  // Build the global row offsets and the row capacities together, so that each
  // target is visited only once before the matrix is allocated.
  m_derivativeFluxResidual_dApertureOffsets.clear();
  stdVector< localIndex > rowCapacities;
  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                      MeshLevel const & mesh,
                                                                      string_array const & regionNames )
  {
    GEOS_UNUSED_VAR( regionNames );
    ElementRegionManager const & elemManager = mesh.getElemManager();

    // These offsets are consumed by the flow sub-solver, which walks its own
    // mesh targets and therefore resolves the discretization level with its own
    // discretization name. The mesh body name is the only part of a target both
    // solvers are guaranteed to agree on, so it is the key; that in turn
    // requires each body to appear exactly once here.
    GEOS_ERROR_IF( m_derivativeFluxResidual_dApertureOffsets.find( meshName ) !=
                   m_derivativeFluxResidual_dApertureOffsets.end(),
                   GEOS_FMT( "{}: mesh body '{}' is targeted at more than one discretization level. The augmented "
                             "Lagrangian contact formulation supports a single level per mesh body.",
                             this->getName(), meshName ) );

    localIndex const rowOffset = rowCapacities.size();
    m_derivativeFluxResidual_dApertureOffsets.get_inserted( meshName ) = rowOffset;

    // The stencil sweeps below index rows by the raw surface-element index, so
    // the contact fracture must be the only face-element region on this target:
    // a second one would alias into its rows. Embedded-surface regions hold a
    // different subregion type and contribute no SurfaceElementStencil here, so
    // they are left alone. The region is required rather than optional because
    // every consumer of this matrix (assembleCouplingTerms,
    // assembleFluidMassResidualDerivativeWrtDisplacement) looks it up
    // unconditionally on every target; skipping a target here would also leave
    // its offset pointing at the next target's rows.
    localIndex numFractureRegions = 0;
    elemManager.forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion const & region )
    {
      if( region.subRegionType() == SurfaceElementRegion::SurfaceSubRegionType::faceElement )
      {
        ++numFractureRegions;
      }
    } );
    GEOS_ERROR_IF_NE_MSG( numFractureRegions, 1,
                          GEOS_FMT( "{}: mesh target '{}' holds {} face-element regions. The augmented Lagrangian "
                                    "contact formulation requires exactly one, named '{}'.",
                                    this->getName(), meshName, numFractureRegions, fractureRegionName ) );
    GEOS_ERROR_IF( !elemManager.hasRegion( fractureRegionName ),
                   GEOS_FMT( "{}: mesh target '{}' does not hold the fracture region '{}' of the contact solver.",
                             this->getName(), meshName, fractureRegionName ) );

    SurfaceElementRegion const & fractureRegion = elemManager.getRegion< SurfaceElementRegion >( fractureRegionName );
    FaceElementSubRegion const & fractureSubRegion = fractureRegion.getUniqueSubRegion< FaceElementSubRegion >();
    rowCapacities.resize( rowOffset + fractureSubRegion.size(), 0 );

    fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
    {
      for( localIndex iconn = 0; iconn < stencil.size(); ++iconn )
      {
        localIndex const numFluxElems = stencil.stencilSize( iconn );
        typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();
        for( localIndex k0 = 0; k0 < numFluxElems; ++k0 )
        {
          localIndex const row = rowOffset + sei[iconn][k0];
          GEOS_ERROR_IF_GE_MSG( row,
                                LvArray::integerConversion< localIndex >( rowCapacities.size() ),
                                "Surface stencil index exceeds the fracture derivative matrix size." );
          rowCapacities[row] += numFluxElems;
        }
      }
    } );
  } );

  std::unique_ptr< CRSMatrix< real64, localIndex > > & derivativeFluxResidual_dAperture = getRefDerivativeFluxResidual_dAperture();
  localIndex const numRows = rowCapacities.size();
  derivativeFluxResidual_dAperture = std::make_unique< CRSMatrix< real64, localIndex > >( numRows, numRows );
  derivativeFluxResidual_dAperture->setName( this->getName() + "/derivativeFluxResidual_dAperture" );
  if( numRows > 0 )
  {
    derivativeFluxResidual_dAperture->resizeFromRowCapacities< parallelHostPolicy >( numRows,
                                                                                     numRows,
                                                                                     rowCapacities.data() );
  }

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                      MeshLevel const & mesh,
                                                                      string_array const & regionNames )
  {
    GEOS_UNUSED_VAR( regionNames );
    localIndex const rowOffset = m_derivativeFluxResidual_dApertureOffsets.at( meshName );

    fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
    {
      for( localIndex iconn = 0; iconn < stencil.size(); ++iconn )
      {
        localIndex const numFluxElems = stencil.stencilSize( iconn );
        typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

        for( localIndex k0 = 0; k0 < numFluxElems; ++k0 )
        {
          localIndex const row = rowOffset + sei[iconn][k0];
          GEOS_ERROR_IF_GE_MSG( row, numRows, "Surface stencil index exceeds the fracture derivative matrix size." );
          for( localIndex k1 = 0; k1 < numFluxElems; ++k1 )
          {
            derivativeFluxResidual_dAperture->insertNonZero( row,
                                                             rowOffset + sei[iconn][k1],
                                                             0.0 );
          }
        }
      }
    } );
  } );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::
addTransmissibilityCouplingNNZ( DomainPartition const & domain,
                                DofManager const & dofManager,
                                arrayView1d< localIndex > const & rowLengths ) const
{
  GEOS_MARK_FUNCTION;

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel const & mesh,
                                                                      string_array const & )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const flowDofKey = dofManager.getKey( m_pressureKey );

    globalIndex const rankOffset = dofManager.rankOffset();

    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( this->flowSolver()->getDiscretizationName() );

    fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
    {
      for( localIndex iconn = 0; iconn < stencil.size(); ++iconn )
      {
        localIndex const numFluxElems = stencil.stencilSize( iconn );
        typename SurfaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
        typename SurfaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
        typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

        FaceElementSubRegion const & elementSubRegion =
          elemManager.getRegion( seri[iconn][0] ).getSubRegion< FaceElementSubRegion >( sesri[iconn][0] );

        ArrayOfArraysView< localIndex const > const elemsToNodes = elementSubRegion.nodeList().toViewConst();

        arrayView1d< globalIndex const > const faceElementDofNumber =
          elementSubRegion.getReference< array1d< globalIndex > >( flowDofKey );

        for( localIndex k0 = 0; k0 < numFluxElems; ++k0 )
        {
          globalIndex const activeFlowDOF = faceElementDofNumber[sei[iconn][k0]];
          globalIndex const rowNumber = activeFlowDOF - rankOffset;

          if( rowNumber >= 0 && rowNumber < rowLengths.size() )
          {
            for( localIndex k1 = 0; k1 < numFluxElems; ++k1 )
            {
              // The coupling with the nodal displacements of the cell itself has already been added by the dofManager
              // so we only add the coupling with the nodal displacements of the neighbors.
              if( k1 != k0 )
              {
                localIndex const numNodesPerElement = elemsToNodes[sei[iconn][k1]].size();
                // Nodal displacement DOFs (3 per node, 2 faces)
                rowLengths[rowNumber] += 3 * numNodesPerElement;
                // Bubble DOFs (3 per face, 2 faces = 6 total)
                rowLengths[rowNumber] += 6;
              }
            }
          }
        }
      }
    } );
  } );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::
addTransmissibilityCouplingPattern( DomainPartition const & domain,
                                    DofManager const & dofManager,
                                    SparsityPatternView< globalIndex > const & pattern ) const
{
  GEOS_MARK_FUNCTION;

  using namespace contact;

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel const & mesh,
                                                                      string_array const & )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    NodeManager const & nodeManager = mesh.getNodeManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
    string const bubbleDofKey = dofManager.getKey( totalBubbleDisplacement::key() );
    string const flowDofKey = dofManager.getKey( m_pressureKey );

    arrayView1d< globalIndex const > const &
    dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
    arrayView1d< globalIndex const > const &
    bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );
    ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();

    // Get the finite volume method used to compute the fluxes
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fvDiscretization = fvManager.getFluxApproximation( this->flowSolver()->getDiscretizationName() );

    SurfaceElementRegion const & fractureRegion =
      elemManager.getRegion< SurfaceElementRegion >( this->solidMechanicsSolver()->getUniqueFractureRegionName() );
    FaceElementSubRegion const & fractureSubRegion =
      fractureRegion.getUniqueSubRegion< FaceElementSubRegion >();

    GEOS_ERROR_IF( !fractureSubRegion.hasWrapper( flow::pressure::key() ),
                   this->getDataContext() << ": The fracture subregion must contain pressure field." );

    arrayView2d< localIndex const > const elem2dToFaces = fractureSubRegion.faceList().toViewConst();

    arrayView1d< globalIndex const > const &
    flowDofNumber = fractureSubRegion.getReference< globalIndex_array >( flowDofKey );

    globalIndex const rankOffset = dofManager.rankOffset();

    fvDiscretization.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
    {
      forAll< serialPolicy >( stencil.size(), [=] ( localIndex const iconn )
      {
        localIndex const numFluxElems = stencil.stencilSize( iconn );

        // A fracture connector has to be an edge shared by two faces
        if( numFluxElems == 2 )
        {
          typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

          // First index: face element. Second index: node
          for( localIndex kf = 0; kf < 2; ++kf )
          {
            // Set row DOF index
            // Note that the 1-kf index is intentional, as this is coupling the pressure of one face cell
            // to the nodes of the adjacent cell
            localIndex const rowIndex = flowDofNumber[sei[iconn][1 - kf]] - rankOffset;

            if( rowIndex >= 0 && rowIndex < pattern.numRows() )
            {
              // Get fracture, face and region/subregion/element indices (for elements on both sides)
              localIndex const fractureIndex = sei[iconn][kf];

              // Get the number of nodes
              localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elem2dToFaces[fractureIndex][0] );

              // Loop over the two sides of each fracture element
              for( localIndex kf1 = 0; kf1 < 2; ++kf1 )
              {
                localIndex const faceIndex = elem2dToFaces[fractureIndex][kf1];

                // Save the list of DOF associated with nodes (displacement)
                for( localIndex a = 0; a < numNodesPerFace; ++a )
                {
                  for( localIndex i = 0; i < 3; ++i )
                  {
                    globalIndex const colIndex = dispDofNumber[faceToNodeMap( faceIndex, a )] + LvArray::integerConversion< globalIndex >( i );
                    pattern.insertNonZero( rowIndex, colIndex );
                  }
                }

                // Save the list of DOF associated with bubble displacement
                for( localIndex i = 0; i < 3; ++i )
                {
                  globalIndex const colIndex = bubbleDofNumber[faceIndex] + LvArray::integerConversion< globalIndex >( i );
                  pattern.insertNonZero( rowIndex, colIndex );
                }
              }
            }
          }
        }
      } );
    } );
  } );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::
addPressureForceCouplingNNZ( DomainPartition const & domain,
                             DofManager const & dofManager,
                             arrayView1d< localIndex > const & rowLengths ) const
{
  GEOS_MARK_FUNCTION;

  using namespace contact;

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel const & mesh,
                                                                      string_array const & )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    NodeManager const & nodeManager = mesh.getNodeManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
    string const bubbleDofKey = dofManager.getKey( totalBubbleDisplacement::key() );

    arrayView1d< globalIndex const > const &
    dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
    arrayView1d< globalIndex const > const &
    bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );
    ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();

    globalIndex const rankOffset = dofManager.rankOffset();

    string const & fractureRegionName = this->solidMechanicsSolver()->getUniqueFractureRegionName();
    SurfaceElementRegion const & fractureRegion =
      elemManager.getRegion< SurfaceElementRegion >( fractureRegionName );
    FaceElementSubRegion const & fractureSubRegion =
      fractureRegion.getUniqueSubRegion< FaceElementSubRegion >();

    arrayView2d< localIndex const > const elem2dToFaces = fractureSubRegion.faceList().toViewConst();

    // For each fracture element, add NNZ for (displacement_row, pressure_col) and (bubble_row, pressure_col)
    forAll< serialPolicy >( fractureSubRegion.size(), [=, &rowLengths] ( localIndex const kfe )
    {
      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elem2dToFaces[kfe][0] );

      // For displacement DOFs: add 1 pressure column per displacement DOF row
      for( localIndex kf = 0; kf < 2; ++kf )
      {
        localIndex const faceIndex = elem2dToFaces[kfe][kf];
        for( localIndex a = 0; a < numNodesPerFace; ++a )
        {
          for( localIndex i = 0; i < 3; ++i )
          {
            globalIndex const rowNumber = dispDofNumber[faceToNodeMap( faceIndex, a )] + i - rankOffset;
            if( rowNumber >= 0 && rowNumber < rowLengths.size() )
            {
              rowLengths[rowNumber] += 1;  // One pressure column
            }
          }
        }
      }

      // For bubble DOFs: add 1 pressure column per bubble DOF row
      for( localIndex kf = 0; kf < 2; ++kf )
      {
        localIndex const faceIndex = elem2dToFaces[kfe][kf];
        for( localIndex i = 0; i < 3; ++i )
        {
          globalIndex const rowNumber = bubbleDofNumber[faceIndex] + i - rankOffset;
          if( rowNumber >= 0 && rowNumber < rowLengths.size() )
          {
            rowLengths[rowNumber] += 1;  // One pressure column
          }
        }
      }
    } );
  } );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::
addPressureForceCouplingPattern( DomainPartition const & domain,
                                 DofManager const & dofManager,
                                 SparsityPatternView< globalIndex > const & pattern ) const
{
  GEOS_MARK_FUNCTION;

  using namespace contact;

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel const & mesh,
                                                                      string_array const & )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    NodeManager const & nodeManager = mesh.getNodeManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
    string const bubbleDofKey = dofManager.getKey( totalBubbleDisplacement::key() );
    string const flowDofKey = dofManager.getKey( m_pressureKey );

    arrayView1d< globalIndex const > const &
    dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
    arrayView1d< globalIndex const > const &
    bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );
    ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();

    globalIndex const rankOffset = dofManager.rankOffset();

    string const & fractureRegionName = this->solidMechanicsSolver()->getUniqueFractureRegionName();
    SurfaceElementRegion const & fractureRegion =
      elemManager.getRegion< SurfaceElementRegion >( fractureRegionName );
    FaceElementSubRegion const & fractureSubRegion =
      fractureRegion.getUniqueSubRegion< FaceElementSubRegion >();

    arrayView2d< localIndex const > const elem2dToFaces = fractureSubRegion.faceList().toViewConst();
    arrayView1d< globalIndex const > const &
    flowDofNumber = fractureSubRegion.getReference< globalIndex_array >( flowDofKey );

    // For each fracture element, add pattern for (displacement_row, pressure_col) and (bubble_row, pressure_col)
    forAll< serialPolicy >( fractureSubRegion.size(), [=] ( localIndex const kfe )
    {
      globalIndex const pressureColIndex = flowDofNumber[kfe];
      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elem2dToFaces[kfe][0] );

      // For displacement DOFs
      for( localIndex kf = 0; kf < 2; ++kf )
      {
        localIndex const faceIndex = elem2dToFaces[kfe][kf];
        for( localIndex a = 0; a < numNodesPerFace; ++a )
        {
          for( localIndex i = 0; i < 3; ++i )
          {
            globalIndex const rowIndex = dispDofNumber[faceToNodeMap( faceIndex, a )] + i - rankOffset;
            if( rowIndex >= 0 && rowIndex < pattern.numRows() )
            {
              pattern.insertNonZero( rowIndex, pressureColIndex );
            }
          }
        }
      }

      // For bubble DOFs
      for( localIndex kf = 0; kf < 2; ++kf )
      {
        localIndex const faceIndex = elem2dToFaces[kfe][kf];
        for( localIndex i = 0; i < 3; ++i )
        {
          globalIndex const rowIndex = bubbleDofNumber[faceIndex] + i - rankOffset;
          if( rowIndex >= 0 && rowIndex < pattern.numRows() )
          {
            pattern.insertNonZero( rowIndex, pressureColIndex );
          }
        }
      }
    } );
  } );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::
assembleForceResidualDerivativeWrtPressure( string const & meshName,
                                            MeshLevel const & mesh,
                                            string_array const & regionNames,
                                            DofManager const & dofManager,
                                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                            arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  using namespace contact;

  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager const & nodeManager = mesh.getNodeManager();

  string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
  string const & bubbleDofKey = dofManager.getKey( totalBubbleDisplacement::key() );
  string const & flowDofKey = dofManager.getKey( m_pressureKey );

  arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
  arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );

  string const & fractureRegionName = this->solidMechanicsSolver()->getUniqueFractureRegionName();

  // Use the same kernel launch pattern as SolidMechanicsAugmentedLagrangianContact::assembleForceResidualPressureContribution
  this->solidMechanicsSolver()->forFiniteElementOnFractureSubRegions( meshName,
                                                                      [&] ( string const &,
                                                                            finiteElement::FiniteElementBase const & subRegionFE,
                                                                            arrayView1d< localIndex const > const & faceElementList )
  {
    // Get pressure DOF number from the fracture subregion
    SurfaceElementRegion const & fractureRegion = mesh.getElemManager().getRegion< SurfaceElementRegion >( fractureRegionName );
    FaceElementSubRegion const & fractureSubRegion = fractureRegion.getUniqueSubRegion< FaceElementSubRegion >();
    arrayView1d< globalIndex const > const pressureDofNumber = fractureSubRegion.getReference< array1d< globalIndex > >( flowDofKey );

    poromechanicsALMKernels::AssembleForceResidualDerivativeWrtPressureFactory
    kernelFactory( dispDofNumber,
                   bubbleDofNumber,
                   dofManager.rankOffset(),
                   localMatrix,
                   localRhs,
                   0.0,  // dt not used
                   faceElementList,
                   pressureDofNumber );

    // Note: const_cast is needed because interfaceBasedKernelApplication takes non-const mesh
    // even though it only modifies the matrix/rhs which are passed separately
    real64 maxResidual = finiteElement::
                           interfaceBasedKernelApplication
                         < parallelDevicePolicy<>,
                           constitutive::NullModel >( const_cast< MeshLevel & >( mesh ),
                                                      fractureRegionName,
                                                      faceElementList,
                                                      subRegionFE,
                                                      "",
                                                      kernelFactory );

    GEOS_UNUSED_VAR( maxResidual );
  } );

  GEOS_UNUSED_VAR( regionNames );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::
assembleFluidMassResidualDerivativeWrtDisplacement( string const & meshName,
                                                    MeshLevel const & mesh,
                                                    string_array const & regionNames,
                                                    DofManager const & dofManager,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( regionNames );

  using namespace contact;

  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager const & nodeManager = mesh.getNodeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();

  // assembleSystem has already brought this matrix to the host after the flux
  // assembly; the traversal below only reads it.
  CRSMatrixView< real64 const, localIndex const > const &
  dFluxResidual_dNormalJump = getDerivativeFluxResidual_dNormalJump().toViewConst();
  auto const derivativeOffsetIt = m_derivativeFluxResidual_dApertureOffsets.find( meshName );
  GEOS_ERROR_IF( derivativeOffsetIt == m_derivativeFluxResidual_dApertureOffsets.end(),
                 GEOS_FMT( "No dR/dAperture row offset is available for mesh body '{}'", meshName ) );
  localIndex const derivativeOffset = derivativeOffsetIt->second;

  string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
  string const & bubbleDofKey = dofManager.getKey( totalBubbleDisplacement::key() );
  string const & presDofKey = dofManager.getKey( m_pressureKey );

  arrayView1d< globalIndex const > const &
  dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
  arrayView1d< globalIndex const > const &
  bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );
  globalIndex const rankOffset = dofManager.rankOffset();

  string const & fractureRegionName = this->solidMechanicsSolver()->getUniqueFractureRegionName();

  // Maximum DOF sizes
  constexpr localIndex maxNumNodesPerFace = m_maxFaceNodes;
  constexpr localIndex maxNumUdofs = 2 * 3 * maxNumNodesPerFace; // 66
  constexpr localIndex numBdofs = 6;

  // Get the fracture subRegion
  SurfaceElementRegion const & fractureRegion = elemManager.getRegion< SurfaceElementRegion >( fractureRegionName );
  FaceElementSubRegion const & subRegion = fractureRegion.getUniqueSubRegion< FaceElementSubRegion >();

  localIndex const numElems = subRegion.size();

  // These arrays are temporary workspaces for one assembly call. Keeping them
  // local avoids retaining potentially large aperture buffers for the solver's
  // lifetime while still giving the device kernels explicit views to capture.
  array2d< real64 > dAperturedU( numElems, maxNumUdofs );
  array2d< real64 > dAperturedB( numElems, numBdofs );
  dAperturedU.zero();
  dAperturedB.zero();
  dAperturedU.move( parallelDeviceMemorySpace, true );
  dAperturedB.move( parallelDeviceMemorySpace, true );
  arrayView2d< real64 > const dAperturedUView = dAperturedU.toView();
  arrayView2d< real64 > const dAperturedBView = dAperturedB.toView();

  // Launch the ComputeApertureDerivatives kernel to fill dAperturedU and dAperturedB
  // This is called for each element type (tri, quad, etc.)
  this->solidMechanicsSolver()->forFiniteElementOnFractureSubRegions( meshName,
                                                                      [&, dAperturedUView, dAperturedBView] ( string const &,
                                                                                                              finiteElement::FiniteElementBase const & subRegionFE,
                                                                                                              arrayView1d< localIndex const > const & faceElementList )
  {
    poromechanicsALMKernels::ComputeApertureDerivativesFactory
    kernelFactory( dispDofNumber,
                   bubbleDofNumber,
                   dofManager.rankOffset(),
                   localMatrix,
                   localRhs,
                   0.0,  // dt not used
                   faceElementList,
                   dAperturedUView,
                   dAperturedBView );

    real64 maxResidual = finiteElement::
                           interfaceBasedKernelApplication
                         < parallelDevicePolicy<>,
                           constitutive::NullModel >( const_cast< MeshLevel & >( mesh ),
                                                      fractureRegionName,
                                                      faceElementList,
                                                      subRegionFE,
                                                      "",
                                                      kernelFactory );

    GEOS_UNUSED_VAR( maxResidual );
  } );

  // Move data to host for serial assembly
  dAperturedU.move( hostMemorySpace );
  dAperturedB.move( hostMemorySpace );

  // Now assemble using the pre-computed derivatives
  string const & fluidName = subRegion.getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() );

  SingleFluidBase const & fluid = this->template getConstitutiveModel< SingleFluidBase >( subRegion, fluidName );
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & density = fluid.density();

  arrayView1d< globalIndex const > const & presDofNumber = subRegion.getReference< array1d< globalIndex > >( presDofKey );

  arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList().toViewConst();

  arrayView1d< integer const > const fractureState = subRegion.getField< contact::fractureState >();

  // Get element area for proper scaling
  // Note: dAperturedU/dB are computed as (1/area) * unitNormal^T * Atu/Atb
  // For accumulation: dR_accum/du = density * unitNormal^T * Atu (no 1/area factor)
  // For flux: dR_flux/du = dR/dAperture * (1/area) * unitNormal^T * Atu (with 1/area factor)
  // So accumulation needs to multiply by area to cancel the 1/area in dAperturedU
  arrayView1d< real64 const > const area = subRegion.getElementArea().toViewConst();

  forAll< serialPolicy >( numElems, [&]( localIndex const kfe )
  {
    localIndex const kf0 = elemsToFaces[kfe][0];
    localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );
    localIndex const numUdofs = 2 * 3 * numNodesPerFace;

    globalIndex nodeDOF[maxNumUdofs];
    globalIndex elemDOF[1];
    elemDOF[0] = presDofNumber[kfe];

    stackArray1d< real64, maxNumUdofs > dRdU( maxNumUdofs );

    bool const isFractureOpen = ( fractureState[kfe] == FractureState::Open );

    // ==== Part 1: Apu - Nodal displacement contribution ====
    // Get DOF indices for displacement
    for( localIndex kf = 0; kf < 2; ++kf )
    {
      for( localIndex a = 0; a < numNodesPerFace; ++a )
      {
        for( localIndex i = 0; i < 3; ++i )
        {
          nodeDOF[kf * 3 * numNodesPerFace + 3 * a + i] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe][kf], a )]
                                                          + LvArray::integerConversion< globalIndex >( i );
        }
      }
    }

    // Accumulation derivative w.r.t. nodal displacement
    // dR_accum/du = density * unitNormal^T * Atu = density * area * dAperturedU
    // (dAperturedU already has 1/area factor, so multiply by area to cancel it)
    if( isFractureOpen )
    {
      for( localIndex j = 0; j < numUdofs; ++j )
      {
        dRdU( j ) = density[kfe][0] * dAperturedU( kfe, j ) * area[kfe];
      }

      localIndex const localRow = LvArray::integerConversion< localIndex >( elemDOF[0] - rankOffset );

      if( localRow >= 0 && localRow < localMatrix.numRows() )
      {
        localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow,
                                                                  nodeDOF,
                                                                  dRdU.data(),
                                                                  numUdofs );
      }
    }

    // Flux derivative w.r.t. nodal displacement
    localIndex const numColumns = dFluxResidual_dNormalJump.numNonZeros( derivativeOffset + kfe );
    arraySlice1d< localIndex const > const & columns = dFluxResidual_dNormalJump.getColumns( derivativeOffset + kfe );
    arraySlice1d< real64 const > const & values = dFluxResidual_dNormalJump.getEntries( derivativeOffset + kfe );

    for( localIndex kfe1 = 0; kfe1 < numColumns; ++kfe1 )
    {
      real64 const dR_dAper = values[kfe1];
      localIndex const kfe2 = columns[kfe1] - derivativeOffset;

      bool const isOpen = ( fractureState[kfe2] == FractureState::Open );
      if( !isOpen && !isFractureOpen )
        continue;

      localIndex const kf0_2 = elemsToFaces[kfe2][0];
      localIndex const numNodesPerFace2 = faceToNodeMap.sizeOfArray( kf0_2 );
      localIndex const numUdofs2 = 2 * 3 * numNodesPerFace2;

      // Get DOF indices for element kfe2
      globalIndex nodeDOF2[maxNumUdofs];
      for( localIndex kf = 0; kf < 2; ++kf )
      {
        for( localIndex a = 0; a < numNodesPerFace2; ++a )
        {
          for( localIndex i = 0; i < 3; ++i )
          {
            nodeDOF2[kf * 3 * numNodesPerFace2 + 3 * a + i] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe2][kf], a )]
                                                              + LvArray::integerConversion< globalIndex >( i );
          }
        }
      }

      // dR_flux/du = dR_flux/dAper * dAper/du (pre-computed for element kfe2)
      stackArray1d< real64, maxNumUdofs > dRdU2( maxNumUdofs );
      for( localIndex j = 0; j < numUdofs2; ++j )
      {
        dRdU2( j ) = dR_dAper * dAperturedU( kfe2, j );
      }

      localIndex const localRow = LvArray::integerConversion< localIndex >( elemDOF[0] - rankOffset );

      if( localRow >= 0 && localRow < localMatrix.numRows() )
      {
        localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow,
                                                                  nodeDOF2,
                                                                  dRdU2.data(),
                                                                  numUdofs2 );
      }
    }

    // ==== Part 2: Apb - Bubble displacement contribution ====
    globalIndex bubbleDOF[numBdofs];

    // Get DOF indices for bubble
    for( localIndex kf = 0; kf < 2; ++kf )
    {
      localIndex const faceIndex = elemsToFaces[kfe][kf];
      for( localIndex i = 0; i < 3; ++i )
      {
        bubbleDOF[kf * 3 + i] = bubbleDofNumber[faceIndex] + LvArray::integerConversion< globalIndex >( i );
      }
    }

    // Accumulation derivative w.r.t. bubble displacement
    // dR_accum/db = density * unitNormal^T * Atb = density * area * dAperturedB
    // (dAperturedB already has 1/area factor, so multiply by area to cancel it)
    if( isFractureOpen )
    {
      real64 dRdB[numBdofs];
      for( localIndex j = 0; j < numBdofs; ++j )
      {
        dRdB[j] = density[kfe][0] * dAperturedB( kfe, j ) * area[kfe];
      }

      localIndex const localRow = LvArray::integerConversion< localIndex >( elemDOF[0] - rankOffset );

      if( localRow >= 0 && localRow < localMatrix.numRows() )
      {
        localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow,
                                                                  bubbleDOF,
                                                                  dRdB,
                                                                  numBdofs );
      }
    }

    // Flux derivative w.r.t. bubble displacement
    for( localIndex kfe1 = 0; kfe1 < numColumns; ++kfe1 )
    {
      real64 const dR_dAper = values[kfe1];
      localIndex const kfe2 = columns[kfe1] - derivativeOffset;

      bool const isOpen = ( fractureState[kfe2] == FractureState::Open );
      if( !isOpen && !isFractureOpen )
        continue;

      // Get DOF indices for bubble of element kfe2
      globalIndex bubbleDOF2[numBdofs];
      for( localIndex kf = 0; kf < 2; ++kf )
      {
        localIndex const faceIndex = elemsToFaces[kfe2][kf];
        for( localIndex i = 0; i < 3; ++i )
        {
          bubbleDOF2[kf * 3 + i] = bubbleDofNumber[faceIndex] + LvArray::integerConversion< globalIndex >( i );
        }
      }

      // dR_flux/db = dR_flux/dAper * dAper/db (pre-computed for element kfe2)
      real64 dRdB2[numBdofs];
      for( localIndex j = 0; j < numBdofs; ++j )
      {
        dRdB2[j] = dR_dAper * dAperturedB( kfe2, j );
      }

      localIndex const localRow = LvArray::integerConversion< localIndex >( elemDOF[0] - rankOffset );

      if( localRow >= 0 && localRow < localMatrix.numRows() )
      {
        localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow,
                                                                  bubbleDOF2,
                                                                  dRdB2,
                                                                  numBdofs );
      }
    }
  } );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::
addMatrixPressureBubbleCouplingNNZ( DomainPartition const & domain,
                                    DofManager const & dofManager,
                                    arrayView1d< localIndex > const & rowLengths ) const
{
  GEOS_MARK_FUNCTION;

  using namespace contact;

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel const & mesh,
                                                                      string_array const & regionNames )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const bubbleDofKey = dofManager.getKey( totalBubbleDisplacement::key() );
    string const flowDofKey = dofManager.getKey( m_pressureKey );
    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );

    globalIndex const rankOffset = dofManager.rankOffset();

    // Loop over matrix cell regions that have bubbles
    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                              [&]( localIndex const, CellElementSubRegion const & subRegion )
    {
      arrayView1d< localIndex const > const bubbleElems = subRegion.bubbleElementsList();
      arrayView2d< localIndex const > const elemsToFaces = subRegion.faceElementsList();
      arrayView1d< globalIndex const > const pressureDofNumber = subRegion.getReference< array1d< globalIndex > >( flowDofKey );

      forAll< serialPolicy >( bubbleElems.size(), [=, &rowLengths]( localIndex const kk )
      {
        localIndex const k = bubbleElems[kk];
        localIndex const faceIndex = elemsToFaces[kk][0];

        // (bubble_row, pressure_col): 1 pressure column for each of the 3 bubble DOFs
        for( localIndex i = 0; i < 3; ++i )
        {
          globalIndex const rowNumber = bubbleDofNumber[faceIndex] + i - rankOffset;
          if( rowNumber >= 0 && rowNumber < rowLengths.size() )
          {
            rowLengths[rowNumber] += 1;  // One pressure DOF from matrix cell
          }
        }

        // (pressure_row, bubble_col): the matrix cell pressure couples to its 3 bubble DOFs (A_pb)
        globalIndex const pRow = pressureDofNumber[k] - rankOffset;
        if( pRow >= 0 && pRow < rowLengths.size() )
        {
          rowLengths[pRow] += 3;  // Three bubble DOFs
        }
      } );
    } );
  } );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::
addMatrixPressureBubbleCouplingPattern( DomainPartition const & domain,
                                        DofManager const & dofManager,
                                        SparsityPatternView< globalIndex > const & pattern ) const
{
  GEOS_MARK_FUNCTION;

  using namespace contact;

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel const & mesh,
                                                                      string_array const & regionNames )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const bubbleDofKey = dofManager.getKey( totalBubbleDisplacement::key() );
    string const flowDofKey = dofManager.getKey( m_pressureKey );

    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );

    globalIndex const rankOffset = dofManager.rankOffset();

    // Loop over matrix cell regions that have bubbles
    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                              [&]( localIndex const, CellElementSubRegion const & subRegion )
    {
      arrayView1d< localIndex const > const bubbleElems = subRegion.bubbleElementsList();
      arrayView2d< localIndex const > const elemsToFaces = subRegion.faceElementsList();
      arrayView1d< globalIndex const > const pressureDofNumber = subRegion.getReference< array1d< globalIndex > >( flowDofKey );

      forAll< serialPolicy >( bubbleElems.size(), [=]( localIndex const kk )
      {
        localIndex const k = bubbleElems[kk];
        localIndex const faceIndex = elemsToFaces[kk][0];
        globalIndex const pressureColIndex = pressureDofNumber[k];

        // (bubble_row, pressure_col) : A_bp
        for( localIndex i = 0; i < 3; ++i )
        {
          globalIndex const rowIndex = bubbleDofNumber[faceIndex] + i - rankOffset;
          if( rowIndex >= 0 && rowIndex < pattern.numRows() )
          {
            pattern.insertNonZero( rowIndex, pressureColIndex );
          }
        }

        // (pressure_row, bubble_col) : A_pb -- transpose location
        globalIndex const pRow = pressureDofNumber[k] - rankOffset;
        if( pRow >= 0 && pRow < pattern.numRows() )
        {
          for( localIndex i = 0; i < 3; ++i )
          {
            pattern.insertNonZero( pRow, bubbleDofNumber[faceIndex] + i );
          }
        }
      } );
    } );
  } );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::
assembleMatrixPressureBubbleContribution( real64 const dt,
                                          DomainPartition & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  using namespace contact;

  string const flowDofKey = dofManager.getKey( m_pressureKey );
  string const mechanicsDiscretizationName = this->solidMechanicsSolver()->getDiscretizationName();

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel & mesh,
                                                                      string_array const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
    string const & bubbleDofKey = dofManager.getKey( totalBubbleDisplacement::key() );
    arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );

    // Identify poromechanics regions (cells with porous solid)
    set< string > poromechanicsRegions;
    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                              [&]( localIndex const regionIndex, CellElementSubRegion const & subRegion )
    {
      if( subRegion.hasWrapper( FlowSolverBase::viewKeyStruct::solidNamesString() ) )
      {
        poromechanicsRegions.insert( regionNames[regionIndex] );
      }
    } );

    string_array poromechanicsRegionNames;
    poromechanicsRegionNames.reserve( poromechanicsRegions.size() );
    for( auto const & region : poromechanicsRegions )
    {
      poromechanicsRegionNames.emplace_back( region );
    }

    // Launch the kernel on matrix cells with bubbles
    poromechanicsMatrixBubbleKernels::MatrixPressureBubbleFactory kernelFactory( dispDofNumber,
                                                                                 bubbleDofNumber,
                                                                                 dofManager.rankOffset(),
                                                                                 localMatrix,
                                                                                 localRhs,
                                                                                 dt,
                                                                                 flowDofKey,
                                                                                 FlowSolverBase::viewKeyStruct::fluidNamesString() );

    real64 maxResidual = finiteElement::regionBasedKernelApplication
                         < parallelDevicePolicy<>,
                           constitutive::PorousSolidBase,
                           CellElementSubRegion >( mesh,
                                                   poromechanicsRegionNames,
                                                   mechanicsDiscretizationName,
                                                   FlowSolverBase::viewKeyStruct::solidNamesString(),
                                                   kernelFactory );

    GEOS_UNUSED_VAR( maxResidual );
  } );
}


template class SinglePhasePoromechanicsConformingFracturesALM<>;
template class SinglePhasePoromechanicsConformingFracturesALM< SinglePhaseReservoirAndWells<> >;

namespace
{
typedef SinglePhasePoromechanicsConformingFracturesALM< SinglePhaseReservoirAndWells<> > SinglePhaseReservoirPoromechanicsConformingFracturesALM;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SinglePhaseReservoirPoromechanicsConformingFracturesALM, string const &, Group * const )
typedef SinglePhasePoromechanicsConformingFracturesALM<> SinglePhasePoromechanicsConformingFracturesALM;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SinglePhasePoromechanicsConformingFracturesALM, string const &, Group * const )
}

} /* namespace geos */
