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
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::setSparsityPattern( DomainPartition & domain,
                                                                                        DofManager & dofManager,
                                                                                        CRSMatrix< real64, globalIndex > & localMatrix,
                                                                                        SparsityPattern< globalIndex > & pattern )
{
  GEOS_MARK_FUNCTION;

  // Recompute fracture face/element geometry and rebuild the ALM contact solver's internal lists.
  // These must happen before assembling the contact-dependent pattern; setSparsityPattern() 
  this->solidMechanicsSolver()->updateFractureGeometry( domain );

  Base::setSparsityPattern(domain,dofManager,localMatrix, pattern);
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
  string const & flowDofKey = dofManager.getKey( this->getFlowDofKey() );

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

  // The mass-balance block of getDerivativeFluxResidual_dNormalJump() below is one row per
  // fracture element, as always. When m_isThermal, a second (advective-only) block for the
  // energy balance is appended after it - see setUpDflux_dApertureMatrix and
  // getDerivativeFluxResidual_dApertureEnergyOffsets(), scattered further down (Part 3).
  using namespace contact;

  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager const & nodeManager = mesh.getNodeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();

  // assembleSystem has already brought this matrix to the host after the flux
  // assembly; the traversal below only reads it.
  CRSMatrixView< real64 const, localIndex const > const &
  dFluxResidual_dNormalJump = this->getDerivativeFluxResidual_dNormalJump().toViewConst();
  auto const derivativeOffsetIt = m_derivativeFluxResidual_dApertureOffsets.find( meshName );
  GEOS_ERROR_IF( derivativeOffsetIt == m_derivativeFluxResidual_dApertureOffsets.end(),
                 GEOS_FMT( "No dR/dAperture row offset is available for mesh body '{}'", meshName ) );
  localIndex const derivativeOffset = derivativeOffsetIt->second;

  string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
  string const & bubbleDofKey = dofManager.getKey( totalBubbleDisplacement::key() );
  string const & presDofKey = dofManager.getKey( this->getFlowDofKey() );

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

  array2d< real64 > dAperturedU( numElems, maxNumUdofs );
  array2d< real64 > dAperturedB( numElems, numBdofs );
  dAperturedU.zero();
  dAperturedB.zero();
  dAperturedU.move( parallelDeviceMemorySpace, true );
  dAperturedB.move( parallelDeviceMemorySpace, true );
  arrayView2d< real64 > const dAperturedUView = dAperturedU.toView();
  arrayView2d< real64 > const dAperturedBView = dAperturedB.toView();

  // Launch the ComputeApertureDerivatives kernel to fill dAperturedU and dAperturedB
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

    // ==== Part 3: Energy-balance flux derivative (advective contribution only) ====
    // Mirrors Parts 1/2 above but reads the energy block appended after all mass rows in
    // dFluxResidual_dNormalJump, reuses the same dAperturedU/dAperturedB chain-rule factors
    // (aperture -> nodal/bubble DOF is equation-agnostic, so the pre-computed derivatives from
    // ComputeApertureDerivativesFactory above apply unchanged), and scatters into the
    // temperature/energy residual row (packed right after pressure) instead of the mass row.
    // The conductive term's aperture sensitivity is not modeled - see
    // setUpDflux_dApertureMatrix and ThermalSinglePhasePoromechanicsConformingFractures.hpp.
    if( this->m_isThermal )
    {
      stdMap< string, localIndex > const & energyOffsets = this->getDerivativeFluxResidual_dApertureEnergyOffsets();
      auto const energyOffsetIt = energyOffsets.find( meshName );
      if( energyOffsetIt != energyOffsets.end() )
      {
        localIndex const energyOffset = energyOffsetIt->second;
        globalIndex elemDOFEnergy[1];
        elemDOFEnergy[0] = presDofNumber[kfe] + 1; // temperature/energy dof, packed right after pressure
        localIndex const localRowEnergy = LvArray::integerConversion< localIndex >( elemDOFEnergy[0] - rankOffset );

        localIndex const numEnergyColumns = dFluxResidual_dNormalJump.numNonZeros( energyOffset + kfe );
        arraySlice1d< localIndex const > const & energyColumns = dFluxResidual_dNormalJump.getColumns( energyOffset + kfe );
        arraySlice1d< real64 const > const & energyValues = dFluxResidual_dNormalJump.getEntries( energyOffset + kfe );

        // Nodal (Apu-energy) contribution
        for( localIndex kfe1 = 0; kfe1 < numEnergyColumns; ++kfe1 )
        {
          real64 const dREnergy_dAper = energyValues[kfe1];
          localIndex const kfe2 = energyColumns[kfe1] - derivativeOffset;

          bool const isOpen = ( fractureState[kfe2] == FractureState::Open );
          if( !isOpen && !isFractureOpen )
            continue;

          localIndex const kf0_2 = elemsToFaces[kfe2][0];
          localIndex const numNodesPerFace2 = faceToNodeMap.sizeOfArray( kf0_2 );
          localIndex const numUdofs2 = 2 * 3 * numNodesPerFace2;

          globalIndex nodeDOF2Energy[maxNumUdofs];
          for( localIndex kf = 0; kf < 2; ++kf )
          {
            for( localIndex a = 0; a < numNodesPerFace2; ++a )
            {
              for( localIndex i = 0; i < 3; ++i )
              {
                nodeDOF2Energy[kf * 3 * numNodesPerFace2 + 3 * a + i] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe2][kf], a )]
                                                                        + LvArray::integerConversion< globalIndex >( i );
              }
            }
          }

          stackArray1d< real64, maxNumUdofs > dRdUEnergy( maxNumUdofs );
          for( localIndex j = 0; j < numUdofs2; ++j )
          {
            dRdUEnergy( j ) = dREnergy_dAper * dAperturedU( kfe2, j );
          }

          if( localRowEnergy >= 0 && localRowEnergy < localMatrix.numRows() )
          {
            localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRowEnergy,
                                                                      nodeDOF2Energy,
                                                                      dRdUEnergy.data(),
                                                                      numUdofs2 );
          }
        }

        // Bubble (Apb-energy) contribution
        for( localIndex kfe1 = 0; kfe1 < numEnergyColumns; ++kfe1 )
        {
          real64 const dREnergy_dAper = energyValues[kfe1];
          localIndex const kfe2 = energyColumns[kfe1] - derivativeOffset;

          bool const isOpen = ( fractureState[kfe2] == FractureState::Open );
          if( !isOpen && !isFractureOpen )
            continue;

          globalIndex bubbleDOF2Energy[numBdofs];
          for( localIndex kf = 0; kf < 2; ++kf )
          {
            localIndex const faceIndex = elemsToFaces[kfe2][kf];
            for( localIndex i = 0; i < 3; ++i )
            {
              bubbleDOF2Energy[kf * 3 + i] = bubbleDofNumber[faceIndex] + LvArray::integerConversion< globalIndex >( i );
            }
          }

          real64 dRdBEnergy[numBdofs];
          for( localIndex j = 0; j < numBdofs; ++j )
          {
            dRdBEnergy[j] = dREnergy_dAper * dAperturedB( kfe2, j );
          }

          if( localRowEnergy >= 0 && localRowEnergy < localMatrix.numRows() )
          {
            localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRowEnergy,
                                                                      bubbleDOF2Energy,
                                                                      dRdBEnergy,
                                                                      numBdofs );
          }
        }
      }
    }
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

  string const flowDofKey = dofManager.getKey( this->getFlowDofKey() );
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
