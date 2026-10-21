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

  Base::setSparsityPattern( domain, dofManager, localMatrix, pattern );
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
  // fracture element. When m_isThermal, a second (advective + conductive) block for the energy
  // balance is appended after it - see setUpDerivativeFluxResidual_dNormalJump and
  // getDerivativeFluxResidual_dNormalJumpEnergyOffset().
  using namespace contact;

  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager const & nodeManager = mesh.getNodeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();

  // assembleSystem has already brought this matrix to the host after the flux
  // assembly; the traversal below only reads it.
  CRSMatrixView< real64 const, localIndex const > const &
  dFluxResidual_dNormalJump = this->getDerivativeFluxResidual_dNormalJump().toViewConst();

  string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
  string const & bubbleDofKey = dofManager.getKey( totalBubbleDisplacement::key() );
  string const & presDofKey = dofManager.getKey( this->getFlowDofKey() );

  arrayView1d< globalIndex const > const &
  dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
  arrayView1d< globalIndex const > const &
  bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );
  globalIndex const rankOffset = dofManager.rankOffset();
  // energy-balance row (temperature DOF) of a fracture element, w.r.t. its mass-balance row (pressure DOF)
  integer const temperatureOffset = this->temperatureDofOffset();

  string const & fractureRegionName = this->solidMechanicsSolver()->getUniqueFractureRegionName();

  // Maximum DOF sizes
  constexpr localIndex maxNumNodesPerFace = m_maxFaceNodes;
  constexpr localIndex maxNumUdofs = 2 * 3 * maxNumNodesPerFace; // 66
  constexpr localIndex numBdofs = 6;

  // Get the fracture subRegion
  SurfaceElementRegion const & fractureRegion = elemManager.getRegion< SurfaceElementRegion >( fractureRegionName );
  FaceElementSubRegion const & subRegion = fractureRegion.getUniqueSubRegion< FaceElementSubRegion >();

  localIndex const numElems = subRegion.size();

  array2d< real64 > dNormalJump_dU( numElems, maxNumUdofs );
  array2d< real64 > dNormalJump_dB( numElems, numBdofs );
  dNormalJump_dU.zero();
  dNormalJump_dB.zero();
  dNormalJump_dU.move( parallelDeviceMemorySpace, true );
  dNormalJump_dB.move( parallelDeviceMemorySpace, true );
  arrayView2d< real64 > const dNormalJump_dUView = dNormalJump_dU.toView();
  arrayView2d< real64 > const dNormalJump_dBView = dNormalJump_dB.toView();

  // Launch the ComputeNormalJumpDerivatives kernel to fill dNormalJump_dU and dNormalJump_dB
  this->solidMechanicsSolver()->forFiniteElementOnFractureSubRegions( meshName,
                                                                      [&, dNormalJump_dUView, dNormalJump_dBView] ( string const &,
                                                                                                              finiteElement::FiniteElementBase const & subRegionFE,
                                                                                                              arrayView1d< localIndex const > const & faceElementList )
  {
    poromechanicsALMKernels::ComputeNormalJumpDerivativesFactory
    kernelFactory( dispDofNumber,
                   bubbleDofNumber,
                   dofManager.rankOffset(),
                   localMatrix,
                   localRhs,
                   0.0,  // dt not used
                   faceElementList,
                   dNormalJump_dUView,
                   dNormalJump_dBView );

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
  dNormalJump_dU.move( hostMemorySpace );
  dNormalJump_dB.move( hostMemorySpace );

  // Now assemble using the pre-computed derivatives
  arrayView1d< globalIndex const > const & presDofNumber = subRegion.getReference< array1d< globalIndex > >( presDofKey );

  arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList().toViewConst();

  // dNormalJump_dU/dB = (1/area) * unitNormal^T * Atu/Atb are the derivatives of the (element-averaged) normal jump.
  // For accumulation: vol = hydraulicAperture * area, so
  //   dR_accum/du = (R/vol) * dVol/dNormalJump * dNormalJump_dU, with dVol/dNormalJump = area * dHydraulicAperture_dNormalJump
  // For flux: dR_flux/du = dR/dNormalJump * dNormalJump_dU (dR_dNormalJump already carries dHydraulicAperture_dNormalJump)
  arrayView1d< real64 const > const area = subRegion.getElementArea().toViewConst();
  arrayView1d< real64 const > const dHydraulicAperture_dNormalJump = subRegion.getField< flow::dHydraulicAperture_dNormalJump >();

  array1d< real64 > massPerVolume( numElems );
  array1d< real64 > energyPerVolume( numElems );
  this->computeAccumulationPerVolume( subRegion, massPerVolume.toView(), energyPerVolume.toView() );

  localIndex const energyOffset = this->getDerivativeFluxResidual_dNormalJumpEnergyOffset();
  GEOS_ERROR_IF( this->m_isThermal && energyOffset < 0,
                 GEOS_FMT( "{}: the energy-balance block of dR/dNormalJump was not built", this->getName() ) );

  // Adds dR_dNormalJump * [ dNormalJump_dU( kfe ), dNormalJump_dB( kfe ) ] to localRow, on the nodal and the bubble
  // displacement columns of fracture element kfe (both faces)
  auto addNormalJumpDerivative = [&]( localIndex const localRow,
                                      localIndex const kfe,
                                      real64 const dR_dNormalJump )
  {
    localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );
    localIndex const numUdofs = 2 * 3 * numNodesPerFace;

    globalIndex nodeDofs[maxNumUdofs];
    real64 dRdU[maxNumUdofs];
    globalIndex bubbleDofs[numBdofs];
    real64 dRdB[numBdofs];
    for( localIndex kf = 0; kf < 2; ++kf )
    {
      localIndex const faceIndex = elemsToFaces[kfe][kf];
      for( localIndex i = 0; i < 3; ++i )
      {
        for( localIndex a = 0; a < numNodesPerFace; ++a )
        {
          localIndex const j = kf * 3 * numNodesPerFace + 3 * a + i;
          nodeDofs[j] = dispDofNumber[faceToNodeMap( faceIndex, a )] + LvArray::integerConversion< globalIndex >( i );
          dRdU[j] = dR_dNormalJump * dNormalJump_dU( kfe, j );
        }
        bubbleDofs[kf * 3 + i] = bubbleDofNumber[faceIndex] + LvArray::integerConversion< globalIndex >( i );
        dRdB[kf * 3 + i] = dR_dNormalJump * dNormalJump_dB( kfe, kf * 3 + i );
      }
    }

    localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow, nodeDofs, dRdU, numUdofs );
    localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow, bubbleDofs, dRdB, numBdofs );
  };

  forAll< serialPolicy >( numElems, [&]( localIndex const kfe )
  {
    // No accumulation coupling when the hydraulic aperture does not follow the jump (e.g. closed Barton-Bandis).
    // Same test as the Lagrange-multiplier solver, on dHydraulicAperture_dNormalJump alone (not scaled by the area)
    bool const hasAccumulationCoupling = !isZero( dHydraulicAperture_dNormalJump[kfe] );
    // vol = hydraulicAperture * area
    real64 const dVol_dNormalJump = area[kfe] * dHydraulicAperture_dNormalJump[kfe];

    // The mass-balance row then, when thermal, the energy-balance row (temperature DOF): each has its own
    // accumulation per unit volume and its own row in dR/dNormalJump
    integer const numRows = this->m_isThermal ? 2 : 1;
    for( integer ir = 0; ir < numRows; ++ir )
    {
      bool const isEnergyRow = ( ir == 1 );

      localIndex const localRow =
        LvArray::integerConversion< localIndex >( presDofNumber[kfe] + ( isEnergyRow ? temperatureOffset : 0 ) - rankOffset );
      if( localRow < 0 || localRow >= localMatrix.numRows() )
      {
        continue;
      }

      // Accumulation: dR_accum/du = (R/vol) * dVol/dNormalJump * dNormalJump_dU, same for the bubble DOFs
      if( hasAccumulationCoupling )
      {
        real64 const accumulationPerVolume = isEnergyRow ? energyPerVolume[kfe] : massPerVolume[kfe];
        addNormalJumpDerivative( localRow, kfe, accumulationPerVolume * dVol_dNormalJump );
      }

      // Flux (mass, or advective + conductive energy): dR_flux/du = dR_flux/dNormalJump( kfe2 ) * dNormalJump_dU( kfe2 )
      // for every element kfe2 of the stencils of kfe. No fractureState gating: dR_flux/dNormalJump already vanishes
      // where the hydraulic aperture does not follow the jump (e.g. closed Barton-Bandis), and is the exact derivative
      // where it does (e.g. HydraulicApertureTable slope at gN < 0)
      localIndex const fluxRow = isEnergyRow ? energyOffset + kfe : kfe;
      localIndex const numColumns = dFluxResidual_dNormalJump.numNonZeros( fluxRow );
      arraySlice1d< localIndex const > const & columns = dFluxResidual_dNormalJump.getColumns( fluxRow );
      arraySlice1d< real64 const > const & values = dFluxResidual_dNormalJump.getEntries( fluxRow );
      for( localIndex icol = 0; icol < numColumns; ++icol )
      {
        addNormalJumpDerivative( localRow, columns[icol], values[icol] );
      }
    }
  } );

  // Fracture-matrix fluxes: their fracture-side half-transmissibility depends on the hydraulic aperture.
  // The matrix cell row gets dFlux/dNormalJump, the fracture row its opposite (both w.r.t. the jump of the fracture element)
  arrayView2d< real64 const > const dFractureMatrixFlux_dNormalJump = this->m_dFractureMatrixFlux_dNormalJump.toViewConst();
  integer const numRows = this->numCoupledFlowDofsPerElem();
  DomainPartition const & domain = this->template getGroupByPath< DomainPartition >( "/Problem/domain" );
  this->forFractureMatrixConnections( domain, mesh, dofManager, [&]( localIndex const iconn,
                                                                     globalIndex const cellDofNumber,
                                                                     localIndex const kfe )
  {
    for( integer ir = 0; ir < numRows; ++ir )
    {
      real64 const dFlux_dNormalJump = dFractureMatrixFlux_dNormalJump( iconn, ir );
      integer const rowOffset = ( ir == 0 ) ? 0 : temperatureOffset;

      localIndex const cellRow = LvArray::integerConversion< localIndex >( cellDofNumber + rowOffset - rankOffset );
      if( cellRow >= 0 && cellRow < localMatrix.numRows() )
      {
        addNormalJumpDerivative( cellRow, kfe, dFlux_dNormalJump );
      }

      localIndex const fractureRow = LvArray::integerConversion< localIndex >( presDofNumber[kfe] + rowOffset - rankOffset );
      if( fractureRow >= 0 && fractureRow < localMatrix.numRows() )
      {
        addNormalJumpDerivative( fractureRow, kfe, -dFlux_dNormalJump );
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

  // The kernel below adds the bubble volumetric strain increment on top of the stored porosity: this is only
  // correct if the bulk kernel has reset the porosity in this same assembly (assembleElementBasedContributions)
  GEOS_ERROR_IF( !this->m_isPorosityReset,
                   GEOS_FMT( "{}: the bubble contribution would be added twice to the porosity "
                             "(the bulk poromechanics kernel did not reset it since the last addition)", this->getName() ) );
  this->m_isPorosityReset = false;

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
                                                                                 FlowSolverBase::viewKeyStruct::fluidNamesString(),
                                                                                 this->m_isThermal,
                                                                                 this->m_performStressInitialization );

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
