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
 * @file PoromechanicsConformingFractures.hpp
 *
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSCONFORMINGFRACTURES_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSCONFORMINGFRACTURES_HPP_

#include "common/logger/Logger.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/multiphysics/SinglePhaseReservoirAndWells.hpp"
#include "physicsSolvers/solidMechanics/contact/SolidMechanicsLagrangeContact.hpp"
#include "physicsSolvers/solidMechanics/contact/ContactFields.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/solidMechanics/contact/ContactFields.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanicsFractures.hpp"
#include "physicsSolvers/multiphysics/PoromechanicsSolver.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/contact/HydraulicApertureBase.hpp"
#include "constitutive/contact/HydraulicApertureRelationSelector.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "common/DataTypes.hpp"
#include "mesh/DomainPartition.hpp"
#include "linearAlgebra/utilities/SparsityPatternUtilities.hpp"

namespace geos
{

template< template< typename, typename > class POROMECHANICS_BASE, typename FLOW_SOLVER = SinglePhaseBase , typename CONTACT_SOLVER = SolidMechanicsLagrangeContact >
class PoromechanicsConformingFractures : public POROMECHANICS_BASE< FLOW_SOLVER, CONTACT_SOLVER >
{
public:
  using Base = POROMECHANICS_BASE< FLOW_SOLVER, CONTACT_SOLVER >;

  PoromechanicsConformingFractures( const string & name,
                                    dataRepository::Group * const parent )
    : Base( name, parent )
  {}

  virtual void setupCoupling( DomainPartition const & domain,
                              DofManager & dofManager ) const override
  {
    /// We need to add 2 coupling terms:
    // 1. Poromechanical coupling in the bulk (p<->disp)
    Base::setupCoupling( domain, dofManager );

    if constexpr (CONTACT_SOLVER::hasContactStabilization) {
        // 2. Pressure - bubble displacement coupling in the fracture
        dofManager.addCoupling( this->getFlowDofKey(),
                          fields::contact::totalBubbleDisplacement::key(),
                          DofManager::Connector::Elem );
    }
    else {
     // 2. Traction - pressure coupling in the fracture
    dofManager.addCoupling( this->getFlowDofKey(),
                            fields::contact::traction::key(),
                            DofManager::Connector::Elem );
    }
  }

  virtual void setSparsityPattern( DomainPartition & domain,
                                   DofManager & dofManager,
                                   CRSMatrix< real64, globalIndex > & localMatrix,
                                   SparsityPattern< globalIndex > & pattern ) override
  {
    // start with the flow solver sparsity pattern (it could be reservoir + wells)
    SparsityPattern< globalIndex > patternOriginal;
    this->flowSolver()->setSparsityPattern( domain, dofManager, localMatrix, patternOriginal );

    SparsityPattern< globalIndex > mechanicsPattern;
    if constexpr (CONTACT_SOLVER::hasContactStabilization) {
    this->solidMechanicsSolver()->setSparsityPattern( domain, dofManager, localMatrix, mechanicsPattern );
    }

    // Get the original row lengths (diagonal blocks only)
    array1d< localIndex > rowLengths( patternOriginal.numRows());
    for( localIndex localRow = 0; localRow < patternOriginal.numRows(); ++localRow )
    {
      rowLengths[localRow] = patternOriginal.numNonZeros( localRow );
      if constexpr ( CONTACT_SOLVER::hasContactStabilization )
      {
        rowLengths[localRow] += mechanicsPattern.numNonZeros( localRow );   // simple sum, see note below
      }
    }

    // Add the number of nonzeros induced by coupling
    //displacement (and opt. bubble) to flow coupling
    addTransmissibilityCouplingNNZ( domain, dofManager, rowLengths.toView());
    if constexpr (CONTACT_SOLVER::hasContactStabilization) {
      //bubble to displacement coupling
      addPressureForceCouplingNNZ( domain, dofManager, rowLengths.toView() );
      addMatrixPressureBubbleCouplingNNZ( domain, dofManager, rowLengths.toView() );//TODO should be brought by CONTACT::STABILIZATION
    }


    // Create a new pattern with enough capacity for coupled matrix
    pattern.resizeFromRowCapacities< parallelHostPolicy >( patternOriginal.numRows(),
                                                           patternOriginal.numColumns(),
                                                           rowLengths.data());

    // Copy the original nonzeros
    appendSparsityPattern( pattern, patternOriginal );
    //ALM appendSparsityPattern( pattern, flowPattern );
    if constexpr (CONTACT_SOLVER::hasContactStabilization) {
       appendSparsityPattern( pattern, mechanicsPattern );
    }

    // Add the nonzeros from coupling
    //displacement (and opt. bubble) to flow coupling
    addTransmissibilityCouplingPattern( domain, dofManager, pattern.toView());
    if constexpr (CONTACT_SOLVER::hasContactStabilization) {
      addPressureForceCouplingPattern( domain, dofManager, pattern.toView() );
      addMatrixPressureBubbleCouplingPattern( domain, dofManager, pattern.toView() );
    }

    setUpDflux_dApertureMatrix( domain );
  }

  //Stabilization specific
  //TODO see refacto with below
void addPressureForceCouplingNNZ( DomainPartition const & domain,
                                  DofManager const & dofManager,
                                  arrayView1d< localIndex > const & rowLengths ) const
{
  GEOS_MARK_FUNCTION;

  integer const numComp = this->flowSolver()->numFluidComponents();

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel const & mesh,
                                                                      string_array const & )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const bubbleDofKey = dofManager.getKey( fields::contact::totalBubbleDisplacement::key() );

    arrayView1d< globalIndex const > const &
    bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );

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
      // For bubble DOFs: add 1 pressure column per bubble DOF row
      for( localIndex kf = 0; kf < 2; ++kf )
      {
        localIndex const faceIndex = elem2dToFaces[kfe][kf];
        for( localIndex i = 0; i < 3; ++i )
        {
          globalIndex const rowNumber = bubbleDofNumber[faceIndex] + i - rankOffset;
          if( rowNumber >= 0 && rowNumber < rowLengths.size() )
          {
            rowLengths[rowNumber] += numComp;  // One pressure column
          }
        }
      }
    } );
  } );
}

void addPressureForceCouplingPattern( DomainPartition const & domain,
                                 DofManager const & dofManager,
                                 SparsityPatternView< globalIndex > const & pattern ) const
{
  GEOS_MARK_FUNCTION;

  integer const numComp = this->flowSolver()->numFluidComponents();

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel const & mesh,
                                                                      string_array const & )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const bubbleDofKey = dofManager.getKey( fields::contact::totalBubbleDisplacement::key() );
    string const flowDofKey = dofManager.getKey( this->getFlowDofKey() );

    arrayView1d< globalIndex const > const &
    bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );

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

      // For bubble DOFs
      for( localIndex kf = 0; kf < 2; ++kf )
      {
        localIndex const faceIndex = elem2dToFaces[kfe][kf];
        for( localIndex i = 0; i < 3; ++i )
        {
          globalIndex const rowIndex = bubbleDofNumber[faceIndex] + i - rankOffset;
          if( rowIndex >= 0 && rowIndex < pattern.numRows() )
          {
            for(integer ic = 0; ic < numComp; ++ic)
              pattern.insertNonZero( rowIndex, pressureColIndex + ic );
          }
        }
      }
    } );
  } );
}

void addMatrixPressureBubbleCouplingNNZ( DomainPartition const & domain,
                                    DofManager const & dofManager,
                                    arrayView1d< localIndex > const & rowLengths ) const
{
  GEOS_MARK_FUNCTION;

  integer const numComp = this->flowSolver()->numFluidComponents();
  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel const & mesh,
                                                                      string_array const & regionNames )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const bubbleDofKey = dofManager.getKey( fields::contact::totalBubbleDisplacement::key() );
    string const flowDofKey = dofManager.getKey( this->getFlowDofKey() );
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
            rowLengths[rowNumber] += numComp;  // One pressure DOF from matrix cell
          }
        }

        // (pressure_row, bubble_col): the matrix cell pressure couples to its 3 bubble DOFs (A_pb)
        globalIndex const pRow = pressureDofNumber[k] - rankOffset;
        if( pRow >= 0 && pRow < rowLengths.size() )
        {
          for( integer ic = 0; ic<numComp; ++ic)
            rowLengths[pRow + ic] += 3;  // Three bubble DOFs
        }
      } );
    } );
  } );
}

void addMatrixPressureBubbleCouplingPattern( DomainPartition const & domain,
                                        DofManager const & dofManager,
                                        SparsityPatternView< globalIndex > const & pattern ) const
{
  GEOS_MARK_FUNCTION;

  integer const numComp = this->flowSolver()->numFluidComponents();

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel const & mesh,
                                                                      string_array const & regionNames )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const bubbleDofKey = dofManager.getKey( fields::contact::totalBubbleDisplacement::key() );
    string const flowDofKey = dofManager.getKey( this->getFlowDofKey() );

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
            for( integer ic = 0; ic < numComp; ++ic )
              pattern.insertNonZero( rowIndex, pressureColIndex + ic);
          }
        }

        // (pressure_row, bubble_col) : A_pb -- transpose location
        globalIndex const pRow = pressureDofNumber[k] - rankOffset;
        if( pRow >= 0 && pRow < pattern.numRows() )
        {
          for( localIndex i = 0; i < 3; ++i )
          {
            for( integer ic = 0; ic < numComp; ++ic )
              pattern.insertNonZero( pRow + ic, bubbleDofNumber[faceIndex] + i );
          }
        }
      } );
    } );
  } );
}


  virtual void assembleSystem( real64 const time_n,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override
  {

    GEOS_MARK_FUNCTION;

    this->solidMechanicsSolver()->synchronizeFractureState( domain );

    // The flux assembly accumulates into this matrix. Clear it before every
    // Newton assembly and make the host copy explicit before the host-side
    // coupling kernels consume it.

    if( !m_derivativeFluxResidual_dAperture )
    {
      setUpDflux_dApertureMatrix( domain );
    }
    m_derivativeFluxResidual_dAperture->move( parallelDeviceMemorySpace, false );
    m_derivativeFluxResidual_dAperture->zero();

    assembleElementBasedContributions( time_n,
                                       dt,
                                       domain,
                                       dofManager,
                                       localMatrix,
                                       localRhs );

    // Assemble fluxes 3D/2D and get dFluidResidualDAperture
    this->flowSolver()->assembleHydrofracFluxTerms( time_n,
                                                    dt,
                                                    domain,
                                                    dofManager,
                                                    localMatrix,
                                                    localRhs,
                                                    getDerivativeFluxResidual_dNormalJump(),
                                                    nullptr );

    m_derivativeFluxResidual_dAperture->move( hostMemorySpace, false );

    // This step must occur after the fluxes are assembled because that's when DerivativeFluxResidual_dAperture is filled.
    assembleCouplingTerms( time_n,
                           dt,
                           domain,
                           dofManager,
                           localMatrix,
                           localRhs );
  if constexpr ( std::is_same_v< FLOW_SOLVER, SinglePhaseReservoirAndWells<> >  )
  {
    this->flowSolver()->wellSolver()->assembleSystem( time_n, dt, domain, dofManager, localMatrix, localRhs );
    this->flowSolver()->assembleCouplingTerms( time_n, dt, domain, dofManager, localMatrix, localRhs );
  }

  }

  virtual void updateState( DomainPartition & domain ) override
  {
    GEOS_MARK_FUNCTION;

    // call base poromechanics update
    Base::updateState( domain );
    // need to call solid mechanics update separately to compute face displacement jump
    this->solidMechanicsSolver()->updateState( domain );

    // remove the contribution of the hydraulic aperture from the stencil weights
    this->flowSolver()->prepareStencilWeights( domain );

    updateHydraulicApertureAndFracturePermeability( domain );

    // update the stencil weights using the updated hydraulic aperture
    this->flowSolver()->updateStencilWeights( domain );
  }

protected:

  /**
   * @Brief add the nnz induced by the flux-aperture coupling
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param rowLenghts the nnz in each row
   */
  void addTransmissibilityCouplingNNZ( DomainPartition const & domain,
                                       DofManager const & dofManager,
                                       arrayView1d< localIndex > const & rowLengths ) const
  {
    GEOS_MARK_FUNCTION;


    integer const numComp = this->flowSolver()->numFluidComponents();

    this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &, //  meshBodyName,
                                                                        MeshLevel const & mesh,
                                                                        string_array const & ) // regionNames
    {
      ElementRegionManager const & elemManager = mesh.getElemManager();

      string const flowDofKey = dofManager.getKey( this->getFlowDofKey() );

      globalIndex const rankOffset = dofManager.rankOffset();

      NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
      FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
      
      //TODO (jafranc) - remove once ALM-bubble is frame as a stab method - tmp runtime is fine as it is tmp
      if(this->solidMechanicsSolver()->hasStabilization())//why is this not done in SolidMech ?
      {

      FluxApproximationBase const & stabilizationMethod = fvManager.getFluxApproximation( this->solidMechanicsSolver()->getStabilizationName() );

      stabilizationMethod.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
      {
        for( localIndex iconn=0; iconn<stencil.size(); ++iconn )
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

          for( localIndex k0=0; k0<numFluxElems; ++k0 )
          {
            globalIndex const activeFlowDOF = faceElementDofNumber[sei[iconn][k0]];
            globalIndex const rowNumber = activeFlowDOF - rankOffset;

            if( rowNumber >= 0 && rowNumber < rowLengths.size() )
            {
              for( localIndex k1=0; k1<numFluxElems; ++k1 )
              {
                // The coupling with the nodal displacements of the cell itself has already been added by the dofManager
                // so we only add the coupling with the nodal displacements of the neighbors.
                if( k1 != k0 )
                {
                  localIndex const numNodesPerElement = elemsToNodes[sei[iconn][k1]].size();
                  for( integer ic = 0; ic < numComp; ic++ )
                  {
                    rowLengths[rowNumber + ic] += 3*numNodesPerElement;
                  }
                }
              }
            }
          }
        }
      } );
      }
      //to decide -- reduce duplication -- can we have both ?
      if constexpr (CONTACT_SOLVER::hasContactStabilization) {

      FluxApproximationBase const & fvMethod = fvManager.getFluxApproximation( this->flowSolver()->getDiscretizationName() );

      fvMethod.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
      {
        for( localIndex iconn=0; iconn<stencil.size(); ++iconn )
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

          for( localIndex k0=0; k0<numFluxElems; ++k0 )
          {
            globalIndex const activeFlowDOF = faceElementDofNumber[sei[iconn][k0]];
            globalIndex const rowNumber = activeFlowDOF - rankOffset;

            if( rowNumber >= 0 && rowNumber < rowLengths.size() )
            {
              for( localIndex k1=0; k1<numFluxElems; ++k1 )
              {
                // The coupling with the nodal displacements of the cell itself has already been added by the dofManager
                // so we only add the coupling with the nodal displacements of the neighbors.
                if( k1 != k0 )
                {
                  for( integer ic = 0; ic < numComp; ic++ )
                  {
                      localIndex const numNodesPerElement = elemsToNodes[sei[iconn][k1]].size();
                      rowLengths[rowNumber + ic] += 3*numNodesPerElement;
                      rowLengths[rowNumber + ic] += 6;
                  }
                }
              }
            }
          }
        }
      } );
      
      }
    } );//end forAll
  }

  /**
   * @brief Set up the Dflux_dApertureMatrix object
   *
   * @param domain
   * @param dofManager
   * @param localMatrix
   */
  void addTransmissibilityCouplingPattern( DomainPartition const & domain,
                                           DofManager const & dofManager,
                                           SparsityPatternView< globalIndex > const & pattern ) const
  {
    GEOS_MARK_FUNCTION;

    this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                        MeshLevel const & mesh,
                                                                        string_array const & )
    {

    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();
    ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();

    addTransmissibilityCouplingPattern( domain, mesh, dofManager, pattern,
        nodeManager.getReference< globalIndex_array >( dofManager.getKey( fields::solidMechanics::totalDisplacement::key() ) ),
        [&faceToNodeMap](localIndex const& faceIndex, localIndex const& a){ return faceToNodeMap(faceIndex,a); },
        [&faceToNodeMap](localIndex const& faceIndex){ return faceToNodeMap.sizeOfArray(faceIndex); } );

    if constexpr (CONTACT_SOLVER::hasContactStabilization)
    {
      addTransmissibilityCouplingPattern( domain, mesh, dofManager, pattern,
          faceManager.getReference< globalIndex_array >( dofManager.getKey( fields::contact::totalBubbleDisplacement::key() ) ),
          [](localIndex const & faceIndex, localIndex const& GEOS_UNUSED_PARAM(a)){ return faceIndex; },
          [](localIndex const & GEOS_UNUSED_PARAM(faceIndex)){ return 1; } );
    }
     
    } );
  }

  template< typename NODE_INDEX_MAP, typename NNODE_PER_FACE >
  void addTransmissibilityCouplingPattern( DomainPartition const & domain,
                                           MeshLevel const & mesh,
                                           DofManager const & dofManager,
                                           SparsityPatternView< globalIndex > const & pattern,
                                           arrayView1d< globalIndex const > const & coupledDisplacementDofNumber,
                                           NODE_INDEX_MAP && dofIndirectionCb,
                                           NNODE_PER_FACE && numNodesPerFace
                                            ) const
  {

      // Get the finite volume method used to compute the stabilization
      NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
      FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
      FluxApproximationBase const & fvDiscretization = fvManager.getFluxApproximation( this->flowSolver()->getDiscretizationName() );

      SurfaceElementRegion const & fractureRegion =
        mesh.getElemManager().getRegion< SurfaceElementRegion >( this->solidMechanicsSolver()->getUniqueFractureRegionName() );
      FaceElementSubRegion const & fractureSubRegion =
        fractureRegion.getUniqueSubRegion< FaceElementSubRegion >();

      GEOS_ERROR_IF( !fractureSubRegion.hasWrapper( fields::flow::pressure::key() ),
                     "The fracture subregion must contain pressure field.", this->getDataContext() );

      arrayView2d< localIndex const > const elem2dToFaces = fractureSubRegion.faceList().toViewConst();

      arrayView1d< globalIndex const > const &
      flowDofNumber = fractureSubRegion.getReference< globalIndex_array >( dofManager.getKey(this->getFlowDofKey()) );
 
        //debug
        // ArrayOfArraysView< localIndex const > const elemsToNodes = fractureSubRegion.nodeList().toViewConst();
      
      globalIndex const rankOffset = dofManager.rankOffset();

      fvDiscretization.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
      {
        forAll< serialPolicy >( stencil.size(), [=,this] ( localIndex const iconn )
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
              localIndex const rowIndex = flowDofNumber[sei[iconn][1-kf]] - rankOffset;

              if( rowIndex >= 0 && rowIndex < pattern.numRows() )
              {

                // Get fracture, face and region/subregion/element indices (for elements on both sides)
                localIndex const fractureIndex = sei[iconn][kf];

                // Loop over the two sides of each fracture element
                for( localIndex kf1 = 0; kf1 < 2; ++kf1 )
                {
                  localIndex const faceIndex = elem2dToFaces[fractureIndex][kf1];

                  // Save the list of DOF associated with nodes
                  for( localIndex a=0; a<numNodesPerFace(elem2dToFaces[fractureIndex][0]); ++a )
                  {
                    for( localIndex i = 0; i < 3; ++i )
                    {
                      globalIndex const colIndex = coupledDisplacementDofNumber[dofIndirectionCb( faceIndex, a )] + LvArray::integerConversion< globalIndex >( i );
                      for( integer ic = 0; ic < this->flowSolver()->numFluidComponents(); ic++ )
                      {
                        pattern.insertNonZero( rowIndex + ic, colIndex );
                      }
                    }
                  }
                }
              }
            }
          }
        } );
      } );

  }


  /**
   * @brief Set up the Dflux_dApertureMatrix object
   *
   * @param domain
   */
  void setUpDflux_dApertureMatrix( DomainPartition & domain )
  {
    integer const numComp = this->flowSolver()->numFluidComponents();
    localIndex numCols = 0.;//number of outerloop pass (not considering innermost component loop)
    
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
    rowCapacities.resize( rowOffset + fractureSubRegion.size() * numComp, 0 );
    numCols += fractureSubRegion.size();

    fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
    {
      for( localIndex iconn = 0; iconn < stencil.size(); ++iconn )
      {
        localIndex const numFluxElems = stencil.stencilSize( iconn );
        typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();
        for( localIndex k0 = 0; k0 < numFluxElems; ++k0 )
        {
          for( integer ic = 0; ic < numComp; ++ic){
            localIndex const row = rowOffset + sei[iconn][k0] * numComp;
            GEOS_ERROR_IF_GE_MSG( row,
                                  LvArray::integerConversion< localIndex >( rowCapacities.size() ),
                                  "Surface stencil index exceeds the fracture derivative matrix size." );
            rowCapacities[ row + ic ] += numFluxElems;
          }
        }
      }
    } );
  } );

  //write real data in structure
  std::unique_ptr< CRSMatrix< real64, localIndex > > & derivativeFluxResidual_dAperture = getRefDerivativeFluxResidual_dAperture();
  localIndex const numRows = rowCapacities.size();
  derivativeFluxResidual_dAperture = std::make_unique< CRSMatrix< real64, localIndex > >( numRows, numCols );
  derivativeFluxResidual_dAperture->setName( this->getName() + "/derivativeFluxResidual_dAperture" );
  if( numRows > 0 )
  {
    derivativeFluxResidual_dAperture->resizeFromRowCapacities< parallelHostPolicy >( numRows,
                                                                                     numCols,
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
          for(integer ic = 0 ; ic<numComp; ++ic){
          localIndex const row = rowOffset + sei[iconn][k0] * numComp + ic;
          GEOS_ERROR_IF_GE_MSG( row, numRows, "Surface stencil index exceeds the fracture derivative matrix size." );
          for( localIndex k1 = 0; k1 < numFluxElems; ++k1 )
          {
            derivativeFluxResidual_dAperture->insertNonZero( row,
                                                             rowOffset/numComp + sei[iconn][k1], //as component-independent indexing
                                                             0.0 );
          }
        }
        }
      }
    } );
  } );
  }

  void assembleElementBasedContributions( real64 const time_n,
                                          real64 const dt,
                                          DomainPartition & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
  {
    GEOS_UNUSED_VAR( time_n, dt );

    /// 3. assemble Force Residual w.r.t. pressure and Flow mass residual w.r.t. displacement

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

    this->solidMechanicsSolver()->assembleContact( time_n, dt, domain, dofManager, localMatrix, localRhs );
  }

  virtual void assembleCouplingTerms( real64 const time_n,
                                      real64 const dt,
                                      DomainPartition const & domain,
                                      DofManager const & dofManager,
                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                      arrayView1d< real64 > const & localRhs ) override
  {
    GEOS_UNUSED_VAR( time_n, dt );
    // These 2 steps need to occur after the fluxes are assembled because that's when DerivativeFluxResidual_dAperture is filled.
    this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                        MeshLevel const & mesh,
                                                                        string_array const & regionNames )
    {
      /// 3. assemble Force Residual w.r.t. pressure and Flow mass residual w.r.t. displacement
      assembleForceResidualDerivativeWrtPressure( meshName, mesh, regionNames, dofManager, localMatrix, localRhs );
      assembleFluidMassResidualDerivativeWrtDisplacement( meshName, mesh, regionNames, dofManager, localMatrix, localRhs );
    } );

    //if hasStabilization via bubble - Apb
    if constexpr (CONTACT_SOLVER::hasContactStabilization)
       assembleMatrixPressureBubbleContribution( dt, const_cast< DomainPartition & >( domain ), dofManager, localMatrix, localRhs );

  }

  virtual void assembleForceResidualDerivativeWrtPressure( string const & GEOS_UNUSED_PARAM(meshName),
                                                   MeshLevel const & mesh,
                                                   string_array const & regionNames,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs ) = 0;

  virtual void assembleFluidMassResidualDerivativeWrtDisplacement( string const& meshName,
                                                                   MeshLevel const & mesh,
                                                                   string_array const & regionNames,
                                                                   DofManager const & dofManager,
                                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                   arrayView1d< real64 > const & localRhs ) = 0;

  /**
   * @Brief assemble the contribution of matrix cell pressure on bubble DOFs
   * with full Jacobian for fully-implicit coupling.
   * @param dt the time step size
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the local system matrix
   * @param localRhs the local system right-hand side vector
   */
  virtual void assembleMatrixPressureBubbleContribution( real64 const GEOS_UNUSED_PARAM(dt),
                                                 DomainPartition & GEOS_UNUSED_PARAM(domain),
                                                 DofManager const & GEOS_UNUSED_PARAM(dofManager),
                                                 CRSMatrixView< real64, globalIndex const > const & GEOS_UNUSED_PARAM(localMatrix),
                                                 arrayView1d< real64 > const & GEOS_UNUSED_PARAM(localRhs) ) 
  { GEOS_WARNING("Should override"); };

  virtual void mapSolutionBetweenSolvers( DomainPartition & domain, integer const solverType ) override
  {
    GEOS_MARK_FUNCTION;

    /// After the solid mechanics solver
    if( solverType == static_cast< integer >( Base::SolverType::SolidMechanics )
        && !this->m_performStressInitialization ) // do not update during poromechanics initialization
    {
      // remove the contribution of the hydraulic aperture from the stencil weights
      this->flowSolver()->prepareStencilWeights( domain );

      updateHydraulicApertureAndFracturePermeability( domain );

      // update the stencil weights using the updated hydraulic aperture
      this->flowSolver()->updateStencilWeights( domain );
    }

    Base::mapSolutionBetweenSolvers( domain, solverType );
  }

  void updateHydraulicApertureAndFracturePermeability( DomainPartition & domain )
  {
    using namespace constitutive;

    this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                        MeshLevel & mesh,
                                                                        string_array const & regionNames )
    {
      ElementRegionManager & elemManager = mesh.getElemManager();

      elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                                [&]( localIndex const,
                                                                     FaceElementSubRegion & subRegion )
      {
        arrayView2d< real64 const > const dispJump           = subRegion.getField< fields::contact::dispJump >();
        arrayView1d< real64 const > const area               = subRegion.getElementArea();
        arrayView1d< real64 const > const volume             = subRegion.getElementVolume();
        arrayView2d< real64 const > const fractureTraction   = subRegion.getField< fields::contact::traction >();
        arrayView1d< real64 const > const pressure           = subRegion.getField< fields::flow::pressure >();
        arrayView1d< real64 const > const oldHydraulicAperture = subRegion.getField< fields::flow::aperture0 >();

        arrayView1d< real64 > const aperture                 = subRegion.getElementAperture();
        arrayView1d< real64 > const hydraulicAperture        = subRegion.getField< fields::flow::hydraulicAperture >();
        arrayView1d< real64 > const deltaVolume              = subRegion.getField< fields::flow::deltaVolume >();
        arrayView1d< integer > const & fractureState   = subRegion.getField< fields::contact::fractureState >();

        string const porousSolidName = subRegion.getReference< string >( FlowSolverBase::viewKeyStruct::solidNamesString() );
        CoupledSolidBase & porousSolid = subRegion.getConstitutiveModel< CoupledSolidBase >( porousSolidName );

        string const & hydraulicApertureRelationName = subRegion.template getReference< string >( viewKeyStruct::hydraulicApertureRelationNameString()  );
        HydraulicApertureBase const & hydraulicApertureModel = this->template getConstitutiveModel< HydraulicApertureBase >( subRegion, hydraulicApertureRelationName );

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
  }

  std::unique_ptr< CRSMatrix< real64, localIndex > > & getRefDerivativeFluxResidual_dAperture()
  {
    return m_derivativeFluxResidual_dAperture;
  }

  CRSMatrixView< real64, localIndex const > getDerivativeFluxResidual_dNormalJump()
  {
    return m_derivativeFluxResidual_dAperture->toViewConstSizes();
  }

  CRSMatrixView< real64 const, localIndex const > getDerivativeFluxResidual_dNormalJump() const
  {
    return m_derivativeFluxResidual_dAperture->toViewConst();
  }

  struct viewKeyStruct : public Base::viewKeyStruct
  {};

  static const localIndex m_maxFaceNodes = 11; // Maximum number of nodes on a contact face

  std::unique_ptr< CRSMatrix< real64, localIndex > > m_derivativeFluxResidual_dAperture;
  stdMap< string, localIndex > m_derivativeFluxResidual_dApertureOffsets;

};

}

#endif //GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSCONFORMINGFRACTURES_HPP_
