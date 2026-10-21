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

template< template< typename, typename > class POROMECHANICS_BASE, typename FLOW_SOLVER = SinglePhaseBase, typename CONTACT_SOLVER = SolidMechanicsLagrangeContact >
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

    if constexpr (CONTACT_SOLVER::hasFaceBubbles) {
      // 2. Pressure - bubble displacement coupling in the fracture
      dofManager.addCoupling( this->getFlowDofKey(),
                              fields::contact::totalBubbleDisplacement::key(),
                              DofManager::Connector::Elem );
    }
    else
    {
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
    if constexpr (CONTACT_SOLVER::hasFaceBubbles) {
      this->solidMechanicsSolver()->setSparsityPattern( domain, dofManager, localMatrix, mechanicsPattern );
    }

    // Get the original row lengths (diagonal blocks only)
    array1d< localIndex > rowLengths( patternOriginal.numRows());
    for( localIndex localRow = 0; localRow < patternOriginal.numRows(); ++localRow )
    {
      rowLengths[localRow] = patternOriginal.numNonZeros( localRow );
      if constexpr ( CONTACT_SOLVER::hasFaceBubbles )
      {
        //compute overlap to have the exact size
        auto const flowCols = patternOriginal.getColumns( localRow );
        auto const mechCols = mechanicsPattern.getColumns( localRow );
        localIndex i = 0, j = 0, numUnion = 0;
        while( i < flowCols.size() && j < mechCols.size() )
        {
          if( flowCols[i] < mechCols[j] ) { ++i; }
          else if( mechCols[j] < flowCols[i] ) { ++j; }
          else { ++i; ++j; }
          ++numUnion;
        }
        rowLengths[localRow] = numUnion + ( flowCols.size() - i ) + ( mechCols.size() - j );
      }
    }

    // Add the number of nonzeros induced by coupling
    //displacement (and opt. bubble) to flow coupling
    addTransmissibilityCouplingNNZ( domain, dofManager, rowLengths.toView());
    if constexpr (CONTACT_SOLVER::hasFaceBubbles) {
      // fracture pressure to bubble, and matrix cell flow to/from bubble couplings
      addPressureForceCouplingNNZ( domain, dofManager, rowLengths.toView() );
      addMatrixTemperaturePressureBubbleCouplingNNZ( domain, dofManager, rowLengths.toView() );
    }


    // Create a new pattern with enough capacity for coupled matrix
    pattern.resizeFromRowCapacities< parallelHostPolicy >( patternOriginal.numRows(),
                                                           patternOriginal.numColumns(),
                                                           rowLengths.data());

    // Copy the original nonzeros
    appendSparsityPattern( pattern, patternOriginal );
    //ALM appendSparsityPattern( pattern, flowPattern );
    if constexpr (CONTACT_SOLVER::hasFaceBubbles) {
      appendSparsityPattern( pattern, mechanicsPattern );
    }

    // Add the nonzeros from coupling
    //displacement (and opt. bubble) to flow coupling
    addTransmissibilityCouplingPattern( domain, dofManager, pattern.toView());
    if constexpr (CONTACT_SOLVER::hasFaceBubbles) {
      addPressureForceCouplingPattern( domain, dofManager, pattern.toView() );
      addMatrixTemperaturePressureBubbleCouplingPattern( domain, dofManager, pattern.toView() );
    }

    setUpDerivativeFluxResidual_dNormalJump( domain );
  }

  // Face bubbles only
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
              rowLengths[rowNumber] += numComp; // One pressure column
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
              for( integer ic = 0; ic < numComp; ++ic )
                pattern.insertNonZero( rowIndex, pressureColIndex + ic );
            }
          }
        }
      } );
    } );
  }

  void addMatrixTemperaturePressureBubbleCouplingNNZ( DomainPartition const & domain,
                                                      DofManager const & dofManager,
                                                      arrayView1d< localIndex > const & rowLengths ) const
  {
    GEOS_MARK_FUNCTION;

    integer const numCoupledFlowDofs = numCoupledFlowDofsPerElem();
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

          // (bubble_row, flow_col): the flow DOFs of the matrix cell (pressure, then temperature when thermal)
          // for each of the 3 bubble DOFs (A_bp, A_bT)
          for( localIndex i = 0; i < 3; ++i )
          {
            globalIndex const rowNumber = bubbleDofNumber[faceIndex] + i - rankOffset;
            if( rowNumber >= 0 && rowNumber < rowLengths.size() )
            {
              rowLengths[rowNumber] += numCoupledFlowDofs;
            }
          }

          // (pressure_row, bubble_col): the matrix cell pressure couples to its 3 bubble DOFs (A_pb),
          // and so does the energy row packed right after the mass rows when thermal (A_Tb)
          globalIndex const pRow = pressureDofNumber[k] - rankOffset;
          if( pRow >= 0 && pRow < rowLengths.size() )
          {
            for( integer ic = 0; ic < numCoupledFlowDofs; ++ic )
              rowLengths[pRow + ic] += 3; // Three bubble DOFs
          }
        } );
      } );
    } );
  }

  void addMatrixTemperaturePressureBubbleCouplingPattern( DomainPartition const & domain,
                                                          DofManager const & dofManager,
                                                          SparsityPatternView< globalIndex > const & pattern ) const
  {
    GEOS_MARK_FUNCTION;

    integer const numCoupledFlowDofs = numCoupledFlowDofsPerElem();

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

          // (bubble_row, flow_col) : A_bp, plus A_bT when thermal
          for( localIndex i = 0; i < 3; ++i )
          {
            globalIndex const rowIndex = bubbleDofNumber[faceIndex] + i - rankOffset;
            if( rowIndex >= 0 && rowIndex < pattern.numRows() )
            {
              for( integer ic = 0; ic < numCoupledFlowDofs; ++ic )
                pattern.insertNonZero( rowIndex, pressureColIndex + ic );
            }
          }

          // (pressure_row, bubble_col) : A_pb -- transpose location, plus (temperature_row, bubble_col) : A_Tb
          globalIndex const pRow = pressureDofNumber[k] - rankOffset;
          if( pRow >= 0 && pRow < pattern.numRows() )
          {
            for( localIndex i = 0; i < 3; ++i )
            {
              for( integer ic = 0; ic < numCoupledFlowDofs; ++ic )
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

    if( !m_derivativeFluxResidual_dNormalJump )
    {
      setUpDerivativeFluxResidual_dNormalJump( domain );
    }
    m_derivativeFluxResidual_dNormalJump->move( parallelDeviceMemorySpace, false );
    m_derivativeFluxResidual_dNormalJump->zero();

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
                                                    m_derivativeFluxResidual_dNormalJumpEnergyOffset );

    m_derivativeFluxResidual_dNormalJump->move( hostMemorySpace, false );

    // This step must occur after the fluxes are assembled because that's when DerivativeFluxResidual_dNormalJump is filled.
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

    // Base::updateState above computed the fracture mass (and energy) with the previous deltaVolume.
    // Refresh them with the volume just updated
    // Note that the order cannot simply be swapped:
    // the flow porosity/permeability update of the surface elements would then overwrite dPerm_dDispJump
    if constexpr ( std::is_base_of_v< SinglePhaseBase, FLOW_SOLVER > || std::is_same_v< FLOW_SOLVER, SinglePhaseReservoirAndWells<> > )
    {
      this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                          MeshLevel & mesh,
                                                                          string_array const & regionNames )
      {
        mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                                            [&]( localIndex const,
                                                                                 FaceElementSubRegion & subRegion )
        {
          this->flowSolver()->updateMass( subRegion );
          if( this->m_isThermal )
          {
            this->flowSolver()->updateEnergy( subRegion );
          }
        } );
      } );
    }
  }

  virtual void resetStateToBeginningOfStep( DomainPartition & domain ) override
  {
    GEOS_MARK_FUNCTION;

    // rewind the flow and mechanics (incl. contact) states
    Base::resetStateToBeginningOfStep( domain );

    // The sequential loop also resets at dtAttempt == 0, before anything was solved: the state is then
    // still consistent with the beginning of step and the (heavy) update can be skipped.
    // The fully implicit loop only resets after a failure (cut or Newton failure), so always update.
    NonlinearSolverParameters const & solverParams = this->getNonlinearSolverParameters();
    bool const isFreshSequentialAttempt =
      solverParams.couplingType() == NonlinearSolverParameters::CouplingType::Sequential &&
      solverParams.m_numTimeStepAttempts == 0;

    if( !isFreshSequentialAttempt )
    {
      // the hydraulic aperture, fracture permeability, stencil weights and fracture mass (and energy)
      // still hold the values of the failed attempt: recompute them from the restored displacement jump
      updateState( domain );
    }
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


    integer const numRowsPerElem = numCoupledFlowDofsPerElem();

    this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &, //  meshBodyName,
                                                                        MeshLevel const & mesh,
                                                                        string_array const & ) // regionNames
    {
      ElementRegionManager const & elemManager = mesh.getElemManager();

      string const flowDofKey = dofManager.getKey( this->getFlowDofKey() );

      globalIndex const rankOffset = dofManager.rankOffset();

      NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
      FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();

      // Lagrange multipliers: the traction stabilization stencil connects the fracture elements
      if constexpr ( !CONTACT_SOLVER::hasFaceBubbles )
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
                    for( integer ic = 0; ic < numRowsPerElem; ic++ )
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
      // Face bubbles (augmented Lagrangian): the flow stencil connects the fracture elements, and the flow rows
      // also couple to the bubble DOFs of the neighbors
      else
      {

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
                    for( integer ic = 0; ic < numRowsPerElem; ic++ )
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
   * @brief Add the nonzeros induced by the flux-normal jump coupling
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param pattern the sparsity pattern
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
                                          [&faceToNodeMap]( localIndex const & faceIndex, localIndex const & a ){ return faceToNodeMap( faceIndex, a ); },
                                          [&faceToNodeMap]( localIndex const & faceIndex ){ return faceToNodeMap.sizeOfArray( faceIndex ); } );

      if constexpr (CONTACT_SOLVER::hasFaceBubbles)
      {
        addTransmissibilityCouplingPattern( domain, mesh, dofManager, pattern,
                                            faceManager.getReference< globalIndex_array >( dofManager.getKey( fields::contact::totalBubbleDisplacement::key() ) ),
                                            []( localIndex const & faceIndex, localIndex const & GEOS_UNUSED_PARAM( a )){ return faceIndex; },
                                            []( localIndex const & GEOS_UNUSED_PARAM( faceIndex )){ return 1; } );
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

    // Get the finite volume method used to compute the fracture fluxes
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
    flowDofNumber = fractureSubRegion.getReference< globalIndex_array >( dofManager.getKey( this->getFlowDofKey()) );

    globalIndex const rankOffset = dofManager.rankOffset();

    // the mass rows and, when thermal, the energy row couple to the same displacements
    integer const numRowsPerElem = numCoupledFlowDofsPerElem();

    fvDiscretization.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
    {
      forAll< serialPolicy >( stencil.size(), [=] ( localIndex const iconn )
      {
        localIndex const numFluxElems = stencil.stencilSize( iconn );
        typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

        // A fracture connector is an edge shared by two faces, or by more at a fracture junction:
        // the flow rows of each element k0 couple to the (nodal/bubble) displacements of every other
        // element k1 of the connector (the element's own displacements come from the dofManager coupling).
        for( localIndex k0 = 0; k0 < numFluxElems; ++k0 )
        {
          localIndex const rowIndex = flowDofNumber[sei[iconn][k0]] - rankOffset;
          if( rowIndex < 0 || rowIndex >= pattern.numRows() )
            continue;

          for( localIndex k1 = 0; k1 < numFluxElems; ++k1 )
          {
            if( k1 == k0 )
              continue;

            {
              localIndex const fractureIndex = sei[iconn][k1];

              // Loop over the two sides of each fracture element
              for( localIndex kf1 = 0; kf1 < 2; ++kf1 )
              {
                localIndex const faceIndex = elem2dToFaces[fractureIndex][kf1];

                // Save the list of DOF associated with nodes
                for( localIndex a=0; a<numNodesPerFace( elem2dToFaces[fractureIndex][0] ); ++a )
                {
                  for( localIndex i = 0; i < 3; ++i )
                  {
                    globalIndex const colIndex = coupledDisplacementDofNumber[dofIndirectionCb( faceIndex, a )] + LvArray::integerConversion< globalIndex >( i );
                    for( integer ic = 0; ic < numRowsPerElem; ic++ )
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
   * @brief Set up the matrix of the derivatives of the flux residual w.r.t. the normal jump
   * @param domain the physical domain object
   */
  void setUpDerivativeFluxResidual_dNormalJump( DomainPartition & domain )
  {
    integer const numComp = this->flowSolver()->numFluidComponents();

    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( this->flowSolver()->getDiscretizationName() );

    string const & fractureRegionName = this->solidMechanicsSolver()->getUniqueFractureRegionName();

    localIndex numMeshTargets = 0;
    this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                        MeshLevel const & mesh,
                                                                        string_array const & regionNames )
    {
      GEOS_UNUSED_VAR( regionNames );
      ElementRegionManager const & elemManager = mesh.getElemManager();

      // The matrix and every consumer of it (the flow solver's fracture flux kernels and the
      // coupling assembly) index it by the raw surface-element index, so a single mesh target
      // is supported: a second one would alias into the rows of the first.
      ++numMeshTargets;
      GEOS_ERROR_IF_GT_MSG( numMeshTargets, 1,
                            GEOS_FMT( "{}: this solver supports a single mesh target; '{}' is the second.",
                                      this->getName(), meshName ) );

      // For the same reason, the contact fracture must be the only face-element region on the
      // target. Embedded-surface regions hold a different subregion type and contribute no
      // SurfaceElementStencil here, so they are left alone.
      localIndex numFractureRegions = 0;
      elemManager.forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion const & region )
      {
        if( region.subRegionType() == SurfaceElementRegion::SurfaceSubRegionType::faceElement )
        {
          ++numFractureRegions;
        }
      } );
      GEOS_ERROR_IF_NE_MSG( numFractureRegions, 1,
                            GEOS_FMT( "{}: mesh target '{}' holds {} face-element regions. The conforming fracture "
                                      "formulation requires exactly one, named '{}'.",
                                      this->getName(), meshName, numFractureRegions, fractureRegionName ) );
      GEOS_ERROR_IF( !elemManager.hasRegion( fractureRegionName ),
                     GEOS_FMT( "{}: mesh target '{}' does not hold the fracture region '{}' of the contact solver.",
                               this->getName(), meshName, fractureRegionName ) );

      SurfaceElementRegion const & fractureRegion = elemManager.getRegion< SurfaceElementRegion >( fractureRegionName );
      FaceElementSubRegion const & fractureSubRegion = fractureRegion.getUniqueSubRegion< FaceElementSubRegion >();

      // Layout: numComp mass-balance rows per fracture element, then (when thermal) one
      // energy-balance row per fracture element appended as a second block, not interleaved.
      // The energy rows carry the advective (enthalpy-carried) and conductive contributions; the
      // latter through the aperture dependence of the stencil weights only (the conductivity itself
      // has no dConductivity/dDispJump constitutive derivative). One column per fracture element.
      localIndex const numElems = fractureSubRegion.size();
      localIndex const numMassRows = numElems * numComp;
      m_derivativeFluxResidual_dNormalJumpEnergyOffset = this->m_isThermal ? numMassRows : -1;
      localIndex const numRows = numMassRows + ( this->m_isThermal ? numElems : 0 );

      array1d< localIndex > rowCapacities( numRows );
      fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
      {
        for( localIndex iconn = 0; iconn < stencil.size(); ++iconn )
        {
          localIndex const numFluxElems = stencil.stencilSize( iconn );
          typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();
          for( localIndex k0 = 0; k0 < numFluxElems; ++k0 )
          {
            GEOS_ERROR_IF_GE_MSG( sei[iconn][k0], numElems,
                                  "Surface stencil index exceeds the fracture derivative matrix size." );
            for( integer ic = 0; ic < numComp; ++ic )
            {
              rowCapacities[ sei[iconn][k0] * numComp + ic ] += numFluxElems;
            }
            if( this->m_isThermal )
            {
              rowCapacities[ numMassRows + sei[iconn][k0] ] += numFluxElems;
            }
          }
        }
      } );

      std::unique_ptr< CRSMatrix< real64, localIndex > > & derivativeFluxResidual_dNormalJump = getRefDerivativeFluxResidual_dNormalJump();
      derivativeFluxResidual_dNormalJump = std::make_unique< CRSMatrix< real64, localIndex > >( numRows, numElems );
      derivativeFluxResidual_dNormalJump->setName( this->getName() + "/derivativeFluxResidual_dNormalJump" );
      if( numRows > 0 )
      {
        derivativeFluxResidual_dNormalJump->resizeFromRowCapacities< parallelHostPolicy >( numRows,
                                                                                         numElems,
                                                                                         rowCapacities.data() );
      }

      fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
      {
        for( localIndex iconn = 0; iconn < stencil.size(); ++iconn )
        {
          localIndex const numFluxElems = stencil.stencilSize( iconn );
          typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();
          for( localIndex k0 = 0; k0 < numFluxElems; ++k0 )
          {
            for( localIndex k1 = 0; k1 < numFluxElems; ++k1 )
            {
              for( integer ic = 0; ic < numComp; ++ic )
              {
                derivativeFluxResidual_dNormalJump->insertNonZero( sei[iconn][k0] * numComp + ic, sei[iconn][k1], 0.0 );
              }
              if( this->m_isThermal )
              {
                derivativeFluxResidual_dNormalJump->insertNonZero( numMassRows + sei[iconn][k0], sei[iconn][k1], 0.0 );
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
    // the bulk kernel has just recomputed the porosity from the nodal strain, discarding any enrichment increment
    m_isPorosityReset = true;

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
    // These 2 steps need to occur after the fluxes are assembled because that's when DerivativeFluxResidual_dNormalJump is filled.
    this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                        MeshLevel const & mesh,
                                                                        string_array const & regionNames )
    {
      /// 3. assemble Force Residual w.r.t. pressure and Flow mass residual w.r.t. displacement
      assembleForceResidualDerivativeWrtPressure( meshName, mesh, regionNames, dofManager, localMatrix, localRhs );
      assembleFluidMassResidualDerivativeWrtDisplacement( meshName, mesh, regionNames, dofManager, localMatrix, localRhs );
    } );

    // face bubbles: matrix cell flow to/from bubble couplings (A_bp, A_bT, A_pb, A_Tb)
    if constexpr (CONTACT_SOLVER::hasFaceBubbles)
      assembleMatrixPressureBubbleContribution( dt, const_cast< DomainPartition & >( domain ), dofManager, localMatrix, localRhs );

  }

  virtual void assembleForceResidualDerivativeWrtPressure( string const & GEOS_UNUSED_PARAM( meshName ),
                                                           MeshLevel const & mesh,
                                                           string_array const & regionNames,
                                                           DofManager const & dofManager,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                           arrayView1d< real64 > const & localRhs ) = 0;

  virtual void assembleFluidMassResidualDerivativeWrtDisplacement( string const & meshName,
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
  virtual void assembleMatrixPressureBubbleContribution( real64 const GEOS_UNUSED_PARAM( dt ),
                                                         DomainPartition & GEOS_UNUSED_PARAM( domain ),
                                                         DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                                         CRSMatrixView< real64, globalIndex const > const & GEOS_UNUSED_PARAM( localMatrix ),
                                                         arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ) )
  { GEOS_WARNING( "Should override" ); };

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
        arrayView1d< real64 > const dHydraulicAperture_dNormalJump = subRegion.getField< fields::flow::dHydraulicAperture_dNormalJump >();
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
                                                dHydraulicAperture_dNormalJump,
                                                fractureTraction,
                                                fractureState );

          } );
        } );
      } );
    } );
  }

  std::unique_ptr< CRSMatrix< real64, localIndex > > & getRefDerivativeFluxResidual_dNormalJump()
  {
    return m_derivativeFluxResidual_dNormalJump;
  }

  CRSMatrixView< real64, localIndex const > getDerivativeFluxResidual_dNormalJump()
  {
    return m_derivativeFluxResidual_dNormalJump->toViewConstSizes();
  }

  CRSMatrixView< real64 const, localIndex const > getDerivativeFluxResidual_dNormalJump() const
  {
    return m_derivativeFluxResidual_dNormalJump->toViewConst();
  }

  /**
   * @brief First row of the energy-balance block appended after all mass-balance rows in
   * getDerivativeFluxResidual_dNormalJump(), or -1 when the flow solver is not thermal.
   */
  localIndex getDerivativeFluxResidual_dNormalJumpEnergyOffset() const
  {
    return m_derivativeFluxResidual_dNormalJumpEnergyOffset;
  }

  /**
   * @brief Offset of the temperature DOF, and of the energy-balance row, from the first flow DOF of an element.
   * The temperature is packed right after the numFluidComponents() mass DOFs (thermal is single-phase only).
   */
  integer temperatureDofOffset() const
  {
    return this->flowSolver()->numFluidComponents();
  }

  /**
   * @brief Number of flow DOFs (and rows) of an element coupled to the displacements: the
   * numFluidComponents() mass ones, then the temperature (energy-balance) one when thermal.
   */
  integer numCoupledFlowDofsPerElem() const
  {
    return this->flowSolver()->numFluidComponents() + ( this->m_isThermal ? 1 : 0 );
  }

  struct viewKeyStruct : public Base::viewKeyStruct
  {};

  static const localIndex m_maxFaceNodes = 11; // Maximum number of nodes on a contact face

  std::unique_ptr< CRSMatrix< real64, localIndex > > m_derivativeFluxResidual_dNormalJump;
  /// First row of the energy-balance block in m_derivativeFluxResidual_dNormalJump (-1 if not thermal)
  localIndex m_derivativeFluxResidual_dNormalJumpEnergyOffset = -1;

  /// True when the bulk poromechanics kernel has reset the porosity and no enrichment (face bubble)
  /// increment has been added to it since: such an increment is added on top of the stored porosity,
  /// so adding it twice without a reset in between would count it twice.
  bool m_isPorosityReset = false;

};

}

#endif //GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSCONFORMINGFRACTURES_HPP_
