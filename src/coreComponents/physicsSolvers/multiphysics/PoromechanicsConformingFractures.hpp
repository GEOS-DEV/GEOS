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

#include "physicsSolvers/solidMechanics/contact/SolidMechanicsLagrangeContact.hpp"
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

template< template< typename, typename > class POROMECHANICS_BASE, typename FLOW_SOLVER >
class PoromechanicsConformingFractures : public POROMECHANICS_BASE< FLOW_SOLVER, SolidMechanicsLagrangeContact >
{
public:
  using Base = POROMECHANICS_BASE< FLOW_SOLVER, SolidMechanicsLagrangeContact >;

  PoromechanicsConformingFractures( const string & name,
                                    dataRepository::Group * const parent )
    : Base( name, parent )
  {}

  virtual void setupCoupling( DomainPartition const & domain,
                              DofManager & dofManager ) const override
  {
    /// We need to add 2 coupling terms:
    // 1. Poromechanical coupling in the bulk
    Base::setupCoupling( domain, dofManager );

    // 2. Traction - pressure coupling in the fracture
    dofManager.addCoupling( this->getFlowDofKey(),
                            fields::contact::traction::key(),
                            DofManager::Connector::Elem );
  }

  virtual void setSparsityPattern( DomainPartition & domain,
                                   DofManager & dofManager,
                                   CRSMatrix< real64, globalIndex > & localMatrix,
                                   SparsityPattern< globalIndex > & pattern ) override
  {
    // start with the flow solver sparsity pattern (it could be reservoir + wells)
    SparsityPattern< globalIndex > patternOriginal;
    this->flowSolver()->setSparsityPattern( domain, dofManager, localMatrix, patternOriginal );

    // Get the original row lengths (diagonal blocks only)
    array1d< localIndex > rowLengths( patternOriginal.numRows());
    for( localIndex localRow = 0; localRow < patternOriginal.numRows(); ++localRow )
    {
      rowLengths[localRow] = patternOriginal.numNonZeros( localRow );
    }

    // Add the number of nonzeros induced by coupling
    addTransmissibilityCouplingNNZ( domain, dofManager, rowLengths.toView());

    // Create a new pattern with enough capacity for coupled matrix
    pattern.resizeFromRowCapacities< parallelHostPolicy >( patternOriginal.numRows(),
                                                           patternOriginal.numColumns(),
                                                           rowLengths.data());

    // Copy the original nonzeros
    appendSparsityPattern( pattern, patternOriginal );

    // Add the nonzeros from coupling
    addTransmissibilityCouplingPattern( domain, dofManager, pattern.toView());

    setUpDflux_dApertureMatrix( domain );
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

    // number of rows a fracture element occupies in the flow block: the mass balance equations, plus
    // the energy balance equation when thermal
    integer const numComp = numFluidComponents() + ( this->m_isThermal ? 1 : 0 );

    this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &, //  meshBodyName,
                                                                        MeshLevel const & mesh,
                                                                        string_array const & ) // regionNames
    {
      ElementRegionManager const & elemManager = mesh.getElemManager();

      string const flowDofKey = dofManager.getKey( this->getFlowDofKey() );

      globalIndex const rankOffset = dofManager.rankOffset();

      NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
      FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
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
    } );
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

    // number of rows a fracture element occupies in the flow block: the mass balance equations, plus
    // the energy balance equation when thermal
    integer const numComp = numFluidComponents() + ( this->m_isThermal ? 1 : 0 );

    this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                        MeshLevel const & mesh,
                                                                        string_array const & )
    {
      FaceManager const & faceManager = mesh.getFaceManager();
      NodeManager const & nodeManager = mesh.getNodeManager();
      ElementRegionManager const & elemManager = mesh.getElemManager();

      string const dispDofKey = dofManager.getKey( fields::solidMechanics::totalDisplacement::key() );
      string const flowDofKey = dofManager.getKey( this->getFlowDofKey() );

      arrayView1d< globalIndex const > const &
      dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
      ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();

      // Get the finite volume method used to compute the stabilization
      NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
      FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
      FluxApproximationBase const & fvDiscretization = fvManager.getFluxApproximation( this->flowSolver()->getDiscretizationName() );

      SurfaceElementRegion const & fractureRegion =
        elemManager.getRegion< SurfaceElementRegion >( this->solidMechanicsSolver()->getUniqueFractureRegionName() );
      FaceElementSubRegion const & fractureSubRegion =
        fractureRegion.getUniqueSubRegion< FaceElementSubRegion >();

      GEOS_ERROR_IF( !fractureSubRegion.hasWrapper( fields::flow::pressure::key() ),
                     "The fracture subregion must contain pressure field.", this->getDataContext() );

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
              localIndex const rowIndex = flowDofNumber[sei[iconn][1-kf]] - rankOffset;

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

                  // Save the list of DOF associated with nodes
                  for( localIndex a=0; a<numNodesPerFace; ++a )
                  {
                    for( localIndex i = 0; i < 3; ++i )
                    {
                      globalIndex const colIndex = dispDofNumber[faceToNodeMap( faceIndex, a )] + LvArray::integerConversion< globalIndex >( i );
                      for( integer ic = 0; ic < numComp; ic++ )
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
    } );
  }

  /**
   * @brief Set up the Dflux_dApertureMatrix object
   *
   * @param domain
   */
  void setUpDflux_dApertureMatrix( DomainPartition & domain )
  {
    // number of rows a fracture element occupies in the flow block: the mass balance equations, plus
    // the energy balance equation when thermal
    integer const numComp = numFluidComponents() + ( this->m_isThermal ? 1 : 0 );
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( this->flowSolver()->getDiscretizationName() );

    localIndex numMeshTargets = 0;
    this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                        MeshLevel const & mesh,
                                                                        string_array const & regionNames )
    {
      std::unique_ptr< CRSMatrix< real64, localIndex > > & derivativeFluxResidual_dAperture = getRefDerivativeFluxResidual_dAperture();

      // The matrix is re-created per target and the flux kernel indexes it by
      // the raw per-target surface element index, so only the last target would
      // survive and the others would write into its rows.
      ++numMeshTargets;
      GEOS_ERROR_IF_GT_MSG( numMeshTargets, 1,
                            GEOS_FMT( "{}: this solver supports a single mesh target; '{}' is the second.",
                                      this->getName(), meshName ) );

      localIndex numRows = 0;
      localIndex numCol = 0;
      {
        // calculate number of fracture elements
        mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                                            [&]( localIndex const, FaceElementSubRegion const & subRegion )
        {
          numRows += subRegion.size();
        } );
        // number of columns (derivatives) = number of fracture elements
        numCol = numRows;
        // number of rows (equations) = number of fracture elements * number of components
        numRows *= numComp;

        derivativeFluxResidual_dAperture = std::make_unique< CRSMatrix< real64, localIndex > >( numRows, numCol );
        derivativeFluxResidual_dAperture->setName( this->getName() + "/derivativeFluxResidual_dAperture" );
      }

      // array1d's sized constructor value-initializes, so no explicit zero().
      array1d< localIndex > rowCapacities( numRows );
      fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
      {
        for( localIndex iconn = 0; iconn < stencil.size(); ++iconn )
        {
          localIndex const numFluxElems = stencil.stencilSize( iconn );
          typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

          for( localIndex k0 = 0; k0 < numFluxElems; ++k0 )
          {
            // The stencil sweep covers every SurfaceElementStencil on the mesh,
            // while numRows only counts the subregions found in regionNames.
            GEOS_ERROR_IF_GE_MSG( sei[iconn][k0] * numComp + numComp - 1, numRows,
                                  "Surface stencil index exceeds the fracture derivative matrix size." );
            for( integer ic = 0; ic < numComp; ic++ )
            {
              rowCapacities[sei[iconn][k0] * numComp + ic] += numFluxElems;
            }
          }
        }
      } );

      if( numRows > 0 )
      {
        derivativeFluxResidual_dAperture->resizeFromRowCapacities< parallelHostPolicy >( numRows,
                                                                                         numCol,
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
            GEOS_ERROR_IF_GE_MSG( sei[iconn][k0] * numComp + numComp - 1, numRows,
                                  "Surface stencil index exceeds the fracture derivative matrix size." );
            for( localIndex k1 = 0; k1 < numFluxElems; ++k1 )
            {
              for( integer ic = 0; ic < numComp; ++ic )
              {
                derivativeFluxResidual_dAperture->insertNonZero( sei[iconn][k0] * numComp + ic,
                                                                 sei[iconn][k1],
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

    this->solidMechanicsSolver()->assembleContact( domain, dofManager, localMatrix, localRhs );
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
    this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                        MeshLevel const & mesh,
                                                                        string_array const & regionNames )
    {
      /// 3. assemble Force Residual w.r.t. pressure and Flow mass residual w.r.t. displacement
      assembleForceResidualDerivativeWrtPressure( mesh, regionNames, dofManager, localMatrix, localRhs );
      assembleFluidMassResidualDerivativeWrtDisplacement( mesh, regionNames, dofManager, localMatrix, localRhs );
    } );
  }

  void assembleForceResidualDerivativeWrtPressure( MeshLevel const & mesh,
                                                   string_array const & regionNames,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs )
  {
    GEOS_MARK_FUNCTION;

    FaceManager const & faceManager = mesh.getFaceManager();
    NodeManager const & nodeManager = mesh.getNodeManager();
    EdgeManager const & edgeManager = mesh.getEdgeManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();
    ArrayOfArraysView< localIndex const > const faceToEdgeMap = faceManager.edgeList().toViewConst();
    arrayView2d< localIndex const > const & edgeToNodeMap = edgeManager.nodeList().toViewConst();
    arrayView2d< real64 const > faceCenters = faceManager.faceCenter();
    arrayView2d< real64 const > const & faceNormal = faceManager.faceNormal();
    arrayView1d< real64 const > faceAreas = faceManager.faceArea();

    string const & dispDofKey = dofManager.getKey( fields::solidMechanics::totalDisplacement::key() );
    string const & flowDofKey = dofManager.getKey( this->getFlowDofKey() );

    arrayView1d< globalIndex const > const &
    dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
    globalIndex const rankOffset = dofManager.rankOffset();

    // Get the coordinates for all nodes
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   FaceElementSubRegion const & subRegion )
    {
      arrayView1d< globalIndex const > const &
      flowDofNumber = subRegion.getReference< globalIndex_array >( flowDofKey );
      arrayView1d< real64 const > const & pressure = subRegion.getReference< array1d< real64 > >( fields::flow::pressure::key() );
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList().toViewConst();

      forAll< serialPolicy >( subRegion.size(), [=, this]( localIndex const kfe )
      {
        localIndex const kf0 = elemsToFaces[kfe][0];
        localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );

        real64 Nbar[3];
        Nbar[ 0 ] = faceNormal[elemsToFaces[kfe][0]][0] - faceNormal[elemsToFaces[kfe][1]][0];
        Nbar[ 1 ] = faceNormal[elemsToFaces[kfe][0]][1] - faceNormal[elemsToFaces[kfe][1]][1];
        Nbar[ 2 ] = faceNormal[elemsToFaces[kfe][0]][2] - faceNormal[elemsToFaces[kfe][1]][2];
        LvArray::tensorOps::normalize< 3 >( Nbar );
        globalIndex rowDOF[3 * m_maxFaceNodes]; // this needs to be changed when dealing with arbitrary element types
        real64 nodeRHS[3 * m_maxFaceNodes];
        stackArray1d< real64, 3 * m_maxFaceNodes > dRdP( 3*m_maxFaceNodes );
        globalIndex colDOF[1];
        colDOF[0] = flowDofNumber[kfe]; // pressure is always first

        for( localIndex kf=0; kf<2; ++kf )
        {
          localIndex const faceIndex = elemsToFaces[kfe][kf];

          // Compute local area contribution for each node
          stackArray1d< real64, FaceManager::maxFaceNodes() > nodalArea;
          this->solidMechanicsSolver()->computeFaceNodalArea( elemsToFaces[kfe][kf],
                                                              nodePosition,
                                                              faceToNodeMap,
                                                              faceToEdgeMap,
                                                              edgeToNodeMap,
                                                              faceCenters,
                                                              faceNormal,
                                                              faceAreas,
                                                              nodalArea );
          for( localIndex a=0; a<numNodesPerFace; ++a )
          {
            real64 const nodalForceMag = -( pressure[kfe] ) * nodalArea[a];
            real64 globalNodalForce[ 3 ];
            LvArray::tensorOps::scaledCopy< 3 >( globalNodalForce, Nbar, nodalForceMag );

            for( localIndex i=0; i<3; ++i )
            {
              rowDOF[3*a+i] = dispDofNumber[faceToNodeMap( faceIndex, a )] + LvArray::integerConversion< globalIndex >( i );
              // Opposite sign w.r.t. theory because of minus sign in stiffness matrix definition (K < 0)
              nodeRHS[3*a+i] = +globalNodalForce[i] * pow( -1, kf );

              // Opposite sign w.r.t. theory because of minus sign in stiffness matrix definition (K < 0)
              dRdP( 3*a+i ) = -nodalArea[a] * Nbar[i] * pow( -1, kf );
            }
          }

          for( localIndex idof = 0; idof < numNodesPerFace * 3; ++idof )
          {
            localIndex const localRow = LvArray::integerConversion< localIndex >( rowDOF[idof] - rankOffset );

            if( localRow >= 0 && localRow < localMatrix.numRows() )
            {
              localMatrix.addToRow< parallelHostAtomic >( localRow,
                                                          colDOF,
                                                          &dRdP[idof],
                                                          1 );
              RAJA::atomicAdd( parallelHostAtomic{}, &localRhs[localRow], nodeRHS[idof] );
            }
          }
        }
      } );
    } );
  }

  virtual void assembleFluidMassResidualDerivativeWrtDisplacement( MeshLevel const & mesh,
                                                                   string_array const & regionNames,
                                                                   DofManager const & dofManager,
                                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                   arrayView1d< real64 > const & localRhs ) = 0;

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

  virtual integer numFluidComponents() const = 0;

  struct viewKeyStruct : public Base::viewKeyStruct
  {};

  static const localIndex m_maxFaceNodes = 11; // Maximum number of nodes on a contact face

  std::unique_ptr< CRSMatrix< real64, localIndex > > m_derivativeFluxResidual_dAperture;

};

}

#endif //GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSCONFORMINGFRACTURES_HPP_
