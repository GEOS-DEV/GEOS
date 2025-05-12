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

#include "finiteVolume/FluxApproximationBase.hpp"

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

  /// We need to add 2 coupling terms:
  // 1. Poromechanical coupling in the bulk
  Base::setupCoupling( domain, dofManager );
}


template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::setupSystem( DomainPartition & domain,
                                                                                 DofManager & dofManager,
                                                                                 CRSMatrix< real64, globalIndex > & localMatrix,
                                                                                 ParallelVector & rhs,
                                                                                 ParallelVector & solution,
                                                                                 bool const setSparsity )
{

  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( setSparsity );

  /// 1. Add all coupling terms handled directly by the DofManager
  //     (Auu, Aup, Apu, App, Apfpf, Apfu, Aupf)
  dofManager.setDomain( domain );
  this->setupDofs( domain, dofManager );
  dofManager.reorderByRank();

  /// 2. Add coupling terms not added by the DofManager.
  //     (Aupf) Add the coupling with the nodal displacements of the neighbors  
  //     due to armonic averaging of the transmissibility
  localIndex const numLocalRows = dofManager.numLocalDofs();

  SparsityPattern< globalIndex > patternOriginal;
  dofManager.setSparsityPattern( patternOriginal );

  // Get the original row lengths (diagonal blocks only)
  array1d< localIndex > rowLengths( patternOriginal.numRows() );
  for( localIndex localRow = 0; localRow < patternOriginal.numRows(); ++localRow )
  {
    rowLengths[localRow] = patternOriginal.numNonZeros( localRow );
  }

  // Add the number of nonzeros induced by coupling
  addTransmissibilityCouplingNNZ( domain, dofManager, rowLengths.toView() );

  // Create a new pattern with enough capacity for coupled matrix
  SparsityPattern< globalIndex > pattern;
  pattern.resizeFromRowCapacities< parallelHostPolicy >( patternOriginal.numRows(),
                                                         patternOriginal.numColumns(),
                                                         rowLengths.data() );

  // Copy the original nonzeros
  for( localIndex localRow = 0; localRow < patternOriginal.numRows(); ++localRow )
  {
    globalIndex const * cols = patternOriginal.getColumns( localRow ).dataIfContiguous();
    pattern.insertNonZeros( localRow, cols, cols + patternOriginal.numNonZeros( localRow ) );
  }

  // Add the nonzeros from coupling
  addTransmissibilityCouplingPattern( domain, dofManager, pattern.toView() );

  localMatrix.setName( this->getName() + "/matrix" );
  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );

  rhs.setName( this->getName() + "/rhs" );
  rhs.create( numLocalRows, MPI_COMM_GEOS );

  solution.setName( this->getName() + "/solution" );
  solution.create( numLocalRows, MPI_COMM_GEOS );

  setUpDflux_dApertureMatrix( domain, dofManager, localMatrix );

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

  this->solidMechanicsSolver()->synchronizeFractureState( domain );

  // Assembly elements-based terms 
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
                                                  getDerivativeFluxResidual_dNormalJump() );

  // This step must occur after the fluxes are assembled because that's when DerivativeFluxResidual_dAperture is filled.
  assembleCouplingTerms( time_n,
                         dt,
                         domain,
                         dofManager,
                         localMatrix,
                         localRhs );
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
  /// 3. assemble Force Residual w.r.t. pressure and Flow mass residual w.r.t. displacement

  Base::assembleElementBasedTerms( time_n, dt, domain, dofManager, localMatrix, localRhs );

  // Flow accumulation for fractures
  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel & mesh,
                                                                      arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          FaceElementSubRegion const & subRegion )
    {
      this->flowSolver()->accumulationAssemblyLaunch( dofManager, subRegion, localMatrix, localRhs );
    } );
  } );

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
  GEOS_UNUSED_VAR( time_n, dt );
  // These 2 steps need to occur after the fluxes are assembled because that's when DerivativeFluxResidual_dAperture is filled.
  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                      MeshLevel const & mesh,
                                                                      arrayView1d< string const > const & regionNames )
  {
    /// 3. assemble Force Residual w.r.t. pressure and Flow mass residual w.r.t. displacement
    assembleForceResidualDerivativeWrtPressure( meshName, mesh, regionNames, dofManager, localMatrix, localRhs );
    assembleFluidMassResidualDerivativeWrtDisplacement( mesh, regionNames, dofManager, localMatrix, localRhs );
  } );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( domain );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::
setUpDflux_dApertureMatrix( DomainPartition & domain,
                            DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                            CRSMatrix< real64, globalIndex > & localMatrix )
{
  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel const & mesh,
                                                                      arrayView1d< string const > const & regionNames )
  {
    std::unique_ptr< CRSMatrix< real64, localIndex > > & derivativeFluxResidual_dAperture = this->getRefDerivativeFluxResidual_dAperture();

    {
      localIndex numRows = 0;
      mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                                          [&]( localIndex const, FaceElementSubRegion const & subRegion )
      {
        numRows += subRegion.size();
      } );

      derivativeFluxResidual_dAperture = std::make_unique< CRSMatrix< real64, localIndex > >( numRows, numRows );
      derivativeFluxResidual_dAperture->setName( this->getName() + "/derivativeFluxResidual_dAperture" );

      derivativeFluxResidual_dAperture->reserveNonZeros( localMatrix.numNonZeros() );
      localIndex maxRowSize = -1;
      for( localIndex row = 0; row < localMatrix.numRows(); ++row )
      {
        localIndex const rowSize = localMatrix.numNonZeros( row );
        maxRowSize = maxRowSize > rowSize ? maxRowSize : rowSize;
      }
      // TODO This is way too much. The With the full system rowSize is not a good estimate for this.
      for( localIndex row = 0; row < numRows; ++row )
      {
        derivativeFluxResidual_dAperture->reserveNonZeros( row, maxRowSize );
      }
    }

    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( this->flowSolver()->getDiscretizationName() );

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
            derivativeFluxResidual_dAperture->insertNonZero( sei[iconn][k0], sei[iconn][k1], 0.0 );
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
                                                                      arrayView1d< string const > const & ) 
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const presDofKey = dofManager.getKey( m_pressureKey );

    globalIndex const rankOffset = dofManager.rankOffset();

    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( this->flowSolver()->getDiscretizationName() );

    fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
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
          elementSubRegion.getReference< array1d< globalIndex > >( presDofKey );

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
                rowLengths[rowNumber] += 3*numNodesPerElement;
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

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel const & mesh,
                                                                      arrayView1d< string const > const & )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    NodeManager const & nodeManager = mesh.getNodeManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
    string const presDofKey = dofManager.getKey( m_pressureKey );

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

    GEOS_ERROR_IF( !fractureSubRegion.hasWrapper( flow::pressure::key() ),
                   this->getDataContext() << ": The fracture subregion must contain pressure field." );

    arrayView2d< localIndex const > const elem2dToFaces = fractureSubRegion.faceList().toViewConst();

    arrayView1d< globalIndex const > const &
    presDofNumber = fractureSubRegion.getReference< globalIndex_array >( presDofKey );

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
            localIndex const rowIndex = presDofNumber[sei[iconn][1-kf]] - rankOffset;

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
                    pattern.insertNonZero( rowIndex, colIndex );
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

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::
assembleForceResidualDerivativeWrtPressure( string const & meshName,
                                            MeshLevel const & mesh,
                                            arrayView1d< string const > const & regionNames,
                                            DofManager const & dofManager,
                                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                            arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;


  NodeManager const & nodeManager = mesh.getNodeManager();
  FaceManager const & faceManager = mesh.getFaceManager();

  string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
  string const & bubbleDofKey = dofManager.getKey( solidMechanics::totalBubbleDisplacement::key() );

  arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
  arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );

  string const & fractureRegionName = this->solidMechanicsSolver()->getUniqueFractureRegionName();

  CRSMatrix< real64, globalIndex > const voidMatrix;
  array1d< real64 > const voidRhs;

  this->solidMechanicsSolver()->forFiniteElementOnFractureSubRegions( meshName, [&] ( string const &,
                                                        finiteElement::FiniteElementBase const & subRegionFE,
                                                        arrayView1d< localIndex const > const & faceElementList )
  {

    GEOS_UNUSED_VAR( subRegionFE, faceElementList, regionNames, localMatrix, localRhs, fractureRegionName );

    //solidMechanicsConformingContactKernels::DispJumpUpdateFactory kernelFactory( dispDofNumber,
    //                                                                             bubbleDofNumber,
    //                                                                             dofManager.rankOffset(),
    //                                                                             voidMatrix.toViewConstSizes(),
    //                                                                             voidRhs.toView(),
    //                                                                             dt,
    //                                                                             faceElementList );

    //real64 maxTraction = finiteElement::
    //                       interfaceBasedKernelApplication
    //                     < parallelDevicePolicy< >,
    //                       constitutive::NullModel >( mesh,
    //                                                  fractureRegionName,
    //                                                  faceElementList,
    //                                                  subRegionFE,
    //                                                  "",
    //                                                  kernelFactory );

    //GEOS_UNUSED_VAR( maxTraction );

  } );

}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::
assembleFluidMassResidualDerivativeWrtDisplacement( MeshLevel const & mesh,
                                                    arrayView1d< string const > const & regionNames,
                                                    DofManager const & dofManager,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ) )
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( mesh, regionNames, dofManager, localMatrix );

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
