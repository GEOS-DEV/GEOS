/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MultiphasePoromechanicsConformingFractures.cpp
 */

#include "MultiphasePoromechanicsConformingFractures.hpp"


namespace geos
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;

template< typename FLOW_SOLVER >
MultiphasePoromechanicsConformingFractures< FLOW_SOLVER >::MultiphasePoromechanicsConformingFractures( const string & name,
                                                                                                       Group * const parent )
  : Base( name, parent )
{
  // TODO: MGR recipe
  // LinearSolverParameters & params = this->m_linearSolverParameters.get();
  // params.mgr.strategy = LinearSolverParameters::MGR::StrategyType::multiphasePoromechanicsConformingFractures;
  // params.mgr.separateComponents = false;
  // params.mgr.displacementFieldName = solidMechanics::totalDisplacement::key();
  // params.dofsPerNode = 3;
}

template< typename FLOW_SOLVER >
void MultiphasePoromechanicsConformingFractures< FLOW_SOLVER >::initializePreSubGroups()
{
  Base::initializePreSubGroups();
  if( this->flowSolver()->isThermal() )
  {
    GEOS_ERROR( "Thermal flow is not yet supported for multiphase poromechanics conforming fractures" );
  }
}

template< typename FLOW_SOLVER >
void MultiphasePoromechanicsConformingFractures< FLOW_SOLVER >::
assembleSystem( real64 const time_n,
                real64 const dt,
                DomainPartition & domain,
                DofManager const & dofManager,
                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  Base::assembleSystem( time_n, dt, domain, dofManager, localMatrix, localRhs );
}

template<>
void MultiphasePoromechanicsConformingFractures< CompositionalMultiphaseReservoirAndWells<> >::
assembleSystem( real64 const time_n,
                real64 const dt,
                DomainPartition & domain,
                DofManager const & dofManager,
                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  Base::assembleSystem( time_n, dt, domain, dofManager, localMatrix, localRhs );

  flowSolver()->wellSolver()->assembleSystem( time_n, dt, domain, dofManager, localMatrix, localRhs );
  flowSolver()->assembleCouplingTerms( time_n, dt, domain, dofManager, localMatrix, localRhs );
}

template< typename FLOW_SOLVER >
void MultiphasePoromechanicsConformingFractures< FLOW_SOLVER >::
assembleFluidMassResidualDerivativeWrtDisplacement( MeshLevel const & mesh,
                                                    string_array const & regionNames,
                                                    DofManager const & dofManager,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ) )
{
  GEOS_MARK_FUNCTION;

  using namespace contact;

  integer const numComp = this->flowSolver()->numFluidComponents();

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

  CRSMatrixView< real64 const, localIndex const > const &
  dFluxResidual_dNormalJump = this->getDerivativeFluxResidual_dNormalJump().toViewConst();

  string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
  string const & flowDofKey = dofManager.getKey( getFlowDofKey() );

  arrayView1d< globalIndex const > const &
  dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
  globalIndex const rankOffset = dofManager.rankOffset();

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                            [&]( localIndex const,
                                                                 FaceElementSubRegion const & subRegion )
  {
    arrayView2d< real64 const, compflow::USD_COMP > const compDens = subRegion.getField< flow::globalCompDensity >();

    arrayView1d< globalIndex const > const & flowDofNumber = subRegion.getReference< array1d< globalIndex > >( flowDofKey );

    arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList().toViewConst();
    arrayView1d< real64 const > const & area = subRegion.getElementArea().toViewConst();

    arrayView1d< integer const > const fractureState = subRegion.getField< contact::fractureState >();

    forAll< serialPolicy >( subRegion.size(), [&]( localIndex const kfe )
    {
      localIndex const kf0 = elemsToFaces[kfe][0], kf1 = elemsToFaces[kfe][1];
      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );
      globalIndex nodeDOF[2*3*m_maxFaceNodes];

      real64 Nbar[3];
      Nbar[ 0 ] = faceNormal[kf0][0] - faceNormal[kf1][0];
      Nbar[ 1 ] = faceNormal[kf0][1] - faceNormal[kf1][1];
      Nbar[ 2 ] = faceNormal[kf0][2] - faceNormal[kf1][2];
      LvArray::tensorOps::normalize< 3 >( Nbar );

      stackArray2d< real64, 2*3*m_maxFaceNodes * MultiFluidBase::MAX_NUM_COMPONENTS > dRdU( MultiFluidBase::MAX_NUM_COMPONENTS, 2*3*m_maxFaceNodes );

      bool const isFractureOpen = ( fractureState[kfe] == FractureState::Open );

      // Accumulation derivative
      if( isFractureOpen )
      {
        for( localIndex kf=0; kf<2; ++kf )
        {
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

          // TODO: move to something like this plus a static method.
          // localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][kf] );
          // stackArray1d<real64, 4> nodalArea( numNodesPerFace );

          for( localIndex a=0; a<numNodesPerFace; ++a )
          {
            real64 const dVolume_dAperture = nodalArea[a];
            for( localIndex i=0; i<3; ++i )
            {
              nodeDOF[ kf*3*numNodesPerFace + 3*a+i ] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe][kf], a )]
                                                        + LvArray::integerConversion< globalIndex >( i );
              real64 const dAper_dU = -pow( -1, kf ) * Nbar[i];
              for( integer ic = 0; ic < numComp; ic++ )
              {
                dRdU[ ic ][ kf*3*numNodesPerFace + 3*a + i ] = dVolume_dAperture * dAper_dU * compDens[kfe][ic]; // assuming poro=1
              }
            }
          }
        }

        localIndex const localRow = LvArray::integerConversion< localIndex >( flowDofNumber[kfe] - rankOffset );
        if( localRow >= 0 && localRow < localMatrix.numRows() )
        {
          integer const numRows = numComp;
          for( integer i = 0; i < numRows; ++i )
          {
            localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow+i,
                                                                      nodeDOF,
                                                                      dRdU[i],
                                                                      2 * 3 * numNodesPerFace );
          }
        }
      }

      // flux derivative
      bool skipAssembly = true;
      for( integer ic = 0; ic < numComp; ic++ )
      {
        localIndex const numColumns = dFluxResidual_dNormalJump.numNonZeros( kfe * numComp + ic );
        arraySlice1d< localIndex const > const & columns = dFluxResidual_dNormalJump.getColumns( kfe * numComp + ic );
        // TODO uncomment and fix build
        arraySlice1d< real64 const > const & values = dFluxResidual_dNormalJump.getEntries( kfe * numComp + ic );

        skipAssembly &= !isFractureOpen;

        for( localIndex kfe1 = 0; kfe1 < numColumns; ++kfe1 )
        {
          real64 const dR_dAper = values[kfe1];
          localIndex const kfe2 = columns[kfe1];

          bool const isOpen = ( fractureState[kfe2] == FractureState::Open );
          skipAssembly &= !isOpen;

          for( localIndex kf = 0; kf < 2; ++kf )
          {
            //TODO: We should avoid allocating LvArrays inside kernel
            stackArray1d< real64, FaceManager::maxFaceNodes() > nodalArea;
            this->solidMechanicsSolver()->computeFaceNodalArea( elemsToFaces[kfe2][kf],
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
              for( localIndex i=0; i<3; ++i )
              {
                nodeDOF[ kf*3*numNodesPerFace + 3*a+i ] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe2][kf], a )]
                                                          + LvArray::integerConversion< globalIndex >( i );
                real64 const dAper_dU = -pow( -1, kf ) * Nbar[i] * ( nodalArea[a] / area[kfe2] );
                dRdU[ ic ][ kf*3*numNodesPerFace + 3*a + i ] = dR_dAper * dAper_dU;
              }
            }
          }
        }

        if( !skipAssembly )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( flowDofNumber[kfe] - rankOffset );
          if( localRow >= 0 && localRow < localMatrix.numRows() )
          {
            integer const numRows = numComp;
            for( integer i = 0; i < numRows; ++i )
            {
              localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow+i,
                                                                        nodeDOF,
                                                                        dRdU[i],
                                                                        2 * 3 * numNodesPerFace );
            }
          }
        }
      }
    } );
  } );
}

template class MultiphasePoromechanicsConformingFractures<>;
// template class MultiphasePoromechanicsConformingFractures< CompositionalMultiphaseReservoirAndWells<> >;

namespace
{
// typedef MultiphasePoromechanicsConformingFractures< CompositionalMultiphaseReservoirAndWells<> > MultiphaseReservoirPoromechanicsConformingFractures;
// REGISTER_CATALOG_ENTRY( PhysicsSolverBase, MultiphaseReservoirPoromechanicsConformingFractures, string const &, Group * const )
typedef MultiphasePoromechanicsConformingFractures<> MultiphasePoromechanicsConformingFractures;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, MultiphasePoromechanicsConformingFractures, string const &, Group * const )
}

} /* namespace geos */
