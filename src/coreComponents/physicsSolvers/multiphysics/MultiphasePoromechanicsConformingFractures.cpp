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
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "codingUtilities/Utilities.hpp"


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
assembleForceResidualDerivativeWrtPressure( string const & GEOS_UNUSED_PARAM( meshName ),
                                            MeshLevel const & mesh,
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

template< typename FLOW_SOLVER >
void MultiphasePoromechanicsConformingFractures< FLOW_SOLVER >::
assembleFluidMassResidualDerivativeWrtDisplacement( string const & GEOS_UNUSED_PARAM( meshName ),
                                                    MeshLevel const & mesh,
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

    // vol = hydraulicAperture * area, so dVol/du = nodalArea * dHydraulicAperture_dNormalJump * dNormalJump/du
    arrayView1d< real64 const > const dHydraulicAperture_dNormalJump = subRegion.getField< flow::dHydraulicAperture_dNormalJump >();

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

      localIndex const localRow = LvArray::integerConversion< localIndex >( flowDofNumber[kfe] - rankOffset );

      // Accumulation derivative: zero when the hydraulic aperture does not follow the jump (e.g. closed Barton-Bandis).
      // No fractureState gating: a closed element may still have dHa/dgN != 0 (e.g. HydraulicApertureTable slope at gN < 0)
      if( !isZero( dHydraulicAperture_dNormalJump[kfe] ) )
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
            real64 const dVolume_dNormalJump = nodalArea[a] * dHydraulicAperture_dNormalJump[kfe];
            for( localIndex i=0; i<3; ++i )
            {
              nodeDOF[ kf*3*numNodesPerFace + 3*a+i ] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe][kf], a )]
                                                        + LvArray::integerConversion< globalIndex >( i );
              real64 const dNormalJump_dU = -pow( -1, kf ) * Nbar[i];
              for( integer ic = 0; ic < numComp; ic++ )
              {
                dRdU[ ic ][ kf*3*numNodesPerFace + 3*a + i ] = dVolume_dNormalJump * dNormalJump_dU * compDens[kfe][ic]; // assuming poro=1
              }
            }
          }
        }

        if( localRow >= 0 && localRow < localMatrix.numRows() )
        {
          for( integer ic = 0; ic < numComp; ++ic )
          {
            localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow+ic,
                                                                      nodeDOF,
                                                                      dRdU[ic],
                                                                      2 * 3 * numNodesPerFace );
          }
        }
      }

      // Flux derivative w.r.t. the nodal displacements of each stencil neighbour kfe2, assembled neighbour by
      // neighbour with the geometry (node count, normal, nodal areas) of kfe2:
      // dR/du = dR/dNormalJump(kfe2) * dNormalJump(kfe2)/du
      // (dR/dNormalJump already vanishes where the hydraulic aperture does not follow the jump)
      for( integer ic = 0; ic < numComp; ic++ )
      {
        localIndex const numColumns = dFluxResidual_dNormalJump.numNonZeros( kfe * numComp + ic );
        arraySlice1d< localIndex const > const & columns = dFluxResidual_dNormalJump.getColumns( kfe * numComp + ic );
        arraySlice1d< real64 const > const & values = dFluxResidual_dNormalJump.getEntries( kfe * numComp + ic );

        for( localIndex kfe1 = 0; kfe1 < numColumns; ++kfe1 )
        {
          real64 const dR_dNormalJump = values[kfe1];
          localIndex const kfe2 = columns[kfe1];

          localIndex const kf0_2 = elemsToFaces[kfe2][0], kf1_2 = elemsToFaces[kfe2][1];
          localIndex const numNodesPerFace2 = faceToNodeMap.sizeOfArray( kf0_2 );

          real64 Nbar2[3];
          Nbar2[ 0 ] = faceNormal[kf0_2][0] - faceNormal[kf1_2][0];
          Nbar2[ 1 ] = faceNormal[kf0_2][1] - faceNormal[kf1_2][1];
          Nbar2[ 2 ] = faceNormal[kf0_2][2] - faceNormal[kf1_2][2];
          LvArray::tensorOps::normalize< 3 >( Nbar2 );

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

            for( localIndex a=0; a<numNodesPerFace2; ++a )
            {
              for( localIndex i=0; i<3; ++i )
              {
                nodeDOF[ kf*3*numNodesPerFace2 + 3*a+i ] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe2][kf], a )]
                                                           + LvArray::integerConversion< globalIndex >( i );
                real64 const dNormalJump_dU = -pow( -1, kf ) * Nbar2[i] * ( nodalArea[a] / area[kfe2] );
                dRdU[ ic ][ kf*3*numNodesPerFace2 + 3*a + i ] = dR_dNormalJump * dNormalJump_dU;
              }
            }
          }

          if( localRow >= 0 && localRow < localMatrix.numRows() )
          {
            localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow+ic,
                                                                      nodeDOF,
                                                                      dRdU[ic],
                                                                      2 * 3 * numNodesPerFace2 );
          }
        }
      }
    } );
  } );
}

template class MultiphasePoromechanicsConformingFractures<>;
template class MultiphasePoromechanicsConformingFractures< CompositionalMultiphaseReservoirAndWells<> >;

namespace
{
typedef MultiphasePoromechanicsConformingFractures< CompositionalMultiphaseReservoirAndWells<> > MultiphaseReservoirPoromechanicsConformingFractures;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, MultiphaseReservoirPoromechanicsConformingFractures, string const &, Group * const )
typedef MultiphasePoromechanicsConformingFractures<> MultiphasePoromechanicsConformingFractures;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, MultiphasePoromechanicsConformingFractures, string const &, Group * const )
}

} /* namespace geos */
