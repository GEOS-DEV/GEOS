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
 * @file SinglePhasePoromechanicsConformingFractures.cpp
 */

#include "SinglePhasePoromechanicsConformingFractures.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;

template< typename FLOW_SOLVER >
SinglePhasePoromechanicsConformingFractures< FLOW_SOLVER >::SinglePhasePoromechanicsConformingFractures( const string & name,
                                                                                                         Group * const parent )
  : Base( name, parent )
{
  Base::template addLogLevel< logInfo::LinearSolverConfiguration >();
}

template<>
void SinglePhasePoromechanicsConformingFractures<>::setMGRStrategy()
{
  LinearSolverParameters & linearSolverParameters = this->m_linearSolverParameters.get();

  if( linearSolverParameters.preconditionerType != LinearSolverParameters::PreconditionerType::mgr )
    return;

  linearSolverParameters.mgr.separateComponents = true;
  linearSolverParameters.dofsPerNode = 3;

  linearSolverParameters.mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanicsConformingFractures;
  GEOS_LOG_LEVEL_RANK_0( logInfo::LinearSolverConfiguration,
                         GEOS_FMT( "{}: MGR strategy set to {}", getName(),
                                   EnumStrings< LinearSolverParameters::MGR::StrategyType >::toString( linearSolverParameters.mgr.strategy )));
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFractures< FLOW_SOLVER >::
assembleFluidMassResidualDerivativeWrtDisplacement( MeshLevel const & mesh,
                                                    string_array const & regionNames,
                                                    DofManager const & dofManager,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ) )
{
  GEOS_MARK_FUNCTION;

  using namespace contact;

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
  string const & presDofKey = dofManager.getKey( getFlowDofKey() );

  arrayView1d< globalIndex const > const &
  dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
  globalIndex const rankOffset = dofManager.rankOffset();

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                            [&]( localIndex const,
                                                                 FaceElementSubRegion const & subRegion )
  {
    string const & fluidName = subRegion.getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() );

    SingleFluidBase const & fluid = this->template getConstitutiveModel< SingleFluidBase >( subRegion, fluidName );
    arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & density = fluid.density();

    arrayView1d< globalIndex const > const & presDofNumber = subRegion.getReference< array1d< globalIndex > >( presDofKey );

    arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList().toViewConst();
    arrayView1d< real64 const > const & area = subRegion.getElementArea().toViewConst();

    arrayView1d< integer const > const & fractureState = subRegion.getField< contact::fractureState >();

    forAll< serialPolicy >( subRegion.size(), [&]( localIndex const kfe )
    {
      localIndex const kf0 = elemsToFaces[kfe][0], kf1 = elemsToFaces[kfe][1];
      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );
      globalIndex nodeDOF[2*3*m_maxFaceNodes];
      globalIndex elemDOF[1];
      elemDOF[0] = presDofNumber[kfe];

      real64 Nbar[3];
      Nbar[ 0 ] = faceNormal[kf0][0] - faceNormal[kf1][0];
      Nbar[ 1 ] = faceNormal[kf0][1] - faceNormal[kf1][1];
      Nbar[ 2 ] = faceNormal[kf0][2] - faceNormal[kf1][2];
      LvArray::tensorOps::normalize< 3 >( Nbar );

      stackArray1d< real64, 2*3*m_maxFaceNodes > dRdU( 2*3*m_maxFaceNodes );

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
            real64 const dAccumulationResidualdAperture = density[kfe][0] * nodalArea[a];
            for( localIndex i=0; i<3; ++i )
            {
              nodeDOF[ kf*3*numNodesPerFace + 3*a+i ] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe][kf], a )]
                                                        + LvArray::integerConversion< globalIndex >( i );
              real64 const dAper_dU = -pow( -1, kf ) * Nbar[i];
              dRdU( kf*3*numNodesPerFace + 3*a+i ) = dAccumulationResidualdAperture * dAper_dU;
            }
          }
        }

        localIndex const localRow = LvArray::integerConversion< localIndex >( elemDOF[0] - rankOffset );

        if( localRow >= 0 && localRow < localMatrix.numRows() )
        {

          localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow,
                                                                    nodeDOF,
                                                                    dRdU.data(),
                                                                    2 * 3 * numNodesPerFace );
        }
      }

      // flux derivative
      bool skipAssembly = true;
      localIndex const numColumns = dFluxResidual_dNormalJump.numNonZeros( kfe );
      arraySlice1d< localIndex const > const & columns = dFluxResidual_dNormalJump.getColumns( kfe );
      arraySlice1d< real64 const > const & values = dFluxResidual_dNormalJump.getEntries( kfe );

      skipAssembly &= !isFractureOpen;

      for( localIndex kfe1=0; kfe1<numColumns; ++kfe1 )
      {
        real64 const dR_dAper = values[kfe1];
        localIndex const kfe2 = columns[kfe1];

        bool const isOpen = ( fractureState[kfe2] == FractureState::Open );
        skipAssembly &= !isOpen;

        for( localIndex kf=0; kf<2; ++kf )
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
              dRdU( kf*3*numNodesPerFace + 3*a+i ) = dR_dAper * dAper_dU;
            }
          }
        }

        if( !skipAssembly )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( elemDOF[0] - rankOffset );

          if( localRow >= 0 && localRow < localMatrix.numRows() )
          {
            localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow,
                                                                      nodeDOF,
                                                                      dRdU.data(),
                                                                      2 * 3 * numNodesPerFace );
          }
        }
      }
    } );
  } );
}

<<<<<<< HEAD
template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFractures< FLOW_SOLVER >::updateState( DomainPartition & domain )
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

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFractures< FLOW_SOLVER >::updateHydraulicApertureAndFracturePermeability( DomainPartition & domain )
{
  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel & mesh,
                                                                      string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   FaceElementSubRegion & subRegion )
    {
      arrayView2d< real64 const > const dispJump           = subRegion.getField< contact::dispJump >();
      arrayView1d< real64 const > const area               = subRegion.getElementArea();
      arrayView1d< real64 const > const volume             = subRegion.getElementVolume();
      arrayView2d< real64 const > const fractureTraction   = subRegion.getField< contact::traction >();
      arrayView1d< real64 const > const pressure           = subRegion.getField< flow::pressure >();
      arrayView1d< real64 const > const oldHydraulicAperture = subRegion.getField< flow::aperture0 >();

      arrayView1d< real64 > const aperture                 = subRegion.getElementAperture();
      arrayView1d< real64 > const hydraulicAperture        = subRegion.getField< flow::hydraulicAperture >();
      arrayView1d< real64 > const deltaVolume              = subRegion.getField< flow::deltaVolume >();

      string const & porousSolidName = subRegion.getReference< string >( FlowSolverBase::viewKeyStruct::solidNamesString() );
      CoupledSolidBase & porousSolid = subRegion.getConstitutiveModel< CoupledSolidBase >( porousSolidName );

      string const & hydraulicApertureRelationName = subRegion.template getReference< string >( viewKeyStruct::hydraulicApertureRelationNameString()  );
      HydraulicApertureBase const & hydraulicApertureModel = this->template getConstitutiveModel< HydraulicApertureBase >( subRegion, hydraulicApertureRelationName );

      constitutiveUpdatePassThru( hydraulicApertureModel, [&] ( auto & castedHydraulicAperture )
      {
        using HydraulicApertureType = TYPEOFREF( castedHydraulicAperture );
        typename HydraulicApertureType::KernelWrapper hydraulicApertureWrapper = castedHydraulicAperture.createKernelWrapper();

        constitutive::ConstitutivePassThru< CompressibleSolidBase >::execute( porousSolid, [=, &subRegion] ( auto & castedPorousSolid )
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
                                              fractureTraction );

        } );
      } );
    } );
  } );
}
=======
>>>>>>> origin/develop

template class SinglePhasePoromechanicsConformingFractures<>;
template class SinglePhasePoromechanicsConformingFractures< SinglePhaseReservoirAndWells<> >;

namespace
{
typedef SinglePhasePoromechanicsConformingFractures< SinglePhaseReservoirAndWells<> > SinglePhaseReservoirPoromechanicsConformingFractures;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SinglePhaseReservoirPoromechanicsConformingFractures, string const &, Group * const )
typedef SinglePhasePoromechanicsConformingFractures<> SinglePhasePoromechanicsConformingFractures;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SinglePhasePoromechanicsConformingFractures, string const &, Group * const )
}

} /* namespace geos */
