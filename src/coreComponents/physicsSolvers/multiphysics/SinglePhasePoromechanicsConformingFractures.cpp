/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhasePoromechanicsConformingFractures.cpp
 */

#include "SinglePhasePoromechanicsConformingFractures.hpp"

#include "constitutive/solid/PorousSolid.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "linearAlgebra/solvers/BlockPreconditioner.hpp"
#include "linearAlgebra/solvers/SeparateComponentPreconditioner.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/multiphysics/SinglePhasePoromechanicsKernel.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"

namespace geosx
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;

SinglePhasePoromechanicsConformingFractures::SinglePhasePoromechanicsConformingFractures( const string & name,
                                                                Group * const parent )
  : Base( name, parent )
{
  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanics;
  m_linearSolverParameters.get().mgr.separateComponents = true;
  m_linearSolverParameters.get().mgr.displacementFieldName = solidMechanics::totalDisplacement::key();
  m_linearSolverParameters.get().dofsPerNode = 3;
}

void SinglePhasePoromechanicsConformingFractures::registerDataOnMesh( Group & meshBodies )
{
  SolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      subRegion.registerWrapper< string >( viewKeyStruct::porousMaterialNamesString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE ).
        setSizedFromParent( 0 );
    } );
  } );
}

void SinglePhasePoromechanicsConformingFractures::setupCoupling( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                                                   DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;

  /// We need to add 2 coupling terms:

  // 1. Poroemechanical coupling in the bulk 
  dofManager.addCoupling( solidMechanics::totalDisplacement::key(),
                          SinglePhaseBase::viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Elem );

  // 2. Traction - pressure coupling in the fracture                        
  dofManager.addCoupling( SinglePhaseBase::viewKeyStruct::elemDofFieldString(),
                          fields::contact::traction::key(),
                          DofManager::Connector::Elem );

  dofManager.printFieldInfo();                        
}

bool SinglePhasePoromechanicsConformingFractures::updateConfiguration( DomainPartition & domain )
{
  return contactSolver()->updateConfiguration( domain );
}

bool SinglePhasePoromechanicsConformingFractures::resetConfigurationToDefault( DomainPartition & domain ) const
{
  return contactSolver()->resetConfigurationToDefault( domain );
}

void SinglePhasePoromechanicsConformingFractures::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elementRegionManager = mesh.getElemManager();
    elementRegionManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                                       [&]( localIndex const,
                                                                            ElementSubRegionBase & subRegion )
    {
      string & porousName = subRegion.getReference< string >( viewKeyStruct::porousMaterialNamesString() );
      porousName = getConstitutiveName< CoupledSolidBase >( subRegion );
      GEOSX_ERROR_IF( porousName.empty(), GEOSX_FMT( "Solid model not found on subregion {}", subRegion.getName() ) );
    } );
  } );
}

void SinglePhasePoromechanicsConformingFractures::setupSystem( DomainPartition & domain,
                                                                 DofManager & dofManager,
                                                                 CRSMatrix< real64, globalIndex > & localMatrix,
                                                                 ParallelVector & rhs,
                                                                 ParallelVector & solution,
                                                                 bool const setSparsity )
{
  GEOSX_MARK_FUNCTION;

  GEOSX_UNUSED_VAR( setSparsity );
  
  /// 1. Add all coupling terms handled directly by the DofManager
  dofManager.setDomain( domain );
  setupDofs( domain, dofManager );
  dofManager.reorderByRank();

  /// 2. Add coupling terms not added by the DofManager.
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
  rhs.create( numLocalRows, MPI_COMM_GEOSX );

  solution.setName( this->getName() + "/solution" );
  solution.create( numLocalRows, MPI_COMM_GEOSX );

  setUpDflux_dApertureMatrix( domain, dofManager, localMatrix );

  if( !m_precond && m_linearSolverParameters.get().solverType != LinearSolverParameters::SolverType::direct )
  {
    createPreconditioner( domain );
  }
}

void SinglePhasePoromechanicsConformingFractures::assembleSystem( real64 const time_n,
                                                                    real64 const dt,
                                                                    DomainPartition & domain,
                                                                    DofManager const & dofManager,
                                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                    arrayView1d< real64 > const & localRhs )
{

   GEOSX_MARK_FUNCTION;

  // Need to synchronize the two iteration counters
  // contactSolver()->getNonlinearSolverParameters().m_numNewtonIterations = m_nonlinearSolverParameters.m_numNewtonIterations;

  contactSolver()->assembleSystem(time_n, dt, domain, dofManager, localMatrix, localRhs );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();

    string const dofKey = dofManager.getKey( fields::solidMechanics::totalDisplacement::key() );
    arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dofKey );

    string const pDofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

    real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

    poromechanicsKernels::SinglePhaseKernelFactory kernelFactory( dispDofNumber,
                                                                  pDofKey,
                                                                  dofManager.rankOffset(),
                                                                  localMatrix,
                                                                  localRhs,
                                                                  gravityVectorData,
                                                                  FlowSolverBase::viewKeyStruct::fluidNamesString() );

    // Cell-based contributions
    finiteElement::
      regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                    constitutive::PorousSolidBase,
                                    CellElementSubRegion >( mesh,
                                                            regionNames,
                                                            contactSolver()->getSolidSolver()->getDiscretizationName(),
                                                            viewKeyStruct::porousMaterialNamesString(),
                                                            kernelFactory );

    updateOpeningForFlow( domain );
    assembleForceResidualDerivativeWrtPressure( domain, dofManager, localMatrix, localRhs );
    assembleFluidMassResidualDerivativeWrtDisplacement( domain, dofManager, localMatrix, localRhs );                                                        

  } );

  // Transmissibility 3D/2D
  flowSolver()->assemblePoroelasticFluxTerms( time_n,
                                              dt,
                                              domain,
                                              dofManager,
                                              localMatrix,
                                              localRhs,
                                              getDerivativeFluxResidual_dAperture() );
}



void SinglePhasePoromechanicsConformingFractures::
setUpDflux_dApertureMatrix( DomainPartition & domain,
                            DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                            CRSMatrix< real64, globalIndex > & localMatrix )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    std::unique_ptr< CRSMatrix< real64, localIndex > > & derivativeFluxResidual_dAperture = this->getRefDerivativeFluxResidual_dAperture();

    {
      localIndex numRows = 0;
                                                                          auto const & elementSubRegion )
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
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_flowSolver->getDiscretizationName() );

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

void SinglePhasePoromechanicsConformingFractures::
  addTransmissibilityCouplingNNZ( DomainPartition const & domain,
                                  DofManager const & dofManager,
                                  arrayView1d< localIndex > const & rowLengths ) const
{
  GEOSX_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &, //  meshBodyName,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & ) // regionNames
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const presDofKey = dofManager.getKey( m_pressureKey );

    globalIndex const rankOffset = dofManager.rankOffset();

    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & stabilizationMethod = fvManager.getFluxApproximation( m_stabilizationName );

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
              // so we only add the coupling with the nodal displacements of the neighbours.
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

void SinglePhasePoromechanicsConformingFractures::
  addTransmissibilityCouplingPattern( DomainPartition const & domain,
                                      DofManager const & dofManager,
                                      SparsityPatternView< globalIndex > const & pattern ) const
{
  GEOSX_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
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
    FluxApproximationBase const & stabilizationMethod = fvManager.getFluxApproximation( m_stabilizationName );

    // Form the SurfaceGenerator, get the fracture name and use it to retrieve the faceMap (from fracture element to face)
    SurfaceGenerator const &
    surfaceGenerator = this->getParent().getGroup< SurfaceGenerator >( "SurfaceGen" );
    SurfaceElementRegion const & fractureRegion = elemManager.getRegion< SurfaceElementRegion >( surfaceGenerator.getFractureRegionName() );
    FaceElementSubRegion const & fractureSubRegion = fractureRegion.getSubRegion< FaceElementSubRegion >( "faceElementSubRegion" );
    GEOSX_ERROR_IF( !fractureSubRegion.hasWrapper( flow::pressure::key() ), "The fracture subregion must contain pressure field." );
    arrayView2d< localIndex const > const faceMap = fractureSubRegion.faceList();
    GEOSX_ERROR_IF( faceMap.size( 1 ) != 2, "A fracture face has to be shared by two cells." );

    arrayView1d< globalIndex const > const &
    presDofNumber = fractureSubRegion.getReference< globalIndex_array >( presDofKey );

    arrayView2d< localIndex const > const & elemsToFaces = fractureSubRegion.faceList();

    stabilizationMethod.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
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
            globalIndex const rowIndex = presDofNumber[sei[iconn][1-kf]];

            // Get fracture, face and region/subregion/element indices (for elements on both sides)
            localIndex fractureIndex = sei[iconn][kf];

            // Get the number of nodes
            localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[fractureIndex][0] );

            // Loop over the two sides of each fracture element
            for( localIndex kf1 = 0; kf1 < 2; ++kf1 )
            {
              localIndex faceIndex = faceMap[fractureIndex][kf1];

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
      } );
    } );
  } );
}

void SinglePhasePoromechanicsConformingFractures::
  assembleForceResidualDerivativeWrtPressure( DomainPartition & domain,
                                              DofManager const & dofManager,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;
//  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    NodeManager & nodeManager = mesh.getNodeManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();
    arrayView2d< real64 const > const & faceNormal = faceManager.faceNormal();

    string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
    string const & presDofKey = dofManager.getKey( m_pressureKey );

    arrayView1d< globalIndex > const &
    dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
    globalIndex const rankOffset = dofManager.rankOffset();

    // Get the coordinates for all nodes
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

    elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const & subRegion )
    {
      if( subRegion.hasWrapper( flow::pressure::key() ) )
      {
        arrayView1d< globalIndex const > const &
        presDofNumber = subRegion.getReference< globalIndex_array >( presDofKey );
        arrayView1d< real64 const > const & pressure = subRegion.getReference< array1d< real64 > >( flow::pressure::key() );
        arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();

        forAll< serialPolicy >( subRegion.size(), [=]( localIndex const kfe )
        {
          localIndex const kf0 = elemsToFaces[kfe][0];
          localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );

          array1d< real64 > Nbar( 3 );
          Nbar[ 0 ] = faceNormal[elemsToFaces[kfe][0]][0] - faceNormal[elemsToFaces[kfe][1]][0];
          Nbar[ 1 ] = faceNormal[elemsToFaces[kfe][0]][1] - faceNormal[elemsToFaces[kfe][1]][1];
          Nbar[ 2 ] = faceNormal[elemsToFaces[kfe][0]][2] - faceNormal[elemsToFaces[kfe][1]][2];
          LvArray::tensorOps::normalize< 3 >( Nbar );

          globalIndex rowDOF[12];
          real64 nodeRHS[12];
          stackArray1d< real64, 12 > dRdP( 3*numNodesPerFace );
          globalIndex colDOF[1];
          colDOF[0] = presDofNumber[kfe];

          for( localIndex kf=0; kf<2; ++kf )
          {
            localIndex const faceIndex = elemsToFaces[kfe][kf];

            for( localIndex a=0; a<numNodesPerFace; ++a )
            {
              // Compute local area contribution for each node
              array1d< real64 > nodalArea;
              contactSolver()->computeFaceNodalArea( nodePosition, faceToNodeMap, elemsToFaces[kfe][kf], nodalArea );

              real64 const nodalForceMag = -( pressure[kfe] ) * nodalArea[a];
              array1d< real64 > globalNodalForce( 3 );
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
      }
    } );
  } );
}

void SinglePhasePoromechanicsConformingFractures::
  assembleFluidMassResidualDerivativeWrtDisplacement( DomainPartition const & domain,
                                                      DofManager const & dofManager,
                                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                      arrayView1d< real64 > const & GEOSX_UNUSED_PARAM( localRhs ) )
{
  GEOSX_MARK_FUNCTION;
//  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    FaceManager const & faceManager = mesh.getFaceManager();
    NodeManager const & nodeManager = mesh.getNodeManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    arrayView2d< real64 const > const & faceNormal = faceManager.faceNormal();
    ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();

    CRSMatrixView< real64 const, localIndex const > const &
    dFluxResidual_dAperture = getDerivativeFluxResidual_dAperture().toViewConst();

    string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
    string const & presDofKey = dofManager.getKey( m_pressureKey );

    arrayView1d< globalIndex const > const &
    dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
    globalIndex const rankOffset = dofManager.rankOffset();

    // Get the coordinates for all nodes
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();
    
    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   FaceElementSubRegion const & subRegion )
    {
      if( subRegion.hasWrapper( flow::pressure::key() ) )
      {
        string const & fluidName = subRegion.getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() );

        SingleFluidBase const & fluid = getConstitutiveModel< SingleFluidBase >( subRegion, fluidName );
        arrayView2d< real64 const > const & density = fluid.density();

        arrayView1d< globalIndex const > const &
        presDofNumber = subRegion.getReference< array1d< globalIndex > >( presDofKey );
        arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();
        arrayView1d< real64 const > const & area = subRegion.getElementArea().toViewConst();

        forAll< serialPolicy >( subRegion.size(), [&]( localIndex const kfe )
        {
          localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );
          globalIndex nodeDOF[24];
          globalIndex elemDOF[1];
          elemDOF[0] = presDofNumber[kfe];

          array1d< real64 > Nbar( 3 );
          Nbar[ 0 ] = faceNormal[elemsToFaces[kfe][0]][0] - faceNormal[elemsToFaces[kfe][1]][0];
          Nbar[ 1 ] = faceNormal[elemsToFaces[kfe][0]][1] - faceNormal[elemsToFaces[kfe][1]][1];
          Nbar[ 2 ] = faceNormal[elemsToFaces[kfe][0]][2] - faceNormal[elemsToFaces[kfe][1]][2];
          LvArray::tensorOps::normalize< 3 >( Nbar );

          stackArray1d< real64, 2*3*4 > dRdU( 2*3*numNodesPerFace );

          // Accumulation derivative
          if( contactSolver()->isElementInOpenState( subRegion, kfe ) )
          {
            for( localIndex kf=0; kf<2; ++kf )
            {
              // Compute local area contribution for each node
              array1d< real64 > nodalArea;
              contactSolver()->computeFaceNodalArea( nodePosition, faceToNodeMap, elemsToFaces[kfe][kf], nodalArea );

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

            localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow,
                                                                      nodeDOF,
                                                                      dRdU.data(),
                                                                      2 * 3 * numNodesPerFace );
          }

          // flux derivative
          bool skipAssembly = true;
          localIndex const numColumns = dFluxResidual_dAperture.numNonZeros( kfe );
          arraySlice1d< localIndex const > const & columns = dFluxResidual_dAperture.getColumns( kfe );
          arraySlice1d< real64 const > const & values = dFluxResidual_dAperture.getEntries( kfe );

          skipAssembly &= !( contactSolver()->isElementInOpenState( subRegion, kfe ) );

          for( localIndex kfe1=0; kfe1<numColumns; ++kfe1 )
          {
            real64 const dR_dAper = values[kfe1];
            localIndex const kfe2 = columns[kfe1];

            skipAssembly &= !( contactSolver()->isElementInOpenState( subRegion, kfe2 ) );

            for( localIndex kf=0; kf<2; ++kf )
            {
              // Compute local area contribution for each node
              array1d< real64 > nodalArea;
              contactSolver()->computeFaceNodalArea( nodePosition, faceToNodeMap, elemsToFaces[kfe2][kf], nodalArea );

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

              localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow,
                                                                        nodeDOF,
                                                                        dRdU.data(),
                                                                        2 * 3 * numNodesPerFace );
            }
          }
        } );
      }
    } );
  } );
}

void SinglePhasePoromechanicsConformingFractures::updateState( DomainPartition & domain )
{
  Base::updateState( domain );
}

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhasePoromechanicsConformingFractures, string const &, Group * const )

} /* namespace geosx */
