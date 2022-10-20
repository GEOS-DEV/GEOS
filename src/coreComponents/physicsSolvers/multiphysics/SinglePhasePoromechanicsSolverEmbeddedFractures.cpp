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
 * @file PoroelasticSolverEmbeddedFractures.cpp
 */

#include "SinglePhasePoromechanicsSolverEmbeddedFractures.hpp"

#include "constitutive/contact/ContactSelector.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "physicsSolvers/contact/SolidMechanicsEFEMKernelsHelper.hpp"
#include "physicsSolvers/contact/SolidMechanicsEmbeddedFractures.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/multiphysics/SinglePhasePoromechanicsEFEMKernel.hpp"
#include "physicsSolvers/multiphysics/SinglePhasePoromechanicsKernel.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsExtrinsicData.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace extrinsicMeshData;

SinglePhasePoromechanicsSolverEmbeddedFractures::SinglePhasePoromechanicsSolverEmbeddedFractures( const std::string & name,
                                                                                                  Group * const parent ):
  SinglePhasePoromechanicsSolver( name, parent ),
  m_fracturesSolverName()
{
  registerWrapper( viewKeyStruct::fracturesSolverNameString(), &m_fracturesSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the fractures solver to use in the fractured poroelastic solver" );

  this->getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::FALSE );

  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanicsEmbeddedFractures;
  m_linearSolverParameters.get().mgr.separateComponents = false;
  m_linearSolverParameters.get().mgr.displacementFieldName = solidMechanics::totalDisplacement::key();
  m_linearSolverParameters.get().dofsPerNode = 3;
}

SinglePhasePoromechanicsSolverEmbeddedFractures::~SinglePhasePoromechanicsSolverEmbeddedFractures()
{}

void SinglePhasePoromechanicsSolverEmbeddedFractures::postProcessInput()
{
  SinglePhasePoromechanicsSolver::postProcessInput();

  m_fracturesSolver  = &this->getParent().getGroup< SolidMechanicsEmbeddedFractures >( m_fracturesSolverName );
}

void SinglePhasePoromechanicsSolverEmbeddedFractures::registerDataOnMesh( dataRepository::Group & meshBodies )
{
  SinglePhasePoromechanicsSolver::registerDataOnMesh( meshBodies );

  using namespace extrinsicMeshData::contact;

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions< EmbeddedSurfaceSubRegion >( regionNames, [&] ( localIndex const,
                                                                                     EmbeddedSurfaceSubRegion & subRegion )
    {
      subRegion.registerExtrinsicData< dTraction_dPressure >( getName() );
    } );
  } );
}

void SinglePhasePoromechanicsSolverEmbeddedFractures::initializePostInitialConditionsPreSubGroups()
{
  updateState( this->getGroupByPath< DomainPartition >( "/Problem/domain" ) );
}

void SinglePhasePoromechanicsSolverEmbeddedFractures::setupDofs( DomainPartition const & domain,
                                                                 DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_fracturesSolver->setupDofs( domain, dofManager );
  flowSolver()->setupDofs( domain, dofManager );

  // Add coupling between displacement and cell pressures
  dofManager.addCoupling( extrinsicMeshData::solidMechanics::totalDisplacement::key(),
                          SinglePhaseBase::viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Elem );

  map< std::pair< string, string >, array1d< string > > meshTargets;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                                MeshLevel const & meshLevel,
                                                                arrayView1d< string const > const & regionNames )
  {
    array1d< string > regions;
    ElementRegionManager const & elementRegionManager = meshLevel.getElemManager();
    elementRegionManager.forElementRegions< SurfaceElementRegion >( regionNames,
                                                                    [&]( localIndex const,
                                                                         SurfaceElementRegion const & region )
    {
      regions.emplace_back( region.getName() );
    } );
    meshTargets[std::make_pair( meshBodyName, meshLevel.getName())] = std::move( regions );
  } );

  dofManager.addCoupling( SinglePhaseBase::viewKeyStruct::elemDofFieldString(),
                          extrinsicMeshData::contact::dispJump::key(),
                          DofManager::Connector::Elem,
                          meshTargets );
}

void SinglePhasePoromechanicsSolverEmbeddedFractures::setupSystem( DomainPartition & domain,
                                                                   DofManager & dofManager,
                                                                   CRSMatrix< real64, globalIndex > & localMatrix,
                                                                   ParallelVector & rhs,
                                                                   ParallelVector & solution,
                                                                   bool const setSparsity )
{
  // Add missing couplings ( matrix pressure with displacement jump and jump - displacement )

  GEOSX_MARK_FUNCTION;

  GEOSX_UNUSED_VAR( setSparsity );

  dofManager.setDomain( domain );
  setupDofs( domain, dofManager );
  dofManager.reorderByRank();

  // Set the sparsity pattern without the Kwu and Kuw blocks.
  SparsityPattern< globalIndex > patternDiag;
  dofManager.setSparsityPattern( patternDiag );

  // Get the original row lengths (diagonal blocks only)
  array1d< localIndex > rowLengths( patternDiag.numRows() );
  for( localIndex localRow = 0; localRow < patternDiag.numRows(); ++localRow )
  {
    rowLengths[localRow] = patternDiag.numNonZeros( localRow );
  }

  // Add the number of nonzeros induced by coupling jump-pm
  addCouplingNumNonzeros( domain, dofManager, rowLengths.toView() );

  // Create a new pattern with enough capacity for coupled matrix
  SparsityPattern< globalIndex > pattern;
  pattern.resizeFromRowCapacities< parallelHostPolicy >( patternDiag.numRows(), patternDiag.numColumns(), rowLengths.data() );

  // Copy the original nonzeros
  for( localIndex localRow = 0; localRow < patternDiag.numRows(); ++localRow )
  {
    globalIndex const * cols = patternDiag.getColumns( localRow ).dataIfContiguous();
    pattern.insertNonZeros( localRow, cols, cols + patternDiag.numNonZeros( localRow ) );
  }

  // Add the nonzeros from coupling
  addCouplingSparsityPattern( domain, dofManager, pattern.toView() );

  // Finally, steal the pattern into a CRS matrix
  localMatrix.setName( this->getName() + "/localMatrix" );
  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );

  rhs.setName( this->getName() + "/rhs" );
  rhs.create( dofManager.numLocalDofs(), MPI_COMM_GEOSX );

  solution.setName( this->getName() + "/solution" );
  solution.create( dofManager.numLocalDofs(), MPI_COMM_GEOSX );
}

void SinglePhasePoromechanicsSolverEmbeddedFractures::addCouplingNumNonzeros( DomainPartition & domain,
                                                                              DofManager & dofManager,
                                                                              arrayView1d< localIndex > const & rowLengths ) const
{
  // 1. Add the number of nonzeros induced by coupling jump-displacement
  m_fracturesSolver->addCouplingNumNonzeros( domain, dofManager, rowLengths );

  // 2. Add the number of nonzeros induced by coupling jump - matrix pressure
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )

  {
    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const jumpDofKey = dofManager.getKey( extrinsicMeshData::contact::dispJump::key() );
    string const pressureDofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

    globalIndex const rankOffset = dofManager.rankOffset();

    elemManager.forElementSubRegions< EmbeddedSurfaceSubRegion >( [&]( EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion )
    {
      localIndex const numEmbeddedElems = embeddedSurfaceSubRegion.size();

      FixedToManyElementRelation const & embeddedSurfacesToCells = embeddedSurfaceSubRegion.getToCellRelation();

      arrayView1d< globalIndex const > const &
      embeddedElementDofNumber = embeddedSurfaceSubRegion.getReference< array1d< globalIndex > >( jumpDofKey );
      arrayView1d< integer const > const & ghostRank = embeddedSurfaceSubRegion.ghostRank();

      for( localIndex k=0; k<numEmbeddedElems; ++k )
      {
        // Get rock matrix element subregion
        CellElementSubRegion const & subRegion =
          elemManager.getRegion( embeddedSurfacesToCells.m_toElementRegion[k][0] ).
            getSubRegion< CellElementSubRegion >( embeddedSurfacesToCells.m_toElementSubRegion[k][0] );

        arrayView1d< globalIndex const > const &
        pressureDofNumber = subRegion.getReference< globalIndex_array >( pressureDofKey );

        localIndex cellElementIndex = embeddedSurfacesToCells.m_toElementIndex[k][0];

        if( ghostRank[k] < 0 )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( embeddedElementDofNumber[k] - rankOffset );
          GEOSX_ASSERT_GE( localRow, 0 );
          GEOSX_ASSERT_GE( rowLengths.size(), localRow + embeddedSurfaceSubRegion.numOfJumpEnrichments()  );

          for( localIndex i=0; i<embeddedSurfaceSubRegion.numOfJumpEnrichments(); ++i )
          {
            rowLengths[localRow + i] += 1;
          }

          localIndex const localPressureRow = LvArray::integerConversion< localIndex >( pressureDofNumber[cellElementIndex] - rankOffset );
          GEOSX_ASSERT_GE( localPressureRow, 0 );
          GEOSX_ASSERT_GE( rowLengths.size(), localPressureRow + embeddedSurfaceSubRegion.numOfJumpEnrichments() );

          rowLengths[ localPressureRow ] += embeddedSurfaceSubRegion.numOfJumpEnrichments();
        }
      }
    } );

    // 3. Add the number of nonzeros induced by coupling jump (aperture) - fracture pressure due to flux term
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( flowSolver()->getDiscretizationName() );

    fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
    {
      for( localIndex iconn=0; iconn<stencil.size(); ++iconn )
      {
        localIndex const numFluxElems = stencil.stencilSize( iconn );
        typename SurfaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
        typename SurfaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
        typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

        EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion =
          elemManager.getRegion( seri[iconn][0] ).getSubRegion< EmbeddedSurfaceSubRegion >( sesri[iconn][0] );

        arrayView1d< globalIndex const > const &
        pressureDofNumber =  embeddedSurfaceSubRegion.getReference< globalIndex_array >( pressureDofKey );

        for( localIndex k0=0; k0<numFluxElems; ++k0 )
        {
          globalIndex const activeFlowDOF = pressureDofNumber[sei[iconn][k0]];
          globalIndex const rowNumber = activeFlowDOF - rankOffset;

          if( rowNumber >= 0 && rowNumber < rowLengths.size() )
          {
            for( localIndex k1=0; k1<numFluxElems; ++k1 )
            {
              // The coupling with the jump of the cell itself has already been added by the dofManager
              // so we only add the coupling with the jumps of the neighbours.
              if( k1 != k0 )
              {
                rowLengths[ rowNumber ] += embeddedSurfaceSubRegion.numOfJumpEnrichments(); // number of jump enrichments.
              }
            }
          }
        }
      }
    } );
  } );
}

void SinglePhasePoromechanicsSolverEmbeddedFractures::addCouplingSparsityPattern( DomainPartition const & domain,
                                                                                  DofManager const & dofManager,
                                                                                  SparsityPatternView< globalIndex > const & pattern ) const
{
  // 1. Add sparsity pattern induced by coupling jump-displacement
  m_fracturesSolver->addCouplingSparsityPattern( domain, dofManager, pattern );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const jumpDofKey = dofManager.getKey( extrinsicMeshData::contact::dispJump::key() );
    string const pressureDofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

    globalIndex const rankOffset = dofManager.rankOffset();

    // 2. Add the sparsity pattern induced by coupling jump - matrix pressure
    elemManager.forElementSubRegions< EmbeddedSurfaceSubRegion >( regionNames, [&]( localIndex const, EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion )
    {
      localIndex const numEmbeddedElems = embeddedSurfaceSubRegion.size();

      FixedToManyElementRelation const & embeddedSurfacesToCells = embeddedSurfaceSubRegion.getToCellRelation();

      arrayView1d< globalIndex const > const &
      jumpDofNumber = embeddedSurfaceSubRegion.getReference< array1d< globalIndex > >( jumpDofKey );
      arrayView1d< integer const > const & ghostRank = embeddedSurfaceSubRegion.ghostRank();

      for( localIndex k=0; k<numEmbeddedElems; ++k )
      {
        // Get rock matrix element subregion
        CellElementSubRegion const & subRegion =
          elemManager.getRegion( embeddedSurfacesToCells.m_toElementRegion[k][0] ).
            getSubRegion< CellElementSubRegion >( embeddedSurfacesToCells.m_toElementSubRegion[k][0] );

        arrayView1d< globalIndex const > const &
        pressureDofNumber = subRegion.getReference< globalIndex_array >( pressureDofKey );

        localIndex cellElementIndex = embeddedSurfacesToCells.m_toElementIndex[k][0];

        if( ghostRank[k] < 0 ) /// TODO is this really necessary?
        {
          localIndex const localJumpRow = LvArray::integerConversion< localIndex >( jumpDofNumber[k] - rankOffset );
          localIndex const localPressureRow = LvArray::integerConversion< localIndex >( pressureDofNumber[cellElementIndex] - rankOffset );

          for( localIndex i=0; i<embeddedSurfaceSubRegion.numOfJumpEnrichments(); ++i )
          {
            if( localJumpRow + i >= 0 && localJumpRow + i < pattern.numRows() )
              pattern.insertNonZero( localJumpRow + i, pressureDofNumber[cellElementIndex] );
            if( localPressureRow >= 0 && localPressureRow < pattern.numRows() )
              pattern.insertNonZero( localPressureRow, jumpDofNumber[k] + i );
          }
        }
      }
    } );

    // 3. Add the sparsity pattern induced by coupling jump (aperture) - fracture pressure due to flux term
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( flowSolver()->getDiscretizationName() );

    fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
    {
      for( localIndex iconn=0; iconn<stencil.size(); ++iconn )
      {
        localIndex const numFluxElems = stencil.stencilSize( iconn );
        typename SurfaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
        typename SurfaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
        typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

        EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion =
          elemManager.getRegion( seri[iconn][0] ).getSubRegion< EmbeddedSurfaceSubRegion >( sesri[iconn][0] );

        arrayView1d< globalIndex const > const &
        pressureDofNumber =  embeddedSurfaceSubRegion.getReference< globalIndex_array >( pressureDofKey );
        arrayView1d< globalIndex const > const &
        jumpDofNumber =  embeddedSurfaceSubRegion.getReference< globalIndex_array >( jumpDofKey );

        for( localIndex k0=0; k0<numFluxElems; ++k0 )
        {
          globalIndex const activeFlowDOF = pressureDofNumber[sei[iconn][k0]];
          globalIndex const rowIndex = activeFlowDOF - rankOffset;

          if( rowIndex >= 0 && rowIndex < pattern.numRows() )
          {
            for( localIndex k1=0; k1<numFluxElems; ++k1 )
            {
              // The coupling with the jump of the cell itself has already been added by the dofManager
              // so we only add the coupling with the jumps of the neighbours.
              if( k1 != k0 )
              {
                for( localIndex i=0; i<embeddedSurfaceSubRegion.numOfJumpEnrichments(); i++ )
                {
                  globalIndex const colIndex = jumpDofNumber[sei[iconn][k1]] + i;
                  pattern.insertNonZero( rowIndex, colIndex );
                }
              }
            }
          }
        }
      }
    } );
  } );

}

void SinglePhasePoromechanicsSolverEmbeddedFractures::assembleSystem( real64 const time_n,
                                                                      real64 const dt,
                                                                      DomainPartition & domain,
                                                                      DofManager const & dofManager,
                                                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                      arrayView1d< real64 > const & localRhs )
{

  GEOSX_MARK_FUNCTION;

  //updateState( domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )

  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();
    SurfaceElementRegion const & region = elemManager.getRegion< SurfaceElementRegion >( m_fracturesSolver->getFractureRegionName() );
    EmbeddedSurfaceSubRegion const & subRegion = region.getSubRegion< EmbeddedSurfaceSubRegion >( 0 );

    string const dofKey = dofManager.getKey( extrinsicMeshData::solidMechanics::totalDisplacement::key() );
    string const jumpDofKey = dofManager.getKey( extrinsicMeshData::contact::dispJump::key() );

    arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dofKey );
    arrayView1d< globalIndex const > const & jumpDofNumber = subRegion.getReference< globalIndex_array >( jumpDofKey );

    string const pDofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

    real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );


    // 1. Cell-based contributions of standard poroelasticity
    poromechanicsKernels::SinglePhaseKernelFactory kernelFactory( dispDofNumber,
                                                                  pDofKey,
                                                                  dofManager.rankOffset(),
                                                                  localMatrix,
                                                                  localRhs,
                                                                  gravityVectorData,
                                                                  FlowSolverBase::viewKeyStruct::fluidNamesString() );

    solidMechanicsSolver()->getMaxForce() =
      finiteElement::
        regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                      constitutive::PorousSolidBase,
                                      CellElementSubRegion >( mesh,
                                                              regionNames,
                                                              solidMechanicsSolver()->getDiscretizationName(),
                                                              viewKeyStruct::porousMaterialNamesString(),
                                                              kernelFactory );

    // 2.  Add EFEM poroelastic contribution
    poromechanicsEFEMKernels::SinglePhaseKernelFactory EFEMkernelFactory( subRegion,
                                                                          dispDofNumber,
                                                                          jumpDofNumber,
                                                                          pDofKey,
                                                                          dofManager.rankOffset(),
                                                                          localMatrix,
                                                                          localRhs,
                                                                          gravityVectorData,
                                                                          FlowSolverBase::viewKeyStruct::fluidNamesString() );

    real64 maxTraction =
      finiteElement::
        regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                      constitutive::PorousSolidBase,
                                      CellElementSubRegion >( mesh,
                                                              regionNames,
                                                              solidMechanicsSolver()->getDiscretizationName(),
                                                              viewKeyStruct::porousMaterialNamesString(),
                                                              EFEMkernelFactory );

    GEOSX_UNUSED_VAR( maxTraction );

    // 3. Assemble poroelastic fluxes and all derivatives
    flowSolver()->assemblePoroelasticFluxTerms( time_n, dt,
                                                domain,
                                                dofManager,
                                                localMatrix,
                                                localRhs,
                                                jumpDofKey );

  } );

}

void SinglePhasePoromechanicsSolverEmbeddedFractures::applyBoundaryConditions( real64 const time_n,
                                                                               real64 const dt,
                                                                               DomainPartition & domain,
                                                                               DofManager const & dofManager,
                                                                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                               arrayView1d< real64 > const & localRhs )
{
  m_fracturesSolver->applyBoundaryConditions( time_n, dt,
                                              domain,
                                              dofManager,
                                              localMatrix,
                                              localRhs );

  flowSolver()->applyBoundaryConditions( time_n, dt,
                                         domain,
                                         dofManager,
                                         localMatrix,
                                         localRhs );
}

void SinglePhasePoromechanicsSolverEmbeddedFractures::implicitStepSetup( real64 const & time_n,
                                                                         real64 const & dt,
                                                                         DomainPartition & domain )
{
  flowSolver()->implicitStepSetup( time_n, dt, domain );
  m_fracturesSolver->implicitStepSetup( time_n, dt, domain );
}

void SinglePhasePoromechanicsSolverEmbeddedFractures::implicitStepComplete( real64 const & time_n,
                                                                            real64 const & dt,
                                                                            DomainPartition & domain )
{
  m_fracturesSolver->implicitStepComplete( time_n, dt, domain );
  flowSolver()->implicitStepComplete( time_n, dt, domain );
}

void SinglePhasePoromechanicsSolverEmbeddedFractures::resetStateToBeginningOfStep( DomainPartition & domain )
{
  flowSolver()->resetStateToBeginningOfStep( domain );
  m_fracturesSolver->resetStateToBeginningOfStep( domain );
}

real64 SinglePhasePoromechanicsSolverEmbeddedFractures::solverStep( real64 const & time_n,
                                                                    real64 const & dt,
                                                                    int const cycleNumber,
                                                                    DomainPartition & domain )
{
  real64 dtReturn = dt;

  /// TODO
  // for (integer outerIter = 0; outerIter < m_maxOuterIter; outerIter++)
  {
    setupSystem( domain,
                 m_dofManager,
                 m_localMatrix,
                 m_rhs,
                 m_solution );

    implicitStepSetup( time_n, dt, domain );

    // Given a fracture state we solve the system
    dtReturn = nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

    implicitStepComplete( time_n, dtReturn, domain );

    // check the fracture state
//    bool fractureStateUnchaged = true; // TODO
//    if ( fractureStateUnChanged )
//    {
//      break;
//    }
  }

  return dtReturn;
}

real64 SinglePhasePoromechanicsSolverEmbeddedFractures::calculateResidualNorm( DomainPartition const & domain,
                                                                               DofManager const & dofManager,
                                                                               arrayView1d< real64 const > const & localRhs )
{
  // compute norm of momentum balance residual equations
  real64 const momentumResidualNorm = m_fracturesSolver->calculateResidualNorm( domain, dofManager, localRhs );

  // compute norm of mass balance residual equations
  real64 const massResidualNorm = flowSolver()->calculateResidualNorm( domain, dofManager, localRhs );

  real64 const residual = sqrt( momentumResidualNorm * momentumResidualNorm + massResidualNorm * massResidualNorm );

  return residual;
}

void SinglePhasePoromechanicsSolverEmbeddedFractures::applySystemSolution( DofManager const & dofManager,
                                                                           arrayView1d< real64 const > const & localSolution,
                                                                           real64 const scalingFactor,
                                                                           DomainPartition & domain )
{
  // update displacement and jump
  m_fracturesSolver->applySystemSolution( dofManager, localSolution, scalingFactor, domain );
  // update pressure field
  flowSolver()->applySystemSolution( dofManager, localSolution, scalingFactor, domain );
}

void SinglePhasePoromechanicsSolverEmbeddedFractures::updateState( DomainPartition & domain )
{
  /// 1. update the reservoir
  SinglePhasePoromechanicsSolver::updateState( domain );

  /// 2. update the fractures
  m_fracturesSolver->updateState( domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< EmbeddedSurfaceSubRegion >( regionNames, [&] ( localIndex const,
                                                                                     auto & subRegion )
    {
      arrayView2d< real64 const > const dispJump =
        subRegion.template getExtrinsicData< extrinsicMeshData::contact::dispJump >();

      arrayView1d< real64 > const aperture = subRegion.getElementAperture();

      arrayView1d< real64 > const hydraulicAperture =
        subRegion.template getExtrinsicData< extrinsicMeshData::flow::hydraulicAperture >();

      arrayView1d< real64 const > const oldHydraulicAperture =
        subRegion.template getExtrinsicData< extrinsicMeshData::flow::aperture0 >();

      arrayView1d< real64 const > const volume = subRegion.getElementVolume();

      arrayView1d< real64 > const deltaVolume =
        subRegion.template getExtrinsicData< extrinsicMeshData::flow::deltaVolume >();

      arrayView1d< real64 const > const area = subRegion.getElementArea().toViewConst();

      arrayView2d< real64 > const & fractureTraction = subRegion.template getExtrinsicData< extrinsicMeshData::contact::traction >();

      arrayView1d< real64 >  const & dTdpf = subRegion.template getExtrinsicData< extrinsicMeshData::contact::dTraction_dPressure >();

      arrayView1d< real64 const > const & pressure =
        subRegion.template getExtrinsicData< extrinsicMeshData::flow::pressure >();

      ContactBase const & contact = getConstitutiveModel< ContactBase >( subRegion, m_fracturesSolver->getContactRelationName() );

      ContactBase::KernelWrapper contactWrapper = contact.createKernelWrapper();

      string const porousSolidName = subRegion.template getReference< string >( FlowSolverBase::viewKeyStruct::solidNamesString() );
      CoupledSolidBase & porousSolid = subRegion.template getConstitutiveModel< CoupledSolidBase >( porousSolidName );

      constitutive::ConstitutivePassThru< CompressibleSolidBase >::execute( porousSolid, [=, &subRegion] ( auto & castedPorousSolid )
      {
        typename TYPEOFREF( castedPorousSolid ) ::KernelWrapper porousMaterialWrapper = castedPorousSolid.createKernelUpdates();

        poromechanicsEFEMKernels::StateUpdateKernel::
          launch< parallelDevicePolicy<> >( subRegion.size(),
                                            contactWrapper,
                                            porousMaterialWrapper,
                                            dispJump,
                                            pressure,
                                            area,
                                            volume,
                                            deltaVolume,
                                            aperture,
                                            oldHydraulicAperture,
                                            hydraulicAperture,
                                            fractureTraction,
                                            dTdpf );

      } );

      // update fracture's permeability and porosity
      flowSolver()->updatePorosityAndPermeability( subRegion );
      // update fluid model
      flowSolver()->updateFluidState( subRegion );

    } );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhasePoromechanicsSolverEmbeddedFractures, std::string const &, Group * const )

} /* namespace geosx */
