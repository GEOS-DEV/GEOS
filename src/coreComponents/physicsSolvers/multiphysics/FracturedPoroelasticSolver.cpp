/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FracturedPoroelasticSolver.cpp
 *
 */


#include "FracturedPoroelasticSolver.hpp"

#include "../solidMechanics/SolidMechanicsPoroElasticKernel.hpp"
#include "common/DataLayouts.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/solid/PoroElastic.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "finiteElement/Kinematics.h"
#include "linearAlgebra/solvers/BlockPreconditioner.hpp"
#include "linearAlgebra/solvers/SeparateComponentPreconditioner.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

FracturedPoroelasticSolver::FracturedPoroelasticSolver( const std::string & name,
                                                        Group * const parent ):
  PoroelasticSolver( name, parent ),
  m_fracturesSolverName()
{
  registerWrapper( viewKeyStruct::fracturesSolverNameString, &m_fracturesSolverName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the fractures solver to use in the fractured poroelastic solver" );
}

FracturedPoroelasticSolver::~FracturedPoroelasticSolver()
{
  // TODO Auto-generated destructor stub
}

void FracturedPoroelasticSolver::PostProcessInput()
{
  PoroelasticSolver::PostProcessInput();

  m_fracturesSolver  = this->getParent()->GetGroup< SolidMechanicsEmbeddedFractures >
  ( m_fracturesSolverName );

  GEOSX_ERROR_IF( m_fracturesSolver == nullptr,
                  "Flow solver not found or invalid type: " << m_fracturesSolverName );
}

void FracturedPoroelasticSolver::SetupDofs( DomainPartition const & domain,
                                            DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_fracturesSolver->SetupDofs( domain, dofManager );
  m_flowSolver->SetupDofs( domain, dofManager );

  // Add coupling between displacement and cell pressures
  dofManager.addCoupling( keys::TotalDisplacement,
                          FlowSolverBase::viewKeyStruct::pressureString,
                          DofManager::Connector::Elem );

  // Add coupling between fracture pressure and displacement jump
  MeshLevel const & meshLevel              = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager const & elemManager = *meshLevel.getElemManager();

  array1d< string > regions;
  elemManager.forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion const & region ) {
    regions.emplace_back( region.getName() );
  } );

  dofManager.addCoupling( FlowSolverBase::viewKeyStruct::pressureString,
                          SolidMechanicsEmbeddedFractures::viewKeyStruct::dispJumpString,
                          DofManager::Connector::Elem,
                          regions );
}

void FracturedPoroelasticSolver::SetupSystem( DomainPartition & domain,
                                              DofManager & dofManager,
                                              CRSMatrix< real64, globalIndex > & localMatrix,
                                              array1d< real64 > & localRhs,
                                              array1d< real64 > & localSolution,
                                              bool const setSparsity )
{
  // Add missing couplings ( matrix pressure with displacement jump and jump - displacement )

  GEOSX_MARK_FUNCTION;

  GEOSX_UNUSED_VAR( setSparsity );

  dofManager.setMesh( domain, 0, 0 );
  SetupDofs( domain, dofManager );
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
  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );
  localRhs.resize( localMatrix.numRows() );
  localSolution.resize( localMatrix.numRows() );

  localMatrix.setName( this->getName() + "/localMatrix" );
  localRhs.setName( this->getName() + "/localRhs" );
  localSolution.setName( this->getName() + "/localSolution" );

  // sets up the object that contains the deravative of the fluxes w.r.t. the aperture ( i.e., the normal jump )
  m_flowSolver->setUpDflux_dApertureMatrix( domain, dofManager, localMatrix );
}

void FracturedPoroelasticSolver::addCouplingNumNonZeros(DomainPartition & domain,
                                                        DofManager & dofManager,
                                                        arrayView1d< localIndex > const & rowLengths ) const
{
  // Add the nonzeros from coupling jump-displacement
  m_fracturesSolver->AddCouplingNumNonzeros( domain, dofManager, rowLengths );

  // Add the number of nonzeros induced by coupling jump-displacement
  MeshLevel const & mesh                   = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager const & nodeManager          = *mesh.getNodeManager();
  ElementRegionManager const & elemManager = *mesh.getElemManager();

  string const jumpDofKey = dofManager.getKey( viewKeyStruct::dispJumpString );
  string const pressureDofKey = dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString );

  globalIndex const rankOffset = dofManager.rankOffset();

  // Coupling jump - matrix pressure
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
      CellElementSubRegion const * const subRegion =
          Group::group_cast< CellElementSubRegion const * const >
      ( elemManager.GetRegion( embeddedSurfacesToCells.m_toElementRegion[k][0] )->
          GetSubRegion( embeddedSurfacesToCells.m_toElementSubRegion[k][0] ));

      arrayView1d< globalIndex const > const &
            pressureDofNumber =  subRegion.getReference< globalIndex_array >( pressureDofKey );

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

  // Coupling jump (aperture) - fracture pressure due to flux term
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_flowSolver->getDiscretization() );

  fluxApprox.forStencils< FaceElementStencil >( mesh, [&]( FaceElementStencil const & stencil )
  {
    for( localIndex iconn=0; iconn<stencil.size(); ++iconn )
    {
      localIndex const numFluxElems = stencil.stencilSize( iconn );
      typename FaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
      typename FaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
      typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

      EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion =
          *elemManager.GetRegion( seri[iconn][0] )->GetSubRegion< EmbeddedSurfaceSubRegion >( sesri[iconn][0] );

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
              rowLengths[ rowNumber ] += 1;
            }
          }
        }
      }
    }
  } );

}

void FracturedPoroelasticSolver::addCouplingSparsityPattern( DomainPartition const & domain,
                                                             DofManager const & dofManager,
                                                             SparsityPatternView< globalIndex > const & pattern ) const
{
  m_fracturesSolver->AddCouplingSparsityPattern( domain, dofManager, pattern );

  MeshLevel const & mesh                   = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager const & nodeManager          = *mesh.getNodeManager();
  ElementRegionManager const & elemManager = *mesh.getElemManager();

  string const jumpDofKey = dofManager.getKey( viewKeyStruct::dispJumpString );
  string const pressureDofKey = dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString );

  globalIndex const rankOffset = dofManager.rankOffset();

  // Coupling jump - matrix pressure
  elemManager.forElementSubRegions< EmbeddedSurfaceSubRegion >( [&]( EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion )
  {
    localIndex const numEmbeddedElems = embeddedSurfaceSubRegion.size();

    FixedToManyElementRelation const & embeddedSurfacesToCells = embeddedSurfaceSubRegion.getToCellRelation();

    arrayView1d< globalIndex const > const &
    jumpDofNumber = embeddedSurfaceSubRegion.getReference< array1d< globalIndex > >( jumpDofKey );
    arrayView1d< integer const > const & ghostRank = embeddedSurfaceSubRegion.ghostRank();

    for( localIndex k=0; k<numEmbeddedElems; ++k )
    {
      // Get rock matrix element subregion
      CellElementSubRegion const * const subRegion =
          Group::group_cast< CellElementSubRegion const * const >
      ( elemManager.GetRegion( embeddedSurfacesToCells.m_toElementRegion[k][0] )->
          GetSubRegion( embeddedSurfacesToCells.m_toElementSubRegion[k][0] ));

      arrayView1d< globalIndex const > const &
            pressureDofNumber =  subRegion.getReference< globalIndex_array >( pressureDofKey );

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
            pattern.insertNonZero( localPressureRow , embeddedElementDofNumber[k] + i );
        }
      }

    }

 } );

  // Coupling fracture pressure - jump (aperture) due to flux term
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_flowSolver->getDiscretization() );

  fluxApprox.forStencils< FaceElementStencil >( mesh, [&]( FaceElementStencil const & stencil )
  {
    for( localIndex iconn=0; iconn<stencil.size(); ++iconn )
    {
      localIndex const numFluxElems = stencil.stencilSize( iconn );
      typename FaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
      typename FaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
      typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

      EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion =
          *elemManager.GetRegion( seri[iconn][0] )->GetSubRegion< EmbeddedSurfaceSubRegion >( sesri[iconn][0] );

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
              globalIndex const colIndex = jumpDofNumber[sei[iconn][k1]];
              pattern.insertNonZero( rowIndex,  colIndex );
            }
          }
        }
      }
    }
  } );

}

void FracturedPoroelasticSolver::AssembleSystem( real64 const time_n,
                                                 real64 const dt,
                                                 DomainPartition & domain,
                                                 DofManager const & dofManager,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs )
{

  // assemble J_SS
//  m_solidSolver->AssembleSystem( time_n, dt,
//                                 domain,
//                                 dofManager,
//                                 localMatrix,
//                                 localRhs );

  m_solidSolver->AssemblyLaunch< constitutive::PoroElasticBase,
                                 SolidMechanicsLagrangianFEMKernels::QuasiStaticPoroElastic >( domain,
                                                                                               dofManager,
                                                                                               localMatrix,
                                                                                               localRhs );

  // assemble J_FF
  m_flowSolver->AssembleSystem( time_n, dt,
                                domain,
                                dofManager,
                                localMatrix,
                                localRhs );

  // assemble J_SF
  AssembleCouplingTerms( domain,
                         dofManager,
                         localMatrix,
                         localRhs );

}

void FracturedPoroelasticSolver::AssembleCouplingTerms( DomainPartition const & domain,
                                                        DofManager const & dofManager,
                                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                        arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

}

void FracturedPoroelasticSolver::ApplyBoundaryConditions( real64 const time_n,
                                                          real64 const dt,
                                                          DomainPartition & domain,
                                                          DofManager const & dofManager,
                                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                          arrayView1d< real64 > const & localRhs )
{
  m_fracturesSolver->ApplyBoundaryConditions( time_n, dt,
                                          domain,
                                          dofManager,
                                          localMatrix,
                                          localRhs );

  m_flowSolver->ApplyBoundaryConditions( time_n, dt,
                                         domain,
                                         dofManager,
                                         localMatrix,
                                         localRhs );
}

void FracturedPoroelasticSolver::ImplicitStepSetup( real64 const & time_n,
                                                    real64 const & dt,
                                                    DomainPartition & domain )
{
  m_flowSolver->ImplicitStepSetup( time_n, dt, domain );
  m_solidSolver->ImplicitStepSetup( time_n, dt, domain );

}

void FracturedPoroelasticSolver::ImplicitStepComplete( real64 const & time_n,
                                                       real64 const & dt,
                                                       DomainPartition & domain )
{
  m_solidSolver->ImplicitStepComplete( time_n, dt, domain );
  m_flowSolver->ImplicitStepComplete( time_n, dt, domain );
}



void FracturedPoroelasticSolver::ResetStateToBeginningOfStep( DomainPartition & domain )
{
  m_flowSolver->ResetStateToBeginningOfStep( domain );
  m_solidSolver->ResetStateToBeginningOfStep( domain );

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegions( mesh, [&] ( localIndex const, ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 const > const & oldTotalMeanStress =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::oldTotalMeanStressString );
    arrayView1d< real64 > const & totalMeanStress =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::totalMeanStressString );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      totalMeanStress[ei] = oldTotalMeanStress[ei];
    } );
  } );
}

real64 FracturedPoroelasticSolver::SolverStep( real64 const & time_n,
                                               real64 const & dt,
                                               int const cycleNumber,
                                               DomainPartition & domain )
{
  real64 dt_return = dt;
  if( m_couplingTypeOption == CouplingTypeOption::SIM_FixedStress )
  {
    GEOSX_ERROR("No sequential coupling for fractured media.");
  }
  else if( m_couplingTypeOption == CouplingTypeOption::FIM )
  {
    SetupSystem( domain,
                 m_dofManager,
                 m_localMatrix,
                 m_localRhs,
                 m_localSolution );

    ImplicitStepSetup( time_n, dt, domain );

    dt_return = NonlinearImplicitStep( time_n, dt, cycleNumber, domain );

    ImplicitStepComplete( time_n, dt_return, domain );
  }
  return dt_return;
}

void FracturedPoroelasticSolver::UpdateDeformationForCoupling( DomainPartition & domain )
{
  PoroelasticSolver::UpdateDeformationForCoupling( domain );

  // update aperture

}

real64 FracturedPoroelasticSolver::CalculateResidualNorm( DomainPartition const & domain,
                                                          DofManager const & dofManager,
                                                          arrayView1d< real64 const > const & localRhs )
{
  // compute norm of momentum balance residual equations
  real64 const momementumResidualNorm = m_solidSolver->CalculateResidualNorm( domain, dofManager, localRhs );

  // compute norm of mass balance residual equations
  real64 const massResidualNorm = m_flowSolver->CalculateResidualNorm( domain, dofManager, localRhs );

  if( getLogLevel() >= 1 && logger::internal::rank==0 )
  {
    char output[200] = {0};
    sprintf( output, "    ( Rsolid, Rfluid ) = ( %4.2e, %4.2e )", momementumResidualNorm, massResidualNorm );
    std::cout << output << std::endl;
  }

  return sqrt( momementumResidualNorm * momementumResidualNorm + massResidualNorm * massResidualNorm );
}

void FracturedPoroelasticSolver::CreatePreconditioner()
{
  if( m_linearSolverParameters.get().preconditionerType == LinearSolverParameters::PreconditionerType::block )
  {
    auto precond = std::make_unique< BlockPreconditioner< LAInterface > >( BlockShapeOption::UpperTriangular,
                                                                           SchurComplementOption::RowsumDiagonalProbing,
                                                                           BlockScalingOption::FrobeniusNorm );

    auto mechPrecond = LAInterface::createPreconditioner( m_solidSolver->getLinearSolverParameters() );
    precond->setupBlock( 0,
                         { { keys::TotalDisplacement, 0, 3 } },
                         std::make_unique< SeparateComponentPreconditioner< LAInterface > >( 3, std::move( mechPrecond ) ) );

    auto flowPrecond = LAInterface::createPreconditioner( m_flowSolver->getLinearSolverParameters() );
    precond->setupBlock( 1,
                         { { SinglePhaseBase::viewKeyStruct::pressureString, 0, 1 } },
                         std::move( flowPrecond ) );

    m_precond = std::move( precond );
  }
  else
  {
    //TODO: Revisit this part such that is coherent across physics solver
    //m_precond = LAInterface::createPreconditioner( m_linearSolverParameters.get() );
  }
}

void FracturedPoroelasticSolver::SolveSystem( DofManager const & dofManager,
                                              ParallelMatrix & matrix,
                                              ParallelVector & rhs,
                                              ParallelVector & solution )
{
  solution.zero();
  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );
}

void FracturedPoroelasticSolver::ApplySystemSolution( DofManager const & dofManager,
                                                      arrayView1d< real64 const > const & localSolution,
                                                      real64 const scalingFactor,
                                                      DomainPartition & domain )
{
  // update displacement field
  m_solidSolver->ApplySystemSolution( dofManager, localSolution, scalingFactor, domain );
  // update pressure field
  m_flowSolver->ApplySystemSolution( dofManager, localSolution, -scalingFactor, domain );
}

REGISTER_CATALOG_ENTRY( SolverBase, FracturedPoroelasticSolver, std::string const &, Group * const )

} /* namespace geosx */
