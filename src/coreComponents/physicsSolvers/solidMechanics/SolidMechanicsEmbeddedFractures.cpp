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

/*
 * SolidMechanicsEmbeddedFractures.cpp
 */

#include "SolidMechanicsEmbeddedFractures.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/ContactRelationBase.hpp"
#include "constitutive/solid/LinearElasticIsotropic.hpp"
#include "finiteElement/elementFormulations/FiniteElementBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"
#include "SolidMechanicsEFEMKernels.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

SolidMechanicsEmbeddedFractures::SolidMechanicsEmbeddedFractures( const string & name,
                                                                  Group * const parent ):
  SolverBase( name, parent ),
  m_solidSolverName(),
  m_fractureRegionName(),
  m_solidSolver( nullptr )
{
  registerWrapper( viewKeyStruct::solidSolverNameString, &m_solidSolverName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the solid mechanics solver in the rock matrix" );

  registerWrapper( viewKeyStruct::fractureRegionNameString, &m_fractureRegionName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the fracture region." );

  registerWrapper( viewKeyStruct::contactRelationNameString, &m_contactRelationName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of contact relation to enforce constraints on fracture boundary." );

  this->getWrapper< string >( viewKeyStruct::discretizationString )->
    setInputFlag( InputFlags::FALSE );

}

SolidMechanicsEmbeddedFractures::~SolidMechanicsEmbeddedFractures()
{
  // TODO Auto-generated destructor stub
}

void SolidMechanicsEmbeddedFractures::registerDataOnMesh( dataRepository::Group * const MeshBodies )
{

  for( auto & mesh : MeshBodies->getSubGroups() )
  {
    MeshLevel * meshLevel = Group::groupCast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    ElementRegionManager * const elemManager = meshLevel->getElemManager();
    {
      elemManager->forElementRegions< SurfaceElementRegion >( [&] ( SurfaceElementRegion & region )
      {
        region.forElementSubRegions< EmbeddedSurfaceSubRegion >( [&]( EmbeddedSurfaceSubRegion & subRegion )
        {
          subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::dispJumpString )->
            setPlotLevel( PlotLevel::LEVEL_0 )->
            reference().resizeDimension< 1 >( 3 );

          subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::deltaDispJumpString )->
            reference().resizeDimension< 1 >( 3 );
        } );
      } );
    }
  }
}


void SolidMechanicsEmbeddedFractures::resetStateToBeginningOfStep( DomainPartition & domain )
{
  m_solidSolver->resetStateToBeginningOfStep( domain );
}

void SolidMechanicsEmbeddedFractures::implicitStepSetup( real64 const & time_n,
                                                         real64 const & dt,
                                                         DomainPartition & domain )
{
  m_solidSolver = this->getParent()->getGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );

  m_solidSolver->implicitStepSetup( time_n, dt, domain );
}

void SolidMechanicsEmbeddedFractures::implicitStepComplete( real64 const & time_n,
                                                            real64 const & dt,
                                                            DomainPartition & domain )
{
  m_solidSolver->implicitStepComplete( time_n, dt, domain );
}

real64 SolidMechanicsEmbeddedFractures::solverStep( real64 const & time_n,
                                                    real64 const & dt,
                                                    int const cycleNumber,
                                                    DomainPartition & domain )
{
  real64 dtReturn = dt;

  implicitStepSetup( time_n,
                     dt,
                     domain );

  setupSystem( domain,
               m_dofManager,
               m_localMatrix,
               m_localRhs,
               m_localSolution );

  // currently the only method is implicit time integration
  dtReturn = this->nonlinearImplicitStep( time_n,
                                          dt,
                                          cycleNumber,
                                          domain );

  // m_solidSolver->updateStress( domain );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  implicitStepComplete( time_n, dtReturn, domain );

  return dtReturn;
}

void SolidMechanicsEmbeddedFractures::setupDofs( DomainPartition const & domain,
                                                 DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->setupDofs( domain, dofManager );

  MeshLevel const & meshLevel              = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager const & elemManager = *meshLevel.getElemManager();

  array1d< string > regions;
  elemManager.forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion const & region ) {
    regions.emplace_back( region.getName() );
  } );

  dofManager.addField( viewKeyStruct::dispJumpString,
                       DofManager::Location::Elem,
                       3,
                       regions );

  dofManager.addCoupling( viewKeyStruct::dispJumpString,
                          viewKeyStruct::dispJumpString,
                          DofManager::Connector::Elem,
                          regions );
}

void SolidMechanicsEmbeddedFractures::setupSystem( DomainPartition & domain,
                                                   DofManager & dofManager,
                                                   CRSMatrix< real64, globalIndex > & localMatrix,
                                                   array1d< real64 > & localRhs,
                                                   array1d< real64 > & localSolution,
                                                   bool const setSparsity )
{
  GEOSX_MARK_FUNCTION;

  GEOSX_UNUSED_VAR( setSparsity );

  dofManager.setMesh( domain, 0, 0 );
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

  // Add the number of nonzeros induced by coupling
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

}

void SolidMechanicsEmbeddedFractures::assembleSystem( real64 const time,
                                                      real64 const dt,
                                                      DomainPartition & domain,
                                                      DofManager const & dofManager,
                                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                      arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  m_solidSolver->assembleSystem( time,
                                 dt,
                                 domain,
                                 dofManager,
                                 localMatrix,
                                 localRhs );

  MeshLevel & mesh = *(domain.getMeshBodies()->getGroup< MeshBody >( 0 )->getMeshLevel( 0 ));

  NodeManager const & nodeManager = *(mesh.getNodeManager());
  ElementRegionManager const & elemManager = *(mesh.getElemManager());
  SurfaceElementRegion const & region = *(elemManager.getRegion< SurfaceElementRegion >( m_fractureRegionName ));
  EmbeddedSurfaceSubRegion const & subRegion = *(region.getSubRegion< EmbeddedSurfaceSubRegion >( 0 ));

  string const dispDofKey = dofManager.getKey( dataRepository::keys::TotalDisplacement );
  string const jumpDofKey = dofManager.getKey( viewKeyStruct::dispJumpString );

  arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
  arrayView1d< globalIndex const > const jumpDofNumber = subRegion.getReference< globalIndex_array >( jumpDofKey );


  real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );


  real64 maxTraction = finiteElement::
                         regionBasedKernelApplication
                       < parallelDevicePolicy< 32 >,
                         constitutive::SolidBase,
                         CellElementSubRegion,
                         SolidMechanicsEFEMKernels::QuasiStatic >( mesh,
                                                                   targetRegionNames(),
                                                                   m_solidSolver->getDiscretizationName(),
                                                                   m_solidSolver->solidMaterialNames(),
                                                                   subRegion,
                                                                   dispDofNumber,
                                                                   jumpDofNumber,
                                                                   dofManager.rankOffset(),
                                                                   localMatrix,
                                                                   localRhs,
                                                                   gravityVectorData );

  GEOSX_UNUSED_VAR( maxTraction );
}

void SolidMechanicsEmbeddedFractures::addCouplingNumNonzeros( DomainPartition & domain,
                                                              DofManager & dofManager,
                                                              arrayView1d< localIndex > const & rowLengths ) const
{
  MeshLevel const & mesh                   = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager const & nodeManager          = *mesh.getNodeManager();
  ElementRegionManager const & elemManager = *mesh.getElemManager();

  string const jumpDofKey = dofManager.getKey( viewKeyStruct::dispJumpString );
  string const dispDofKey = dofManager.getKey( keys::TotalDisplacement );

  arrayView1d< globalIndex const > const &
  dispDofNumber =  nodeManager.getReference< globalIndex_array >( dispDofKey );

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
      CellElementSubRegion const * const subRegion =
        Group::groupCast< CellElementSubRegion const * const >
          ( elemManager.getRegion( embeddedSurfacesToCells.m_toElementRegion[k][0] )->
            getSubRegion( embeddedSurfacesToCells.m_toElementSubRegion[k][0] ));

      localIndex cellElementIndex = embeddedSurfacesToCells.m_toElementIndex[k][0];

      if( ghostRank[k] < 0 )
      {
        localIndex const localRow = LvArray::integerConversion< localIndex >( embeddedElementDofNumber[k] - rankOffset );
        GEOSX_ASSERT_GE( localRow, 0 );
        GEOSX_ASSERT_GE( rowLengths.size(), localRow + embeddedSurfaceSubRegion.numOfJumpEnrichments()  );

        for( localIndex i=0; i<embeddedSurfaceSubRegion.numOfJumpEnrichments(); ++i )
        {
          rowLengths[localRow + i] += 3*subRegion->numNodesPerElement();
        }

        for( localIndex a=0; a<subRegion->numNodesPerElement(); ++a )
        {
          const localIndex & node = subRegion->nodeList( cellElementIndex, a );
          localIndex const localDispRow = LvArray::integerConversion< localIndex >( dispDofNumber[node] - rankOffset );
          GEOSX_ASSERT_GE( localDispRow, 0 );
          GEOSX_ASSERT_GE( rowLengths.size(), localDispRow + 3*subRegion->numNodesPerElement() );

          for( int d=0; d<3; ++d )
          {
            rowLengths[localDispRow + d] += embeddedSurfaceSubRegion.numOfJumpEnrichments();
          }
        }
      }
    }
  } );
}

void SolidMechanicsEmbeddedFractures::addCouplingSparsityPattern( DomainPartition const & domain,
                                                                  DofManager const & dofManager,
                                                                  SparsityPatternView< globalIndex > const & pattern ) const
{
  MeshLevel const & mesh                   = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager const & nodeManager          = *mesh.getNodeManager();
  ElementRegionManager const & elemManager = *mesh.getElemManager();

  string const jumpDofKey = dofManager.getKey( viewKeyStruct::dispJumpString );
  string const dispDofKey = dofManager.getKey( keys::TotalDisplacement );

  arrayView1d< globalIndex const > const &
  dispDofNumber =  nodeManager.getReference< globalIndex_array >( dispDofKey );

  globalIndex const rankOffset = dofManager.rankOffset();

  static constexpr int maxNumDispDof = 3 * 8; // this is hard-coded for now.

  elemManager.forElementSubRegions< EmbeddedSurfaceSubRegion >( [&]( EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion )
  {

    FixedToManyElementRelation const & embeddedSurfacesToCells = embeddedSurfaceSubRegion.getToCellRelation();

    arrayView1d< globalIndex const > const &
    embeddedElementDofNumber = embeddedSurfaceSubRegion.getReference< array1d< globalIndex > >( jumpDofKey );

    // Insert the entries corresponding to jump-disp coupling
    // This will fill K_wu, and K_uw
    for( localIndex k=0; k<embeddedSurfaceSubRegion.size(); ++k )
    {
      CellBlock const * const elemSubRegion = Group::groupCast< CellBlock const * const >( elemManager.getRegion( embeddedSurfacesToCells.m_toElementRegion[k][0] )->
                                                                                             getSubRegion( embeddedSurfacesToCells.m_toElementSubRegion[k][0] ));


      localIndex cellElementIndex = embeddedSurfacesToCells.m_toElementIndex[k][0];

      // working arrays
      stackArray1d< globalIndex, maxNumDispDof > eqnRowIndicesDisp ( 3*elemSubRegion->numNodesPerElement() );
      stackArray1d< globalIndex, 3 > eqnRowIndicesJump( embeddedSurfaceSubRegion.numOfJumpEnrichments() );
      stackArray1d< globalIndex, maxNumDispDof > dofColIndicesDisp ( 3*elemSubRegion->numNodesPerElement() );
      stackArray1d< globalIndex, 3 > dofColIndicesJump( embeddedSurfaceSubRegion.numOfJumpEnrichments() );


      for( localIndex idof = 0; idof < embeddedSurfaceSubRegion.numOfJumpEnrichments(); ++idof )
      {
        eqnRowIndicesJump[idof] = embeddedElementDofNumber[k] + idof - rankOffset;
        dofColIndicesJump[idof] = embeddedElementDofNumber[k] + idof;
      }

      for( localIndex a=0; a<elemSubRegion->numNodesPerElement(); ++a )
      {
        const localIndex & node = elemSubRegion->nodeList( cellElementIndex, a );
        for( localIndex idof = 0; idof < 3; ++idof )
        {
          eqnRowIndicesDisp[3*a + idof] = dispDofNumber[node] + idof - rankOffset;
          dofColIndicesDisp[3*a + idof] = dispDofNumber[node] + idof;
        }
      }

      for( localIndex i = 0; i < eqnRowIndicesDisp.size(); ++i )
      {
        if( eqnRowIndicesDisp[i] >= 0 && eqnRowIndicesDisp[i] < pattern.numRows() )
        {
          for( localIndex j = 0; j < dofColIndicesJump.size(); ++j )
          {
            pattern.insertNonZero( eqnRowIndicesDisp[i], dofColIndicesJump[j] );
          }
        }
      }


      for( localIndex i = 0; i < eqnRowIndicesJump.size(); ++i )
      {
        if( eqnRowIndicesJump[i] >= 0 && eqnRowIndicesJump[i] < pattern.numRows() )
        {
          for( localIndex j=0; j < dofColIndicesDisp.size(); ++j )
          {
            pattern.insertNonZero( eqnRowIndicesJump[i], dofColIndicesDisp[j] );
          }
        }
      }
    }

  } );
}

void SolidMechanicsEmbeddedFractures::applyBoundaryConditions( real64 const time,
                                                               real64 const dt,
                                                               DomainPartition & domain,
                                                               DofManager const & dofManager,
                                                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                               arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  m_solidSolver->applyBoundaryConditions( time,
                                          dt,
                                          domain,
                                          dofManager,
                                          localMatrix,
                                          localRhs );
}

real64 SolidMechanicsEmbeddedFractures::calculateResidualNorm( DomainPartition const & domain,
                                                               DofManager const & dofManager,
                                                               arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  // Matrix residual
  real64 const solidResidualNorm = m_solidSolver->calculateResidualNorm( domain,
                                                                         dofManager,
                                                                         localRhs );
  // Fracture residual
  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  string const jumpDofKey = dofManager.getKey( viewKeyStruct::dispJumpString );

  globalIndex const rankOffset = dofManager.rankOffset();

  RAJA::ReduceSum< parallelDeviceReduce, real64 > localSum( 0.0 );

  // globalResidualNorm[0]: the sum of all the local sum(rhs^2).
  // globalResidualNorm[1]: max of max force of each rank. Basically max force globally
  real64 globalResidualNorm[2] = {0, 0};

  forTargetSubRegions< EmbeddedSurfaceSubRegion >( mesh, [&]( localIndex const targetIndex,
                                                              EmbeddedSurfaceSubRegion const & subRegion )
  {
    GEOSX_UNUSED_VAR( targetIndex );

    arrayView1d< globalIndex const > const &
    dofNumber = subRegion.getReference< array1d< globalIndex > >( jumpDofKey );
    arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();

    forAll< parallelDevicePolicy<> >( subRegion.size(),
                                      [localRhs, localSum, dofNumber, rankOffset, ghostRank] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      if( ghostRank[k] < 0 )
      {
        localIndex const localRow = LvArray::integerConversion< localIndex >( dofNumber[k] - rankOffset );
        for( localIndex i = 0; i < 3; ++i )
        {
          localSum += localRhs[localRow + i] * localRhs[localRow + i];
        }
      }
    } );

    real64 const localResidualNorm[2] = { localSum.get(), m_solidSolver->getMaxForce() };


    int const rank     = MpiWrapper::commRank( MPI_COMM_GEOSX );
    int const numRanks = MpiWrapper::commSize( MPI_COMM_GEOSX );
    array1d< real64 > globalValues( numRanks * 2 );

    // Everything is done on rank 0
    MpiWrapper::gather( localResidualNorm,
                        2,
                        globalValues.data(),
                        2,
                        0,
                        MPI_COMM_GEOSX );

    if( rank==0 )
    {
      for( int r=0; r<numRanks; ++r )
      {
        // sum/max across all ranks
        globalResidualNorm[0] += globalValues[r*2];
        globalResidualNorm[1] = std::max( globalResidualNorm[1], globalValues[r*2+1] );
      }
    }

    MpiWrapper::bcast( globalResidualNorm, 2, 0, MPI_COMM_GEOSX );
  } );

  real64 const fractureResidualNorm = sqrt( globalResidualNorm[0] )/(globalResidualNorm[1]+1);  // the + 1 is for the first
  // time-step when maxForce = 0;

  if( getLogLevel() >= 1 && logger::internal::rank==0 )
  {
    char output[200] = {0};
    sprintf( output,
             "( RFracture ) = (%4.2e) ; ",
             fractureResidualNorm );
    std::cout<<output;
  }

  return sqrt( solidResidualNorm * solidResidualNorm + fractureResidualNorm * fractureResidualNorm );
}

void SolidMechanicsEmbeddedFractures::applySystemSolution( DofManager const & dofManager,
                                                           arrayView1d< real64 const > const & localSolution,
                                                           real64 const scalingFactor,
                                                           DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  m_solidSolver->applySystemSolution( dofManager,
                                      localSolution,
                                      scalingFactor,
                                      domain );

  dofManager.addVectorToField( localSolution, viewKeyStruct::dispJumpString, viewKeyStruct::deltaDispJumpString, -scalingFactor );

  dofManager.addVectorToField( localSolution, viewKeyStruct::dispJumpString, viewKeyStruct::dispJumpString, -scalingFactor );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::dispJumpString ) );
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaDispJumpString ) );

  CommunicationTools::synchronizeFields( fieldNames,
                                         domain.getMeshBody( 0 )->getMeshLevel( 0 ),
                                         domain.getNeighbors(),
                                         true );

}

REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsEmbeddedFractures, string const &, Group * const )
} /* namespace geosx */
