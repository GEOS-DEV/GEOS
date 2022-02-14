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

/*
 * SolidMechanicsEmbeddedFractures.cpp
 */

#include "SolidMechanicsEmbeddedFractures.hpp"

#include "SolidMechanicsEFEMKernels.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/ContactSelector.hpp"
#include "constitutive/solid/ElasticIsotropic.hpp"
#include "finiteElement/elementFormulations/FiniteElementBase.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"
#include "mesh/DomainPartition.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

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
  registerWrapper( viewKeyStruct::solidSolverNameString(), &m_solidSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the solid mechanics solver in the rock matrix" );

  registerWrapper( viewKeyStruct::fractureRegionNameString(), &m_fractureRegionName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the fracture region." );

  registerWrapper( viewKeyStruct::contactRelationNameString(), &m_contactRelationName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of contact relation to enforce constraints on fracture boundary." );

  this->getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::FALSE );

}

SolidMechanicsEmbeddedFractures::~SolidMechanicsEmbeddedFractures()
{
  // TODO Auto-generated destructor stub
}

void SolidMechanicsEmbeddedFractures::postProcessInput()
{
  m_solidSolver = &this->getParent().getGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );
}

void SolidMechanicsEmbeddedFractures::registerDataOnMesh( dataRepository::Group & meshBodies )
{
  forMeshTargets( meshBodies, [&] ( string const &,
                                    MeshLevel & mesh,
                                    arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions< EmbeddedSurfaceSubRegion >( regionNames, [&] ( localIndex const, EmbeddedSurfaceSubRegion & subRegion )
    {
      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::dispJumpString() ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        reference().resizeDimension< 1 >( 3 );

      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::deltaDispJumpString() ).
        reference().resizeDimension< 1 >( 3 );

      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::oldDispJumpString() ).
        reference().resizeDimension< 1 >( 3 );

      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::fractureTractionString() ).
        reference().resizeDimension< 1 >( 3 );

      subRegion.registerWrapper< array3d< real64 > >( viewKeyStruct::dTraction_dJumpString() ).
        reference().resizeDimension< 1, 2 >( 3, 3 );
    } );
  } );
}

void SolidMechanicsEmbeddedFractures::initializePostInitialConditionsPreSubGroups()
{
  updateState( this->getGroupByPath< DomainPartition >( "/Problem/domain" ) );
}


void SolidMechanicsEmbeddedFractures::resetStateToBeginningOfStep( DomainPartition & domain )
{
  m_solidSolver->resetStateToBeginningOfStep( domain );

  updateState( domain );
}

void SolidMechanicsEmbeddedFractures::implicitStepSetup( real64 const & time_n,
                                                         real64 const & dt,
                                                         DomainPartition & domain )
{
  m_solidSolver->implicitStepSetup( time_n, dt, domain );
}

void SolidMechanicsEmbeddedFractures::implicitStepComplete( real64 const & time_n,
                                                            real64 const & dt,
                                                            DomainPartition & domain )
{
  m_solidSolver->implicitStepComplete( time_n, dt, domain );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();
    SurfaceElementRegion & region = elemManager.getRegion< SurfaceElementRegion >( m_fractureRegionName );
    EmbeddedSurfaceSubRegion & subRegion = region.getSubRegion< EmbeddedSurfaceSubRegion >( 0 );

    arrayView2d< real64 > oldDispJump = subRegion.getReference< array2d< real64 > >( viewKeyStruct::oldDispJumpString() );
    arrayView2d< real64 const > const dispJump = subRegion.getReference< array2d< real64 > >( viewKeyStruct::dispJumpString() );

    forAll< parallelDevicePolicy<> >( subRegion.size(),
                                      [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      LvArray::tensorOps::copy< 3 >( oldDispJump[k], dispJump[k] );
    } );
  } );

      // Laura print
      arrayView1d< real64 const > const & faceArea = subRegion.getElementArea().toViewConst();
      real64 maxNorm   = 0.0;
      real64 Jmean = 0.0;
      real64 area_slip = 0.0;
      forAll< parallelHostPolicy >( subRegion.size(), [&] ( localIndex const kfe )
      {
        real64 const normD = sqrt(dispJump[kfe][1]*dispJump[kfe][1]+dispJump[kfe][2]*dispJump[kfe][2]);
        if( normD > maxNorm )
        {
          maxNorm = normD;
        }
        real64 const area = faceArea[kfe];
        if( std::abs( normD ) > 1.e-4 )
        {
          area_slip += area;
          Jmean += normD*area;
        }
      } );
      if( std::abs( area_slip ) > 0 )
      {
        Jmean /= area_slip;
      }
      GEOSX_LOG_RANK_0( GEOSX_FMT( " max disp jump {:15.6e}", maxNorm ) );
      GEOSX_LOG_RANK_0( GEOSX_FMT( "mean disp jump {:15.6e}", Jmean ) );
      // end Laura

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
               m_rhs,
               m_solution );

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

  map< string, array1d< string > > meshTargets;
  forMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
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
    meshTargets[meshBodyName] = std::move( regions );
  } );

  dofManager.addField( viewKeyStruct::dispJumpString(),
                       DofManager::Location::Elem,
                       3,
                       meshTargets );

  dofManager.addCoupling( viewKeyStruct::dispJumpString(),
                          viewKeyStruct::dispJumpString(),
                          DofManager::Connector::Elem );
}

void SolidMechanicsEmbeddedFractures::setupSystem( DomainPartition & domain,
                                                   DofManager & dofManager,
                                                   CRSMatrix< real64, globalIndex > & localMatrix,
                                                   ParallelVector & rhs,
                                                   ParallelVector & solution,
                                                   bool const setSparsity )
{
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
  localMatrix.setName( this->getName() + "/localMatrix" );

  rhs.setName( this->getName() + "/rhs" );
  rhs.create( dofManager.numLocalDofs(), MPI_COMM_GEOSX );

  solution.setName( this->getName() + "/solution" );
  solution.create( dofManager.numLocalDofs(), MPI_COMM_GEOSX );
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

  // If specified as a b.c. apply traction
  applyTractionBC( time, dt, domain );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();
    SurfaceElementRegion const & region = elemManager.getRegion< SurfaceElementRegion >( m_fractureRegionName );
    EmbeddedSurfaceSubRegion const & subRegion = region.getSubRegion< EmbeddedSurfaceSubRegion >( 0 );

    string const dispDofKey = dofManager.getKey( dataRepository::keys::TotalDisplacement );
    string const jumpDofKey = dofManager.getKey( viewKeyStruct::dispJumpString() );

    arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
    arrayView1d< globalIndex const > const jumpDofNumber = subRegion.getReference< globalIndex_array >( jumpDofKey );

    real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

    SolidMechanicsEFEMKernels::QuasiStaticFactory kernelFactory( subRegion,
                                                                 dispDofNumber,
                                                                 jumpDofNumber,
                                                                 dofManager.rankOffset(),
                                                                 localMatrix,
                                                                 localRhs,
                                                                 gravityVectorData );

    real64 maxTraction = finiteElement::
                           regionBasedKernelApplication
                         < parallelDevicePolicy< 32 >,
                           constitutive::SolidBase,
                           CellElementSubRegion >( mesh,
                                                   regionNames,
                                                   m_solidSolver->getDiscretizationName(),
                                                   SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString(),
                                                   kernelFactory );

    GEOSX_UNUSED_VAR( maxTraction );

  } );
}

void SolidMechanicsEmbeddedFractures::addCouplingNumNonzeros( DomainPartition & domain,
                                                              DofManager & dofManager,
                                                              arrayView1d< localIndex > const & rowLengths ) const
{
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel const & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager const & nodeManager          = mesh.getNodeManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const jumpDofKey = dofManager.getKey( viewKeyStruct::dispJumpString() );
    string const dispDofKey = dofManager.getKey( keys::TotalDisplacement );

    arrayView1d< globalIndex const > const &
    dispDofNumber =  nodeManager.getReference< globalIndex_array >( dispDofKey );

    globalIndex const rankOffset = dofManager.rankOffset();

    SurfaceElementRegion const & fractureRegion = elemManager.getRegion< SurfaceElementRegion >( getFractureRegionName() );

    EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion = fractureRegion.getSubRegion< EmbeddedSurfaceSubRegion >( 0 );

    arrayView1d< globalIndex const > const &
    jumpDofNumber = embeddedSurfaceSubRegion.getReference< array1d< globalIndex > >( jumpDofKey );

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const, CellElementSubRegion const & cellElementSubRegion )
    {

      SortedArrayView< localIndex const > const fracturedElements = cellElementSubRegion.fracturedElementsList();

      ArrayOfArraysView< localIndex const > const cellsToEmbeddedSurfaces = cellElementSubRegion.embeddedSurfacesList().toViewConst();

      localIndex const numDispDof = 3*cellElementSubRegion.numNodesPerElement();

      for( localIndex ei=0; ei<fracturedElements.size(); ++ei )
      {
        localIndex const cellIndex = fracturedElements[ei];

        localIndex k = cellsToEmbeddedSurfaces[cellIndex][0];
        localIndex const localRow = LvArray::integerConversion< localIndex >( jumpDofNumber[k] - rankOffset );
        if( localRow >= 0 && localRow < rowLengths.size() )
        {
          for( localIndex i=0; i<3; ++i )
          {
            rowLengths[localRow + i] += numDispDof;
          }
        }

        for( localIndex a=0; a<cellElementSubRegion.numNodesPerElement(); ++a )
        {
          const localIndex & node = cellElementSubRegion.nodeList( cellIndex, a );
          localIndex const localDispRow = LvArray::integerConversion< localIndex >( dispDofNumber[node] - rankOffset );

          if( localDispRow >= 0 && localDispRow < rowLengths.size() )
          {
            for( int d=0; d<3; ++d )
            {
              rowLengths[localDispRow + d] += 3;
            }
          }
        }
      }

    } );
  } );
}

void SolidMechanicsEmbeddedFractures::addCouplingSparsityPattern( DomainPartition const & domain,
                                                                  DofManager const & dofManager,
                                                                  SparsityPatternView< globalIndex > const & pattern ) const
{
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel const & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager const & nodeManager          = mesh.getNodeManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const jumpDofKey = dofManager.getKey( viewKeyStruct::dispJumpString() );
    string const dispDofKey = dofManager.getKey( keys::TotalDisplacement );

    arrayView1d< globalIndex const > const &
    dispDofNumber =  nodeManager.getReference< globalIndex_array >( dispDofKey );

    globalIndex const rankOffset = dofManager.rankOffset();

    SurfaceElementRegion const & fractureRegion = elemManager.getRegion< SurfaceElementRegion >( getFractureRegionName() );

    EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion = fractureRegion.getSubRegion< EmbeddedSurfaceSubRegion >( 0 );

    arrayView1d< globalIndex const > const &
    jumpDofNumber = embeddedSurfaceSubRegion.getReference< array1d< globalIndex > >( jumpDofKey );

    static constexpr int maxNumDispDof = 3 * 8;

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                CellElementSubRegion const & cellElementSubRegion )
    {

      SortedArrayView< localIndex const > const fracturedElements = cellElementSubRegion.fracturedElementsList();

      ArrayOfArraysView< localIndex const > const cellsToEmbeddedSurfaces = cellElementSubRegion.embeddedSurfacesList().toViewConst();

      localIndex const numDispDof = 3*cellElementSubRegion.numNodesPerElement();

      for( localIndex ei=0; ei<fracturedElements.size(); ++ei )
      {
        localIndex const cellIndex = fracturedElements[ei];
        localIndex const k = cellsToEmbeddedSurfaces[cellIndex][0];

        // working arrays
        stackArray1d< globalIndex, maxNumDispDof > eqnRowIndicesDisp ( numDispDof );
        stackArray1d< globalIndex, 3 > eqnRowIndicesJump( 3 );
        stackArray1d< globalIndex, maxNumDispDof > dofColIndicesDisp ( numDispDof );
        stackArray1d< globalIndex, 3 > dofColIndicesJump( 3 );

        for( localIndex idof = 0; idof < 3; ++idof )
        {
          eqnRowIndicesJump[idof] = jumpDofNumber[k] + idof - rankOffset;
          dofColIndicesJump[idof] = jumpDofNumber[k] + idof;
        }

        for( localIndex a=0; a<cellElementSubRegion.numNodesPerElement(); ++a )
        {
          const localIndex & node = cellElementSubRegion.nodeList( cellIndex, a );
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

void SolidMechanicsEmbeddedFractures::applyTractionBC( real64 const time_n,
                                                       real64 const dt,
                                                       DomainPartition & domain )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  fsManager.apply( time_n+ dt,
                   domain,
                   "ElementRegions",
                   viewKeyStruct::fractureTractionString(),
                   [&] ( FieldSpecificationBase const & fs,
                         string const &,
                         SortedArrayView< localIndex const > const & targetSet,
                         Group & subRegion,
                         string const & )
  {
    fs.applyFieldValue< FieldSpecificationEqual, parallelHostPolicy >( targetSet,
                                                                       time_n+dt,
                                                                       subRegion,
                                                                       viewKeyStruct::fractureTractionString() );
  } );

}

real64 SolidMechanicsEmbeddedFractures::calculateResidualNorm( DomainPartition const & domain,
                                                               DofManager const & dofManager,
                                                               arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  // Matrix residual
  real64 const solidResidualNorm = m_solidSolver->calculateResidualNorm( domain, dofManager, localRhs );

  RAJA::ReduceSum< parallelDeviceReduce, real64 > localSum( 0.0 );

  // globalResidualNorm[0]: the sum of all the local sum(rhs^2).
  // globalResidualNorm[1]: max of max force of each rank. Basically max force globally
  real64 globalResidualNorm[2] = {0, 0};

  // Fracture residual
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel const & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    string const jumpDofKey = dofManager.getKey( viewKeyStruct::dispJumpString() );

    globalIndex const rankOffset = dofManager.rankOffset();

    mesh.getElemManager().forElementSubRegions< EmbeddedSurfaceSubRegion >( regionNames, [&]( localIndex const,
                                                                                              EmbeddedSurfaceSubRegion const & subRegion )
    {
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
    std::cout << GEOSX_FMT( "( RFracture ) = ( {:4.2e} ) ; ", fractureResidualNorm );
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

  dofManager.addVectorToField( localSolution, viewKeyStruct::dispJumpString(), viewKeyStruct::deltaDispJumpString(), -scalingFactor );

  dofManager.addVectorToField( localSolution, viewKeyStruct::dispJumpString(), viewKeyStruct::dispJumpString(), -scalingFactor );


  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & )
  {
    std::map< string, string_array > fieldNames;
    fieldNames["elems"].emplace_back( string( viewKeyStruct::dispJumpString() ) );
    fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaDispJumpString() ) );
    CommunicationTools::getInstance().synchronizeFields( fieldNames,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );
}

void SolidMechanicsEmbeddedFractures::updateState( DomainPartition & domain )
{

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();


    elemManager.forElementSubRegions< EmbeddedSurfaceSubRegion >( regionNames, [&]( localIndex const,
                                                                                    EmbeddedSurfaceSubRegion & subRegion )
    {
      ContactBase const & contact = getConstitutiveModel< ContactBase >( subRegion, m_contactRelationName );

      arrayView2d< real64 const > const & jump  =
        subRegion.getReference< array2d< real64 > >( viewKeyStruct::dispJumpString() );

      arrayView2d< real64 const > const & oldJump  =
        subRegion.getReference< array2d< real64 > >( viewKeyStruct::oldDispJumpString() );

      arrayView2d< real64 > const & fractureTraction =
        subRegion.getReference< array2d< real64 > >( viewKeyStruct::fractureTractionString() );

      arrayView3d< real64 > const & dFractureTraction_dJump =
        subRegion.getReference< array3d< real64 > >( viewKeyStruct::dTraction_dJumpString() );

      constitutiveUpdatePassThru( contact, [&] ( auto & castedContact )
      {
        using ContactType = TYPEOFREF( castedContact );
        typename ContactType::KernelWrapper contactWrapper = castedContact.createKernelWrapper();

        SolidMechanicsEFEMKernels::StateUpdateKernel::
          launch< parallelDevicePolicy<> >( subRegion.size(),
                                            contactWrapper,
                                            oldJump,
                                            jump,
                                            fractureTraction,
                                            dFractureTraction_dJump );
      } );
    } );
  } );
}


REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsEmbeddedFractures, string const &, Group * const )
} /* namespace geosx */
