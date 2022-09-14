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

#include "common/TimingMacros.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
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
#include "physicsSolvers/contact/SolidMechanicsEFEMKernels.hpp"
#include "physicsSolvers/contact/SolidMechanicsEFEMStaticCondensationKernels.hpp"
#include "physicsSolvers/contact/SolidMechanicsEFEMJumpUpdateKernels.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

SolidMechanicsEmbeddedFractures::SolidMechanicsEmbeddedFractures( const string & name,
                                                                  Group * const parent ):
  ContactSolverBase( name, parent )
{
  registerWrapper( viewKeyStruct::useStaticCondensationString(), &m_useStaticCondensation ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Defines whether to use static condensation or not." );

  this->getWrapper< string >( SolverBase::viewKeyStruct::discretizationString() ).
    setInputFlag( dataRepository::InputFlags::FALSE );


}

SolidMechanicsEmbeddedFractures::~SolidMechanicsEmbeddedFractures()
{
  // TODO Auto-generated destructor stub
}

void SolidMechanicsEmbeddedFractures::postProcessInput()
{
  m_solidSolver = &this->getParent().getGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );

  LinearSolverParameters & linParams = m_linearSolverParameters.get();
  linParams.dofsPerNode = 3;

  if( m_useStaticCondensation )
  {
    linParams.isSymmetric = true;
    linParams.amg.separateComponents = true;
  }
  else
  {
    linParams.mgr.strategy = LinearSolverParameters::MGR::StrategyType::solidMechanicsEmbeddedFractures;
    linParams.mgr.separateComponents = true;
    linParams.mgr.displacementFieldName = keys::TotalDisplacement;
  }
}

void SolidMechanicsEmbeddedFractures::registerDataOnMesh( dataRepository::Group & meshBodies )
{
  ContactSolverBase::registerDataOnMesh( meshBodies );

  using namespace extrinsicMeshData::contact;

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions< EmbeddedSurfaceSubRegion >( regionNames, [&] ( localIndex const, EmbeddedSurfaceSubRegion & subRegion )
    {
      subRegion.registerExtrinsicData< dTraction_dJump >( getName() ).
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

  // reset displacementJump
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< EmbeddedSurfaceSubRegion >( regionNames, [&]( localIndex const,
                                                                                    EmbeddedSurfaceSubRegion & subRegion )
    {
      arrayView2d< real64 > const & jump  =
        subRegion.getExtrinsicData< extrinsicMeshData::contact::dispJump >();

      arrayView2d< real64 const > const & oldJump  =
        subRegion.getExtrinsicData< extrinsicMeshData::contact::oldDispJump >();

      arrayView2d< real64 > const & deltaJump  =
        subRegion.getExtrinsicData< extrinsicMeshData::contact::deltaDispJump >();


      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const kfe )
      {
        for( localIndex i = 0; i < 3; ++i )
        {
          jump( kfe, i ) = oldJump( kfe, i );
          deltaJump( kfe, i ) = 0.0;
        }
      } );
    } );
  } );

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

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();
    SurfaceElementRegion & region = elemManager.getRegion< SurfaceElementRegion >( m_fractureRegionName );
    EmbeddedSurfaceSubRegion & subRegion = region.getSubRegion< EmbeddedSurfaceSubRegion >( 0 );

    arrayView2d< real64 > oldDispJump = subRegion.getExtrinsicData< extrinsicMeshData::contact::oldDispJump >();
    arrayView2d< real64 const > const dispJump = subRegion.getExtrinsicData< extrinsicMeshData::contact::dispJump >();

    forAll< parallelDevicePolicy<> >( subRegion.size(),
                                      [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      LvArray::tensorOps::copy< 3 >( oldDispJump[k], dispJump[k] );
    } );

    // Laura print
    arrayView1d< real64 const > const & faceArea = subRegion.getElementArea().toViewConst();
    real64 maxNorm   = 0.0;
    real64 Jmean = 0.0;
    real64 area_slip = 0.0;
    forAll< parallelHostPolicy >( subRegion.size(), [&] ( localIndex const kfe )
    {
      real64 const normD = sqrt( dispJump[kfe][1]*dispJump[kfe][1]+dispJump[kfe][2]*dispJump[kfe][2] );
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

  } );
}

void SolidMechanicsEmbeddedFractures::setupDofs( DomainPartition const & domain,
                                                 DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->setupDofs( domain, dofManager );

  if( !m_useStaticCondensation )
  {
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

    dofManager.addField( extrinsicMeshData::contact::dispJump::key(),
                         FieldLocation::Elem,
                         3,
                         meshTargets );

    dofManager.addCoupling( extrinsicMeshData::contact::dispJump::key(),
                            extrinsicMeshData::contact::dispJump::key(),
                            DofManager::Connector::Elem );
  }
}
void SolidMechanicsEmbeddedFractures::setupSystem( DomainPartition & domain,
                                                   DofManager & dofManager,
                                                   CRSMatrix< real64, globalIndex > & localMatrix,
                                                   ParallelVector & rhs,
                                                   ParallelVector & solution,
                                                   bool const setSparsity )
{
  GEOSX_MARK_FUNCTION;

  if( !m_useStaticCondensation )
  {

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
  else
  {
    m_solidSolver->setupSystem( domain, dofManager, localMatrix, rhs, solution, setSparsity );
  }
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

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    ElementRegionManager & elemManager = mesh.getElemManager();
    SurfaceElementRegion & region = elemManager.getRegion< SurfaceElementRegion >( m_fractureRegionName );
    EmbeddedSurfaceSubRegion & subRegion = region.getSubRegion< EmbeddedSurfaceSubRegion >( 0 );

    string const dispDofKey = dofManager.getKey( dataRepository::keys::TotalDisplacement );

    arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );

    real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

    if( !m_useStaticCondensation )
    {
      string const jumpDofKey = dofManager.getKey( extrinsicMeshData::contact::dispJump::key() );
      arrayView1d< globalIndex const > const jumpDofNumber = subRegion.getReference< globalIndex_array >( jumpDofKey );

      solidMechanicsEFEMKernels::EFEMFactory kernelFactory( subRegion,
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

    }
    else
    {
      solidMechanicsEFEMKernels::EFEMStaticCondensationFactory kernelFactory( subRegion,
                                                                              dispDofNumber,
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


    }

  } );
}

void SolidMechanicsEmbeddedFractures::addCouplingNumNonzeros( DomainPartition & domain,
                                                              DofManager & dofManager,
                                                              arrayView1d< localIndex > const & rowLengths ) const
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager const & nodeManager          = mesh.getNodeManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const jumpDofKey = dofManager.getKey( extrinsicMeshData::contact::dispJump::key() );
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
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager const & nodeManager          = mesh.getNodeManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const jumpDofKey = dofManager.getKey( extrinsicMeshData::contact::dispJump::key() );
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

void SolidMechanicsEmbeddedFractures::applyTractionBC( real64 const time_n,
                                                       real64 const dt,
                                                       DomainPartition & domain )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & )
  {

    fsManager.apply< ElementSubRegionBase >( time_n+ dt,
                                             mesh,
                                             extrinsicMeshData::contact::traction::key(),
                                             [&] ( FieldSpecificationBase const & fs,
                                                   string const &,
                                                   SortedArrayView< localIndex const > const & targetSet,
                                                   ElementSubRegionBase & subRegion,
                                                   string const & )
    {
      fs.applyFieldValue< FieldSpecificationEqual, parallelHostPolicy >( targetSet,
                                                                         time_n+dt,
                                                                         subRegion,
                                                                         extrinsicMeshData::contact::traction::key() );
    } );
  } );
}

real64 SolidMechanicsEmbeddedFractures::calculateResidualNorm( DomainPartition const & domain,
                                                               DofManager const & dofManager,
                                                               arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  // Matrix residual
  real64 const solidResidualNorm = m_solidSolver->calculateResidualNorm( domain, dofManager, localRhs );

  if( !m_useStaticCondensation )
  {

    string const jumpDofKey = dofManager.getKey( extrinsicMeshData::contact::dispJump::key() );

    globalIndex const rankOffset = dofManager.rankOffset();

    RAJA::ReduceSum< parallelDeviceReduce, real64 > localSum( 0.0 );

    // globalResidualNorm[0]: the sum of all the local sum(rhs^2).
    // globalResidualNorm[1]: max of max force of each rank. Basically max force globally
    real64 globalResidualNorm[2] = {0, 0};

    // Fracture residual
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                  MeshLevel const & mesh,
                                                                  arrayView1d< string const > const & regionNames )
    {
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
  else
  {
    return solidResidualNorm;
  }
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

  if( !m_useStaticCondensation )
  {
    dofManager.addVectorToField( localSolution, extrinsicMeshData::contact::dispJump::key(), extrinsicMeshData::contact::deltaDispJump::key(), scalingFactor );

    dofManager.addVectorToField( localSolution, extrinsicMeshData::contact::dispJump::key(), extrinsicMeshData::contact::dispJump::key(), scalingFactor );
  }
  else
  {
    updateJump( dofManager, domain );
  }

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    FieldIdentifiers fieldsToBeSync;

    fieldsToBeSync.addElementFields( { extrinsicMeshData::contact::dispJump::key(),
                                       extrinsicMeshData::contact::deltaDispJump::key() },
                                     { getFractureRegionName() } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );
}

void SolidMechanicsEmbeddedFractures::updateJump( DofManager const & dofManager,
                                                  DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    ElementRegionManager & elemManager = mesh.getElemManager();
    SurfaceElementRegion & region = elemManager.getRegion< SurfaceElementRegion >( m_fractureRegionName );
    EmbeddedSurfaceSubRegion & subRegion = region.getSubRegion< EmbeddedSurfaceSubRegion >( 0 );

    string const dispDofKey = dofManager.getKey( dataRepository::keys::TotalDisplacement );

    arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );

    real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

    CRSMatrix< real64, globalIndex >  voidMatrix;
    array1d< real64 > voidRhs;

    solidMechanicsEFEMKernels::EFEMJumpUpdateFactory kernelFactory( subRegion,
                                                                    dispDofNumber,
                                                                    dofManager.rankOffset(),
                                                                    voidMatrix.toViewConstSizes(),
                                                                    voidRhs.toView(),
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

void SolidMechanicsEmbeddedFractures::updateState( DomainPartition & domain )
{

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();


    elemManager.forElementSubRegions< EmbeddedSurfaceSubRegion >( regionNames, [&]( localIndex const,
                                                                                    EmbeddedSurfaceSubRegion & subRegion )
    {
      ContactBase const & contact = getConstitutiveModel< ContactBase >( subRegion, m_contactRelationName );

      arrayView2d< real64 const > const & jump = subRegion.getExtrinsicData< extrinsicMeshData::contact::dispJump >();

      arrayView2d< real64 const > const & oldJump = subRegion.getExtrinsicData< extrinsicMeshData::contact::oldDispJump >();

      arrayView2d< real64 > const & fractureTraction = subRegion.getExtrinsicData< extrinsicMeshData::contact::traction >();

      arrayView3d< real64 > const & dFractureTraction_dJump = subRegion.getExtrinsicData< extrinsicMeshData::contact::dTraction_dJump >();

      arrayView1d< integer const > const & fractureState = subRegion.getExtrinsicData< extrinsicMeshData::contact::fractureState >();

      constitutiveUpdatePassThru( contact, [&] ( auto & castedContact )
      {
        using ContactType = TYPEOFREF( castedContact );
        typename ContactType::KernelWrapper contactWrapper = castedContact.createKernelWrapper();

        solidMechanicsEFEMKernels::StateUpdateKernel::
          launch< parallelDevicePolicy<> >( subRegion.size(),
                                            contactWrapper,
                                            oldJump,
                                            jump,
                                            fractureTraction,
                                            dFractureTraction_dJump,
                                            fractureState );
      } );
    } );
  } );
}

bool SolidMechanicsEmbeddedFractures::updateConfiguration( DomainPartition & domain )
{
  int hasConfigurationConverged = true;

  using namespace extrinsicMeshData::contact;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();


    // We want to update the configuration (fracture state in this case) and check if it has changed. If it hasn't changed
    // then we know the configuartion loop has converged and we can return true.
    elemManager.forElementSubRegions< EmbeddedSurfaceSubRegion >( regionNames, [&]( localIndex const,
                                                                                    EmbeddedSurfaceSubRegion & subRegion )
    {
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView2d< real64 const > const & dispJump = subRegion.getExtrinsicData< extrinsicMeshData::contact::dispJump >();
      arrayView2d< real64 const > const & traction = subRegion.getExtrinsicData< extrinsicMeshData::contact::traction >();
      arrayView1d< integer > const & fractureState = subRegion.getExtrinsicData< extrinsicMeshData::contact::fractureState >();

      ContactBase const & contact = getConstitutiveModel< ContactBase >( subRegion, m_contactRelationName );

      constitutiveUpdatePassThru( contact, [&] ( auto & castedContact )
      {
        using ContactType = TYPEOFREF( castedContact );
        typename ContactType::KernelWrapper contactWrapper = castedContact.createKernelWrapper();

        RAJA::ReduceMin< parallelHostReduce, integer > checkActiveSetSub( 1 );

        forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
        {
          if( ghostRank[kfe] < 0 )
          {
            integer const originalFractureState = fractureState[kfe];
            contactWrapper.updateFractureState( kfe, dispJump[kfe], traction[kfe], fractureState[kfe] );
            checkActiveSetSub.min( compareFractureStates( originalFractureState, fractureState[kfe] ) );
          }
        } );
        hasConfigurationConverged &= checkActiveSetSub.get();
      } );
    } );
  } );
  // Need to synchronize the fracture state due to the use will be made of in AssemblyStabilization
  synchronizeFractureState( domain );

  // Compute if globally the fracture state has changed
  int hasConfigurationConvergedGlobally;
  MpiWrapper::allReduce( &hasConfigurationConverged,
                         &hasConfigurationConvergedGlobally,
                         1,
                         MPI_LAND,
                         MPI_COMM_GEOSX );

  // for this solver it makes sense to reset the state.
  // if( !hasConfigurationConvergedGlobally )
  //   resetStateToBeginningOfStep( domain );

  return hasConfigurationConvergedGlobally;
}

REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsEmbeddedFractures, string const &, Group * const )
} /* namespace geosx */
