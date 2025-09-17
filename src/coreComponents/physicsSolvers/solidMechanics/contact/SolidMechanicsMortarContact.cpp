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

/*
 * SolidMechanicsMortarContact.cpp
 */

#include "mesh/DomainPartition.hpp"
#include "SolidMechanicsMortarContact.hpp"


#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsConformingContactKernelsBase.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsALMKernels.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsALMKernelsBase.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsALMSimultaneousKernels.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsDisplacementJumpUpdateKernels.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsContactFaceBubbleKernels.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsMortarContactKernels.hpp"
#include "physicsSolvers/solidMechanics/contact/LogLevelsInfo.hpp"
#include "physicsSolvers/solidMechanics/contact/ContactFields.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "constitutive/contact/FrictionSelector.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "denseLinearAlgebra/denseLASolvers.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;

SolidMechanicsMortarContact::SolidMechanicsMortarContact( const string & name,
                                                          Group * const parent ):
  ContactSolverBase( name, parent )
{

  m_faceTypeToMortarFiniteElements[ElementShape::Quadrilateral] = std::make_unique< finiteElement::H1_QuadrilateralFace_Lagrange1_GaussLegendre2 >();
  m_faceTypeToMortarFiniteElements[ElementShape::Triangle] = std::make_unique< finiteElement::H1_TriangleFace_Lagrange1_Gauss4 >();

  registerWrapper( viewKeyStruct::masterString(), &m_masterName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the master surface" );

  registerWrapper( viewKeyStruct::slaveString(), &m_slaveName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the slave surface" );

}

SolidMechanicsMortarContact::~SolidMechanicsMortarContact()
{}

void SolidMechanicsMortarContact::registerDataOnMesh( dataRepository::Group & meshBodies )
{
  GEOS_MARK_FUNCTION;
  ContactSolverBase::registerDataOnMesh( meshBodies );

}

void SolidMechanicsMortarContact::registerMortarDataOnMesh( )
{

  MeshLevel & mesh = *m_mortarSide.at( MortarSide::Slave ).mesh;
  FaceElementSubRegion & subRegion =
    m_mortarSide.at( MortarSide::Slave ).surface->getUniqueSubRegion< FaceElementSubRegion >();

  FaceManager & faceManager = mesh.getFaceManager();

  // Register the total bubble displacement
  faceManager.registerField< contact::totalBubbleDisplacement >( this->getName() ).
    reference().resizeDimension< 1 >( 3 );

  // Register the incremental bubble displacement
  faceManager.registerField< contact::incrementalBubbleDisplacement >( this->getName() ).
    reference().resizeDimension< 1 >( 3 );

  // Register the rotation matrix
  subRegion.registerField< contact::rotationMatrix >( this->getName() ).
    reference().resizeDimension< 1, 2 >( 3, 3 );

   subRegion.registerField< contact::deltaTraction >( getName() ).
    reference().resizeDimension< 1 >( 3 );

   subRegion.registerField< contact::targetIncrementalJump >( getName() ).
    reference().resizeDimension< 1 >( 3 );

}


void SolidMechanicsMortarContact::setupDofs( DomainPartition const & domain,
                                             DofManager & dofManager) const
{
  GEOS_MARK_FUNCTION;

  SolidMechanicsLagrangianFEM::setupDofs( domain, dofManager );

  map< std::pair< string, string >, string_array > meshTargets;
  
  // bubble and tractions are defined on the slave side only
  string_array regions;

  MeshLevel & meshSlave = *m_mortarSide.at(MortarSide::Slave).mesh;
  ElementRegionManager const & elementRegionManager = meshSlave.getElemManager();
  elementRegionManager.forElementRegions< SurfaceElementRegion >([&]( SurfaceElementRegion const & region )
    {
      regions.emplace_back( region.getName() );
    } );

  meshTargets[std::make_pair( m_meshSlaveName, meshSlave.getName())] = std::move( regions );


  dofManager.addField( contact::traction::key(),
                       FieldLocation::Elem,
                       3,
                       meshTargets );
  
  dofManager.addField( contact::totalBubbleDisplacement::key(),
                       FieldLocation::Face,
                       3,
                       meshTargets );


  // Add coupling between bubble
  dofManager.addCoupling( contact::totalBubbleDisplacement::key(),
                          contact::totalBubbleDisplacement::key(),
                          DofManager::Connector::Elem );

  dofManager.addCoupling( contact::traction::key(),
                          contact::traction::key(),
                          DofManager::Connector::Elem );

  dofManager.addCoupling( solidMechanics::totalDisplacement::key(),
                          contact::traction::key(),
                          DofManager::Connector::Elem,
                          meshTargets );

  dofManager.addCoupling( contact::totalBubbleDisplacement::key(),
                          contact::traction::key(),
                          DofManager::Connector::Elem,
                          meshTargets );

  //dofManager.addCoupling( solidMechanics::totalDisplacement::key(),
  //                        contact::totalBubbleDisplacement::key(),
  //                        DofManager::Connector::Elem,
  //                        meshTargets );
}


void SolidMechanicsMortarContact::setupSystem( DomainPartition & domain,
                                                            DofManager & dofManager,
                                                            CRSMatrix< real64, globalIndex > & localMatrix,
                                                            ParallelVector & rhs,
                                                            ParallelVector & solution,
                                                            bool const setSparsity )
{

  ElementShape shapes[2] = { ElementShape::Triangle, ElementShape::Quadrilateral };

  setMortarSurfaces(domain);

  registerMortarDataOnMesh();

  createFaceTypeListMortar(MortarSide::Master);
  createFaceTypeListMortar(MortarSide::Slave);

  std::cout << "Number of master quad cells: " <<  m_faceTypeToElementList[MortarSide::Master][ElementShape::Quadrilateral].size() << std::endl;
  std::cout << "Number of master tri cells: " <<  m_faceTypeToElementList[MortarSide::Master][ElementShape::Triangle].size() << std::endl;
  std::cout << "Number of slave quad cells: " <<  m_faceTypeToElementList[MortarSide::Slave][ElementShape::Quadrilateral].size() << std::endl;
  std::cout << "Number of slave tri cells: " <<  m_faceTypeToElementList[MortarSide::Slave][ElementShape::Triangle].size() << std::endl;

  // create list of bubbles for the slave side
  this->createBubbleCellList();

  // perform rough screen to find potential connections between master and slave faces
  connectivityMapType connectivityMap;
  getConnectivityMap( connectivityMap );

  // preprocess mortar integration data
  computeMortarInterpolation( connectivityMap );

  // setup dofs for the problem
  dofManager.setDomain( domain );
  setupDofs( domain, dofManager);
  dofManager.reorderByRank();

  // Set the sparsity pattern without the Abu and Aub blocks.
  SparsityPattern< globalIndex > patternDiag;
  dofManager.setSparsityPattern( patternDiag );

  // Manually define the sparsity pattern induced by the slave traction <-> master displacements coupling
  // Get the original row lengths (diagonal blocks only)
  array1d< localIndex > rowLengths( patternDiag.numRows() );
  for( localIndex localRow = 0; localRow < patternDiag.numRows(); ++localRow )
  {
    rowLengths[localRow] = patternDiag.numNonZeros( localRow );
  }

  // Add the number of nonzeros induced by coupling
  this->addBubbleCouplingNumNonzeros( dofManager, rowLengths.toView() );

  for (auto slaveShape : shapes)
  {
    for (auto masterShape : shapes)
    {
      this->addMortarCouplingNumNonzeros(dofManager, 
                                         slaveShape,
                                         masterShape,
                                         connectivityMap,
                                         rowLengths.toView());
    }
  }  

  
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
  this->addBubbleCouplingSparsityPattern( dofManager, pattern.toView() );

  for (auto slaveShape : shapes)
  {
    for (auto masterShape : shapes)
    {
      this->addMortarCouplingSparsityPattern( dofManager, 
                                              slaveShape,
                                              masterShape,
                                              pattern.toView() );
    }
  }

  // // Finally, steal the pattern into a CRS matrix
  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );
  localMatrix.setName( this->getName() + "/localMatrix" );

  rhs.setName( this->getName() + "/rhs" );
  rhs.create( dofManager.numLocalDofs(), MPI_COMM_GEOS );

  solution.setName( this->getName() + "/solution" );
  solution.create( dofManager.numLocalDofs(), MPI_COMM_GEOS );


  //Write the matrix pattern to a file for debugging
  ParallelMatrix parallel_matrix;
  parallel_matrix.create( localMatrix.toViewConst(), dofManager.numLocalDofs(), MPI_COMM_GEOS );
  parallel_matrix.write("mortar_sparsity_new.mtx");

  computeRotationMatrices( ); 

  //GEOS_ERROR("Interrupt");
  
  //////////////////////////////////////////////////////////////////////////////////
  
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( setSparsity );

}

void SolidMechanicsMortarContact::implicitStepSetup( real64 const & time_n,
                                                     real64 const & dt,
                                                     DomainPartition & domain )
{

  GEOS_MARK_FUNCTION;
  SolidMechanicsLagrangianFEM::implicitStepSetup( time_n, dt, domain );

 

}

void SolidMechanicsMortarContact::assembleSystem( real64 const time,
                                                               real64 const dt,
                                                               DomainPartition & domain,
                                                               DofManager const & dofManager,
                                                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                               arrayView1d< real64 > const & localRhs )
{

  GEOS_MARK_FUNCTION;
  synchronizeFractureState( domain );

  ParallelMatrix parallel_matrix_0;
  parallel_matrix_0.create( localMatrix.toViewConst(), dofManager.numLocalDofs(), MPI_COMM_GEOS );
  parallel_matrix_0.write("matrix_preMech.mtx");

  SolidMechanicsLagrangianFEM::assembleSystem( time,
                                               dt,
                                               domain,
                                               dofManager,
                                               localMatrix,
                                               localRhs );

  ParallelMatrix parallel_matrix_1;
  parallel_matrix_1.create( localMatrix.toViewConst(), dofManager.numLocalDofs(), MPI_COMM_GEOS );
  parallel_matrix_1.write("matrix_postMech.mtx");

  //GEOS_ERROR("Interrupt");

  assembleBubbles( dt, domain, dofManager, localMatrix, localRhs );

  ParallelMatrix parallel_matrix_2;
  parallel_matrix_2.create( localMatrix.toViewConst(), dofManager.numLocalDofs(), MPI_COMM_GEOS );
  parallel_matrix_2.write("matrix_postBubble.mtx");

  std::cout << "Assembling mortar matrices" << std::endl;

  assembleMortar( dt, dofManager, localMatrix, localRhs );

  ParallelMatrix parallel_matrix_3;
  parallel_matrix_3.create( localMatrix.toViewConst(), dofManager.numLocalDofs(), MPI_COMM_GEOS );
  parallel_matrix_3.write("matrix_postMortar.mtx");

  //GEOS_ERROR("Mortar solver is not implemented yet. ");

}


void SolidMechanicsMortarContact::assembleBubbles( real64 const dt,
                                                   DomainPartition & domain,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs )
{
  
  // Loop for assembling contributes of bubble elements (Abb, Abu, Aub)
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {

    if (meshName != m_meshSlaveName)
    {
      return;
    }
    else
    {
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
    string const & bubbleDofKey = dofManager.getKey( contact::totalBubbleDisplacement::key() );

    arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );

    real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );


    solidMechanicsConformingContactKernels::FaceBubbleFactory kernelFactory( dispDofNumber,
                                                                             bubbleDofNumber,
                                                                             dofManager.rankOffset(),
                                                                             localMatrix,
                                                                             localRhs,
                                                                             dt,
                                                                             gravityVectorData );

    real64 maxTraction = finiteElement::
                           regionBasedKernelApplication
                         < parallelDevicePolicy< >,
                           ElasticIsotropic,
                           CellElementSubRegion >( mesh,
                                                   regionNames,
                                                   getDiscretizationName(),
                                                   SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString(),
                                                   kernelFactory );

    GEOS_UNUSED_VAR( maxTraction );


    }

  } );

  // assembling interface bubble contribution: Abt, Atb, inactive Abb
  assembleMortarBubbles<ElementShape::Triangle>( dofManager, localMatrix, localRhs );
  assembleMortarBubbles<ElementShape::Quadrilateral>( dofManager, localMatrix, localRhs );

}


template< ElementShape shape>
void SolidMechanicsMortarContact::assembleMortarBubbles( DofManager const & dofManager,
                                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                         arrayView1d< real64 > const & localRhs )
  {

    arrayView1d< localIndex const > elementList = m_faceTypeToElementList.at( MortarSide::Slave ).at( shape );
    
    if (elementList.size() == 0) return;
    
    auto const & feType = getFE< shape >();
    localIndex constexpr numQuadraturePoints = feType.numQuadraturePoints;
    localIndex constexpr numNodesPerElem = feType.numNodes;

    MeshLevel & mesh = *m_mortarSide.at( MortarSide::Slave ).mesh;
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();
    FaceElementSubRegion const & surfRegion = m_mortarSide.at(MortarSide::Slave).surface->getUniqueSubRegion< FaceElementSubRegion >();
    arrayView2d< localIndex const > const elemsToFaces = surfRegion.faceList().toViewConst();

    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const coords =  nodeManager.referencePosition();
    ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  
    globalIndex const rankOffset = dofManager.rankOffset();
    string const bubbleDofKey = dofManager.getKey( contact::totalBubbleDisplacement::key() );
    string const tractionDofKey = dofManager.getKey( contact::traction::key() );

    arrayView2d< real64 const > const traction = surfRegion.getField< fields::contact::traction >().toViewConst();
    arrayView2d< real64 const > const bubbleDisplacement = faceManager.getField< fields::contact::totalBubbleDisplacement >().toViewConst();
    arrayView3d< real64 const > const rotationMatrix = surfRegion.getField< fields::contact::rotationMatrix >().toViewConst();

    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );
    arrayView1d< globalIndex const > const tractionDofNumber = surfRegion.getReference< globalIndex_array >( tractionDofKey );

    forAll< parallelHostPolicy >( elementList.size(), [=] ( localIndex const kk )
    {
      localIndex const k = elementList[kk];
      localIndex const kf0 = elemsToFaces[k][0];
      localIndex const kf1 = elemsToFaces[k][1];

      globalIndex bDof1[3], bDof2[3], tDof[3]; 
      real64 bLoc[3], tLoc[3];
      for (localIndex i=0; i<3; ++i)
      {
        bDof1[i] = bubbleDofNumber[kf0]+i;
        bDof2[i] = bubbleDofNumber[kf1]+i;
        tDof[i] = tractionDofNumber[k]+i;
        bLoc[i] = bubbleDisplacement[kf0][i];
        tLoc[i] = traction[k][i];
      }

      real64 X[numNodesPerElem][3];
      for (localIndex a=0; a<numNodesPerElem; ++a)
      {
        localIndex const nodeIndex = faceToNodeMap[kf0][a];
        X[a][0] = coords[nodeIndex][0];
        X[a][1] = coords[nodeIndex][1];
        X[a][2] = coords[nodeIndex][2];
      }

      // compute the scalar abt contribution
      real64 abt = 0.0;
      real64 N[1];
      for (localIndex q = 0; q < numQuadraturePoints; ++q)
      {
         feType.calcBubbleN(q, N);
         real64 detJ = feType.transformedQuadratureWeight(q, X);
         abt += N[0] * detJ;
      }

      real64 R[3][3];
      real64 RtAbt[3][3];
      real64 localAbt[3][3] = {{ }};
      real64 localAtb[3][3] = {{ }};
      real64 rhsB[3], rhsT[3];

      for( int j=0; j<3; ++j )
      {
        localAtb[j][j] = abt;
        for( int i=0; i<3; ++i )
        {
          R[ i ][ j ] = rotationMatrix( k, i, j );
        }
      }

      LvArray::tensorOps::Rij_eq_AkiBkj< 3, 3, 3 >( RtAbt, R, localAtb );
      LvArray::tensorOps::copy< 3, 3 >( localAtb, RtAbt );
      LvArray::tensorOps::transpose< 3, 3 >(localAbt, localAtb);
      LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( rhsB, localAbt, tLoc );
      LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( rhsT, localAtb, bLoc );

      // assemble matrix and rhs contribution of Abt
      for( localIndex i = 0; i < 3; ++i )
      {
        localIndex const localRow = LvArray::integerConversion< localIndex >( bDof1[i] - rankOffset );
        if( localRow >= 0 && localRow < localMatrix.numRows() )
        {
          localMatrix.addToRow< parallelHostAtomic >( localRow,
                                                      tDof,
                                                      localAbt[i],
                                                      3 );
          RAJA::atomicAdd( parallelHostAtomic{}, &localRhs[localRow], rhsB[i] );
        }
      }

      // assemble matrix and rhs contribution of Atb
      for( localIndex i = 0; i < 3; ++i )
      {
        localIndex const localRow = LvArray::integerConversion< localIndex >( tDof[i] - rankOffset );
        if( localRow >= 0 && localRow < localMatrix.numRows() )
        {
          localMatrix.addToRow< parallelHostAtomic >( localRow,
                                                      bDof1,
                                                      localAtb[i],
                                                      3 );
          RAJA::atomicAdd( parallelHostAtomic{}, &localRhs[localRow], rhsT[i] );
        }
      }

      real64 diagVal[1];
      globalIndex colDof[1];
      diagVal[0] = 1.0;

      // assemble diagonal contribution of inactive bubbles (dirichlet bcs)
      for( localIndex i = 0; i < 3; ++i )
      {
        localIndex const localRow = LvArray::integerConversion< localIndex >( bDof2[i] - rankOffset );
        if( localRow >= 0 && localRow < localMatrix.numRows() )
        {
          colDof[0] = bDof2[i];
          localMatrix.addToRow< parallelHostAtomic >( localRow,
                                                      colDof,
                                                      diagVal,
                                                      1 );
        }
      }


    });

  }


void SolidMechanicsMortarContact::assembleMortar( real64 const dt,
                                                  DofManager const & dofManager,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs )
{
  // Call the templated version of assembleMortar with the appropriate element shapes
  assembleMortar< ElementShape::Triangle, ElementShape::Triangle >( dt, dofManager, localMatrix, localRhs );
  assembleMortar< ElementShape::Triangle, ElementShape::Quadrilateral >( dt, dofManager, localMatrix, localRhs );
  assembleMortar< ElementShape::Quadrilateral, ElementShape::Triangle >( dt, dofManager, localMatrix, localRhs );
  assembleMortar< ElementShape::Quadrilateral, ElementShape::Quadrilateral >( dt, dofManager, localMatrix, localRhs );
}

template< ElementShape shpS, ElementShape shpM>
void SolidMechanicsMortarContact::assembleMortar( real64 const dt,
                                                  DofManager const & dofManager,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs )
{

  if (m_triCellsDet.find(std::make_pair(shpS, shpM)) == m_triCellsDet.end())
  {
    // no connections for the given element shapes
    std::cout << "No connection found for the slave element" << std::endl;
    return;
  }

  MeshLevel & meshSlave = *m_mortarSide.at(MortarSide::Slave).mesh;
  FaceElementSubRegion & subRegionSlave = 
    m_mortarSide.at(MortarSide::Slave).surface->getUniqueSubRegion< FaceElementSubRegion >();
  NodeManager const & nodeManagerSlave = meshSlave.getNodeManager();
  //FaceManager const & faceManagerSlave = meshSlave.getFaceManager();

  MeshLevel & meshMaster = *m_mortarSide.at(MortarSide::Master).mesh;
  FaceElementSubRegion const & subRegionMaster =
    m_mortarSide.at(MortarSide::Master).surface->getUniqueSubRegion< FaceElementSubRegion >();
  NodeManager const & nodeManagerMaster = meshMaster.getNodeManager();
  //FaceManager const & faceManagerMaster = meshMaster.getFaceManager();

  string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
  string const & tractionDofKey = dofManager.getKey( contact::traction::key() );

  arrayView1d< globalIndex const > const dispDofNumberSlave = nodeManagerSlave.getReference< globalIndex_array >( dispDofKey );
  arrayView1d< globalIndex const > const dispDofNumberMaster = nodeManagerMaster.getReference< globalIndex_array >( dispDofKey );

  finiteElement::FiniteElementBase const & feSlave = *(m_faceTypeToMortarFiniteElements.at( shpS ));
  finiteElement::FiniteElementBase const & feMaster = *(m_faceTypeToMortarFiniteElements.at( shpM ));

  arrayView1d< real64 const> detPlus = m_triCellsDet.at(std::make_pair(shpS, shpM)).toView();
  // arrayView1d< real64 > detMinus =  m_triCellsDet.at(std::make_pair(shpS, shpM)).toView();
  localIndex const nTri = detPlus.size();
  array1d< real64 > detMinus( nTri );

  for (localIndex i=0; i<nTri; ++i)
  {
    detMinus[i] = -detPlus[i];
  }

  arrayView1d< localIndex const > list1 =  m_triCells.at(MortarSide::Slave).at(std::make_pair(shpS, shpM)).toView();
  arrayView1d< localIndex const > list2 =  m_triCells.at(MortarSide::Master).at(std::make_pair(shpS, shpM)).toView();
  arrayView3d< real64 const > gpLocalCoordsSlave = m_gpLocalCoords.at(MortarSide::Slave).at(std::make_pair(shpS, shpM)).toView();
  arrayView3d< real64 const > gpLocalCoordsMaster = m_gpLocalCoords.at(MortarSide::Master).at(std::make_pair(shpS, shpM)).toView();

  solidMechanicsMortarContactKernels::MortarFactory
  kernelFactorySlave( dispDofNumberSlave,
                      dofManager.rankOffset(),
                      localMatrix,
                      localRhs,
                      dt,
                      subRegionSlave,
                      list1,
                      list1,
                      detPlus,
                      gpLocalCoordsSlave,
                      tractionDofKey );

  std::cout << "Slave Kernel launching" << std::endl;
  real64 maxTraction1 = solidMechanicsMortarContactKernels::
                        interfaceBasedMortarKernelApplication
                        < parallelDevicePolicy< >,
                        FrictionBase >( meshSlave,
                                        subRegionSlave,
                                        nTri,
                                        feSlave,
                                        viewKeyStruct::frictionLawNameString(),
                                        kernelFactorySlave );


  solidMechanicsMortarContactKernels::MortarFactory
  kernelFactoryMaster( dispDofNumberMaster,
                       dofManager.rankOffset(),
                       localMatrix,
                       localRhs,
                       dt,
                       subRegionMaster,
                       list1,
                       list2,
                       detMinus.toView(),
                       gpLocalCoordsMaster,
                       tractionDofKey );

  std::cout << "Master Kernel launching" << std::endl;
  real64 maxTraction2 = solidMechanicsMortarContactKernels::
                        interfaceBasedMortarKernelApplication
                        < parallelDevicePolicy< >,
                        FrictionBase >( meshMaster,
                                        subRegionSlave,
                                        nTri,
                                        feMaster,
                                        viewKeyStruct::frictionLawNameString(),
                                        kernelFactoryMaster );

  GEOS_UNUSED_VAR(maxTraction1);
  GEOS_UNUSED_VAR(maxTraction2);

}


void SolidMechanicsMortarContact::implicitStepComplete( real64 const & time_n,
                                                        real64 const & dt,
                                                        DomainPartition & domain )
{

  SolidMechanicsLagrangianFEM::implicitStepComplete( time_n, dt, domain );

  FaceElementSubRegion & subRegion = 
    m_mortarSide.at( MortarSide::Slave ).surface->getUniqueSubRegion< FaceElementSubRegion >();

  arrayView2d< real64 > const deltaTraction  = subRegion.getField< contact::deltaTraction >();
  // arrayView2d< real64 > const deltaDispJump    = subRegion.getField< contact::deltaDispJump >();
  // arrayView2d< real64 const > const dispJump = subRegion.getField< contact::dispJump >();
  // arrayView2d< real64 > const oldDispJump    = subRegion.getField< contact::oldDispJump >()

  forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
  {
    //LvArray::tensorOps::fill< 3 >( deltaDispJump[kfe], 0.0 );
    LvArray::tensorOps::fill< 3 >( deltaTraction[kfe], 0.0 );
    //LvArray::tensorOps::copy< 3 >( oldDispJump[kfe], dispJump[kfe] );
  } );

}

real64 SolidMechanicsMortarContact::calculateResidualNorm( real64 const & time,
                                                           real64 const & dt,
                                                           DomainPartition const & domain,
                                                           DofManager const & dofManager,
                                                           arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  real64 const solidResidual = SolidMechanicsLagrangianFEM::calculateResidualNorm( time, dt, domain, dofManager, localRhs );

  real64 const contactResidual = calculateContactResidualNorm( domain, dofManager, localRhs );

  return sqrt( solidResidual * solidResidual + contactResidual * contactResidual );
}

real64 SolidMechanicsMortarContact::calculateContactResidualNorm( DomainPartition const & domain,
                                                                  DofManager const & dofManager,
                                                                  arrayView1d< real64 const > const & localRhs )
{

  GEOS_UNUSED_VAR( domain );

  string const & dofKey = dofManager.getKey( contact::traction::key() );
  globalIndex const rankOffset = dofManager.rankOffset();

  real64 stickResidual = 0.0;

   FaceElementSubRegion const & subRegion = 
    m_mortarSide.at( MortarSide::Slave ).surface->getUniqueSubRegion< FaceElementSubRegion >();

  arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
  arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
  arrayView1d< real64 const > const & area = subRegion.getElementArea();
  RAJA::ReduceSum< parallelHostReduce, real64 > stickSum( 0.0 );
  forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const k )
  {
    if( ghostRank[k] < 0 )
    {
      localIndex const localRow = LvArray::integerConversion< localIndex >( dofNumber[k] - rankOffset );
      for( localIndex dim = 0; dim < 3; ++dim )
      {
        real64 const norm = localRhs[localRow + dim] / area[k];
        stickSum += norm * norm;
      }
    }
  } );
  stickResidual += stickSum.get();

  stickResidual = MpiWrapper::sum( stickResidual );
  stickResidual = sqrt( stickResidual );

  if( getLogLevel() >= 1 && logger::internal::rank==0 )
  {
    std::cout << GEOS_FMT( "        ( Rt  ) = ( {:15.6e}  )", stickResidual );
  }

  return sqrt( stickResidual * stickResidual );
}

void SolidMechanicsMortarContact::applySystemSolution( DofManager const & dofManager,
                                                       arrayView1d< real64 const > const & localSolution,
                                                       real64 const scalingFactor,
                                                       real64 const dt,
                                                       DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  SolidMechanicsLagrangianFEM::applySystemSolution( dofManager, localSolution, scalingFactor, dt, domain );

  dofManager.addVectorToField( localSolution,
                               contact::traction::key(),
                               contact::deltaTraction::key(),
                               scalingFactor );

  dofManager.addVectorToField( localSolution,
                               contact::traction::key(),
                               contact::traction::key(),
                               scalingFactor );

  dofManager.addVectorToField( localSolution,
                               contact::totalBubbleDisplacement::key(),
                               contact::totalBubbleDisplacement::key(),
                               scalingFactor );

  dofManager.addVectorToField( localSolution,
                               contact::totalBubbleDisplacement::key(),
                               contact::incrementalBubbleDisplacement::key(),
                               scalingFactor );

  MeshLevel & mesh = *m_mortarSide.at( MortarSide::Slave ).mesh;

  FieldIdentifiers fieldsToBeSync;

  fieldsToBeSync.addFields( FieldLocation::Face,
                            { contact::incrementalBubbleDisplacement::key(),
                              contact::totalBubbleDisplacement::key() } );
  fieldsToBeSync.addElementFields( { contact::traction::key(),
                                     contact::deltaTraction::key(),
                                     contact::dispJump::key() },
                                   { getUniqueFractureRegionName() } );

  CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                       mesh,
                                                       domain.getNeighbors(),
                                                       true );

}



void SolidMechanicsMortarContact::addBubbleCouplingNumNonzeros( DofManager & dofManager,
                                                                arrayView1d< localIndex > const & rowLengths ) const
{

  MeshLevel & mesh = *m_mortarSide.at( MortarSide::Slave ).mesh;

  ElementRegionManager const & elemManagerSlave = mesh.getElemManager();
  NodeManager const & nodeManagerSlave = mesh.getNodeManager();
  FaceManager const & faceManagerSlave = mesh.getFaceManager();

  globalIndex const rankOffset = dofManager.rankOffset();

  string const bubbleDofKey = dofManager.getKey( contact::totalBubbleDisplacement::key() );
  string const dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );

  arrayView1d< globalIndex const > const bubbleDofNumber = faceManagerSlave.getReference< globalIndex_array >( bubbleDofKey );
  arrayView1d< globalIndex const > const dispSlaveDofNumber =  nodeManagerSlave.getReference< globalIndex_array >( dispDofKey );

  // add coupling between bubble and displacement in the slave side 
  elemManagerSlave.forElementSubRegions< CellElementSubRegion >([&]( const CellElementSubRegion & cellElementSubRegion )
  {

    arrayView1d< localIndex const > const bubbleElemsList = cellElementSubRegion.bubbleElementsList();
    arrayView2d< localIndex const > const faceElemsList = cellElementSubRegion.faceElementsList();

    localIndex const numDispDof = 3*cellElementSubRegion.numNodesPerElement();

    for( localIndex bi=0; bi<bubbleElemsList.size(); ++bi )
    {
      localIndex const cellIndex = bubbleElemsList[bi];
      localIndex const k = faceElemsList[bi][0];

      localIndex const localRow = LvArray::integerConversion< localIndex >( bubbleDofNumber[k] - rankOffset );

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
        localIndex const localDispRow = LvArray::integerConversion< localIndex >( dispSlaveDofNumber[node] - rankOffset );

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

}


void SolidMechanicsMortarContact::addMortarCouplingNumNonzeros( DofManager & dofManager,
                                                                ElementShape const & slaveShape,
                                                                ElementShape const & masterShape,
                                                                connectivityMapType const & connectivityMap,
                                                                arrayView1d< localIndex > const & rowLengths ) const
{

  FaceElementSubRegion const & surfRegionSlave = 
    m_mortarSide.at( MortarSide::Slave ).surface->getUniqueSubRegion< FaceElementSubRegion >();

  MeshLevel & meshMaster = *m_mortarSide.at( MortarSide::Master ).mesh;
  FaceElementSubRegion const & surfRegionMaster = 
    m_mortarSide.at( MortarSide::Master ).surface->getUniqueSubRegion< FaceElementSubRegion >();
  NodeManager const & nodeManagerMaster = meshMaster.getNodeManager();
  FaceManager const & faceManagerMaster = meshMaster.getFaceManager();
  ArrayOfArraysView< localIndex const > const faceToNodeMapMaster = faceManagerMaster.nodeList().toViewConst();
  arrayView2d< localIndex const > const elemsToFacesMaster = surfRegionMaster.faceList().toViewConst();

  ArrayOfArraysView< localIndex const > const connections = connectivityMap.at({slaveShape, masterShape}).toViewConst();
  arrayView1d< localIndex const > const masterElementList = m_faceTypeToElementList.at(MortarSide::Master).at(masterShape).toViewConst();
  arrayView1d< localIndex const > const slaveElementList = m_faceTypeToElementList.at(MortarSide::Slave).at(slaveShape).toViewConst();

  if (slaveElementList.size() == 0 || masterElementList.size() == 0)
  {
    return;
  }

  globalIndex const rankOffset = dofManager.rankOffset();

  string const dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
  string const tractionDofKey = dofManager.getKey( contact::traction::key() );
  arrayView1d< globalIndex const > const dispMasterDofNumber =  nodeManagerMaster.getReference< globalIndex_array >( dispDofKey );
  arrayView1d< globalIndex const > const tractionDofNumber = surfRegionSlave.getReference< globalIndex_array >( tractionDofKey );


  // add coupling between tractions and master displacements
  for( localIndex im=0; im<masterElementList.size(); ++im )
  {
    localIndex numDispDof;

    localIndex const kMaster = masterElementList[im];
    localIndex const kfaceMaster = elemsToFacesMaster[kMaster][0];
    localIndex const numNodesPerFace = faceToNodeMapMaster.sizeOfArray( kfaceMaster );

    numDispDof = 3*numNodesPerFace;

    localIndex const nSlave = connections.sizeOfArray(im);

  
    for( localIndex a=0; a<numNodesPerFace; ++a )
    {
      const localIndex & node = faceToNodeMapMaster( kfaceMaster, a );
      localIndex const localDispRow = LvArray::integerConversion< localIndex >( dispMasterDofNumber[node] - rankOffset );

      if( localDispRow >= 0 && localDispRow < rowLengths.size() )
      {
        for( int d=0; d<3; ++d )
        {
          rowLengths[localDispRow + d] += 3 * nSlave;
        }
      }
    }

    // loop over connected slave elements
    for (localIndex is=0; is<nSlave; ++is)
    {

      localIndex kSlave = slaveElementList(connections(im,is));
      localIndex const localRow = LvArray::integerConversion< localIndex >( tractionDofNumber[kSlave] - rankOffset );
      if( localRow >= 0 && localRow < rowLengths.size() )
      {
        for( localIndex i=0; i<3; ++i )
        {
          rowLengths[localRow + i] += numDispDof;
        }
      }
    }

  }
}


void SolidMechanicsMortarContact::addBubbleCouplingSparsityPattern( DofManager const & dofManager,
                                                                    SparsityPatternView< globalIndex > const & pattern ) const
{

  MeshLevel & mesh = *m_mortarSide.at( MortarSide::Slave ).mesh;
  ElementRegionManager const & elemManagerSlave = mesh.getElemManager();
  NodeManager const & nodeManagerSlave = mesh.getNodeManager();
  FaceManager const & faceManagerSlave = mesh.getFaceManager();

  FaceElementSubRegion const & surfRegionSlave = m_mortarSide.at( MortarSide::Slave ).surface->getUniqueSubRegion< FaceElementSubRegion >();
  arrayView2d< localIndex const > const elemsToFace = surfRegionSlave.faceList().toViewConst();


  globalIndex const rankOffset = dofManager.rankOffset();

  string const bubbleDofKey = dofManager.getKey( contact::totalBubbleDisplacement::key() );
  string const dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
  arrayView1d< globalIndex const > const bubbleDofNumber = faceManagerSlave.getReference< globalIndex_array >( bubbleDofKey );
  arrayView1d< globalIndex const > const dispDofNumber =  nodeManagerSlave.getReference< globalIndex_array >( dispDofKey );

  // for (localIndex i=0; i<elemsToFace.size(0); ++i)
  // {
  //   localIndex f = elemsToFace[i][1];
  //   std::cout << "Bubble dof " << i << " : " << bubbleDofNumber[f] << std::endl;
  // }

  static constexpr int maxNumDispDof = 3 * 8;

  elemManagerSlave.forElementSubRegions< CellElementSubRegion >( [&]( const CellElementSubRegion & cellElementSubRegion )
  {
    arrayView1d< localIndex const > const bubbleElemsList = cellElementSubRegion.bubbleElementsList();
    arrayView2d< localIndex const > const faceElemsList = cellElementSubRegion.faceElementsList();
    localIndex const numDispDof = 3*cellElementSubRegion.numNodesPerElement();
    for( localIndex bi=0; bi<bubbleElemsList.size(); ++bi )
    {
      localIndex const cellIndex = bubbleElemsList[bi];
      localIndex const k = faceElemsList[bi][0];

      // working arrays
      stackArray1d< globalIndex, maxNumDispDof > eqnRowIndicesDisp ( numDispDof );
      stackArray1d< globalIndex, 3 > eqnRowIndicesBubble( 3 );
      stackArray1d< globalIndex, maxNumDispDof > dofColIndicesDisp ( numDispDof );
      stackArray1d< globalIndex, 3 > dofColIndicesBubble( 3 );
      for( localIndex idof = 0; idof < 3; ++idof )
      {
        eqnRowIndicesBubble[idof] = bubbleDofNumber[k] + idof - rankOffset;
        dofColIndicesBubble[idof] = bubbleDofNumber[k] + idof;
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
          for( localIndex j = 0; j < dofColIndicesBubble.size(); ++j )
          {
            pattern.insertNonZero( eqnRowIndicesDisp[i], dofColIndicesBubble[j] );
            //std::cout << "Displacement: " << eqnRowIndicesDisp[i] << " Bubble: " << dofColIndicesBubble[j] << std::endl;
          }
        }
      }
      for( localIndex i = 0; i < eqnRowIndicesBubble.size(); ++i )
      {
        if( eqnRowIndicesBubble[i] >= 0 && eqnRowIndicesBubble[i] < pattern.numRows() )
        {
          for( localIndex j=0; j < dofColIndicesDisp.size(); ++j )
          {
            pattern.insertNonZero( eqnRowIndicesBubble[i], dofColIndicesDisp[j] );
            //std::cout << "Bubble: " << eqnRowIndicesBubble[i] << " Displacement: " << dofColIndicesDisp[j] << std::endl;
          }
        }
      }
    }
  } );

}


void SolidMechanicsMortarContact::addMortarCouplingSparsityPattern( DofManager & dofManager,
                                                                    ElementShape const & slaveShape,
                                                                    ElementShape const & masterShape,
                                                                    SparsityPatternView< globalIndex > const & pattern ) const
{

  FaceElementSubRegion const & surfRegionSlave = m_mortarSide.at( MortarSide::Slave ).surface->getUniqueSubRegion< FaceElementSubRegion >();

  MeshLevel & meshMaster = *m_mortarSide.at( MortarSide::Master ).mesh;
  FaceElementSubRegion const & surfRegionMaster = m_mortarSide.at( MortarSide::Master ).surface->getUniqueSubRegion< FaceElementSubRegion >();
  NodeManager const & nodeManagerMaster = meshMaster.getNodeManager();
  FaceManager const & faceManagerMaster = meshMaster.getFaceManager();
  ArrayOfArraysView< localIndex const > const faceToNodeMapMaster = faceManagerMaster.nodeList().toViewConst();
  arrayView2d< localIndex const > const elemsToFacesMaster = surfRegionMaster.faceList().toViewConst();

  //ArrayOfArraysView< localIndex const > const connections = connectivityMap.at({slaveShape, masterShape}).toViewConst();
  //arrayView1d< localIndex const > const masterElementList = m_faceTypeToElementList.at(MortarSide::Master).at(masterShape).toViewConst();
  //arrayView1d< localIndex const > const slaveElementList = m_faceTypeToElementList.at(MortarSide::Slave).at(slaveShape).toViewConst();

  if (m_triCellsDet.find(std::make_pair(slaveShape, masterShape)) == m_triCellsDet.end())
  {
    // no connections for the given element shapes
    return;
  }

  arrayView1d< localIndex const > const masterElementList = 
    m_triCells.at(MortarSide::Master).at(std::make_pair(slaveShape, masterShape)).toViewConst();

  arrayView1d< localIndex const > const slaveElementList = 
    m_triCells.at(MortarSide::Slave).at(std::make_pair(slaveShape, masterShape)).toViewConst();


  globalIndex const rankOffset = dofManager.rankOffset();

  string const dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
  string const tractionDofKey = dofManager.getKey( contact::traction::key() );
  arrayView1d< globalIndex const > const dispMasterDofNumber =  nodeManagerMaster.getReference< globalIndex_array >( dispDofKey );
  arrayView1d< globalIndex const > const tractionDofNumber = surfRegionSlave.getReference< globalIndex_array >( tractionDofKey );

  // populate sparsity pattern for coupling between tractions and master elements
  static constexpr int maxNumDispFaceDof = 3 * 4;

  for( localIndex k=0; k<masterElementList.size(); ++k )
  {
    localIndex kElMaster = masterElementList[k];
    localIndex kElSlave = slaveElementList[k];

    localIndex const kfaceMaster = elemsToFacesMaster[kElMaster][0];
    localIndex const numNodesPerFace = faceToNodeMapMaster.sizeOfArray( kfaceMaster );
    localIndex const numDispDof = 3*numNodesPerFace;
    // working arrays
    stackArray1d< globalIndex, maxNumDispFaceDof > eqnRowIndicesDisp ( numDispDof );
    stackArray1d< globalIndex, 3 > eqnRowIndicesTraction( 3 );
    stackArray1d< globalIndex, maxNumDispFaceDof > dofColIndicesDisp( numDispDof );
    stackArray1d< globalIndex, 3 > dofColIndicesTraction( 3 );
    for( localIndex idof = 0; idof < 3; ++idof )
    {
      eqnRowIndicesTraction[idof] = tractionDofNumber[kElSlave] + idof - rankOffset;
      dofColIndicesTraction[idof] = tractionDofNumber[kElSlave] + idof;
    }
    for( localIndex a=0; a<numNodesPerFace; ++a )
    {
      const localIndex & node = faceToNodeMapMaster( kfaceMaster, a );
      //std::cout << "Displacement dofs for face " << kfaceMaster << ", node " << node << std::endl;
      for( localIndex idof = 0; idof < 3; ++idof )
      {
        eqnRowIndicesDisp[3*a + idof] = dispMasterDofNumber[node] + idof - rankOffset;
        dofColIndicesDisp[3*a + idof] = dispMasterDofNumber[node] + idof;
      }
    }

    for( localIndex i = 0; i < eqnRowIndicesDisp.size(); ++i )
    {
      if( eqnRowIndicesDisp[i] >= 0 && eqnRowIndicesDisp[i] < pattern.numRows() )
      {
        for( localIndex j = 0; j < dofColIndicesTraction.size(); ++j )
        {
          pattern.insertNonZero( eqnRowIndicesDisp[i], dofColIndicesTraction[j] );
          //std::cout << "Displacement: " << eqnRowIndicesDisp[i] << " Traction: " << dofColIndicesTraction[j] << std::endl;
        }
      }
    }
    for( localIndex i = 0; i < eqnRowIndicesTraction.size(); ++i )
    {
      if( eqnRowIndicesTraction[i] >= 0 && eqnRowIndicesTraction[i] < pattern.numRows() )
      {
        for( localIndex j=0; j < dofColIndicesDisp.size(); ++j )
        {
          pattern.insertNonZero( eqnRowIndicesTraction[i], dofColIndicesDisp[j] );
          //std::cout << "Traction: " << eqnRowIndicesTraction[i] << " Displacement: " << dofColIndicesDisp[j] << std::endl;
        }
      }
    }
  }
  

}


void SolidMechanicsMortarContact::computeRotationMatrices( )
{

  MeshLevel & mesh = *m_mortarSide.at(MortarSide::Slave).mesh;
  FaceElementSubRegion & subRegion = 
    m_mortarSide.at(MortarSide::Slave).surface->getUniqueSubRegion< FaceElementSubRegion >();

  FaceManager & faceManager = mesh.getFaceManager();
  //ElementRegionManager & elemManager = mesh.getElemManager();

  arrayView2d< real64 const > const faceNormal = faceManager.faceNormal();
  arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList().toViewConst();

  // remember, view modify the original object
  arrayView3d< real64 > const rotationMatrix =
    subRegion.getField< contact::rotationMatrix >().toView();

  arrayView2d< real64 > const unitNormal   = subRegion.getNormalVector();
  arrayView2d< real64 > const unitTangent1 = subRegion.getTangentVector1();
  arrayView2d< real64 > const unitTangent2 = subRegion.getTangentVector2();

   // Compute rotation matrices
  solidMechanicsConformingContactKernels::ComputeRotationMatricesKernel::launch< parallelDevicePolicy<> >( subRegion.size(),
                                                                                                           faceNormal,
                                                                                                           elemsToFaces,
                                                                                                           rotationMatrix,
                                                                                                           unitNormal,
                                                                                                           unitTangent1,
                                                                                                           unitTangent2 );
}



void SolidMechanicsMortarContact::updateState( DomainPartition & domain )
{
  GEOS_UNUSED_VAR( domain );
}



void SolidMechanicsMortarContact::createFaceTypeListMortar( MortarSide side)
{
    // Generate lists containing elements of various face types
    MeshLevel & mesh = *m_mortarSide.at(side).mesh;
    FaceElementSubRegion const & surfRegion = m_mortarSide.at(side).surface->getUniqueSubRegion< FaceElementSubRegion >();

    FaceManager const & faceManager = mesh.getFaceManager();
    //ElementRegionManager const & elemManager = mesh.getElemManager();
    ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

    array1d< localIndex > keys( surfRegion.size());
    array1d< localIndex > vals( surfRegion.size());
    array1d< localIndex > quadList;
    array1d< localIndex > triList;
    RAJA::ReduceSum< ReducePolicy< parallelDevicePolicy<> >, localIndex > nTri_r( 0 );
    RAJA::ReduceSum< ReducePolicy< parallelDevicePolicy<> >, localIndex > nQuad_r( 0 );

    arrayView1d< localIndex > const keys_v = keys.toView();
    arrayView1d< localIndex > const vals_v = vals.toView();
    // Determine the size of the lists and generate the vector keys and vals for parallel indexing into lists.
    // (With RAJA, parallelizing this operation seems the most viable approach.)
    forAll< parallelDevicePolicy<> >( surfRegion.size(),
                                      [ = ] GEOS_HOST_DEVICE ( localIndex const kfe )
    {

      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kfe );
      if( numNodesPerFace == 3 )
      {
        keys_v[kfe]=0;
        vals_v[kfe]=kfe;
        nTri_r += 1;
        GEOS_ERROR( "SolidMechanicsMortarContact:: triangular face type not yet available" );
      }
      else if( numNodesPerFace == 4 )
      {
        keys_v[kfe]=1;
        vals_v[kfe]=kfe;
        nQuad_r += 1;
      }
      else
      {
        GEOS_ERROR( "SolidMechanicsMortarContact:: invalid face type" );
      }
    } );

    localIndex nQuad = static_cast< localIndex >(nQuad_r.get());
    localIndex nTri = static_cast< localIndex >(nTri_r.get());

    // Sort vals according to keys to ensure that
    // elements of the same type are adjacent in the vals list.
    // This arrangement allows for efficient copying into the container
    // by leveraging parallelism.
    RAJA::sort_pairs< parallelDevicePolicy<> >( keys_v, vals_v );

    quadList.resize( nQuad );
    triList.resize( nTri );
    arrayView1d< localIndex > const quadList_v = quadList.toView();
    arrayView1d< localIndex > const triList_v = triList.toView();

    forAll< parallelDevicePolicy<> >( nTri, [ = ] GEOS_HOST_DEVICE ( localIndex const kfe )
    {
      triList_v[kfe] = vals_v[kfe];
    } );

    forAll< parallelDevicePolicy<> >( nQuad, [ = ] GEOS_HOST_DEVICE ( localIndex const kfe )
    {
      quadList_v[kfe] = vals_v[nTri+kfe];
    } );

    std::map< string, array1d< localIndex > > faceTypeList;

    m_faceTypeToElementList[side][ElementShape::Quadrilateral] =  quadList;
    m_faceTypeToElementList[side][ElementShape::Triangle] =  triList;

}


void SolidMechanicsMortarContact::createBubbleCellList( ) const
{
    MeshLevel & mesh = *m_mortarSide.at(MortarSide::Slave).mesh;
    FaceElementSubRegion const & surfRegion = m_mortarSide.at(MortarSide::Slave).surface->getUniqueSubRegion< FaceElementSubRegion >();

    ElementRegionManager & elemManager = mesh.getElemManager();

    // Array to store face indexes
    array1d< localIndex > tmpSpace( 2*surfRegion.size());
    SortedArray< localIndex > faceIdList;

    arrayView1d< localIndex > const tmpSpace_v = tmpSpace.toView();
    // Store indexes of faces in the temporany array.
    {
      arrayView2d< localIndex const > const elemsToFaces = surfRegion.faceList().toViewConst();

      forAll< parallelDevicePolicy<> >( surfRegion.size(), [ = ] GEOS_HOST_DEVICE ( localIndex const kfe )
      {

        localIndex const kf0 = elemsToFaces[kfe][0], kf1 = elemsToFaces[kfe][1];
        tmpSpace_v[2*kfe] = kf0, tmpSpace_v[2*kfe+1] = kf1;

      } );
    }

    // Sort indexes to enable efficient searching using binary search.
    RAJA::stable_sort< parallelDevicePolicy<> >( tmpSpace_v );
    faceIdList.insert( tmpSpace_v.begin(), tmpSpace_v.end());

    // Search for bubble element on each CellElementSubRegion and
    // store element indexes, global and local face indexes.
    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & cellElementSubRegion )
    {

      arrayView2d< localIndex const > const elemsToFaces = cellElementSubRegion.faceList().toViewConst();

      RAJA::ReduceSum< ReducePolicy< parallelDevicePolicy<> >, localIndex > nBubElems_r( 0 );

      localIndex const n_max = cellElementSubRegion.size() * elemsToFaces.size( 1 );
      array1d< localIndex > keys( n_max );
      array1d< localIndex > perms( n_max );
      array1d< localIndex > vals( n_max );
      array1d< localIndex > localFaceIds( n_max );

      arrayView1d< localIndex > const keys_v = keys.toView();
      arrayView1d< localIndex > const perms_v = perms.toView();
      arrayView1d< localIndex > const vals_v = vals.toView();
      arrayView1d< localIndex > const localFaceIds_v = localFaceIds.toView();
      SortedArrayView< localIndex const > const faceIdList_v = faceIdList.toViewConst();

      forAll< parallelDevicePolicy<> >( cellElementSubRegion.size(),
                                        [ = ]
                                        GEOS_HOST_DEVICE ( localIndex const kfe )
      {
        for( int i=0; i < elemsToFaces.size( 1 ); ++i )
        {
          perms_v[kfe*elemsToFaces.size( 1 )+i] = kfe*elemsToFaces.size( 1 )+i;
          if( faceIdList_v.contains( elemsToFaces[kfe][i] ))
          {
            keys_v[kfe*elemsToFaces.size( 1 )+i] = 0;
            vals_v[kfe*elemsToFaces.size( 1 )+i] = kfe;
            localFaceIds_v[kfe*elemsToFaces.size( 1 )+i] = i;
            nBubElems_r += 1;
          }
          else
          {
            keys_v[kfe*elemsToFaces.size( 1 )+i] = 1;
            vals_v[kfe*elemsToFaces.size( 1 )+i] = -1;
            localFaceIds_v[kfe*elemsToFaces.size( 1 )+i] = -1;
          }
        }
      } );

      // Sort perms according to keys to ensure that bubble elements are adjacent
      // and occupy the first positions of the list.
      // This arrangement allows for efficient copying into the container
      // by leveraging parallelism.
      localIndex nBubElems = static_cast< localIndex >(nBubElems_r.get());
      RAJA::sort_pairs< parallelDevicePolicy<> >( keys_v, perms_v );

      array1d< localIndex > bubbleElemsList;
      bubbleElemsList.resize( nBubElems );

      arrayView1d< localIndex > const bubbleElemsList_v = bubbleElemsList.toView();

      forAll< parallelDevicePolicy<> >( n_max, [ = ] GEOS_HOST_DEVICE ( localIndex const k )
      {
        keys_v[k] = vals_v[perms_v[k]];
      } );

      forAll< parallelDevicePolicy<> >( nBubElems, [ = ] GEOS_HOST_DEVICE ( localIndex const k )
      {
        bubbleElemsList_v[k] = keys_v[k];
      } );
      cellElementSubRegion.setBubbleElementsList( bubbleElemsList.toViewConst());

      forAll< parallelDevicePolicy<> >( n_max, [ = ] GEOS_HOST_DEVICE ( localIndex const k )
      {
        keys_v[k] = localFaceIds_v[perms_v[k]];
      } );

      array2d< localIndex > faceElemsList;
      faceElemsList.resize( nBubElems, 2 );

      arrayView2d< localIndex > const faceElemsList_v = faceElemsList.toView();

      forAll< parallelDevicePolicy<> >( nBubElems,
                                        [ = ]
                                        GEOS_HOST_DEVICE ( localIndex const k )
      {
        localIndex const kfe =  bubbleElemsList_v[k];
        faceElemsList_v[k][0] = elemsToFaces[kfe][keys_v[k]];
        faceElemsList_v[k][1] = keys_v[k];
      } );
      cellElementSubRegion.setFaceElementsList( faceElemsList.toViewConst());

    } );

}

void SolidMechanicsMortarContact::setMortarSurfaces( DomainPartition & domain)
{

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {

    GEOS_UNUSED_VAR(meshName);
    ElementRegionManager  & elemManager = mesh.getElemManager();
    elemManager.forElementRegions< SurfaceElementRegion >([&,this]( SurfaceElementRegion & region )
    {
      // TO DO: assign different subregions for triangles and quadrilaterals?
      // assign surface to master or slave depending on object path (again, naive)
      string surfacePath = region.getPath();
      std::cout << surfacePath << std::endl;
      std::cout << mesh.getName() << std::endl;
      // check if region is master or slave and populate member maps
      if (surfacePath.find(m_slaveName) != std::string::npos)
      {
        MortarSurface surfaceSlave;
        surfaceSlave.mesh = &mesh;
        surfaceSlave.surface = &region;
        m_mortarSide[MortarSide::Slave] = surfaceSlave;
        m_meshSlaveName = meshName;
      }
      else if (surfacePath.find(m_masterName) != std::string::npos)
      {
        MortarSurface surfaceMaster;
        surfaceMaster.mesh = &mesh;
        surfaceMaster.surface = &region;
        m_mortarSide[MortarSide::Master] = surfaceMaster;
      }
    });
  } );

  // debug log surfaces path
  string pathMaster;
  string pathMeshMaster;
  pathMaster = m_mortarSide[MortarSide::Master].surface->getPath();
  pathMeshMaster = m_mortarSide[MortarSide::Master].mesh->getPath();
  std::cout << "Path of master surface: " << pathMaster << std::endl;
  std::cout << "Path of master mesh level: " << pathMeshMaster << std::endl; 
  string pathSlave;
  string pathMeshSlave;
  pathSlave = m_mortarSide[MortarSide::Slave].surface->getPath();
  pathMeshSlave = m_mortarSide[MortarSide::Slave].mesh->getPath();
  std::cout << "Path of slave surface: " << pathSlave << std::endl;
  std::cout << "Path of slave mesh level: " << pathMeshSlave << std::endl;

}

void SolidMechanicsMortarContact::computeMortarInterpolation ( 
  connectivityMapType & connectivityMap)
{
  computeMortarInterpolation< ElementShape::Triangle, ElementShape::Triangle>( 
    connectivityMap[{ElementShape::Triangle, ElementShape::Triangle}]);

  computeMortarInterpolation< ElementShape::Triangle, ElementShape::Quadrilateral>( 
    connectivityMap[{ElementShape::Triangle, ElementShape::Quadrilateral}]);

  computeMortarInterpolation< ElementShape::Quadrilateral, ElementShape::Triangle>( 
    connectivityMap[{ElementShape::Quadrilateral, ElementShape::Triangle}]);

  computeMortarInterpolation< ElementShape::Quadrilateral, ElementShape::Quadrilateral>( 
    connectivityMap[{ElementShape::Quadrilateral, ElementShape::Quadrilateral}]);
}

template< ElementShape slaveShape, ElementShape masterShape >
void SolidMechanicsMortarContact::computeMortarInterpolation( ArrayOfArrays<localIndex> const & connections )
{

  arrayView1d< localIndex const> faceListMaster = m_faceTypeToElementList.at(MortarSide::Master).at(masterShape).toViewConst();
  arrayView1d< localIndex const> faceListSlave = m_faceTypeToElementList.at(MortarSide::Slave).at(slaveShape).toViewConst();

  localIndex numbConnections;
  numbConnections = 0;
  for (localIndex i = 0; i < connections.size(); ++i)
  {
    numbConnections += connections.sizeOfArray(i);
  }
  std::cout << "Found " << numbConnections << " connections." << std::endl;

  // allocate mortar quantities
  auto const & slaveFE = getFE< slaveShape >();
  auto const & masterFE = getFE< masterShape >();
  localIndex constexpr nNodeMaster = masterFE.numNodes;
  localIndex constexpr nNodeSlave = slaveFE.numNodes;
  // maximum number of expected subtriangles for each mortar pair
  localIndex constexpr maxSubTri = nNodeMaster + nNodeSlave - 2;
  localIndex  nSubTriangles = maxSubTri * numbConnections;

  array1d<real64> numbTrianglesInPairs(numbConnections);
  array2d<localIndex> cellPairs(nSubTriangles,2);
  array1d<real64> subTriDeterminants(nSubTriangles);
  array3d<real64> localCoordsSlave(nSubTriangles, nGPtri, 2);
  array3d<real64> localCoordsMaster(nSubTriangles, nGPtri, 2);
  // get views
  arrayView1d<real64> numbTrianglesInPairs_v = numbTrianglesInPairs.toView();
  arrayView2d<localIndex> cellPairs_v = cellPairs.toView();
  arrayView1d<real64> subTriDeterminants_v = subTriDeterminants.toView();
  arrayView3d<real64> localCoordsSlave_v = localCoordsSlave.toView();
  arrayView3d<real64> localCoordsMaster_v = localCoordsMaster.toView();

  if ( numbConnections==0 )
  {
    return; 
  } 

  // get mesh manager objects for both sides
  MeshLevel & meshMaster = *m_mortarSide.at(MortarSide::Master).mesh;
  FaceElementSubRegion const & surfMaster = m_mortarSide.at(MortarSide::Master).surface->getUniqueSubRegion< FaceElementSubRegion >();
  FaceManager const & faceManagerMaster = meshMaster.getFaceManager();
  NodeManager const & nodeManagerMaster = meshMaster.getNodeManager();
  arrayView2d< localIndex const > const elemsToFacesMaster = surfMaster.faceList().toViewConst();
  ArrayOfArraysView< localIndex const > const & faceToNodeMapMaster = faceManagerMaster.nodeList().toViewConst();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const coordsMaster =  nodeManagerMaster.referencePosition();

  MeshLevel & meshSlave = *m_mortarSide.at(MortarSide::Slave).mesh;
  FaceElementSubRegion const & surfSlave = m_mortarSide.at(MortarSide::Slave).surface->getUniqueSubRegion< FaceElementSubRegion >();
  FaceManager const & faceManagerSlave = meshSlave.getFaceManager();
  NodeManager const & nodeManagerSlave = meshSlave.getNodeManager();
  arrayView2d< localIndex const > const elemsToFacesSlave = surfSlave.faceList().toViewConst();
  ArrayOfArraysView< localIndex const > const & faceToNodeMapSlave = faceManagerSlave.nodeList().toViewConst();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const coordsSlave =  nodeManagerSlave.referencePosition();

  localIndex k = 0;
  localIndex totTri = 0;

  for (localIndex im = 0; im < faceListMaster.size(); ++im)
  {
    localIndex masterFaceId = faceListMaster[im];

    for (localIndex is = 0; is < connections.sizeOfArray(im); ++is)
    {
      localIndex slaveFaceId = faceListSlave[connections(im, is)];

      localIndex kfeM = elemsToFacesMaster[masterFaceId][0];
      localIndex kfeS = elemsToFacesSlave[slaveFaceId][0];

      arraySlice1d< localIndex const > nodesMaster = faceToNodeMapMaster[kfeM];
      arraySlice1d< localIndex const > nodesSlave = faceToNodeMapSlave[kfeS];

      localIndex nTriPair = processMortarPair < slaveShape, masterShape > ( slaveFaceId, 
                                                                            masterFaceId,
                                                                            nodesSlave,
                                                                            nodesMaster,
                                                                            coordsSlave,
                                                                            coordsMaster,
                                                                            cellPairs_v,
                                                                            subTriDeterminants_v,
                                                                            localCoordsSlave_v,
                                                                            localCoordsMaster_v,
                                                                            k
                                                                          );

      numbTrianglesInPairs[k] = nTriPair;
      totTri += nTriPair;
      ++k;
    }
  }

  // allocate and populate member maps with segment based preprocessed info
  array1d<localIndex> triCellsListSlave(totTri);
  array1d<localIndex> triCellsListMaster(totTri);
  array1d<real64> triCellsDetList(totTri);
  array3d<real64> gpLocalCoordsSlave(totTri, nGPtri, 2);
  array3d<real64> gpLocalCoordsMaster(totTri, nGPtri, 2);

  localIndex kp = 0;
  for (localIndex i=0; i<numbConnections; ++i)
  {
    for (localIndex j=0; j<numbTrianglesInPairs[i]; ++j)
    {
      triCellsListSlave(kp) = cellPairs_v(i*maxSubTri+j,0);
      triCellsListMaster(kp) = cellPairs_v(i*maxSubTri+j,1);
      triCellsDetList(kp) = subTriDeterminants_v(i*maxSubTri+j);
      for (localIndex q=0; q < nGPtri; ++q)
      {
        gpLocalCoordsSlave(kp,q,0) = localCoordsSlave_v(i*maxSubTri+j,q,0);
        gpLocalCoordsSlave(kp,q,1) = localCoordsSlave_v(i*maxSubTri+j,q,1);
        gpLocalCoordsMaster(kp,q,0) = localCoordsMaster_v(i*maxSubTri+j,q,0);
        gpLocalCoordsMaster(kp,q,1) = localCoordsMaster_v(i*maxSubTri+j,q,1);
      }
      ++kp;
    }
  }

  // Log the segment based information
  std::cout << "======================================================================================================================" << std::endl;
  std::cout << "                                     SEGMENT BASED MORTAR PREPROCESSING                                               " << std::endl;
  std::cout << "======================================================================================================================" << std::endl;
  std::cout << std::setw(10) << "Tri ID"
            << std::setw(15) << "Slave Face"
            << std::setw(15) << "Master Face"
            << std::setw(20) << "Determinant" << std::endl;
  std::cout << "----------------------------------------------------------------------------------------------------------------------" << std::endl;

  for (localIndex i = 0; i < totTri; ++i)
  {
    std::cout << std::setw(10) << i
              << std::setw(15) << triCellsListSlave(i)
              << std::setw(15) << triCellsListMaster(i)
              << std::setw(20) << std::scientific << triCellsDetList(i) << std::endl;

    std::cout << "  GP Coords Slave: ";
    for (localIndex q = 0; q < nGPtri; ++q)
    {
      std::cout << "(" << std::fixed << std::setprecision(4) << gpLocalCoordsSlave(i, q, 0) << ", "
                << std::fixed << std::setprecision(4) << gpLocalCoordsSlave(i, q, 1) << ") ";
    }
    std::cout << std::endl;

    std::cout << "  GP Coords Master: ";
    for (localIndex q = 0; q < nGPtri; ++q)
    {
      std::cout << "(" << std::fixed << std::setprecision(4) << gpLocalCoordsMaster(i, q, 0) << ", "
                << std::fixed << std::setprecision(4) << gpLocalCoordsMaster(i, q, 1) << ") ";
    }
    std::cout << std::endl;
    std::cout << "----------------------------------------------------------------------------------------------------------------------" << std::endl;
  }

  m_triCells[MortarSide::Slave][{slaveShape, masterShape}] = triCellsListSlave;
  m_triCells[MortarSide::Master][{slaveShape, masterShape}] = triCellsListMaster;
  m_triCellsDet[{slaveShape, masterShape}] = triCellsDetList;
  m_gpLocalCoords[MortarSide::Slave][{slaveShape, masterShape}] = gpLocalCoordsSlave;
  m_gpLocalCoords[MortarSide::Master][{slaveShape, masterShape}] = gpLocalCoordsMaster;
}


template<ElementShape slaveShape, ElementShape masterShape>
localIndex SolidMechanicsMortarContact::processMortarPair( localIndex const slaveFaceId, 
                                                           localIndex const masterFaceId,
                                                           arraySlice1d<localIndex const> const & nodesSlave,
                                                           arraySlice1d<localIndex const> const & nodesMaster,
                                                           arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & coordsSlave,
                                                           arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & coordsMaster,
                                                           arrayView2d<localIndex> & cellPairs,
                                                           arrayView1d<real64> & subTriDeterminants,
                                                           arrayView3d<real64> & localCoordsSlave,
                                                           arrayView3d<real64> & localCoordsMaster,
                                                           localIndex const & kPair)
{
  
  std::cout << "Processing mortar pair slave/master: " << slaveFaceId << " - " << masterFaceId << std::endl;

  auto const & slaveFE = getFE< slaveShape >();
  auto const & masterFE = getFE< masterShape >();

  // compile time knowledge ?
  constexpr static localIndex numNodeMaster = masterFE.numNodes;
  constexpr static localIndex numNodeSlave = slaveFE.numNodes;

  // STEPS FOR THE SEGMENT BASED SCHEME
  // 1. find auxiliary plane on the slave element (use computational geometry)
  real64 planeCentroid[3];
  real64 planeNormal[3];
  real64 area = computationalGeometry::centroid_3DPolygon( nodesSlave,
                                                           coordsSlave,
                                                           planeCentroid,
                                                           planeNormal);

  GEOS_UNUSED_VAR(area);

  // 2. project slave and master nodes on the auxiliary plane (use computational geometry)
  array2d<real64> projSlave(numNodeSlave, 2);
  array2d<real64> projMaster(numNodeMaster, 2);

  for (localIndex i=0; i<numNodeSlave; ++i)
  {
    // find the 3D intersection point
    real64 coord3d[3] = {coordsSlave(nodesSlave[i],0), coordsSlave(nodesSlave[i],1), coordsSlave(nodesSlave[i],2)};
    real64 projected2d[2];
    projectPointInPlane< MortarSide::Slave >(coord3d, planeNormal, planeCentroid, projected2d);
    projSlave[i][0] = projected2d[0];
    projSlave[i][1] = projected2d[1];
  }

  for (localIndex i=0; i<numNodeMaster; ++i)
  {
    // find the 3D intersection point
    real64 coord3d[3] = {coordsMaster(nodesMaster[i],0), coordsMaster(nodesMaster[i],1), coordsMaster(nodesMaster[i],2)};
    real64 projected2d[2];
    projectPointInPlane< MortarSide::Master >(coord3d, planeNormal, planeCentroid, projected2d);
    projMaster[i][0] = projected2d[0];
    projMaster[i][1] = projected2d[1];
  }

  // std::cout << "Coordinates of projected slave element:" << std::endl;
  // for ( localIndex i=0; i<projSlave.size(0); ++i)
  // {
  //  std::cout << "(" << projSlave(i,0) << ", " << projSlave(i,1) << ")" << std::endl;
  // }
  // std::cout << "Coordinates of projected master element:" << std::endl;
  // for ( localIndex i=0; i<projMaster.size(0); ++i)
  // {
  //  std::cout << "(" << projMaster(i,0) << ", " << projMaster(i,1) << ")" << std::endl;
  // }

  // 3. intersect slave and master element and return the resulting intersecting polygon
  array2d<real64> clippedPoly(numNodeSlave,2);
  LvArray::tensorOps::copy<numNodeSlave,2>(clippedPoly, projSlave);
  polygonClipping<numNodeSlave,numNodeMaster>(clippedPoly, projMaster);

  localIndex clipSize = clippedPoly.size(0);
  //std::cout << "Found intersection with " << clipSize << " edges." << std::endl;
  arrayView2d<real64 const> const clipPoly = clippedPoly.toViewConst();


  // 4. project gauss point for each local triangle into master and slave side
  for (localIndex i = 1; i < clipSize - 1; ++i)
  {
    real64 coordsTri[3][2];
    // get coordinates of sub-triangle
    coordsTri[0][0] = clipPoly(0, 0);
    coordsTri[0][1] = clipPoly(0, 1);
    coordsTri[1][0] = clipPoly(i, 0);
    coordsTri[1][1] = clipPoly(i, 1);
    coordsTri[2][0] = clipPoly(i+1, 0);
    coordsTri[2][1] = clipPoly(i+1, 1);

    real64 const subTriDeterminant = (coordsTri[1][0] - coordsTri[0][0]) * (coordsTri[2][1] - coordsTri[0][1]) -
                                     (coordsTri[2][0] - coordsTri[0][0]) * (coordsTri[1][1] - coordsTri[0][1]);

    if ( subTriDeterminant < 1e-10)
    {
      GEOS_ERROR("Found degenerate triangle");
    }

    real64 xiSlave[nGPtri][2] = {{}};  
    real64 xiMaster[nGPtri][2] = {{}};  

    projectGP< slaveShape >(coordsTri, projSlave.toViewConst(), xiSlave);

    projectGP< masterShape >(coordsTri, projMaster.toViewConst(), xiMaster);

    // 5. populate output lists
    localIndex maxTri = numNodeSlave+numNodeMaster-2;
    cellPairs(kPair*maxTri+i-1, 0) = slaveFaceId;
    cellPairs(kPair*maxTri+i-1, 1) = masterFaceId;

    subTriDeterminants[kPair*maxTri+i-1] = subTriDeterminant;

    for (localIndex j = 0; j < nGPtri; ++j)
    {
      localCoordsSlave(kPair*maxTri+i-1, j, 0) = xiSlave[j][0];
      localCoordsSlave(kPair*maxTri+i-1, j, 1) = xiSlave[j][1];

      localCoordsMaster(kPair*maxTri+i-1, j, 0) = xiMaster[j][0];
      localCoordsMaster(kPair*maxTri+i-1, j, 1) = xiMaster[j][1];
    }

    std::cout << std::endl;

  }

  return clipSize - 2; // return number of triangles for the mortar pair

}


template< MortarSide side > 
void SolidMechanicsMortarContact::projectPointInPlane( real64 const (& coord3d)[3],
                                                       real64 const (& normal)[3],
                                                       real64 const (& origin)[3],
                                                       real64 (& proj2d)[2])
{

  real64 proj3d[3];
  computationalGeometry::LinePlaneIntersection(normal,
                                               coord3d,
                                               normal,
                                               origin,
                                               proj3d);

  // find coordinates in the local plane
  real64 tmp[3];

  if (std::fabs(normal[0]) < 0.9)
  {
    tmp[0] = 1.0;
    tmp[1] = 0.0;
    tmp[2] = 0.0;
  }
  else
  {
    tmp[0] = 0.0;
    tmp[1] = 1.0;
    tmp[2] = 0.0;
  }
  real64 u[3];
  real64 v[3];

  LvArray::tensorOps::crossProduct( u, normal, tmp );
  LvArray::tensorOps::crossProduct( v, u, normal );
  LvArray::tensorOps::normalize<3>(u);
  LvArray::tensorOps::normalize<3>(v);

  if constexpr( side == MortarSide::Slave)
  {
    // flip v to ensure the slave polygon is in CCW order in 2D.
    LvArray::tensorOps::scale<3>(v, -1);
  }

  LvArray::tensorOps::subtract<3>(proj3d, origin);

  proj2d[0] = LvArray::tensorOps::AiBi<3>(u, proj3d);
  proj2d[1] = LvArray::tensorOps::AiBi<3>(v,proj3d);
}

template< localIndex sizePoly, localIndex sizeClipper>
void SolidMechanicsMortarContact::polygonClipping(array2d<real64>& poly,
                                                  array2d<real64> & clipPoly)
{
  // input polygons must be in CCW order!
  for (localIndex i = 0; i < sizeClipper; ++i)
  {
    localIndex k = (i + 1) % sizeClipper;
    clip( poly, clipPoly(i,0), clipPoly(i,1), clipPoly(k,0), clipPoly(k,1));
  }

  for (localIndex i = 0; i < poly.size(0); ++i)
  {
    std::cout << "Clipped polygon: (" << poly(i,0) << ", " << poly(i,1) << ")" << std::endl;
  }

}


// intersect lines [(x1,y1),(x2,y2)] with [(x3,y3),(x4,y4)]
void SolidMechanicsMortarContact::intersect(real64 x1, real64 y1,real64 x2, real64 y2,
                                       real64 x3, real64 y3,real64 x4, real64 y4,
                                       real64 & xInt, real64 & yInt)
{
  real64 numX = (x1*y2 - y1*x2) * (x3-x4) -
              (x1-x2) * (x3*y4 - y3*x4);
  real64 denX = (x1-x2) * (y3-y4) - (y1-y2) * (x3-x4);

  real64 numY = (x1*y2 - y1*x2) * (y3-y4) -
              (y1-y2) * (x3*y4 - y3*x4);
  real64 denY = (x1-x2) * (y3-y4) - (y1-y2) * (x3-x4);

  if (std::fabs(denX) < 1e-10 || std::fabs(denY) < 1e-10)
  {
    GEOS_ERROR("Lines are parallel");
  }
  else
  {
    xInt = numX / denX;
    yInt = numY / denY;
  }

}

template< ElementShape shape >
void SolidMechanicsMortarContact::projectGP( real64 const (& coordsTri)[3][2],
                                             arrayView2d<real64 const> const & coordsElem,
                                             real64 (& xi)[nGPtri][2])
{

  constexpr localIndex numNodes = (shape == ElementShape::Triangle) ? 3 : 4;

  auto const & FE = getFE< shape >();

  // newton params
  localIndex itMax = 3;
  real64 tol = 1e-9;

  for (localIndex i = 0; i < nGPtri; ++i)
  {
    real64 xiProj[2] = {0.0, 0.0};

    real64 Ntri[3];
    real64 xiq[2] = {qCoords[i][0], qCoords[i][1]};
 
    feTriangleCell::calcN( xiq, Ntri );
    real64 coordGP[2];
    LvArray::tensorOps::Ri_eq_AjiBj<2, 3>(coordGP, coordsTri, Ntri);  

    //std::cout << "Coordinate of gp: ( " << coordGP[0] << " , " << coordGP[1] << ")" << std::endl;
 
    real64 Nq[numNodes];
    FE.calcN(xiProj, Nq);
    permuteN<numNodes>(Nq);

    real64 rhs[2] = {0.0, 0.0};
    LvArray::tensorOps::Ri_eq_AjiBj<2, numNodes>(rhs, coordsElem, Nq);
    LvArray::tensorOps::subtract<2>(rhs, coordGP);
 
    localIndex iter = 0;
 
    while (LvArray::tensorOps::l2Norm<2>(rhs) > tol && iter < itMax)
    {
      iter = iter + 1;
      real64 dN[2][numNodes]= {{}};
      calcGradN< numNodes >(xiProj, dN);

      real64 Jt[2][2] = {{}};
      real64 J[2][2] = {{}};
      LvArray::tensorOps::Rij_eq_AikBkj<2, 2, numNodes>(Jt, dN, coordsElem);
      LvArray::tensorOps::transpose<2,2>(J,Jt);

      real64 dxi[2];

      bool success = denseLinearAlgebra::details::solveTwoByTwoSystem( J, rhs, dxi);

      if (success)
      {
        xiProj[0] -= dxi[0];
        xiProj[1] -= dxi[1];
      }
      else
      {
        GEOS_ERROR("Failed to solve linear system in GP projection algorithm");
      }

      FE.calcN(xiProj, Nq);
      permuteN<numNodes>(Nq);

      //std::cout << std::endl;
      LvArray::tensorOps::Ri_eq_AjiBj<2, numNodes>(rhs, coordsElem, Nq);
      LvArray::tensorOps::subtract<2>(rhs, coordGP);
      //std::cout << "rhs: " << rhs[0] << ", " << rhs[1] << std::endl;
    }

    if ( iter == itMax)
    {
      GEOS_ERROR("SolidMechanicsMortarContact::projectGP - Newton raphson not converged");
    }

    if ( !checkInFE< shape >(xiProj[0], xiProj[1]))
    {
      GEOS_ERROR("Projected GP is outside the reference element");
    }

    xi[i][0] = xiProj[0];
    xi[i][1] = xiProj[1];

  }
}         

template<ElementShape shape>
bool SolidMechanicsMortarContact::checkInFE(real64 xi0, real64 xi1)
{
    real64 tol = 1e-7;
    bool isInRange = false;

    if constexpr (shape == ElementShape::Triangle)
    {
        // barycentric coords must be >= 0 and xi1+xi2 <= 1
        isInRange = (xi0 >= -tol && xi1 >= -tol &&
                     xi0 <= 1.0 + tol && xi1 <= 1.0 + tol &&
                     xi0 + xi1 <= 1.0 + tol);
    }
    else if constexpr (shape == ElementShape::Quadrilateral)
    {
        // in square [-1,1] x [-1,1]
        isInRange = (xi0 >= -1.0 - tol && xi0 <= 1.0 + tol &&
                     xi1 >= -1.0 - tol && xi1 <= 1.0 + tol);
    }

    return isInRange;
}


  // clip polygon against clipping line [(xc1,yc1),(xc2,yc2)]
void SolidMechanicsMortarContact::clip( array2d<real64> & poly,
                                        real64 xc1, real64 yc1, real64 xc2, real64 yc2)
{
  constexpr localIndex maxPoints = 8;

  real64 clippedPoly[maxPoints][2];
  localIndex clippedSize = 0;

  localIndex polySize = poly.size(0);
  //std::cout << "Clipping polygon with " << polySize << " points." << std::endl;

  for (localIndex i = 0; i < polySize; ++i)
  {
    localIndex k = (i+1) % polySize;
    real64 x1 = poly(i,0);
    real64 y1 = poly(i,1);
    real64 x2 = poly(k,0);
    real64 y2 = poly(k,1);

    // positions of points w.r.t clipper line
    real64 i_pos = (xc2 - xc1) * (y1 - yc1) - (yc2 - yc1) * (x1 - xc1);
    real64 k_pos = (xc2 - xc1) * (y2 - yc1) - (yc2 - yc1) * (x2 - xc1);

    real64 xInt = 0.0;
    real64 yInt = 0.0;

    //std::cout << i_pos << ", " << k_pos << std::endl;

    // case 1: both points inside
    if (i_pos >= 0 &&  k_pos >= 0)
    {
      // add only second point
      clippedPoly[clippedSize][0] = x2;
      clippedPoly[clippedSize][1] = y2;
      ++clippedSize;
    }
    // case 2: only first point is outside
    else if (i_pos < 0 && k_pos >= 0)
    {
      // compute intersection point
      intersect(x1, y1, x2, y2, xc1, yc1, xc2, yc2, xInt, yInt);

      // add intersection point and second point
      clippedPoly[clippedSize][0] = xInt;
      clippedPoly[clippedSize][1] = yInt;
      ++clippedSize;
      clippedPoly[clippedSize][0] = x2;
      clippedPoly[clippedSize][1] = y2;
      ++clippedSize;
    }
    // case 3: only second point is outside
    else if (i_pos >= 0 && k_pos < 0)
    {
      // add intersection point
      intersect(x1, y1, x2, y2, xc1, yc1, xc2, yc2, xInt, yInt);
      clippedPoly[clippedSize][0] = xInt;
      clippedPoly[clippedSize][1] = yInt;
      ++clippedSize;      
    }
    // case 4: both points are outside, no points are added
  }

  // Update the original polygon with the clipped points
  poly.resize(clippedSize, 2);
  for (localIndex i = 0; i < clippedSize; ++i)
  {
    poly(i,0) = clippedPoly[i][0];
    poly(i,1) = clippedPoly[i][1];
  }
}

template<localIndex numNodeElement>
void SolidMechanicsMortarContact::calcGradN( real64 const (& xi)[2], real64 (& dN)[2][numNodeElement] )
{
  // provisional method here, to be moved in Element classes
  // for the gradients permutation has already been applied!
  if constexpr ( numNodeElement == 3)
  {
    dN[0][1] = -1.0;
    dN[0][1] = 1.0;
    dN[0][2] = 0.0;
    dN[1][0] = -1.0;
    dN[1][1] = 0.0;
    dN[1][2] = 1.0;
  }
  else if constexpr ( numNodeElement == 4)
  {
    dN[0][0] = -0.25 * (1.0 - xi[1]);
    dN[0][1] =  0.25 * (1.0 - xi[1]);
    dN[0][2] =  0.25 * (1.0 + xi[1]);
    dN[0][3] = -0.25 * (1.0 + xi[1]);
    dN[1][0] = -0.25 * (1.0 - xi[0]);
    dN[1][1] = -0.25 * (1.0 + xi[0]);
    dN[1][2] =  0.25 * (1.0 + xi[0]);
    dN[1][3] =  0.25 * (1.0 - xi[0]);
  }
  else
  {
    GEOS_ERROR("SolidMechanicsMortarContact:: invalid number of nodes for computing gradient of basis functions");
  }
}

template<localIndex numNodes>
void SolidMechanicsMortarContact::permuteN(real64 (& N)[numNodes])
{
  if constexpr (numNodes == 3)
  {
    // permutation not need for triangles
    return;
  }
  else if  constexpr (numNodes == 4)
  {
    real64 Ntmp[numNodes];
    LvArray::tensorOps::copy<numNodes>(Ntmp, N);
    localIndex permutation[4];
    permutation[0] = 0;
    permutation[1] = 1;
    permutation[2] = 3;
    permutation[3] = 2;

    for (localIndex i=0; i<numNodes; ++i)
    {
      N[permutation[i]] = Ntmp[i];
    }
  }
}


void SolidMechanicsMortarContact::getConnectivityMap( connectivityMapType & connectivityMap )
{
  ElementShape shapes[2] = { ElementShape::Triangle, ElementShape::Quadrilateral };
  for (auto slaveShape : shapes)
  {
    for (auto masterShape : shapes)
    {
      // loop over all master-slave element shape pairs
      ArrayOfArrays<localIndex> connections;
      getMortarConnections(slaveShape, masterShape, connections);
      connectivityMap[{slaveShape, masterShape}] = connections;
    }
  }  
}


void SolidMechanicsMortarContact::getMortarConnections( ElementShape slaveShape, 
                                                      ElementShape masterShape, 
                                                      ArrayOfArrays<localIndex> & connections )
{

  // using smart pointers for the trees
  std::unique_ptr<TreeNodeMortar> treeMaster = std::make_unique<TreeNodeMortar>();
  std::unique_ptr<TreeNodeMortar> treeSlave = std::make_unique<TreeNodeMortar>();

  // get list of surfaces
  array1d<localIndex> surfMaster = m_faceTypeToElementList.at(MortarSide::Master).at(masterShape);
  array1d<localIndex> surfSlave = m_faceTypeToElementList.at(MortarSide::Slave).at(slaveShape);

  if (surfMaster.size() == 0 || surfSlave.size() == 0)
  {
    return;
  }

  
  array1d<localIndex> surfRootMaster(surfMaster.size());
  array1d<localIndex> surfRootSlave(surfSlave.size());
  for (localIndex i=0; i<surfMaster.size(); ++i)
  {
    surfRootMaster(i) = i;
  }

  for (localIndex i=0; i<surfSlave.size(); ++i)
  {
    surfRootSlave(i) = i;
  }

  // create binary trees recursively
  std::cout << "Populating master tree" << std::endl;
  treeMaster->createNode( *m_mortarSide.at(MortarSide::Master).mesh,
                           m_mortarSide.at(MortarSide::Master).surface->getUniqueSubRegion< FaceElementSubRegion >(),
                           surfMaster,
                           surfRootMaster);
  std::cout << "Populating slave tree" << std::endl;
  treeSlave->createNode( *m_mortarSide.at(MortarSide::Slave).mesh,
                          m_mortarSide.at(MortarSide::Slave).surface->getUniqueSubRegion< FaceElementSubRegion >(),
                          surfSlave,
                          surfRootSlave);

  // initialize connectivity map
  connections.resize(surfRootMaster.size());

  // perform contact search and populate connectivity map
  contactSearch(treeSlave, treeMaster, connections);

  //localIndex numConnections = 0;

  for (localIndex i=0; i<surfRootMaster.size(); ++i)
  {
    localIndex N = connections.sizeOfArray(i);
    //numConnections += connectivityMap.sizeOfArray(i);

    std::cout << "MASTER ELEMENT: " << i << std::endl;
    std::cout << "SLAVE ELEMENT:"; 
    for (localIndex j=0; j<N; ++j)
    {  
      std::cout << " " << connections[i][j];
    }
    std::cout << std::endl;
    std::cout << "__________________________________________" << std::endl;
  }

 // return numConnections;
}

void SolidMechanicsMortarContact::contactSearch(std::unique_ptr<TreeNodeMortar> const & nodeSlave, 
                                                std::unique_ptr<TreeNodeMortar> const & nodeMaster, 
                                                ArrayOfArrays<localIndex> & connections)
{
  if (!checkIntersection(nodeSlave,nodeMaster)) return;
  
  if ((nodeSlave->isLeaf) && (nodeMaster->isLeaf))
  {
    //std::cout << "found intersection between slave " << nodeSlave->leafId << " and master " << nodeMaster->leafId << std::endl;
    connections.emplaceBack(nodeMaster->leafId,nodeSlave->leafId);
    return;
  }

  if (!nodeSlave->isLeaf)
  {
    contactSearch(nodeSlave->left, nodeMaster, connections);
    contactSearch(nodeSlave->right, nodeMaster, connections);
  }
  else
  {
    contactSearch(nodeSlave, nodeMaster->left, connections);
    contactSearch(nodeSlave, nodeMaster->right, connections);
  }
}

bool SolidMechanicsMortarContact::checkIntersection(std::unique_ptr<TreeNodeMortar> const & nodeSlave, std::unique_ptr<TreeNodeMortar> const & nodeMaster)
{
  double* pS = nodeSlave->polytop;
  double* pM = nodeMaster->polytop;

  for (localIndex k=0; k<9; ++k)
  {
    if (pS[2*k]>=pM[2*k+1]) return false;
    if (pM[2*k]>=pS[2*k+1]) return false;
  }

  // if pass all checks, the two nodes intersect
  return true;
}



void TreeNodeMortar::createNode( MeshLevel const & mesh, 
                                 FaceElementSubRegion const & surf, 
                                 arrayView1d<localIndex> & surfId, 
                                 array1d<localIndex> & surfList ) 
{
  
  //GEOS_UNUSED_VAR(surfList);
  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager const & nodeManager = mesh.getNodeManager();
  arrayView2d< localIndex const > const elemsToFaces = surf.faceList().toViewConst();
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();
  arrayView2d< double const > const surfCenter = faceManager.faceCenter().toViewConst();
  localIndex nSurf = surfList.size();
  arrayView2d<double const> const coords =  nodeManager.referencePosition();

  double const pP[3][9] = POLYTOP_PRIMITIVES;
  
  // compute primitive of polytopal bounding box for input list of surfaces
  for (localIndex i=0; i<nSurf; ++i)
  {
    localIndex kface = elemsToFaces[surfId[surfList[i]]][0]; 
    localIndex nNodes = faceToNodeMap.sizeOfArray( kface );

    for (localIndex j=0; j<nNodes; ++j)
    { 
      localIndex id = faceToNodeMap[kface][j];

      for (localIndex k = 0; k<9; ++k)
      {
        real64 P = coords[id][0]*pP[0][k] + 
                   coords[id][1]*pP[1][k] +
                   coords[id][2]*pP[2][k];
        // store minimum polytop primitive
        if (P < this->polytop[2*k]) this->polytop[2*k] = P;     
        // store maximum polytop primitive 
        if (P > this->polytop[2*k+1]) this->polytop[2*k+1] = P; 
      }
    }
  }

  // get direction of polytop maximum size
  localIndex dir = 0;
  real64 diff = -1;
  for (localIndex k = 0; k<9; ++k)
  {
    // add a small tolerance to avoid degenerate boxes 
    real64 d = polytop[2*k+1] - polytop[2*k] + 1e-3;

    // expand the polytop for safety
    polytop[2*k] -= BOUNDING_BOX_EXPANSION*d; 
    polytop[2*k+1] += BOUNDING_BOX_EXPANSION*d; 

    if (d > diff) 
    {
      dir = k;
      diff = d;
    } 
  }

  if (surfList.size() == 1) // leaf node
  { 
    this -> isLeaf = true;
    this -> leafId = surfList[0];
    return;
  }

  // split the surface list into left and right child
  std::vector<real64> surfacePrimitive(nSurf);
  for (localIndex i=0; i<nSurf; ++i)
  {
    // compute the polytopal primitives of the surface centroids
    localIndex kface = elemsToFaces[surfList[i]][0];
    surfacePrimitive[i] = surfCenter[kface][0]*pP[0][dir] + 
                           surfCenter[kface][1]*pP[1][dir] +
                           surfCenter[kface][2]*pP[2][dir];
  }

  // get median m of surface centroid primitives
  std::vector<real64> surfSort = surfacePrimitive;
  std::sort(surfSort.begin(), surfSort.end());
  real64 m = (nSurf % 2 == 1) ? surfSort[nSurf / 2] :  (surfSort[nSurf / 2 - 1] + surfSort[nSurf / 2]) / 2.0;
  
  array1d<localIndex> surfLeft;
  array1d<localIndex> surfRight;

  // divide left and right set of surface based on primitive value w.r.t median
  for (localIndex i=0; i<nSurf; ++i)
  {
    (surfacePrimitive[i] <= m) ? surfLeft.emplace_back(surfList[i]) : surfRight.emplace_back(surfList[i]);
  }

  /*
  for (localIndex i = 0; i < surfLeft.size(); ++i)
  {
    std::cout << "Left child" << std::endl;
    std::cout << surfLeft[i] << std::endl;
  }

    for (localIndex i = 0; i < surfRight.size(); ++i)
  {
    std::cout << "Right child" << std::endl;
    std::cout << surfRight[i] << std::endl;
  }
  std::cout << "_____________________________________________" << std::endl;
  */

  //  create left and right child recursively
  this->left = std::make_unique<TreeNodeMortar>();
  this->left->createNode(mesh,surf,surfId,surfLeft);
  this->right = std::make_unique<TreeNodeMortar>();
  this->right->createNode(mesh,surf,surfId,surfRight);
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SolidMechanicsMortarContact, string const &, Group * const )
}
 /* namespace geos */
