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
#include "physicsSolvers/solidMechanics/contact/LogLevelsInfo.hpp"
#include "physicsSolvers/solidMechanics/contact/ContactFields.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "constitutive/contact/FrictionSelector.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;

SolidMechanicsMortarContact::SolidMechanicsMortarContact( const string & name,
                                                          Group * const parent ):
  ContactSolverBase( name, parent )
{
  // higher order quadrature for element-based mortar integration
  m_faceTypeToMortarFiniteElements[ElementShape::Quadrilateral] = std::make_unique< finiteElement::H1_QuadrilateralFace_Lagrange1_GaussLegendre6 >();
  m_faceTypeToMortarFiniteElements[ElementShape::Triangle] = std::make_unique< finiteElement::H1_TriangleFace_Lagrange1_Gauss6 >();

  m_faceTypeToFiniteElements["Quadrilateral"] = std::make_unique< finiteElement::H1_QuadrilateralFace_Lagrange1_GaussLegendre6 >();
  m_faceTypeToFiniteElements["Triangle"] = std::make_unique< finiteElement::H1_TriangleFace_Lagrange1_Gauss6 >();

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

void SolidMechanicsMortarContact::setupDofs( DomainPartition const & domain,
                                            DofManager & dofManager,
                                            string const & nameSlave  ) const
{
  GEOS_MARK_FUNCTION;

  SolidMechanicsLagrangianFEM::setupDofs( domain, dofManager );

  map< std::pair< string, string >, string_array > meshTargets;
  
  // bubble and tractions leave on the slave side of the surface
  string_array regions;
  ElementRegionManager const & elementRegionManager = m_meshSlave->getElemManager();
  elementRegionManager.forElementRegions< SurfaceElementRegion >([&]( SurfaceElementRegion const & region )
    {
      regions.emplace_back( region.getName() );
    } );

  meshTargets[std::make_pair( nameSlave, m_meshSlave->getName())] = std::move( regions );


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

  setMortarSurfaces(domain);

  //string const meshMasterName;
  string meshSlaveName;

   //////////////////////////////////////////////////////////////////////////////////
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    // naive way to assign master and slave mesh levels
    if (meshName == "meshMaster"){
      this->m_meshMaster = &mesh;
    }

    if (meshName == "meshSlave"){
      this->m_meshSlave = &mesh;
      meshSlaveName = meshName;
    }
    
    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< FaceElementSubRegion >([&](const FaceElementSubRegion & subRegion )
    {
      // TO DO: assign different subregions for triangles and quadrilaterals?

      // assign surface to master or slave depending on object path (again, naive)
      string surfacePath = subRegion.getPath();
      std::cout << surfacePath << std::endl;
      if (surfacePath.find("surfaceSlave") != std::string::npos){
        GEOS_ERROR_IF(m_surfaceSlave != nullptr,"Slave surface has already been assigned.");
        this->m_surfaceSlave = &subRegion;
      }

      if (surfacePath.find("surfaceMaster") != std::string::npos){
        GEOS_ERROR_IF(m_surfaceMaster != nullptr,"Master surface has already been assigned.");
        this->m_surfaceMaster = &subRegion;
      }
    });

  } );

  createFaceTypeListMortar(MortarSide::Master);
  createFaceTypeListMortar(MortarSide::Slave);

  std::cout << "Number of master quad cells: " <<  m_faceTypeToElementList[MortarSide::Master][ElementShape::Quadrilateral].size() << std::endl;
  std::cout << "Number of master tri cells: " <<  m_faceTypeToElementList[MortarSide::Master][ElementShape::Triangle].size() << std::endl;
  std::cout << "Number of slave quad cells: " <<  m_faceTypeToElementList[MortarSide::Slave][ElementShape::Quadrilateral].size() << std::endl;
  std::cout << "Number of slave tri cells: " <<  m_faceTypeToElementList[MortarSide::Slave][ElementShape::Triangle].size() << std::endl;

  std::cout << "Processing master/slave connectivity" << std::endl;

  this->getConnectivityMap();

  m_faceTypesToFaceElementsMaster = this->createFaceTypeList(*m_meshMaster, *m_surfaceMaster, "master");
  m_faceTypesToFaceElementsSlave = this->createFaceTypeList(*m_meshSlave, *m_surfaceSlave, "slave");

  // create list of bubbles for the slave side
  this->createBubbleCellList();

  
  dofManager.setDomain( domain );
  setupDofs( domain, dofManager, meshSlaveName);
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
  this->addCouplingNumNonzeros( dofManager, rowLengths.toView() );

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
  this->addCouplingSparsityPattern( dofManager, pattern.toView() );

  // Finally, steal the pattern into a CRS matrix
  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );
  localMatrix.setName( this->getName() + "/localMatrix" );

  rhs.setName( this->getName() + "/rhs" );
  rhs.create( dofManager.numLocalDofs(), MPI_COMM_GEOS );

  solution.setName( this->getName() + "/solution" );
  solution.create( dofManager.numLocalDofs(), MPI_COMM_GEOS );


  /* Write the matrix pattern to a file for debugging
  {
    ParallelMatrix parallel_matrix;
    parallel_matrix.create( localMatrix.toViewConst(), dofManager.numLocalDofs(), MPI_COMM_GEOS );
    parallel_matrix.write("mortar_sparsity_3.mtx");
  }
  */


  // compute mortar interpolation: populate maps from gauss points to master basis functions
  //this -> computeMortarInterpolation();


  forMortarSurfaces(MortarInterpolation{});

  //computeRotationMatrices( domain );
  
  //////////////////////////////////////////////////////////////////////////////////
  GEOS_ERROR("Mortar solver is not implemented yet. ");

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

  SolidMechanicsLagrangianFEM::assembleSystem( time,
                                               dt,
                                               domain,
                                               dofManager,
                                               localMatrix,
                                               localRhs );

  //ParallelMatrix parallel_matrix;
  //parallel_matrix.create( localMatrix.toViewConst(), dofManager.numLocalDofs(), MPI_COMM_GEOS );
  //parallel_matrix.write("mech.mtx");

}

void SolidMechanicsMortarContact::implicitStepComplete( real64 const & time_n,
                                                                     real64 const & dt,
                                                                     DomainPartition & domain )
{

  SolidMechanicsLagrangianFEM::implicitStepComplete( time_n, dt, domain );

}

real64 SolidMechanicsMortarContact::calculateResidualNorm( real64 const & time,
                                                                        real64 const & dt,
                                                                        DomainPartition const & domain,
                                                                        DofManager const & dofManager,
                                                                        arrayView1d< real64 const > const & localRhs )
{

  GEOS_MARK_FUNCTION;
  real64 const solidResidualNorm = SolidMechanicsLagrangianFEM::calculateResidualNorm( time, dt, domain, dofManager, localRhs );

  return solidResidualNorm;

}

void SolidMechanicsMortarContact::applySystemSolution( DofManager const & dofManager,
                                                                    arrayView1d< real64 const > const & localSolution,
                                                                    real64 const scalingFactor,
                                                                    real64 const dt,
                                                                    DomainPartition & domain )
{

  GEOS_MARK_FUNCTION;
  SolidMechanicsLagrangianFEM::applySystemSolution( dofManager,
                                                    localSolution,
                                                    scalingFactor,
                                                    dt,
                                                    domain );

}



void SolidMechanicsMortarContact::addCouplingNumNonzeros( DofManager & dofManager,
                                                          arrayView1d< localIndex > const & rowLengths ) const
{

  ElementRegionManager const & elemManagerSlave = m_meshSlave->getElemManager();
  NodeManager const & nodeManagerSlave = m_meshSlave->getNodeManager();
  FaceManager const & faceManagerSlave = m_meshSlave->getFaceManager();

  FaceElementSubRegion const & subRegionSlave = *m_surfaceSlave;
  FaceElementSubRegion const & subRegionMaster = *m_surfaceMaster;

  //ArrayOfArraysView< localIndex const > const faceToNodeMapSlave = faceManagerSlave.nodeList().toViewConst();


  //ElementRegionManager const & elemManagerMaster = m_meshMaster->getElemManager();
  NodeManager const & nodeManagerMaster = m_meshMaster->getNodeManager();
  FaceManager const & faceManagerMaster = m_meshMaster->getFaceManager();

  ArrayOfArraysView< localIndex const > const faceToNodeMapMaster = faceManagerMaster.nodeList().toViewConst();
  arrayView2d< localIndex const > const elemsToFacesMaster = subRegionMaster.faceList().toViewConst();

  globalIndex const rankOffset = dofManager.rankOffset();

  string const bubbleDofKey = dofManager.getKey( contact::totalBubbleDisplacement::key() );
  string const dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
  string const tractionDofKey = dofManager.getKey( contact::traction::key() );

  arrayView1d< globalIndex const > const bubbleDofNumber = faceManagerSlave.getReference< globalIndex_array >( bubbleDofKey );
  arrayView1d< globalIndex const > const dispSlaveDofNumber =  nodeManagerSlave.getReference< globalIndex_array >( dispDofKey );
  arrayView1d< globalIndex const > const dispMasterDofNumber =  nodeManagerMaster.getReference< globalIndex_array >( dispDofKey );
  arrayView1d< globalIndex const > const tractionDofNumber = subRegionSlave.getReference< globalIndex_array >( tractionDofKey );

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

  // add coupling between tractions and master displacements

  for( localIndex kfe=0; kfe<subRegionSlave.size(); ++kfe )
  {

    localIndex const nMaster = m_connectivityMapSlave.sizeOfArray(kfe);
    localIndex numDispDof = 0;

    // loop over connected master elements
    for (int i=0; i<nMaster; ++i)
    {
      localIndex kmaster = m_connectivityMapSlave(kfe,i);
      localIndex const kface = elemsToFacesMaster[kmaster][0];
      localIndex const numNodesPerFace = faceToNodeMapMaster.sizeOfArray( kface );

      for( localIndex a=0; a<numNodesPerFace; ++a )
      {
         const localIndex & node = faceToNodeMapMaster( kface, a );
        localIndex const localDispRow = LvArray::integerConversion< localIndex >( dispMasterDofNumber[node] - rankOffset );

        if( localDispRow >= 0 && localDispRow < rowLengths.size() )
        {
          for( int d=0; d<3; ++d )
          {
            rowLengths[localDispRow + d] += 3;
          }
        }
      }

      numDispDof += 3*numNodesPerFace;
    }

    localIndex const localRow = LvArray::integerConversion< localIndex >( tractionDofNumber[kfe] - rankOffset );

    if( localRow >= 0 && localRow < rowLengths.size() )
    {
      for( localIndex i=0; i<3; ++i )
      {
        rowLengths[localRow + i] += numDispDof;
      }
    }
  }
}


void SolidMechanicsMortarContact::addCouplingSparsityPattern( DofManager const & dofManager,
                                                              SparsityPatternView< globalIndex > const & pattern ) const
{

  std::cout<<"Creating sparsity pattern" << std::endl;
    
  ElementRegionManager const & elemManagerSlave = m_meshSlave->getElemManager();
  NodeManager const & nodeManagerSlave = m_meshSlave->getNodeManager();
  FaceManager const & faceManagerSlave = m_meshSlave->getFaceManager();

  FaceElementSubRegion const & subRegionSlave = *m_surfaceSlave;
  FaceElementSubRegion const & subRegionMaster = *m_surfaceMaster;


  //ArrayOfArraysView< localIndex const > const faceToNodeMapSlave = faceManagerSlave.nodeList().toViewConst();

  //ElementRegionManager const & elemManagerMaster = m_meshMaster->getElemManager();
  NodeManager const & nodeManagerMaster = m_meshMaster->getNodeManager();
  FaceManager const & faceManagerMaster = m_meshMaster->getFaceManager();

  ArrayOfArraysView< localIndex const > const faceToNodeMapMaster = faceManagerMaster.nodeList().toViewConst();
    arrayView2d< localIndex const > const elemsToFacesMaster = subRegionMaster.faceList().toViewConst();

  globalIndex const rankOffset = dofManager.rankOffset();

  string const bubbleDofKey = dofManager.getKey( contact::totalBubbleDisplacement::key() );
  string const dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
  string const tractionDofKey = dofManager.getKey( contact::traction::key() );

  arrayView1d< globalIndex const > const bubbleDofNumber = faceManagerSlave.getReference< globalIndex_array >( bubbleDofKey );
  arrayView1d< globalIndex const > const dispSlaveDofNumber =  nodeManagerSlave.getReference< globalIndex_array >( dispDofKey );
  arrayView1d< globalIndex const > const dispMasterDofNumber =  nodeManagerMaster.getReference< globalIndex_array >( dispDofKey );
  arrayView1d< globalIndex const > const tractionDofNumber = subRegionSlave.getReference< globalIndex_array >( tractionDofKey );

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
          eqnRowIndicesDisp[3*a + idof] = dispSlaveDofNumber[node] + idof - rankOffset;
          dofColIndicesDisp[3*a + idof] = dispSlaveDofNumber[node] + idof;
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

  // populate sparsity pattern for coupling between tractions and master elements
  
  static constexpr int maxNumDispFaceDof = 3 * 4;

  for( localIndex kfe=0; kfe<subRegionSlave.size(); ++kfe )
  {
    localIndex const nMaster = m_connectivityMapSlave.sizeOfArray(kfe);

    for (int m=0; m<nMaster; ++m)
    {
      localIndex kmaster = m_connectivityMapSlave(kfe,m);
      localIndex const kface = elemsToFacesMaster[kmaster][0];
      localIndex const numNodesPerFace = faceToNodeMapMaster.sizeOfArray( kface );
      localIndex const numDispDof = 3*numNodesPerFace;

      // working arrays
      stackArray1d< globalIndex, maxNumDispFaceDof > eqnRowIndicesDisp ( numDispDof );
      stackArray1d< globalIndex, 3 > eqnRowIndicesTraction( 3 );
      stackArray1d< globalIndex, maxNumDispFaceDof > dofColIndicesDisp( numDispDof );
      stackArray1d< globalIndex, 3 > dofColIndicesTraction( 3 );

      for( localIndex idof = 0; idof < 3; ++idof )
      {
        eqnRowIndicesTraction[idof] = tractionDofNumber[kfe] + idof - rankOffset;
        dofColIndicesTraction[idof] = tractionDofNumber[kfe] + idof;
      }

      for( localIndex a=0; a<numNodesPerFace; ++a )
      {
        const localIndex & node = faceToNodeMapMaster( kface, a );
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

  /*
  for( localIndex kfe=0; kfe<subRegion.size(); ++kfe )
  {
    localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kfe );
    localIndex const numDispDof = 3*numNodesPerFace;
    for( int k=0; k<2; ++k )
    {
      localIndex const kf = elemsToFaces[kfe][k];
      localIndex const kf_other = elemsToFaces[kfe][(1+k)%2];
      // working arrays
      stackArray1d< globalIndex, maxNumDispFaceDof > eqnRowIndicesDisp ( numDispDof );
      stackArray1d< globalIndex, 3 > eqnRowIndicesBubble( 3 );
      stackArray1d< globalIndex, maxNumDispFaceDof > dofColIndicesDisp ( numDispDof );
      stackArray1d< globalIndex, 3 > dofColIndicesBubble( 3 );
      for( localIndex idof = 0; idof < 3; ++idof )
      {
        eqnRowIndicesBubble[idof] = bubbleDofNumber[kf] + idof - rankOffset;
        dofColIndicesBubble[idof] = bubbleDofNumber[kf] + idof;
      }
      for( localIndex a=0; a<numNodesPerFace; ++a )
      {
        const localIndex & node = faceToNodeMap( kf_other, a );
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
          }
        }
      }
    }
  }
  */

}





void SolidMechanicsMortarContact::updateState( DomainPartition & domain )
{
  GEOS_UNUSED_VAR( domain );
}



SolidMechanicsMortarContact::FaceTypeMap
SolidMechanicsMortarContact::createFaceTypeList( MeshLevel const & mesh,  
                                                 FaceElementSubRegion const & surfRegion,
                                                 string meshName)
{

  // Generate lists containing elements of various face types

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

    SolidMechanicsMortarContact::FaceTypeMap faceTypeList;

    faceTypeList[meshName]["Quadrilateral"] =  quadList;
    //faceTypeList[meshName]["Triangle"] =  triList;

    return faceTypeList;

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
      // // naive way to assign master and slave mesh levels
      // if (meshName == "meshMaster"){
      //   this->m_meshMaster = &mesh;
      // }

      // if (meshName == "meshSlave"){
      //   this->m_meshSlave = &mesh;
      //   meshSlaveName = meshName;
      // }
      GEOS_UNUSED_VAR(meshName);

      ElementRegionManager const & elemManager = mesh.getElemManager();

      elemManager.forElementRegions< SurfaceElementRegion >([&,this](const SurfaceElementRegion & region )
      {
        // TO DO: assign different subregions for triangles and quadrilaterals?

        // assign surface to master or slave depending on object path (again, naive)
        string surfacePath = region.getPath();
        std::cout << surfacePath << std::endl;

        // check if region is master or slave and populate member maps
        if (surfacePath.find(m_slaveName) != std::string::npos)
        {
          MortarSurface surfaceSlave;
          surfaceSlave.mesh = &mesh;
          surfaceSlave.surface = &region;
          m_mortarSide[MortarSide::Slave] = surfaceSlave;
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

void SolidMechanicsMortarContact::computeMortarInterpolation()
{

  FaceManager const & faceManagerMaster = m_meshMaster->getFaceManager();
  NodeManager const & nodeManagerMaster = m_meshMaster->getNodeManager();
  arrayView2d< localIndex const > const elemsToFacesMaster = m_surfaceMaster->faceList().toViewConst();
  ArrayOfArraysView< localIndex const > const & faceToNodeMapMaster = faceManagerMaster.nodeList().toViewConst();
  arrayView2d<double const> const nodeCoordsMaster =  nodeManagerMaster.referencePosition();

  FaceManager const & faceManagerSlave = m_meshSlave->getFaceManager();
  NodeManager const & nodeManagerSlave = m_meshSlave->getNodeManager();
  arrayView2d< double const > const slaveCenters = faceManagerSlave.faceCenter().toViewConst();
  arrayView2d< localIndex const > const elemsToFacesSlave = m_surfaceSlave->faceList().toViewConst();
  ArrayOfArraysView< localIndex const > const & faceToNodeMapSlave = faceManagerSlave.nodeList().toViewConst();
  arrayView2d<double const> const nodeCoordsSlave =  nodeManagerSlave.referencePosition();


  std::cout<< "Interpolating master basis functions" << std::endl;
  localIndex nInt = 5; // number of interpolation points in each direction of the reference element

  // initialize gp maps assuming only quadrilaterals for the moment (will be replaced by a loop over face types)
  finiteElement::FiniteElementBase const & subRegionSlaveFE = *(m_faceTypeToFiniteElements.at( "Quadrilateral" ));
  localIndex nGP = subRegionSlaveFE.getNumQuadraturePoints();
  std::cout << "Number of quadrature points (slave): " << nGP << std::endl;
  localIndex nSlaveElements = m_faceTypesToFaceElementsSlave["slave"]["Quadrilateral"].size();
  m_gpToMasterId["Quadrilateral"] = array2d< localIndex >( nSlaveElements, nGP );
  arrayView2d< localIndex > gpToMasterId_v = m_gpToMasterId["Quadrilateral"].toView();

  std:: cout << "size of gp to master map: " << gpToMasterId_v.size(0) << " x " << gpToMasterId_v.size(1) << std::endl;

  // prepare list of local coordinates, based on element type quad/triangle.
  for( const auto & [finiteElementName, faceElementList] : m_faceTypesToFaceElementsMaster.at("master") )
  {
    array2d< real64 > localCoordsInterpolation;

    this->getLocalInterpolationPoints( nInt , finiteElementName, localCoordsInterpolation );

    localIndex nIntPts = localCoordsInterpolation.size(0);

    arrayView1d< localIndex const > const faceElemList = faceElementList.toViewConst();

    //finiteElement::FiniteElementBase const & subRegionMasterFE = *(m_faceTypeToFiniteElements.at( finiteElementName ));

    // must be dispatched at compile time!
    constexpr localIndex numNodeperElement = finiteElement::H1_QuadrilateralFace_Lagrange1_GaussLegendre6::numNodes;

    int permutation[numNodeperElement];
    finiteElement::H1_QuadrilateralFace_Lagrange1_GaussLegendre6::getPermutation( permutation );

    m_gpToMasterBasis["Quadrilateral"] = ArrayOfArrays< real64 >(nGP*nSlaveElements, numNodeperElement);
    ArrayOfArraysView< real64 > gpToMasterBasis_v = m_gpToMasterBasis["Quadrilateral"].toView();

    for ( localIndex iFE = 0; iFE < faceElemList.size(); ++iFE )
    {
      std::cout << std::endl << "MASTER ELEMENT " << iFE << std::endl << std::endl;
      // get interpolation points in global coordinates
      array2d<real64> realCoordsInterpolation(nIntPts, 3);
      localIndex const kfe = faceElemList[iFE];
      localIndex kface = elemsToFacesMaster[kfe][0];

      real64 N[numNodeperElement]; // shape functions at a given point
      array2d<real64> Nm(nIntPts, numNodeperElement+1); // shape functions at interpolation points

      for ( localIndex i=0; i<nIntPts; ++i )
      {
        // get shape functions at interpolation point
        real64 localCoords[2] = {localCoordsInterpolation(i,0), localCoordsInterpolation(i,1)};
        //std::cout << localCoordsInterpolation(i,0) << " " << localCoordsInterpolation(i,1) << std::endl;
        finiteElement::H1_QuadrilateralFace_Lagrange1_GaussLegendre6::calcN( localCoords, N );
        for ( localIndex a=0; a<numNodeperElement; ++a )
        {
          localIndex const node = faceToNodeMapMaster(kface,a);
          for ( localIndex d=0; d<3; ++d )
          {
            realCoordsInterpolation[i][d] += N[permutation[a]]*nodeCoordsMaster(node,d);
          }

          Nm[i][a] = N[permutation[a]];
        }
        Nm[i][numNodeperElement] = 1.0; // unit function for RBF rescaling
        //std::cout << "IntPt " << i << " Coords " << realCoordsInterpolation[i][0] << " " << realCoordsInterpolation[i][1] << " " << realCoordsInterpolation[i][2] << std::endl;
        //std::cout << "IntPt " << i << " Basis function " << Nm[i][0] << " " << Nm[i][1] << " " << Nm[i][2] << " " << Nm[i][3] << std::endl;
      }

      array2d<real64> RBFmatrix(nIntPts, nIntPts);
      real64 radius = this->computeRBFmatrix( realCoordsInterpolation.toView(), RBFmatrix.toView(), nIntPts );

      // compute RBF weights
      array2d<real64> weightsRBF(nIntPts, numNodeperElement+1);

      this -> computeRBFweights( RBFmatrix.toView(), Nm.toView(), weightsRBF.toView(), nIntPts, numNodeperElement );
      
      // get connected slave elements
      localIndex numConnectedSlaveElements = m_connectivityMapMaster.sizeOfArray(kfe);

      ArrayOfArraysView< localIndex >  m_connectivityMapMaster_v = m_connectivityMapMaster.toView();

      real64 gpCoord[3];

      for (localIndex j=0; j<numConnectedSlaveElements; ++j)
      {
        std::cout <<  m_connectivityMapMaster.capacityOfArray(kfe) << std::endl;
        localIndex kslave = m_connectivityMapMaster(kfe,j);
        std::cout << m_connectivityMapMaster_v.getOffsets()[kfe] << std::endl;
        localIndex kfaceSlave = elemsToFacesSlave[kslave][0];
        //std::cout << "SLAVE ELEMENT: " << kslave << " CENTROID " << slaveCenters(kfaceSlave,0) << " " << slaveCenters(kfaceSlave,1) << " " << slaveCenters(kfaceSlave,2) << std::endl;
        //arrayView1d< localIndex > masterList =  gpToMasterId_v.slice( kslave, RAJA::ALL );

        // loop over gauss points
        for (localIndex q=0; q<nGP; ++q)
        {
          gpCoord[0] = 0.0;
          gpCoord[1] = 0.0;
          gpCoord[2] = 0.0;
          //std::cout << std::endl;
          //std::cout << " SLAVE " << kslave << " GP " << q << " ";
          if (gpToMasterId_v(kslave,q) != 0) continue; // already computed

          real64 Nslave[numNodeperElement];
          finiteElement::H1_QuadrilateralFace_Lagrange1_GaussLegendre6::calcN( q, Nslave );

          for (localIndex a=0; a<numNodeperElement; ++a)
          {
            //std::cout << Nslave[permutation[a]] << " ";
            localIndex const node = faceToNodeMapSlave(kfaceSlave,a);
            for ( localIndex d=0; d<3; ++d )
            {
              gpCoord[d] += Nslave[permutation[a]]*nodeCoordsSlave(node,d);
            }
            //std::cout << nodeCoordsSlave(node,0) << " " << nodeCoordsSlave(node,1) << " " << nodeCoordsSlave(node,2) << std::endl;
          }

          //std::cout << gpCoord[0] << " " << gpCoord[1] << " " << gpCoord[2] << std::endl;

          real64 f[numNodeperElement] = {};
          real64 f1 = 0.0;

          // interpolate basis functions for all nodes
          for (localIndex i=0; i < nIntPts; ++i)
          {
            real64 intCoords[3] = {realCoordsInterpolation[i][0], realCoordsInterpolation[i][1], realCoordsInterpolation[i][2]}; 
            real64 rbfEntry = computeRBF( gpCoord, intCoords, radius );
            f1 += rbfEntry * weightsRBF[i][numNodeperElement];
            for (localIndex a=0; a<numNodeperElement; ++a)
            {
              f[a] += computeRBF( gpCoord, intCoords, radius ) * weightsRBF[i][a];
            }
          }

          bool isInSupport = true;

          // support detection
          for (localIndex a=0; a<numNodeperElement; ++a)
          {
            f[a] /= f1; // rescaling RBF
            if (f[a] < -1e-4 || f[a] > 1+1.1e-4) 
            {
              isInSupport = false;
              //std::cout << "SLAVE " << kslave << " - GP " << q << " - " << gpCoord[0] << " " << gpCoord[1] << " " << gpCoord[2] << " || "<<  " ND " << std::endl;
              break; // gp not in support of master element
            }
          }

          if (isInSupport)
          {
            //std::cout << "SLAVE " << kslave << " - GP " << q << " " << gpCoord[0] << " " << gpCoord[1] << " " << gpCoord[2] << " ||";
            for (localIndex a=0; a<numNodeperElement; ++a)
            {
            gpToMasterBasis_v[nGP*kslave+q][a] = f[a];
            std::cout << " " << f[a];
            //masterList[q] = kfe;
            }
            std::cout << std::endl;
          }
        }
      }
    }
  }

} 


void SolidMechanicsMortarContact::getLocalInterpolationPoints( localIndex nInt, string const & finiteElementName, array2d< real64 > & localCoordsMaster )
{
  if (finiteElementName == "Quadrilateral")
  {
    localCoordsMaster.resize(nInt*nInt,2);
    arrayView2d< real64 > const localCoordsMaster_v = localCoordsMaster.toView();

    localIndex k=0;
    for (localIndex i=0; i<nInt; ++i)
    {
      real64 xi = -1.0 + 2.0*i/(nInt-1);
      for (localIndex j=0; j<nInt; ++j)
      {
        real64 eta = -1.0 + 2.0*j/(nInt-1);
        localCoordsMaster_v(k,0) = eta;
        localCoordsMaster_v(k,1) = xi;
        ++k;
      }
    }
  }
  else if (finiteElementName == "Triangle")
  {
    GEOS_ERROR( "SolidMechanicsMortarContact:: triangular face type not yet available" );
  }
  else
  {
    GEOS_ERROR( "SolidMechanicsMortarContact:: invalid face type" );
  }
}

real64 SolidMechanicsMortarContact::computeRBFmatrix( arrayView2d<real64> realCoordsInterpolation, arrayView2d<real64> RBFmatrix, localIndex nIntPts)
{

  array2d<real64> distanceMatrix(nIntPts, nIntPts);
  real64 radius = 0.0;                    // RBF radius (maximum distance between points)


  for (localIndex i=0; i<nIntPts; ++i)
  {
    for (localIndex j=0; j<nIntPts; ++j)
    {
      real64 d = sqrt( (realCoordsInterpolation[i][0]-realCoordsInterpolation[j][0])*(realCoordsInterpolation[i][0]-realCoordsInterpolation[j][0]) +
                (realCoordsInterpolation[i][1]-realCoordsInterpolation[j][1])*(realCoordsInterpolation[i][1]-realCoordsInterpolation[j][1]) +
                (realCoordsInterpolation[i][2]-realCoordsInterpolation[j][2])*(realCoordsInterpolation[i][2]-realCoordsInterpolation[j][2]) );
      
      if (d>radius) radius = d;
      distanceMatrix[i][j] = d;
    }
  }

  for (localIndex i=0; i<nIntPts; ++i)
  {
    for (localIndex j=0; j<nIntPts; ++j)
    {
      // gaussian spline RBF 
      RBFmatrix[i][j] = computeRBF( distanceMatrix[i][j], radius );
    }
  }

  return radius;
}

real64 SolidMechanicsMortarContact::computeRBF( real64 const d, real64 const radius )
{
  return exp(-d*d/(radius*radius));
}

real64 SolidMechanicsMortarContact::computeRBF( real64 const (&pt1)[3], real64 const (&pt2)[3], real64 const radius )
{
  real64 d = sqrt( (pt1[0]-pt2[0])*(pt1[0]-pt2[0]) +
                   (pt1[1]-pt2[1])*(pt1[1]-pt2[1]) +
                   (pt1[2]-pt2[2])*(pt1[2]-pt2[2]) );

  return exp(-d*d/(radius*radius));
}





void SolidMechanicsMortarContact::computeRBFweights( arrayView2d<real64> M, arrayView2d<real64> Nm, arrayView2d<real64> wRBF, localIndex nIntPts, localIndex numNodeperElement )
{

  std::cout << "Computing RBF weights.: FACTORIZATION" << std::endl;

  // Choleskij factorization
  for (localIndex k=0; k<nIntPts; ++k)
  {
    M[k][k] = sqrt(M[k][k]);
    for (localIndex i=k+1; i<nIntPts; ++i)
    {
      M[i][k] = M[i][k]/M[k][k];
    }
    for (localIndex j=k+1; j<nIntPts; ++j)
    {
      for (localIndex i=j; i<nIntPts; ++i)
      {
        M[i][j] = M[i][j] - M[i][k]*M[j][k];
      }
    }
  }

  array1d<real64> y(nIntPts);
  array1d<real64> x(nIntPts);

  for (localIndex col=0; col<numNodeperElement+1; ++col)
  {
    // forward substitution to solve Ly = b
    for (localIndex i=0; i<nIntPts; ++i)
    {
      y[i] = Nm[i][col];
      for (localIndex j=0; j<i; ++j)
      {
        y[i] -= M[i][j]*y[j];
      }
      y[i] = y[i]/M[i][i];
    }


    // backward substitution to solve L^Tx = y
    for (localIndex i=nIntPts-1; i>=0; --i)
    {
      x[i] = y[i];
      for (localIndex j=i+1; j<nIntPts; ++j)
      {
        x[i] -= M[j][i]*x[j];
      }
      x[i] = x[i]/M[i][i];
      wRBF[i][col] = x[i];
    }
  }

  std::cout << "RBF weights computed successfully." << std::endl;

}

void SolidMechanicsMortarContact::getConnectivityMap()
{
  // using smart pointers for the trees
  std::unique_ptr<TreeNodeMortar> treeMaster = std::make_unique<TreeNodeMortar>();
  std::unique_ptr<TreeNodeMortar> treeSlave = std::make_unique<TreeNodeMortar>();

  // progressive list of surfaces to create tree roots
  array1d<localIndex> surfRootMaster(m_surfaceMaster->size()) ;
  array1d<localIndex> surfRootSlave(m_surfaceSlave->size()) ;

  for (localIndex i=0; i < m_surfaceMaster->size(); ++i)
  {
    surfRootMaster(i) = i;
  }

  for (localIndex i=0; i < m_surfaceSlave->size(); ++i)
  {
    surfRootSlave(i) = i;
  }

  // create binary trees recursively
  std::cout << "Populating master tree" << std::endl;
  treeMaster->createNode(*m_meshMaster,*m_surfaceMaster,surfRootMaster);
  std::cout << "Populating slave tree" << std::endl;
  treeSlave->createNode(*m_meshSlave,*m_surfaceSlave,surfRootSlave);

  // initialize connectivity map
  localIndex nSurfSlave = m_surfaceSlave->size();
  localIndex nSurfMaster = m_surfaceMaster->size();
  m_connectivityMapSlave.resize(nSurfSlave);
  m_connectivityMapMaster.resize(nSurfMaster);

  // perform contact search and populate connectivity map
  contactSearch(treeSlave,treeMaster);

    for (localIndex i=0; i<nSurfMaster; ++i)
  {
    localIndex N = m_connectivityMapMaster.sizeOfArray(i);
    std::cout << "MASTER ELEMENT: " << i << std::endl;
    std::cout << "SLAVE ELEMENT:"; 
    for (localIndex j=0; j<N; ++j)
    {  
      std::cout << " " << m_connectivityMapMaster[i][j];
    }
    std::cout << std::endl;
    std::cout << "__________________________________________" << std::endl;
  }
}

void SolidMechanicsMortarContact::contactSearch(std::unique_ptr<TreeNodeMortar> const & nodeSlave, std::unique_ptr<TreeNodeMortar> const & nodeMaster)
{
  if (!checkIntersection(nodeSlave,nodeMaster)) return;
  
  if ((nodeSlave->isLeaf) && (nodeMaster->isLeaf))
  {
    //std::cout << "found intersection between slave " << nodeSlave->leafId << " and master " << nodeMaster->leafId << std::endl;
    m_connectivityMapSlave.emplaceBack(nodeSlave->leafId,nodeMaster->leafId);
    m_connectivityMapMaster.emplaceBack(nodeMaster->leafId,nodeSlave->leafId);
    return;
  }

  if (!nodeSlave->isLeaf)
  {
    contactSearch(nodeSlave->left,nodeMaster);
    contactSearch(nodeSlave->right,nodeMaster);
  }
  else
  {
    contactSearch(nodeSlave,nodeMaster->left);
    contactSearch(nodeSlave,nodeMaster->right);
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



void TreeNodeMortar::createNode(MeshLevel const & mesh, 
                FaceElementSubRegion const & surf, 
                array1d<localIndex> & surfList) 
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
    localIndex kface = elemsToFaces[surfList[i]][0];
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
  this->left->createNode(mesh,surf,surfLeft);
  this->right = std::make_unique<TreeNodeMortar>();
  this->right->createNode(mesh,surf,surfRight);
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SolidMechanicsMortarContact, string const &, Group * const )
}
 /* namespace geos */
