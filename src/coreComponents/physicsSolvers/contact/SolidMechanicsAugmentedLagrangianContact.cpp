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
 * SolidMechanicsAugmentedLagrangianContact.cpp
 */

#include "SolidMechanicsAugmentedLagrangianContact.hpp"

#include "physicsSolvers/contact/SolidMechanicsALMKernels.hpp"
#include "physicsSolvers/contact/SolidMechanicsALMJumpUpdateKernels.hpp"
#include "physicsSolvers/contact/SolidMechanicsALMBubbleKernels.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;

SolidMechanicsAugmentedLagrangianContact::SolidMechanicsAugmentedLagrangianContact( const string & name,
                                                                                    Group * const parent ):
  ContactSolverBase( name, parent )
{
  std::cout << "SolidMechanicsAugmentedLagrangianContact" << std::endl;

  //m_fe = 
  m_faceTypeToFiniteElements["Quadrilateral"] =  std::make_unique<finiteElement::H1_QuadrilateralFace_Lagrange1_GaussLegendre2>();
  m_faceTypeToFiniteElements["Triangle"] =  std::make_unique<finiteElement::H1_TriangleFace_Lagrange1_Gauss1>();
  
  // TODO Implement the constructor 

  // Set the default linear solver parameters
  //LinearSolverParameters & linParams = m_linearSolverParameters.get();
  //linParams.dofsPerNode = 3;
  //linParams.isSymmetric = true;
  //linParams.amg.separateComponents = true;
}

SolidMechanicsAugmentedLagrangianContact::~SolidMechanicsAugmentedLagrangianContact()
{
  std::cout << "~SolidMechanicsAugmentedLagrangianContact" << std::endl;
  // TODO Auto-generated destructor stub
}

void SolidMechanicsAugmentedLagrangianContact::registerDataOnMesh( dataRepository::Group & meshBodies )
{

  std::cout << "registerDataOnMesh" << std::endl;
  
  ContactSolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & meshLevel,
                                                    arrayView1d< string const > const & )
  {
    FaceManager & faceManager = meshLevel.getFaceManager();

    std::cout << "getName: " << getName() << std::endl;
    faceManager.registerField< solidMechanics::totalBubbleDisplacement >( getName() ).
      reference().resizeDimension< 1 >( 3 );
  });
  
  using namespace fields::contact;

  forFractureRegionOnMeshTargets( meshBodies, [&] ( SurfaceElementRegion & fractureRegion )
  {
    fractureRegion.forElementSubRegions< SurfaceElementSubRegion >( [&]( SurfaceElementSubRegion & subRegion )
    {

      // Register the rotation matrix
      subRegion.registerField< rotationMatrix >( this->getName() ).
        reference().resizeDimension< 1, 2 >( 3, 3 );
      //subRegion.registerWrapper< array3d< real64 > >( viewKeyStruct::rotationMatrixString() ).
      //  setPlotLevel( PlotLevel::NOPLOT ).
      //  setRegisteringObjects( this->getName()).
      //  setDescription( "An array that holds the rotation matrices on the fracture." ).
      //  reference().resizeDimension< 1, 2 >( 3, 3 );

      // Register the traction field
      subRegion.registerField< traction >( this->getName() ).
        reference().resizeDimension< 1 >( 3 );

      // Register the displacement jump
      subRegion.registerField< dispJump >( this->getName() ).
        reference().resizeDimension< 1 >( 3 );

      // Register the displacement jump
      subRegion.registerField< deltaDispJump >( this->getName() ).
        reference().resizeDimension< 1 >( 3 );

      // Register the displacement jump old
      subRegion.registerField< oldDispJump >( this->getName() ).
        reference().resizeDimension< 1 >( 3 );

      // Register the displacement jump
      subRegion.registerField< penalty >( this->getName() ).
        reference().resizeDimension< 1 >( 2 );

    } );
  } );
  
}

void SolidMechanicsAugmentedLagrangianContact::initializePreSubGroups()
{
  std::cout << "initializePreSubGroups" << std::endl;
  ContactSolverBase::initializePreSubGroups();
  std::cout << "end initializePreSubGroups" << std::endl;

}

//void SolidMechanicsAugmentedLagrangianContact::initializePostInitialConditionsPreSubGroups()
//{
//  std::cout << "initializePostInitialConditionsPreSubGroups" << std::endl;
//  SolidMechanicsLagrangianFEM::initializePostInitialConditionsPreSubGroups();

  //array1d< localIndex > quadList;
  //array1d< localIndex > triList;
  //quadList.resize(10);
  //quadList[0] = 10;
  //quadList[1] = 9;
  //std::cout << quadList[0] << std::endl;
  //this->m_faceTypesToFaceElements["Quadrilateral"] =  quadList;
  //this->m_faceTypesToFaceElements["Triangle"] =  quadList;


/*
  auto & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forFractureRegionOnMeshTargets( domain.getMeshBodies(), [&] ( SurfaceElementRegion & fractureRegion )
  {
    fractureRegion.forElementSubRegions< SurfaceElementSubRegion >( [&]( SurfaceElementSubRegion & subRegion )
    {
      string const & contactRelationName = subRegion.template getReference< string >( viewKeyStruct::contactRelationNameString() );
      std::cout << contactRelationName << std::endl;
      std::cout << subRegion.size() << std::endl;

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] ( localIndex const kfe )
      {
        std::cout << kfe << std::endl;
      });
    });
  });
  std::cout << "End initializePostInitialConditionsPreSubGroups" << std::endl;
  */
  /*
  SolidMechanicsLagrangianFEM::initializePostInitialConditionsPreSubGroups();
  this->updateState( this->getGroupByPath< DomainPartition >( "/Problem/domain" ) );
  */
//}

void SolidMechanicsAugmentedLagrangianContact::setupDofs( DomainPartition const & domain,
                                                          DofManager & dofManager ) const
{
  
  std::cout << "setupDofs" << std::endl;
  GEOS_UNUSED_VAR( domain, dofManager );
  
  GEOS_MARK_FUNCTION;
  SolidMechanicsLagrangianFEM::setupDofs( domain, dofManager );

  map< std::pair< string, string >, array1d< string > > meshTargets;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                                MeshLevel const & meshLevel,
                                                                arrayView1d< string const > const & )
  {
    array1d< string > regions;
    regions.emplace_back( getUniqueFractureRegionName() );
    meshTargets[std::make_pair( meshBodyName, meshLevel.getName())] = std::move( regions );
  } );

  dofManager.addField( solidMechanics::totalBubbleDisplacement::key(),
                       FieldLocation::Face,
                       3,
                       meshTargets );
                       //getMeshTargets());

  dofManager.addCoupling( solidMechanics::totalBubbleDisplacement::key(),
                          solidMechanics::totalBubbleDisplacement::key(),
                          DofManager::Connector::Elem);

  //dofManager.addCoupling( solidMechanics::totalDisplacement::key(),
  //                        solidMechanics::totalBubbleDisplacement::key(),
  //                        DofManager::Connector::Elem);

}

void SolidMechanicsAugmentedLagrangianContact::setupSystem( DomainPartition & domain,
                                                            DofManager & dofManager,
                                                            CRSMatrix< real64, globalIndex > & localMatrix,
                                                            ParallelVector & rhs,
                                                            ParallelVector & solution,
                                                            bool const setSparsity )
{

  //GEOS_UNUSED_VAR( dofManager, localMatrix, rhs, solution, setSparsity );

  std::cout << "setupSystem" << std::endl;
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( setSparsity );
  //SolidMechanicsLagrangianFEM::setupSystem( domain, dofManager, localMatrix, rhs, solution, true );
  //SolverBase::setupSystem( domain, dofManager, localMatrix, rhs, solution, true ); // "true" is to force setSparsity
  
  //this->m_faceTypesToFaceElements["Quadrilateral"] = quadList;
  //this->m_faceTypesToFaceElements["Quadrilateral"].toView()
  //arrayView1d< localIndex const > const quadList = this->m_faceTypesToFaceElements["Quadrilateral"].toView();
  //std::cout << quadList[0] << std::endl;
  //std::cout << quadList[0] << " " << quadList[1] << " " << quadList[2] << std::endl;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const regionNames )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();
    ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

    SurfaceElementRegion const & region = elemManager.getRegion< SurfaceElementRegion >( getUniqueFractureRegionName() );
    FaceElementSubRegion const & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();


    array1d< localIndex > tmpSpace(2*subRegion.size());
    SortedArray< localIndex > faceIdList;
    {
    ArrayOfArraysView< localIndex const > const elemsToFaces = subRegion.faceList().toViewConst();

    array1d< localIndex > keys(subRegion.size());
    array1d< localIndex > vals(subRegion.size());
    array1d< localIndex > quadList;
    array1d< localIndex > triList;
    RAJA::ReduceSum< ReducePolicy< parallelDevicePolicy<> >, localIndex > nTri_r( 0 );
    RAJA::ReduceSum< ReducePolicy< parallelDevicePolicy<> >, localIndex > nQuad_r( 0 );

    forAll< parallelDevicePolicy<> > ( subRegion.size(), [&] GEOS_HOST_DEVICE ( localIndex const kfe )
    {

      localIndex const kf0 = elemsToFaces[kfe][0], kf1 = elemsToFaces[kfe][1];
      tmpSpace[2*kfe] = kf0, tmpSpace[2*kfe+1] = kf1;

      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kfe );
      if (numNodesPerFace == 3)
      {
        keys[kfe]=0;
        vals[kfe]=kfe;
        nTri_r+=1;
      }
      else if (numNodesPerFace == 4) 
      {
        keys[kfe]=1;
        vals[kfe]=kfe;
        nQuad_r+=1;
      }
      else 
      {
        GEOS_ERROR( "SolidMechanicsAugmentedLagrangianContact:: invalid face type" );
      }
    });
    localIndex nQuad = static_cast<localIndex>(nQuad_r.get());
    localIndex nTri = static_cast<localIndex>(nTri_r.get());
    RAJA::sort_pairs< parallelDevicePolicy<> >(keys, vals);

    std::cout << nQuad << " " << nTri << std::endl;
    quadList.resize(nQuad);
    triList.resize(nTri);

    forAll< parallelDevicePolicy<> > ( nTri, [&] GEOS_HOST_DEVICE ( localIndex const kfe )
    {
      quadList[kfe] = vals[nTri+kfe];
    });

    forAll< parallelDevicePolicy<> > ( nQuad, [&] GEOS_HOST_DEVICE ( localIndex const kfe )
    {
      quadList[kfe] = vals[nTri+kfe];
    });

    this->m_faceTypesToFaceElements[meshName]["Quadrilateral"] =  quadList;
    this->m_faceTypesToFaceElements[meshName]["Triangle"] =  triList;
    }

    RAJA::stable_sort< parallelDevicePolicy<> >(tmpSpace);
    faceIdList.insert( tmpSpace.begin(), tmpSpace.end());

    //elemManager.forElementSubRegionsComplete< CellElementSubRegion >( [&]( localIndex const, 
    //                                                                       localIndex const, 
    //                                                                       ElementRegionBase &, 
    //                                                                       CellElementSubRegion & subRegion1 )

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const, CellElementSubRegion & cellElementSubRegion )
    {

      arrayView2d< localIndex const  > const elemsToFaces = cellElementSubRegion.faceList().toViewConst();
      std::cout << "ElemToFaces Size: " << elemsToFaces.size(0) << " " << elemsToFaces.size(1) << std::endl;

      RAJA::ReduceSum< ReducePolicy< parallelDevicePolicy<> >, localIndex > nBubElems_r( 0 );

      localIndex const n_max = cellElementSubRegion.size() * elemsToFaces.size(1);
      array1d< localIndex > keys(n_max);
      array1d< localIndex > perms(n_max);
      array1d< localIndex > vals(n_max);
      array1d< localIndex > localFaceIds(n_max);

      forAll< parallelDevicePolicy<> >( cellElementSubRegion.size(), [&] GEOS_HOST_DEVICE ( localIndex const kfe )
      {
        //std::cout << "Elem: " <<  kfe << std::endl;
        for (int i=0; i < elemsToFaces.size(1); ++i) 
        {
          //std::cout << "      " <<  elemsToFaces[kfe][i] << " ";
          perms[kfe*elemsToFaces.size(1)+i] = kfe*elemsToFaces.size(1)+i;
          if (faceIdList.contains(elemsToFaces[kfe][i]))
          {
            keys[kfe*elemsToFaces.size(1)+i] = 0;
            vals[kfe*elemsToFaces.size(1)+i] = kfe;
            localFaceIds[kfe*elemsToFaces.size(1)+i] = i;
            nBubElems_r += 1;
            //std::cout << "elem - faceId - locFaceId: " << kfe << " " << elemsToFaces[kfe][i] << " " << i << std::endl;
          }
          else 
          {
            keys[kfe*elemsToFaces.size(1)+i] = 1;
            vals[kfe*elemsToFaces.size(1)+i] = -1;
            localFaceIds[kfe*elemsToFaces.size(1)+i] = -1;
          }
        }
        //std::cout << std::endl;
      });

      localIndex nBubElems = static_cast<localIndex>(nBubElems_r.get());
      RAJA::sort_pairs< parallelDevicePolicy<> >(keys, perms);
      std::cout << "# bubbles: " << nBubElems << std::endl;
      //for (int i=0; i<keys.size(); ++i)
      //{
      //  std::cout << keys[i] << " " << perms[i] << std::endl;
      //}
      array1d< localIndex > bubbleElemsList;
      bubbleElemsList.resize(nBubElems);

      forAll< parallelDevicePolicy<> >( n_max, [&] GEOS_HOST_DEVICE ( localIndex const k )
      {
        keys[k] = vals[perms[k]];
      });
      //RAJA::stable_sort< parallelDevicePolicy<> >(vals);
      //bubbleElemsList.insert( keys.begin(), keys.begin() + nBubElems);
      forAll< parallelDevicePolicy<> >( nBubElems, [&] GEOS_HOST_DEVICE ( localIndex const k )
      {
        bubbleElemsList[k] = keys[k];
      });
      cellElementSubRegion.setBubbleElementsList(bubbleElemsList.toViewConst());

      forAll< parallelDevicePolicy<> >( n_max, [&] GEOS_HOST_DEVICE ( localIndex const k )
      {
        keys[k] = localFaceIds[perms[k]];
      });

      array2d< localIndex > faceElemsList;
      faceElemsList.resize(nBubElems,2);
      forAll< parallelDevicePolicy<> >( nBubElems, [&] GEOS_HOST_DEVICE ( localIndex const k )
      {
        localIndex const kfe =  bubbleElemsList[k];
        //if (faceIdList.contains(elemsToFaces[kfe][i]))
        //{
        faceElemsList[k][0] = elemsToFaces[kfe][keys[k]];
        faceElemsList[k][1] = keys[k];
        //}
      });
      cellElementSubRegion.setFaceElementsList(faceElemsList.toViewConst());

    });

  });

  //SolverBase::setupSystem( domain, dofManager, localMatrix, rhs, solution, true ); // "true" is to force setSparsity

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
  this->addCouplingNumNonzeros( domain, dofManager, rowLengths.toView() );

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
  this->addCouplingSparsityPattern( domain, dofManager, pattern.toView() );

  // Finally, steal the pattern into a CRS matrix
  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );
  localMatrix.setName( this->getName() + "/localMatrix" );

  rhs.setName( this->getName() + "/rhs" );
  rhs.create( dofManager.numLocalDofs(), MPI_COMM_GEOSX );

  solution.setName( this->getName() + "/solution" );
  solution.create( dofManager.numLocalDofs(), MPI_COMM_GEOSX );

}

void SolidMechanicsAugmentedLagrangianContact::addCouplingNumNonzeros( DomainPartition & domain,
                                                                       DofManager & dofManager,
                                                                       arrayView1d< localIndex > const & rowLengths ) const
{

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();
    NodeManager const &          nodeManager = mesh.getNodeManager();
    FaceManager const &          faceManager = mesh.getFaceManager();

    ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

    globalIndex const rankOffset = dofManager.rankOffset();

    string const bubbleDofKey = dofManager.getKey( solidMechanics::totalBubbleDisplacement::key() );
    string const dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );

    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );
    arrayView1d< globalIndex const > const dispDofNumber =  nodeManager.getReference< globalIndex_array >( dispDofKey );

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const, CellElementSubRegion const & cellElementSubRegion )
    {

      arrayView1d< localIndex const > const bubbleElemsList = cellElementSubRegion.bubbleElementsList();
      arrayView2d< localIndex const > const faceElemsList = cellElementSubRegion.faceElementsList();

      localIndex const numDispDof = 3*cellElementSubRegion.numNodesPerElement();

      //for( localIndex bi=0; bi<bubbleElemsList.size(); ++bi )
      forAll< parallelDevicePolicy<> >( bubbleElemsList.size(), [&] GEOS_HOST_DEVICE ( localIndex const bi )
      {
        localIndex const cellIndex = bubbleElemsList[bi];
        localIndex const k = faceElemsList[bi][0];

        localIndex const localRow = LvArray::integerConversion< localIndex >( bubbleDofNumber[k] - rankOffset );

        if( localRow >= 0 && localRow < rowLengths.size() )
        {
          for( localIndex i=0; i<3; ++i )
          {
            //rowLengths[localRow + i] += numDispDof;
            RAJA::atomicAdd<parallelDeviceAtomic>(&rowLengths[localRow + i], numDispDof);
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
              //rowLengths[localDispRow + d] += 3;
              RAJA::atomicAdd<parallelDeviceAtomic>( &rowLengths[localDispRow + d], 3);
            }
          }
        }
      });

    });

    SurfaceElementRegion const & region = elemManager.getRegion< SurfaceElementRegion >( getUniqueFractureRegionName() );
    FaceElementSubRegion const & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();
    ArrayOfArraysView< localIndex const > const  elemsToFaces = subRegion.faceList().toViewConst();

    forAll< parallelDevicePolicy<> > ( subRegion.size(), [&] GEOS_HOST_DEVICE ( localIndex const kfe )
    {
      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kfe );
      localIndex const numDispDof = 3*numNodesPerFace;

      for (int k=0; k<2; ++k)
      {
        localIndex const kf = elemsToFaces[kfe][k];
     
        localIndex const localRow = LvArray::integerConversion< localIndex >( bubbleDofNumber[kf] - rankOffset );
     
        if( localRow >= 0 && localRow < rowLengths.size() )
        {
          for( localIndex i=0; i<3; ++i )
          {
            //rowLengths[localRow + i] += numDispDof;
            RAJA::atomicAdd<parallelDeviceAtomic>(&rowLengths[localRow + i], numDispDof);
          }
        }

        for( localIndex a=0; a<numNodesPerFace; ++a )
        {
          const localIndex & node = faceToNodeMap( kf, a );
          localIndex const localDispRow = LvArray::integerConversion< localIndex >( dispDofNumber[node] - rankOffset );

          if( localDispRow >= 0 && localDispRow < rowLengths.size() )
          {
            for( int d=0; d<3; ++d )
            {
              //rowLengths[localDispRow + d] += 3;
              RAJA::atomicAdd<parallelDeviceAtomic>( &rowLengths[localDispRow + d], 3);
            }
          }
        }
      }

    });

  });
}

void SolidMechanicsAugmentedLagrangianContact::addCouplingSparsityPattern( DomainPartition const & domain,
                                                                           DofManager const & dofManager,
                                                                           SparsityPatternView< globalIndex > const & pattern ) const
{

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();
    NodeManager const &          nodeManager = mesh.getNodeManager();
    FaceManager const &          faceManager = mesh.getFaceManager();

    globalIndex const rankOffset = dofManager.rankOffset();

    string const bubbleDofKey = dofManager.getKey( solidMechanics::totalBubbleDisplacement::key() );
    string const dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );

    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );
    arrayView1d< globalIndex const > const dispDofNumber =  nodeManager.getReference< globalIndex_array >( dispDofKey );

    static constexpr int maxNumDispDof = 3 * 8;

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const, CellElementSubRegion const & cellElementSubRegion )
    {

      arrayView1d< localIndex const > const bubbleElemsList = cellElementSubRegion.bubbleElementsList();
      arrayView2d< localIndex const > const faceElemsList = cellElementSubRegion.faceElementsList();

      localIndex const numDispDof = 3*cellElementSubRegion.numNodesPerElement();

      // Can't I update the pattern atomically?
      //forAll< parallelDevicePolicy<> >( bubbleElemsList.size(), [&] GEOS_HOST_DEVICE ( localIndex const bi )
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

    });

    SurfaceElementRegion const & region = elemManager.getRegion< SurfaceElementRegion >( getUniqueFractureRegionName() );
    FaceElementSubRegion const & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();
    ArrayOfArraysView< localIndex const > const  elemsToFaces = subRegion.faceList().toViewConst();
    ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

    static constexpr int maxNumDispFaceDof = 3 * 4;
    // Can't I update the pattern atomically?
    //forAll< parallelDevicePolicy<> > ( subRegion.size(), [&] GEOS_HOST_DEVICE ( localIndex const kfe )
    for( localIndex kfe=0; kfe<subRegion.size(); ++kfe )
    {

      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kfe );
      localIndex const numDispDof = 3*numNodesPerFace;

      for (int k=0; k<2; ++k)
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
  });

}

void SolidMechanicsAugmentedLagrangianContact::implicitStepSetup( real64 const & time_n,
                                                                  real64 const & dt,
                                                                  DomainPartition & domain )
{
  /*
  computeRotationMatrices( domain );
  computeTolerances( domain );
  computeFaceDisplacementJump( domain );

  */
  SolidMechanicsLagrangianFEM::implicitStepSetup( time_n, dt, domain );

  std::cout << "implicitStepSetup" << std::endl;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {

    FaceManager & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    SurfaceElementRegion & region = elemManager.getRegion< SurfaceElementRegion >( getUniqueFractureRegionName() );
    FaceElementSubRegion & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();

    arrayView2d< real64 const > const & faceNormal = faceManager.faceNormal();
    ArrayOfArraysView< localIndex const > const  elemsToFaces = subRegion.faceList().toViewConst();

    //std::cout << "Before Rotation" << std::endl;
    arrayView3d< real64 > const 
      rotationMatrix = subRegion.getField< fields::contact::rotationMatrix >().toView();
      //rotationMatrix = subRegion.getReference< array3d< real64 > >( viewKeyStruct::rotationMatrixString() );
      

    //std::cout << "Before Map" << std::endl;
    //std::map< string, 
    //          array1d< localIndex > > const & faceTypesToFaceElements = m_faceTypesToFaceElements.at(meshName); 


    //for (const auto& [finiteElementName, faceElementList] : faceTypesToFaceElements)
    //for (auto it = faceTypesToFaceElements.begin(); it != faceTypesToFaceElements.end(); ++it) 
    //{

    solidMechanicsALMKernels::ComputeRotationMatricesKernel::
      launch< parallelDevicePolicy<> >( subRegion.size(),
                                        faceNormal,
                                        elemsToFaces,
                                        rotationMatrix );
    //}

    arrayView2d< real64 > const 
      penalty = subRegion.getField< fields::contact::penalty >().toView();

    forAll< parallelDevicePolicy<> >( subRegion.size(), [&] GEOS_HOST_DEVICE ( localIndex const k )
    {
      penalty[k] [0] = 1.e+9;
      penalty[k] [1] = 1.e+9;
    });

  } );
}

void SolidMechanicsAugmentedLagrangianContact::assembleSystem( real64 const time,
                                                               real64 const dt,
                                                               DomainPartition & domain,
                                                               DofManager const & dofManager,
                                                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                               arrayView1d< real64 > const & localRhs )
{

  //GEOS_UNUSED_VAR( time, dt, domain, dofManager, localMatrix, localRhs );
  GEOS_UNUSED_VAR( time);
  synchronizeFractureState( domain );
  std::cout << "assembleSystem" << std::endl;
  //ParallelMatrix parallel_matrix;
  //parallel_matrix.create( localMatrix.toViewConst(), dofManager.numLocalDofs(), MPI_COMM_GEOSX );
  //parallel_matrix.write("mech0.mtx");
  //abort();

  GEOS_MARK_FUNCTION;

  SolidMechanicsLagrangianFEM::assembleSystem( time,
                                               dt,
                                               domain,
                                               dofManager,
                                               localMatrix,
                                               localRhs );

  //ParallelMatrix parallel_matrix;
  //parallel_matrix.create( localMatrix.toViewConst(), dofManager.numLocalDofs(), MPI_COMM_GEOSX );
  //parallel_matrix.write("mech.mtx");
  //abort();

  // If specified as a b.c. apply traction
  //applyTractionBC( time, dt, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )

  {

    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
    string const & bubbleDofKey = dofManager.getKey( solidMechanics::totalBubbleDisplacement::key() );

    arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );

    string const & fractureRegionName = this->getUniqueFractureRegionName();

    //std::map< string, 
    //          array1d< localIndex > > const & faceTypesToFaceElements = m_faceTypesToFaceElements.at(meshName); 


    //for (const auto& [finiteElementName, faceElementList] : faceTypesToFaceElements)
    //for (auto it = faceTypesToFaceElements.begin(); it != faceTypesToFaceElements.end(); ++it) 
    forFiniteElementOnFractureSubRegions( meshName, [&] (string const & finiteElementName,
                                                        arrayView1d< localIndex const > const & faceElementList )
    {

      //string const & finiteElementName = it->first;
      //array1d< localIndex > const &  faceElementList = faceTypesToFaceElements.at(finiteElementName);  

      finiteElement::FiniteElementBase & subRegionFE = *(m_faceTypeToFiniteElements[finiteElementName]);

      //arrayView1d< localIndex const > const faceElemList = faceElementList.toViewConst();

      solidMechanicsALMKernels::ALMFactory kernelFactory( dispDofNumber,
                                                          bubbleDofNumber,
                                                          dofManager.rankOffset(),
                                                          localMatrix,
                                                          localRhs,
                                                          dt,
                                                          faceElementList );

      //GEOS_UNUSED_VAR(kernelFactory, fractureRegionName, subRegionFE);
      real64 maxTraction = finiteElement::
                             interfaceBasedKernelApplication
                           < parallelDevicePolicy< >,
                             constitutive::NullModel >( mesh,
                                                        fractureRegionName,
                                                        faceElementList,
                                                        subRegionFE,
                                                        "",
                                                        kernelFactory );

      //solidMechanicsALMKernels::ALMJumpUpdateFactory kernelFactory1( dispDofNumber,
      //                                                               dofManager.rankOffset(),
      //                                                               localMatrix,
      //                                                               localRhs,
      //                                                               dt,
      //                                                               faceElementList );

      //real64 maxTraction1 = finiteElement::
      //                       interfaceBasedKernelApplication
      //                     < parallelDevicePolicy< >,
      //                       constitutive::NullModel >( mesh,
      //                                                  fractureRegionName,
      //                                                  faceElementList,
      //                                                  subRegionFE,
      //                                                  "",
      //                                                  kernelFactory1 );

      GEOS_UNUSED_VAR( maxTraction );

    } );

    /*ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const, CellElementSubRegion const & subRegion1 )
    {
      arrayView1d< localIndex const > const bubElems = subRegion1.bubbleElementsList();
      arrayView2d< localIndex const > const faceElems = subRegion1.faceElementsList();
      std::cout << bubElems.size() << " " << faceElems.size(0) << " " << faceElems.size(1) << std::endl;
      for (int i=0; i<bubElems.size(); ++i)
      {
        std::cout << bubElems[i] << " " << faceElems[i][0] << " " << faceElems[i][1] << std::endl;
      }

    });
    abort();
    */
  });

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();
    //ElementRegionManager & elemManager = mesh.getElemManager();

    string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
    string const & bubbleDofKey = dofManager.getKey( solidMechanics::totalBubbleDisplacement::key() );

    arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );

    real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );


    solidMechanicsALMKernels::ALMBubbleFactory kernelFactory( dispDofNumber,
                                                              bubbleDofNumber,
                                                              dofManager.rankOffset(),
                                                              localMatrix,
                                                              localRhs,
                                                              dt,
                                                              gravityVectorData );

    real64 maxTraction = finiteElement::
                           regionBasedKernelApplication
                         < parallelDevicePolicy< >,
                           constitutive::ElasticIsotropic,
                           CellElementSubRegion >( mesh,
                                                  regionNames,
                                                   getDiscretizationName(),
                                                   SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString(),
                                                   kernelFactory );

    GEOS_UNUSED_VAR( maxTraction );
    //GEOS_UNUSED_VAR( regionNames );

  } );
  
  //ParallelMatrix parallel_matrix_1;
  //parallel_matrix_1.create( localMatrix.toViewConst(), dofManager.numLocalDofs(), MPI_COMM_GEOSX );
  //parallel_matrix_1.write("mech_bubble.mtx");
  //abort();
  //std::ofstream ofile;
  //ofile.open("amech.rhs");

  //for (int i=0; i<localRhs.size(); ++i)
  //{
  //  ofile << localRhs[i] << std::endl;
  //}

  /*
  GEOS_MARK_FUNCTION;

  SolidMechanicsLagrangianFEM::assembleSystem( time,
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
    SurfaceElementRegion & region = elemManager.getRegion< SurfaceElementRegion >( getUniqueFractureRegionName() );
    FaceElementSubRegion & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();

    string const dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );

    arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );

    real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

    solidMechanicsALMKernels::ALMStaticCondensationFactory kernelFactory( subRegion,
                                                                          dispDofNumber,
                                                                          dofManager.rankOffset(),
                                                                          localMatrix,
                                                                          localRhs,
                                                                          dt,
                                                                          gravityVectorData );
    real64 maxTraction = finiteElement::
                           regionBasedKernelApplication
                         < parallelDevicePolicy< >,
                           constitutive::SolidBase,
                           CellElementSubRegion >( mesh,
                                                   regionNames,
                                                   getDiscretizationName(),
                                                   SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString(),
                                                   kernelFactory );

    GEOS_UNUSED_VAR( maxTraction );
  } );
  */
  //abort();

}

void SolidMechanicsAugmentedLagrangianContact::implicitStepComplete( real64 const & time_n,
                                                                     real64 const & dt,
                                                                     DomainPartition & domain )
{

  GEOS_UNUSED_VAR( time_n, dt, domain);

  std::cout << "implicitStepComplete" << std::endl;
  SolidMechanicsLagrangianFEM::implicitStepComplete( time_n, dt, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      arrayView2d< real64 const > const & dispJump = subRegion.getField< contact::dispJump >();
      arrayView2d< real64 > const & oldDispJump = subRegion.getField< contact::oldDispJump >();
      //arrayView1d< integer const > const & fractureState = subRegion.getField< contact::fractureState >();
      //arrayView1d< integer > const & oldFractureState = subRegion.getField< contact::oldFractureState >();

      forAll< parallelHostPolicy >( subRegion.size(), [&] ( localIndex const kfe )
      {
        for( localIndex i = 0; i < 3; ++i )
        {
          oldDispJump[kfe][i] = dispJump[kfe][i];
        }
        //oldFractureState[kfe] = fractureState[kfe];
      } );
    } );

  });

  /*
  SolidMechanicsLagrangianFEM::implicitStepComplete( time_n, dt, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      arrayView2d< real64 > const & deltaTraction = subRegion.getField< contact::deltaTraction >();
      arrayView2d< real64 const > const & dispJump = subRegion.getField< contact::dispJump >();
      arrayView2d< real64 > const & oldDispJump = subRegion.getField< contact::oldDispJump >();
      arrayView1d< integer const > const & fractureState = subRegion.getField< contact::fractureState >();
      arrayView1d< integer > const & oldFractureState = subRegion.getField< contact::oldFractureState >();

      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        for( localIndex i = 0; i < 3; ++i )
        {
          deltaTraction[kfe][i] = 0.0;
          oldDispJump[kfe][i] = dispJump[kfe][i];
        }
        oldFractureState[kfe] = fractureState[kfe];
      } );
    } );

    // Need a synchronization of deltaTraction as will be used in AssembleStabilization
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( { contact::deltaTraction::key() },
                                     { getUniqueFractureRegionName() } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );

  } );
  */
}

real64 SolidMechanicsAugmentedLagrangianContact::calculateResidualNorm( real64 const & time,
                                                                        real64 const & dt,
                                                                        DomainPartition const & domain,
                                                                        DofManager const & dofManager,
                                                                        arrayView1d< real64 const > const & localRhs )
{

  GEOS_MARK_FUNCTION;

  std::cout << "calculateResidualNorm" << std::endl;

  //std::ofstream ofile;
  //ofile.open("amech_crn.rhs");

  //for (int i=0; i<localRhs.size(); ++i)
  //{
  //  ofile << localRhs[i] << std::endl;
  //}

  //std::cout << viewKeyStruct::targetNodesString() << std::endl;
  //real64 normrhs = 0.0;
  //for (int i=0; i<localRhs.size(); ++i)
  //{
  //  normrhs += localRhs[i]*localRhs[i];
  //}
  //std::cout << normrhs << std::endl;

  //real64 normrhs = 0.0;
  //for (localIndex i=0; i<localRhs.size(); i++ )
  //{
  //  std::cout << normrhs << " " << localRhs[i] << std::endl;
  //  normrhs += localRhs[i] * localRhs[i];
  //}

  //std::cout << normrhs << std::endl;

  // Matrix residual
  real64 const solidResidualNorm = SolidMechanicsLagrangianFEM::calculateResidualNorm( time, dt, domain, dofManager, localRhs );

  string const bubbleDofKey = dofManager.getKey( solidMechanics::totalBubbleDisplacement::key() );

  globalIndex const rankOffset = dofManager.rankOffset();

  RAJA::ReduceSum< parallelDeviceReduce, real64 > localSum( 0.0 );

  // globalResidualNorm[0]: the sum of all the local sum(rhs^2).
  // globalResidualNorm[1]: max of max force of each rank. Basically max force globally
  real64 globalResidualNorm[2] = {0, 0};

  // Bubble residual
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    SurfaceElementRegion const & region = elemManager.getRegion< SurfaceElementRegion >( getUniqueFractureRegionName() );
    FaceElementSubRegion const & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();

    arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();

    ArrayOfArraysView< localIndex const > const elemsToFaces = subRegion.faceList().toViewConst();

    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );

    forAll< parallelDevicePolicy<> > ( subRegion.size(), [ elemsToFaces, localRhs, localSum, bubbleDofNumber, rankOffset, ghostRank] GEOS_HOST_DEVICE ( localIndex const kfe )
    {

      for (int kk=0; kk<2; ++kk)
      {
        localIndex const k = elemsToFaces[kfe][kk];
        if( ghostRank[k] < 0 )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( bubbleDofNumber[k] - rankOffset );
          for( localIndex i = 0; i < 3; ++i )
          {
            localSum += localRhs[localRow + i] * localRhs[localRow + i];
          }
        }
      }

    });
    std::cout << "localSum: " << localSum.get() << std::endl;
    real64 const localResidualNorm[2] = { localSum.get(), SolidMechanicsLagrangianFEM::getMaxForce() };



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
  });

  real64 const bubbleResidualNorm = sqrt( globalResidualNorm[0] )/(globalResidualNorm[1]+1);  // the + 1 is for the first
    // time-step when maxForce = 0;

  if( getLogLevel() >= 1 && logger::internal::rank==0 )
  {
    std::cout << GEOS_FMT( "        ( RBubbleDisp ) = ( {:4.2e} )", bubbleResidualNorm );
  }

    return sqrt( solidResidualNorm * solidResidualNorm + bubbleResidualNorm * bubbleResidualNorm );
  //abort();

}

void SolidMechanicsAugmentedLagrangianContact::applySystemSolution( DofManager const & dofManager,
                                                                  arrayView1d< real64 const > const & localSolution,
                                                                  real64 const scalingFactor,
                                                                  real64 const dt,
                                                                  DomainPartition & domain )
{

  std::cout << "applySystemSolution" << std::endl;
  //abort();

  GEOS_MARK_FUNCTION;

  SolidMechanicsLagrangianFEM::applySystemSolution( dofManager,
                                                    localSolution,
                                                    scalingFactor,
                                                    dt,
                                                    domain );

  dofManager.addVectorToField( localSolution,
                               solidMechanics::totalBubbleDisplacement::key(),
                               solidMechanics::totalBubbleDisplacement::key(),
                               scalingFactor );

  //dofManager.addVectorToField( localSolution,
  //                             solidMechanics::totalBubbleDisplacement::key(),
  //                             solidMechanics::incrementalBubbleDisplacement::key(),
  //                             scalingFactor );

 // dofManager.addVectorToField( localSolution, contact::dispJump::key(), contact::deltaDispJump::key(), scalingFactor );


  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )

  {

    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
    string const & bubbleDofKey = dofManager.getKey( solidMechanics::totalBubbleDisplacement::key() );

    arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );

    string const & fractureRegionName = this->getUniqueFractureRegionName();

    CRSMatrix< real64, globalIndex > const voidMatrix;
    array1d< real64 > const voidRhs;

    forFiniteElementOnFractureSubRegions( meshName, [&] (string const & finiteElementName,
                                                        arrayView1d< localIndex const > const & faceElementList )
    {

      finiteElement::FiniteElementBase & subRegionFE = *(m_faceTypeToFiniteElements[finiteElementName]);

      solidMechanicsALMKernels::ALMJumpUpdateFactory kernelFactory( dispDofNumber,
                                                                    bubbleDofNumber,
                                                                    dofManager.rankOffset(),
                                                                    voidMatrix.toViewConstSizes(),
                                                                    voidRhs.toView(),
                                                                    dt,
                                                                    faceElementList );

      real64 maxTraction = finiteElement::
                             interfaceBasedKernelApplication
                           < parallelDevicePolicy< >,
                             constitutive::NullModel >( mesh,
                                                        fractureRegionName,
                                                        faceElementList,
                                                        subRegionFE,
                                                        "",
                                                        kernelFactory );

    GEOS_UNUSED_VAR( maxTraction );

    } );
  });


/*
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    FieldIdentifiers fieldsToBeSync;

    fieldsToBeSync.addElementFields( { contact::dispJump::key(),
                                       contact::deltaDispJump::key() },
                                     { getUniqueFractureRegionName() } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );
  */

/*
  // EFEM /////////////////////////
  updateJump( dofManager, dt, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    FieldIdentifiers fieldsToBeSync;

    fieldsToBeSync.addElementFields( { contact::dispJump::key(),
                                       contact::deltaDispJump::key() },
                                     { getUniqueFractureRegionName() } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );
  // ContactMechanics /////////////////////////
  dofManager.addVectorToField( localSolution,
                               contact::traction::key(),
                               contact::deltaTraction::key(),
                               scalingFactor );

  dofManager.addVectorToField( localSolution,
                               contact::traction::key(),
                               contact::traction::key(),
                               scalingFactor );

  // fractureStateString is synchronized in UpdateFractureState
  // oldFractureStateString and oldDispJumpString used locally only

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    FieldIdentifiers fieldsToBeSync;

    fieldsToBeSync.addElementFields( { contact::traction::key(),
                                       contact::deltaTraction::key(),
                                       contact::dispJump::key() },
                                     { getUniqueFractureRegionName() } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );
  */
}

//void SolidMechanicsEmbeddedFractures::resetStateToBeginningOfStep( DomainPartition & domain )
//{
  /*
  SolidMechanicsLagrangianFEM::resetStateToBeginningOfStep( domain );

  // reset displacementJump
  forFractureRegionOnMeshTargets( domain.getMeshBodies(), [&] ( SurfaceElementRegion & fractureRegion )
  {
    fractureRegion.forElementSubRegions< EmbeddedSurfaceSubRegion >( [&]( EmbeddedSurfaceSubRegion & subRegion )
    {
      arrayView2d< real64 > const & jump  =
        subRegion.getField< contact::dispJump >();

      arrayView2d< real64 const > const & oldJump  =
        subRegion.getField< contact::oldDispJump >();

      arrayView2d< real64 > const & deltaJump  =
        subRegion.getField< contact::deltaDispJump >();


      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const kfe )
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
  */
//}


// Contact //////////////////////////////////
void SolidMechanicsAugmentedLagrangianContact::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( domain );

  std::cout << "updateState" << std::endl;
//  /*
//  GEOS_MARK_FUNCTION;
//
//  computeFaceDisplacementJump( domain );
//  */
}

bool SolidMechanicsAugmentedLagrangianContact::updateConfiguration( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( domain );

  std::cout << "updateConfiguration" << std::endl;

  real64 tol[4];
  tol[0] = 1e-5;
  tol[1] = 1e-5;
  tol[2] = 5.e-2; 
  tol[3] = 1e0;

  int hasConfigurationConverged = true;
  int condConv = true;
  array2d< real64 > traction_new;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                FaceElementSubRegion & subRegion )
    {

      //arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView2d< real64 const > const  traction = subRegion.getField< contact::traction >();
      arrayView2d< real64 const > const  dispJump = subRegion.getField< contact::dispJump >();
      //Only to test quickly
      //arrayView2d< real64 const > const  deltaDispJump = subRegion.getField< contact::deltaDispJump >();
      arrayView2d< real64 const > const  deltaDispJump = subRegion.getField< contact::dispJump >();
      arrayView2d< real64 const > const  penalty = subRegion.getField< contact::penalty >();

      std::cout << traction.size() << " " << subRegion.size() << std::endl;
      std::ptrdiff_t const sizes[ 2 ] = {subRegion.size(), 3};
      traction_new.resize(2, sizes );

      forAll< parallelDevicePolicy<> >( subRegion.size(), [&] ( localIndex const kfe )
      {
        real64 eps_N = penalty[kfe][0];
        real64 eps_T = penalty[kfe][1];
        traction_new[kfe][0] = traction[kfe][0] + eps_N * dispJump[kfe][0];
        traction_new[kfe][1] = traction[kfe][1] + eps_T * deltaDispJump[kfe][1];
        traction_new[kfe][2] = traction[kfe][2] + eps_T * deltaDispJump[kfe][2];
      });
      std::cout << "Computed traction_new" << std::endl;


      forAll< parallelDevicePolicy<> >( subRegion.size(), [&] ( localIndex const kfe )
      {
        //if( ghostRank[kfe] < 0 )
        //{
          if( traction_new[kfe][0] >= tol[3] )
          {
            std::cout << "Traction_N: " << traction_new[kfe][0] << std::endl;
            std::cout << "Disp:" << std::endl;
            std::cout << kfe << " " << dispJump[kfe][0] << " " 
                      << deltaDispJump[kfe][0] << " " << deltaDispJump[kfe][1] << std::endl;
            std::cout << "Open" << std::endl;
            abort();
          }
          else 
          {
            real64 deltaDisp = sqrt(pow(deltaDispJump[kfe][1],2) + pow(deltaDispJump[kfe][2],2));
            std::cout << kfe << " " << std::abs(dispJump[kfe][0]) << " " << deltaDisp << std::endl;
            if( std::abs(dispJump[kfe][0]) > tol[0] )
            {
              std::cout << "Stick and g > tol1 => compenetration" << std::endl;
              condConv = false;
            }
            if( deltaDisp > tol[1] )
            {
              std::cout << "Stick and dg > tol2" << std::endl;
              condConv = false;
            }
          }
        //}
      });
      
    });
  });

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
 
    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                FaceElementSubRegion & subRegion )
    {
      
      arrayView2d< real64 > const  traction = subRegion.getField< contact::traction >();
      forAll< parallelDevicePolicy<> >( subRegion.size(), [&] ( localIndex const kfe )
      {
         traction[kfe][0] = traction_new(kfe, 0);
         traction[kfe][1] = traction_new(kfe, 1);
         traction[kfe][2] = traction_new(kfe, 2);
      });
    });
  });

  if (!condConv)
  {
    hasConfigurationConverged = false;
  }

  /*forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                FaceElementSubRegion & subRegion )
    {
      string const & contactRelationName = subRegion.template getReference< string >( viewKeyStruct::contactRelationNameString() );
      ContactBase const & contact = getConstitutiveModel< ContactBase >( subRegion, contactRelationName );

      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView2d< real64 const > const & traction = subRegion.getField< contact::traction >();
      arrayView2d< real64 const > const & dispJump = subRegion.getField< contact::dispJump >();
      arrayView1d< integer > const & fractureState = subRegion.getField< contact::fractureState >();

      arrayView1d< real64 const > const & normalTractionTolerance =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalTractionToleranceString() );
      arrayView1d< real64 const > const & normalDisplacementTolerance =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalDisplacementToleranceString() );

      RAJA::ReduceMin< parallelHostReduce, integer > checkActiveSetSub( 1 );

      constitutiveUpdatePassThru( contact, [&] ( auto & castedContact )
      {
        using ContactType = TYPEOFREF( castedContact );
        typename ContactType::KernelWrapper contactWrapper = castedContact.createKernelWrapper();

        forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
        {
          if( ghostRank[kfe] < 0 )
          {
            integer const originalFractureState = fractureState[kfe];
            if( originalFractureState == contact::FractureState::Open )
            {
              if( dispJump[kfe][0] > -normalDisplacementTolerance[kfe] )
              {
                fractureState[kfe] = contact::FractureState::Open;
              }
              else
              {
                fractureState[kfe] = contact::FractureState::Stick;
              }
            }
            else if( traction[kfe][0] > normalTractionTolerance[kfe] )
            {
              fractureState[kfe] = contact::FractureState::Open;
            }
            else
            {
              real64 currentTau = sqrt( traction[kfe][1]*traction[kfe][1] + traction[kfe][2]*traction[kfe][2] );

              real64 dLimitTangentialTractionNorm_dTraction = 0.0;
              real64 const limitTau =
                contactWrapper.computeLimitTangentialTractionNorm( traction[kfe][0],
                                                                   dLimitTangentialTractionNorm_dTraction );

              if( originalFractureState == contact::FractureState::Stick && currentTau >= limitTau )
              {
                currentTau *= (1.0 - m_slidingCheckTolerance);
              }
              else if( originalFractureState != contact::FractureState::Stick && currentTau <= limitTau )
              {
                currentTau *= (1.0 + m_slidingCheckTolerance);
              }
              if( currentTau > limitTau )
              {
                if( originalFractureState == contact::FractureState::Stick )
                {
                  fractureState[kfe] = contact::FractureState::NewSlip;
                }
                else
                {
                  fractureState[kfe] = contact::FractureState::Slip;
                }
              }
              else
              {
                fractureState[kfe] = contact::FractureState::Stick;
              }
            }
            checkActiveSetSub.min( compareFractureStates( originalFractureState, fractureState[kfe] ) );
          }
        } );
      } );

      hasConfigurationConverged &= checkActiveSetSub.get();
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

  return hasConfigurationConvergedGlobally;
  */
  return hasConfigurationConverged;
}

// EFEM //////////////////////////////////
//void SolidMechanicsEmbeddedFractures::updateState( DomainPartition & domain )
//{
//  /*
//  GEOS_MARK_FUNCTION;
//
//  forFractureRegionOnMeshTargets( domain.getMeshBodies(), [&] ( SurfaceElementRegion & fractureRegion )
//  {
//    fractureRegion.forElementSubRegions< SurfaceElementSubRegion >( [&]( SurfaceElementSubRegion & subRegion )
//    {
//      string const & contactRelationName = subRegion.template getReference< string >( viewKeyStruct::contactRelationNameString() );
//      ContactBase const & contact = getConstitutiveModel< ContactBase >( subRegion, contactRelationName );
//
//      arrayView2d< real64 const > const & jump = subRegion.getField< contact::dispJump >();
//
//      arrayView2d< real64 const > const & oldJump = subRegion.getField< contact::oldDispJump >();
//
//      arrayView2d< real64 > const & fractureTraction = subRegion.getField< contact::traction >();
//
//      arrayView3d< real64 > const & dFractureTraction_dJump = subRegion.getField< contact::dTraction_dJump >();
//
//      arrayView1d< integer const > const & fractureState = subRegion.getField< contact::fractureState >();
//
//      constitutiveUpdatePassThru( contact, [&] ( auto & castedContact )
//      {
//        using ContactType = TYPEOFREF( castedContact );
//        typename ContactType::KernelWrapper contactWrapper = castedContact.createKernelWrapper();
//
//        solidMechanicsEFEMKernels::StateUpdateKernel::
//          launch< parallelDevicePolicy<> >( subRegion.size(),
//                                            contactWrapper,
//                                            oldJump,
//                                            jump,
//                                            fractureTraction,
//                                            dFractureTraction_dJump,
//                                            fractureState );
//      } );
//    } );
//  } );
//  */
//}

REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsAugmentedLagrangianContact, string const &, Group * const )
} /* namespace geos */
