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

      // Register the delta traction field
      subRegion.registerField< deltaTraction >( this->getName() ).
        reference().resizeDimension< 1 >( 3 );

    } );
  } );
  
}

void SolidMechanicsAugmentedLagrangianContact::initializePostInitialConditionsPreSubGroups()
{
  std::cout << "initializePostInitialConditionsPreSubGroups" << std::endl;

  //array1d< localIndex > quadList;
  //array1d< localIndex > triList;
  //quadList.resize(10);
  //quadList[0] = 10;
  //quadList[1] = 9;
  //std::cout << quadList[0] << std::endl;
  //this->m_faceTypesToFaceElements["Quadrilateral"] =  quadList;
  //this->m_faceTypesToFaceElements["Triangle"] =  quadList;

  SolidMechanicsLagrangianFEM::initializePostInitialConditionsPreSubGroups();
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
  /*
  SolidMechanicsLagrangianFEM::initializePostInitialConditionsPreSubGroups();
  this->updateState( this->getGroupByPath< DomainPartition >( "/Problem/domain" ) );
  */
}

void SolidMechanicsAugmentedLagrangianContact::setupDofs( DomainPartition const & domain,
                                                          DofManager & dofManager ) const
{
  
  std::cout << "setupDofs" << std::endl;
  //GEOS_UNUSED_VAR( domain, dofManager );
  
  GEOS_MARK_FUNCTION;
  SolidMechanicsLagrangianFEM::setupDofs( domain, dofManager );
  
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
  SolidMechanicsLagrangianFEM::setupSystem( domain, dofManager, localMatrix, rhs, solution, setSparsity );
  
  //this->m_faceTypesToFaceElements["Quadrilateral"] = quadList;
  //this->m_faceTypesToFaceElements["Quadrilateral"].toView()
  //arrayView1d< localIndex const > const quadList = this->m_faceTypesToFaceElements["Quadrilateral"].toView();
  //std::cout << quadList[0] << std::endl;
  //std::cout << quadList[0] << " " << quadList[1] << " " << quadList[2] << std::endl;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    const FaceManager* const faceManager = &(mesh.getFaceManager());
    const ElementRegionManager*  const elemManager = &(mesh.getElemManager());
    ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager->nodeList().toViewConst();

    const SurfaceElementRegion* const region = &(elemManager->getRegion< SurfaceElementRegion >( getUniqueFractureRegionName() ));
    const FaceElementSubRegion* const subRegion = &(region->getUniqueSubRegion< FaceElementSubRegion >());

    array1d< localIndex > keys(subRegion->size());
    array1d< localIndex > vals(subRegion->size());
    array1d< localIndex > quadList;
    array1d< localIndex > triList;
    RAJA::ReduceSum< ReducePolicy< parallelDevicePolicy<> >, localIndex > nTri_r( 0 );
    RAJA::ReduceSum< ReducePolicy< parallelDevicePolicy<> >, localIndex > nQuad_r( 0 );

    forAll< parallelDevicePolicy<> > ( subRegion->size(), [&] GEOS_HOST_DEVICE ( localIndex const kfe )
    {
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

  });
}


void SolidMechanicsAugmentedLagrangianContact::implicitStepSetup( real64 const & time_n,
                                                                  real64 const & dt,
                                                                  DomainPartition & domain )
{
  GEOS_UNUSED_VAR( time_n, dt);

  /*
  computeRotationMatrices( domain );
  computeTolerances( domain );
  computeFaceDisplacementJump( domain );

  SolidMechanicsLagrangianFEM::implicitStepSetup( time_n, dt, domain );
  */

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

  std::cout << "assembleSystem" << std::endl;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )

  {

    NodeManager const & nodeManager = mesh.getNodeManager();
    string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
    arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );

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
                                                          dofManager.rankOffset(),
                                                          localMatrix,
                                                          localRhs,
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

  abort();
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

}

void SolidMechanicsAugmentedLagrangianContact::implicitStepComplete( real64 const & time_n,
                                                                     real64 const & dt,
                                                                     DomainPartition & domain )
{

  GEOS_UNUSED_VAR( time_n, dt, domain);
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

  // Matrix residual
  real64 const solidResidualNorm = SolidMechanicsLagrangianFEM::calculateResidualNorm( time, dt, domain, dofManager, localRhs );
  

  return solidResidualNorm;
}

void SolidMechanicsAugmentedLagrangianContact::applySystemSolution( DofManager const & dofManager,
                                                                  arrayView1d< real64 const > const & localSolution,
                                                                  real64 const scalingFactor,
                                                                  real64 const dt,
                                                                  DomainPartition & domain )
{

  GEOS_UNUSED_VAR( dofManager, localSolution, scalingFactor, dt, domain);
  /*
  GEOS_MARK_FUNCTION;

  SolidMechanicsLagrangianFEM::applySystemSolution( dofManager,
                                                    localSolution,
                                                    scalingFactor,
                                                    dt,
                                                    domain );

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
//void SolidMechanicsLagrangeContact::updateState( DomainPartition & domain )
//{
//  /*
//  GEOS_MARK_FUNCTION;
//
//  computeFaceDisplacementJump( domain );
//  */
//}

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