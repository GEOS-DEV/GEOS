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
 * @file MultiResolutionHFSolver.cpp
 *
 */


#include "MultiResolutionHFSolver.hpp"

#include "common/GEOS_RAJA_Interface.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/ContactSelector.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "constitutive/fluid/SingleFluidFields.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "finiteElement/Kinematics.h"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MapMeshLevels.hpp"
#include "mesh/MeshManager.hpp"
#include "mesh/generators/InternalMeshGenerator.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"
#include "physicsSolvers/contact/SolidMechanicsEFEMKernelsHelper.hpp"
#include "physicsSolvers/contact/SolidMechanicsEmbeddedFractures.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/surfaceGeneration/SurfaceGenerator.hpp"
#include "physicsSolvers/surfaceGeneration/EmbeddedSurfaceGenerator.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include <unistd.h>



namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

GEOSX_HOST_DEVICE inline
void coarseToFineStructuredElemMap(localIndex const coarseElemIndex,
                                   localIndex const GEOSX_UNUSED_PARAM(coarseNx), 
                                   localIndex const coarseNy,
                                   localIndex const coarseNz,
                                   localIndex const fineRx,
                                   localIndex const fineRy,
                                   localIndex const fineRz,
                                   localIndex triplet[3])   
{
    localIndex X = coarseElemIndex/(coarseNy*coarseNz);
    localIndex Y = (coarseElemIndex - X*(coarseNy*coarseNz)) / coarseNz;
    localIndex Z = (coarseElemIndex - X*(coarseNy*coarseNz) - Y * coarseNz);
    //this only returns the triplet for the first element 
    //to generate all, the equation is the following
    //0<= u < fineRx
    //0<= v < fineRy
    //0<= w < fineRz
    //fineIndex = (X*fineRx + u)*fineRy*coarseNy*fineRz*coarseNz 
    // + (Y*fineRy + v)*fineRz*coarseNz + (Z*fineRz + w)
    //vary u,v,w
    triplet[0] = X*fineRx;
    triplet[1] = Y*fineRy;
    triplet[2] = Z*fineRz;
    //R1Tensor fineRefElem = {X*fineRx, Y*fineRy, Z*fineRz};
} 

GEOSX_HOST_DEVICE inline
localIndex fineToCoarseStructuredElemMap(localIndex const fineElemIndex,
                                               localIndex const GEOSX_UNUSED_PARAM(coarseNx), 
                                               localIndex const coarseNy,
                                               localIndex const coarseNz,
                                               localIndex const fineRx,
                                               localIndex const fineRy,
                                               localIndex const fineRz)   
{
    localIndex coarseElemIndex = 0;
    localIndex x = fineElemIndex/(fineRy*coarseNy*fineRz*coarseNz);
    localIndex y = (fineElemIndex - x * (fineRy*coarseNy*fineRz*coarseNz)) / (fineRz*coarseNz);
    localIndex z = fineElemIndex - x * (fineRy*coarseNy*fineRz*coarseNz) - y * (fineRz*coarseNz);
    coarseElemIndex = (x/fineRx) * coarseNy * coarseNz + (y/fineRy) * coarseNz + (z/fineRz); 
    return coarseElemIndex;    
}                                               

MultiResolutionHFSolver::MultiResolutionHFSolver( const string & name,
                                                  Group * const parent ):
  SolverBase( name, parent ),
  m_baseSolverName(),
  m_patchSolverName(),
  m_nodeFixDamage(),
  m_nodeFixDisp(),
  m_fixedDispList(),
  m_maxNumResolves( 10 ),
  m_baseCrackFront(),
  m_addedFractureElements(false)
{
  registerWrapper( viewKeyStruct::baseSolverNameString(), &m_baseSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription(
    "Name of the EFEM (SolidMechanicsEmbeddedFractures) solver to be used as the base solver in the MultiResolution scheme" );

  registerWrapper( viewKeyStruct::patchSolverNameString(), &m_patchSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription(
    "Name of the Phase-Field (PhaseFieldFracture) solver to be used as patch solver in the MultiResolution scheme" );

  registerWrapper( viewKeyStruct::maxNumResolvesString(), &m_maxNumResolves ).
    setApplyDefaultValue( 10 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value to indicate how many resolves may be executed to perform surface generation after the execution of base and patch scale solvers. " );

  registerWrapper( viewKeyStruct::initialBaseTipString(), &m_baseTip ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Position of the initial crack tip in the EFEM mesh." );

  registerWrapper( viewKeyStruct::initialTipElementIndexString(), &m_baseTipElementIndex ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Index of base element ahead of the crack tip." );

}

void MultiResolutionHFSolver::registerDataOnMesh( dataRepository::Group & meshBodies )
{
  //GEOSX_UNUSED_VAR( meshBodies );
  //the patch elem manager
  // forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
  //                                                   MeshLevel & meshLevel,
  //                                                   arrayView1d< string const > const & )
  // {

    //ElementRegionManager & patchElemManager = meshLevel.getElemManager();
    ElementRegionManager & baseElemManager = meshBodies.getGroup<MeshBody>(0).getBaseDiscretization().getElemManager();

    baseElemManager.forElementSubRegions< CellElementSubRegion,
                                          FaceElementSubRegion >( [&] ( auto & elementSubRegion )
    {                                  
    //register mapped elem index baseToPatch
    //register mapped elem index patchToBase
    //register fake pressures for testing purposes
      elementSubRegion.template registerWrapper< array1d< localIndex > >( "frontIndicator" ).
        setPlotLevel( PlotLevel::LEVEL_1 ).
        setDescription( "indicator of the crack front elements" );
    });


    ElementRegionManager & patchElemManager = meshBodies.getGroup<MeshBody>(1).getBaseDiscretization().getElemManager();

    patchElemManager.forElementSubRegions< CellElementSubRegion,
                                          FaceElementSubRegion >( [&] ( auto & elementSubRegion )
    {                                  
    //register mapped elem index baseToPatch
    //register mapped elem index patchToBase
    //register fake pressures for testing purposes
      elementSubRegion.template registerWrapper< array1d< localIndex > >( "coarseToFineMap" ).
        setPlotLevel( PlotLevel::LEVEL_1 ).
        setDescription( "mapping elem indices from coarse to fine" );
      elementSubRegion.template registerWrapper< array1d< localIndex > >( "fineToCoarseMap" ).
        setPlotLevel( PlotLevel::LEVEL_1 ).
        setDescription( "mapping elem indices from fine to coarse" );
    });

  // } );
}

MultiResolutionHFSolver::~MultiResolutionHFSolver()
{
  // TODO Auto-generated destructor stub
}

void MultiResolutionHFSolver::implicitStepSetup( real64 const & time_n,
                                                 real64 const & dt,
                                                 DomainPartition & domain )
{
  m_baseSolver->implicitStepSetup( time_n, dt, domain );
  m_patchSolver->implicitStepSetup( time_n, dt, domain );
}

void MultiResolutionHFSolver::postProcessInput()
{
  SolverBase::postProcessInput();
  m_baseSolver = &this->getParent().getGroup< SolidMechanicsEmbeddedFractures >( m_baseSolverName );
  m_patchSolver = &this->getParent().getGroup< PhaseFieldFractureSolver >( m_patchSolverName );
}

void MultiResolutionHFSolver::initializePostInitialConditionsPreSubGroups()
{}

void MultiResolutionHFSolver::resetStateToBeginningOfStep( DomainPartition & domain )
{
  m_baseSolver->resetStateToBeginningOfStep( domain );
  m_patchSolver->resetStateToBeginningOfStep( domain );
}

real64 MultiResolutionHFSolver::solverStep( real64 const & time_n,
                                            real64 const & dt,
                                            int const cycleNumber,
                                            DomainPartition & domain )
{
  real64 dtReturn = dt;
  dtReturn = splitOperatorStep( time_n, dt, cycleNumber, domain );
  return dtReturn;
}

// this function will loop over all subdomain nodes and check their distance to the prescribed discrete crack. If the
// distance is smaller than 1 element size (subdomain), we set the damage in this node to be fixed at 1.
void MultiResolutionHFSolver::setInitialCrackDamageBCs( DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                                                        CRSMatrixView< real64, globalIndex const > const & GEOSX_UNUSED_PARAM( localMatrix ),
                                                        MeshLevel const & patch,
                                                        MeshLevel const & GEOSX_UNUSED_PARAM(base) )
{
  GEOSX_MARK_FUNCTION;
  //ElementRegionManager const & baseElemManager = base.getElemManager();
  ElementRegionManager const & baseElemManager = patch.getElemManager();//TODO: wrong name
  baseElemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & cellElementSubRegion )
  {
    m_nodeFixDamage.zero();
    //for this parte, it might be a good idea to use the initial bounded circles to cut the fine mesh as well, and help 
    //initialize the damage
    SortedArrayView< localIndex const > const fracturedElements = cellElementSubRegion.fracturedElementsList();
    m_nodeFixDamage.resize( cellElementSubRegion.numNodesPerElement()*fracturedElements.size() );
    localIndex count = 0;
    for( localIndex a : fracturedElements )
    {
      GEOSX_LOG_LEVEL_RANK_0( 3, "patch fractured elem: " << a << "\n" ); 
      //get all nodes of fracturedElements(a)
      for( localIndex b=0; b < cellElementSubRegion.numNodesPerElement(); b++ )
      {
        //TODO: if meshes arent identical, this should be nodeList(elem_mapped_to_patch(a),b) or check (2)
        localIndex c = cellElementSubRegion.nodeList( a, b );
        //append c = node b of element a
        if( std::find( m_nodeFixDamage.begin(), m_nodeFixDamage.end(), c ) == m_nodeFixDamage.end() ) //avoid repetition
        {
          //TODO: use node_mapped_to_patch(cellElementSubRegion.nodeList(a,b))
          m_nodeFixDamage( count ) = cellElementSubRegion.nodeList( a, b );
          ++count;
        }
      }
    }
    m_nodeFixDamage.resize( count );
  } );
}

void MultiResolutionHFSolver::findNewlyDamagedElements(SortedArray<localIndex> GEOSX_UNUSED_PARAM(toCutElems), 
                                                       MeshLevel const & GEOSX_UNUSED_PARAM(base), 
                                                       MeshLevel const & GEOSX_UNUSED_PARAM(patch))
{
  //loop over all coarse elements - meshLevel base
  // baseElementManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & cellElementSubRegion )
  // {
    //discard elements already fractured
    //loop over all fine elements associated with base element
      //test if there is a plane of fully damaged elements inside with d_elemental > thr
      //if yes, add base element to cut list
  // });
}


// this function will read the patch solution and locate the crack tip to uptade the crack geometry in the base solver
void MultiResolutionHFSolver::findPhaseFieldTip( R1Tensor & tip,
                                                 MeshLevel const & patch )
{
  GEOSX_MARK_FUNCTION;
  //reference point must be prescribed in base coordinate system (usually injection source)
  R1Tensor m_referencePoint = {-1.0, 0.0, 0.0}; //add this as input file parameter
  real64 threshold = 0.95;
  //get mpi communicator
  MPI_Comm const & comm = MPI_COMM_GEOSX;
  ElementRegionManager const & patchElemManager = patch.getElemManager();
  //this is risky, not sure FE_TYPE will come from patch
  //each rank has these arrays
  real64 rankMaxDist = 1e-20;
  R1Tensor rankFarthestCenter = {0.0, 0.0, 0.0};
  //loop over all elements in patch and compute elemental average damage
  //observe that patch coordinates are relative, so, reference point should be added
  //for elements in patch
  patchElemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & cellElementSubRegion )
  {
    localIndex numQuadraturePointsPerElem = 0;
    finiteElement::FiniteElementBase const & fe = cellElementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
    finiteElement::FiniteElementDispatchHandler< ALL_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
    {
      using FE_TYPE = TYPEOFREF( finiteElement );
      numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
    } );
    string const & damageModelName = cellElementSubRegion.getReference< string >( PhaseFieldDamageFEM::viewKeyStruct::solidModelNamesString());
    //TODO: this call may need to constitutive pass-thru loop to be generalized to multiple damage types
    const constitutive::Damage< ElasticIsotropic > & damageUpdates = cellElementSubRegion.getConstitutiveModel< Damage< ElasticIsotropic > >( damageModelName );
    arrayView2d< const real64 > allElemCenters = cellElementSubRegion.getElementCenter();
    const arrayView2d< real64 const > qp_damage = damageUpdates.getDamage();
    forAll< serialPolicy >( cellElementSubRegion.size(), [&] ( localIndex const k )
    {
      //compute elemental averaged damage
      real64 average_d = 0;
      R1Tensor elemCenter;
      elemCenter[0] = allElemCenters[k][0];
      elemCenter[1] = allElemCenters[k][1];
      elemCenter[2] = allElemCenters[k][2]; //this is trying to create a R1Tensor from a array2d< real64 >
      for( localIndex q=0; q<numQuadraturePointsPerElem; q++ )
      {
        //get damage at quadrature point i
        average_d = average_d + qp_damage( k, q )/numQuadraturePointsPerElem;
      }
      //if elemental damage > 0.95 (or another threshold)
      if( average_d > threshold )
      {
        //check if this element is farther than current farthest
        R1Tensor elemVec = LVARRAY_TENSOROPS_INIT_LOCAL_3( m_referencePoint );
        LvArray::tensorOps::subtract< 3 >( elemVec, elemCenter );
        real64 dist = LvArray::tensorOps::l2Norm< 3 >( elemVec );
        if( dist > rankMaxDist )
        {
          rankMaxDist = dist;
          rankFarthestCenter = elemCenter;
        }
      }
    } );
  } );

  real64 globalMax = MpiWrapper::max< real64 >( rankMaxDist, comm );
  if( std::abs( rankMaxDist-globalMax )<1e-12 )
  {
    tip = rankFarthestCenter;
  }

}

//The function prepares a list of dofs and u values that will be used by the patch solver to set the boundary conditions
void MultiResolutionHFSolver::prepareSubProblemBCs( MeshLevel const & base,
                                                    MeshLevel & patch )
{
  GEOSX_MARK_FUNCTION;
  m_nodeFixDisp.zero();
  m_fixedDispList.zero();
  // get list of nodes on the boundary of the patch
  FaceManager const & patchFaceManager = patch.getFaceManager();
  NodeManager & patchNodeManager = patch.getNodeManager();
  //this function finds the nodes in the boundary assuming that the domain is 2D extruded in the z direction
  patchNodeManager.setDomain2DBoundaryObjects( patchFaceManager );
  arrayView1d< integer const > patch2DBoundaryIndicator = patchNodeManager.getDomain2DBoundaryIndicator();
  arrayView1d< real64 const > const patchDamage = patchNodeManager.getReference< array1d< real64 > >( "Damage" );
  NodeManager const & baseNodeManager = base.getNodeManager();
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const baseDisp = baseNodeManager.getField< fields::solidMechanics::totalDisplacement >();
  real64 damage_threshold = 0.3;
  localIndex count=0;
  for( localIndex nodeIndex=0; nodeIndex<patchNodeManager.size(); ++nodeIndex )
  {
    if( patch2DBoundaryIndicator[nodeIndex]==1 )
    {
      if( patchDamage[nodeIndex] < damage_threshold )
      {
        // NOTE: there needs to be a translation between patch and base mesh for the indices and values/weights.

        //append node nodeIndex
        //USE REGISTERED ARRAY HERE
        m_nodeFixDisp.resize( count+1 );
        m_fixedDispList.resize( count+1, 3 );
        m_nodeFixDisp( count ) = nodeIndex;
        localIndex const numBaseNodes = 1;  //TODO: this wont be 1 in other cases, m_nodeMapIndices.sizeOfArray( a );
        for( localIndex b=0; b<numBaseNodes; ++b )
        {
          //write displacements from the base domain to become a boundary condition in the patch domain
          //TODO: a should be replaced by mapped(a) if meshes arent identical
          // mapped(a) is the node in base that has the same physical coordinates of a in patch
          m_fixedDispList( count, 0 ) = baseDisp( nodeIndex, 0 );
          m_fixedDispList( count, 1 ) = baseDisp( nodeIndex, 1 );
          m_fixedDispList( count, 2 ) = baseDisp( nodeIndex, 2 );
        }

        ++count;
      }
    }
  }

}

void MultiResolutionHFSolver::testElemMappingPatchToBase( MeshLevel & GEOSX_UNUSED_PARAM(base),
                                                          MeshLevel & patch )
{
  //for every elem in base
  //call coarseToFineStructuredElemMap
  //loop over all associated patch elements and write the base elem number
  //get patch elem manager
  ElementRegionManager & patchElemManager = patch.getElemManager();

  //get accessor to elemental field
  ElementRegionManager::ElementViewAccessor< arrayView1d< localIndex > > const fineToCoarseMap =
  patchElemManager.constructViewAccessor< array1d< localIndex >, arrayView1d< localIndex > >( "fineToCoarseMap" );
  //fineToCoarseMap[0][0][0] = fineToCoarseStructuredElemMap(0,5,5,5,3,3,3);
  patchElemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & cellElementSubRegion )
  {

    forAll< serialPolicy >( cellElementSubRegion.size(), [&] ( localIndex const k )
    {
      fineToCoarseMap[0][0][k] = fineToCoarseStructuredElemMap(k,5,5,5,3,3,3);
    } );
  } );  
}

void MultiResolutionHFSolver::testElemMappingBaseToPatch( MeshLevel & base,
                                                          MeshLevel & patch )
{
  //for every elem in patch
  //call fineToCoarseStructuredElemMap
  //write index of the associated base elem to all patch elems
  //this should be exactly the same as the function above
  ElementRegionManager & coarseElemManager = base.getElemManager();
  ElementRegionManager & patchElemManager = patch.getElemManager();

  //get accessor to elemental field
  ElementRegionManager::ElementViewAccessor< arrayView1d< localIndex > > const coarseToFineMap =
  patchElemManager.constructViewAccessor< array1d< localIndex >, arrayView1d< localIndex > >( "coarseToFineMap" );
  localIndex rx=3;
  localIndex ry=3;
  localIndex rz=3;
  //R1Tensor refElem = coarseToFineStructuredElemMap(0,5,5,5,rx,ry,rz);
  coarseElemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & cellElementSubRegion )
  {
    forAll< serialPolicy >( cellElementSubRegion.size(), [&] ( localIndex const K )
    {
      localIndex triplet[3];
      coarseToFineStructuredElemMap(K,5,5,5,rx,ry,rz,triplet);
      GEOSX_LOG_LEVEL_RANK_0( 3, "Base elem index: " << K << " triplet: ( "<<triplet[0]<<", "<<triplet[1]<<", "<<triplet[2]<<" )\n" );
      for(int i=0; i<rx; i++){
        for(int j=0; j<ry; j++){
          for(int k=0; k<rz; k++){
            localIndex fineK = (triplet[0] + i)*ry*5*rz*5 + (triplet[1] + j)*rz*5 + (triplet[2] + k);
            if(fineK > 0){
                GEOSX_LOG_LEVEL_RANK_0( 3, "fineK: " << fineK << " from base elem: "<<K<<"\n" );
            }
            coarseToFineMap[0][0][fineK] = K;
          }        
        }
      }
    } );
  } );

}


real64 MultiResolutionHFSolver::utilGetElemAverageDamage(globalIndex const patchElemNum,
                                                         MeshLevel const & patch)
{

    ElementRegionManager const & elemManager = patch.getElemManager();
    // Hard-coded region and subRegion
    ElementRegionBase const & elementRegion = elemManager.getRegion( 0 );

    CellElementSubRegion const & subRegion = elementRegion.getSubRegion< CellElementSubRegion >( 0 );

    string const & damageModelName = subRegion.getReference< string >( PhaseFieldDamageFEM::viewKeyStruct::solidModelNamesString());
    //TODO: this call may need to constitutive pass-thru loop to be generalized to multiple damage types

    const constitutive::Damage< ElasticIsotropic > & damageUpdates = subRegion.getConstitutiveModel< Damage< ElasticIsotropic > >( damageModelName );

    //arrayView2d< const real64 > allElemCenters = subRegion.getElementCenter();

    const arrayView2d< real64 const > qp_damage = damageUpdates.getDamage();

    /////////////LONG PATH TO GET numQuadPointPerElem  
    localIndex numQuadraturePointsPerElem = 0;
    finiteElement::FiniteElementBase const & fe = subRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
    finiteElement::FiniteElementDispatchHandler< ALL_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
    {
      using FE_TYPE = TYPEOFREF( finiteElement );
      numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
    } );
    ///////////////

    //compute elemental averaged damage
    real64 average_d = 0;
    for( localIndex q=0; q<numQuadraturePointsPerElem; q++ )
    {
      //get damage at quadrature point i
      //std::cout<<"this rank is: "<<MpiWrapper::commRank( MPI_COMM_GEOSX )<<" patchElemNum is: "<<patchElemNum<<std::endl;
      average_d = average_d + qp_damage( subRegion.globalToLocalMap().at(patchElemNum), q )/8.0;
    }
    //GEOSX_LOG_LEVEL( 2, "after damage call\n" );
    //std::cout<<"this rank is: "<<MpiWrapper::commRank( MPI_COMM_GEOSX )<<" is exiting the damageFunction"<<std::endl;
    return average_d;    
}


void MultiResolutionHFSolver::initializeCrackFront( MeshLevel & base )
{
  ElementRegionManager & baseElemManager = base.getElemManager();
  FaceManager & baseFaceManager = base.getFaceManager();
  auto faceToElemList = baseFaceManager.elementList();
  auto faceNormals = baseFaceManager.faceNormal();
  
  
  //get accessor to elemental field
  ElementRegionManager::ElementViewAccessor< arrayView1d< localIndex > > const frontIndicator =
  baseElemManager.constructViewAccessor< array1d< localIndex >, arrayView1d< localIndex > >( "frontIndicator" );
  //GEOSX_LOG_LEVEL_RANK_0( 1, "num regions: "<<baseElemManager.numRegions()<<"\n" );

  baseElemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & cellElementSubRegion )
  {
    m_baseCrackFront.clear();
    //auto elemToFaceList = cellElementSubRegion.faceList();
    array2d< localIndex > const & elemToFaceList = cellElementSubRegion.faceList();
    //arrayView1d< integer const > const ghostRank = cellElementSubRegion.ghostRank();
    SortedArrayView< localIndex const > const fracturedElements = cellElementSubRegion.fracturedElementsList();
    //TODO: can we make this a forall?
    //PARALLEL: Do ghost elements count as fracturedElements? - YES
    for(auto && fracElem:fracturedElements)
    {
      //getElemFaces
      //PARALLEL: assume each rank owns all associated faces
      for(auto && face:elemToFaceList[fracElem])
      {
        if(pow(faceNormals[face][1],2) + pow(faceNormals[face][2],2) < 1e-6) //remove faces with normal pointing in x direction
        {
            continue; //these faces dont contain fractures
        }
        //PARALLEL: are all elem here owned by current node? some may be ghosts 
        //if ghost is in fracturedList, then this may contain elements from another rank
        //if ghost are not in fracturedList, this may contain ghost, but not elems from other rank <--- assumed this
        for(auto && elem:faceToElemList[face])
        {
          //PARALLEL: what does faceToElemList return if elem is ghost??
          //PARALLEL ISSUE: If elem is ghost, it might be fractured, but contains() will return it is not, so, we cant add ghosts here
          //on the other had, if this possible front element is only on the front because of this neighbor in the current rank, it will
          //not be added to the front from another rank - TODO: global fracturedElements check? 
          if(!fracturedElements.contains(elem) && elem > -1)// && ghostRank[elem] < 0)
          {
            m_baseCrackFront.insert(elem);//Maybe I should insert the ghosts here and them synchronize later using global indices?
            frontIndicator[0][0][elem] = 1;
          }
        }
      }
    }
  } );

}                                                   

void MultiResolutionHFSolver::cutDamagedElements( MeshLevel & base,
                                                  MeshLevel const & patch )
{

  EmbeddedSurfaceGenerator &
  efemGenerator = this->getParent().getGroup< EmbeddedSurfaceGenerator >( "SurfaceGenerator" ); //this is hard coded
  
  ElementRegionManager & baseElemManager = base.getElemManager();

  // Hard-coded region and subRegion
  ElementRegionBase const & elementRegion = baseElemManager.getRegion( 0 );
  CellElementSubRegion const & subRegion = elementRegion.getSubRegion< CellElementSubRegion >( 0 );

  //ghosts
  arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
  
  //FaceManager const & baseFaceManager = base.getFaceManager();
  //array2d< localIndex > const & elemToFaceList = subRegion.faceList();
  //auto faceToElemList = baseFaceManager.elementList();
  //auto faceNormals = baseFaceManager.faceNormal();
  // Get domain
  MeshManager & meshManager = this->getGroupByPath< MeshManager >( "/Problem/Mesh");
  integer Nx = meshManager.getGroup<InternalMeshGenerator>(0).getNx();
  integer Ny = meshManager.getGroup<InternalMeshGenerator>(0).getNy();
  integer Nz = meshManager.getGroup<InternalMeshGenerator>(0).getNz();
  integer nx = meshManager.getGroup<InternalMeshGenerator>(1).getNx();
  integer ny = meshManager.getGroup<InternalMeshGenerator>(1).getNy();
  integer nz = meshManager.getGroup<InternalMeshGenerator>(1).getNz();

  localIndex rx=nx/Nx;
  localIndex ry=ny/Ny;
  localIndex rz=nz/Nz;
  SortedArrayView< localIndex const > const fracturedElements = subRegion.fracturedElementsList();
  //get accessor to elemental field
  ElementRegionManager::ElementViewAccessor< arrayView1d< localIndex > > const frontIndicator =
  baseElemManager.constructViewAccessor< array1d< localIndex >, arrayView1d< localIndex > >( "frontIndicator" );

  SortedArray<localIndex> baseFrontCopy = m_baseCrackFront;

  for(auto && elem:baseFrontCopy)
  {
    // if(ghostRank[elem]>=0){
    //   continue;
    // }
    //PARALLEL: split the work to the right rank
    localIndex triplet[3];
    integer fracCount = 0;
    GEOSX_LOG_LEVEL( 3, "Entering coarseMap\n" );
    coarseToFineStructuredElemMap(subRegion.localToGlobalMap()[elem],Nx,Ny,Nz,rx,ry,rz,triplet);
    GEOSX_LOG_LEVEL( 3, "Base elem index: " << elem << " triplet: ( "<<triplet[0]<<", "<<triplet[1]<<", "<<triplet[2]<<" )\n" );
    for(int i=0; i<rx; i++){
      for(int j=0; j<ry; j++){
        for(int k=0; k<rz; k++){
          globalIndex fineElem = (triplet[0] + i)*ry*Ny*rz*Nz + (triplet[1] + j)*rz*Nz + (triplet[2] + k);
          if(fineElem > 0){
              GEOSX_LOG_LEVEL( 3, "fineElem: " << fineElem << " from base elem: "<<elem<<"\n" );
          }
          //get averageDamage in fineElem
          real64 averageDamage;
          if(ghostRank[elem]>=0){
             averageDamage = 0.0;
          }
          else{
             averageDamage = utilGetElemAverageDamage(fineElem, patch);
          }
          if (averageDamage > 0.9)
          {
            fracCount++;  
          }
        }        
      }
    }

    if (fracCount >= ry*rz || fracturedElements.contains(elem))
    {
      efemGenerator.insertToCut(elem);
    }
  }

  m_addedFractureElements = efemGenerator.propagationStep3D();
  initializeCrackFront(base);

  //SortedArray<localIndex> const & toCutList = efemGenerator.getCutList();
  // for (auto elem:toCutList)
  // {
  //   //add non-fractured neighbors to fron
  //   for(auto && face:elemToFaceList[elem])
  //   {
  //     if(pow(faceNormals[face][1],2) + pow(faceNormals[face][2],2) < 1e-6) //remove faces with normal pointing in x direction
  //     {
  //       continue; //these faces dont contain fractures
  //     }
  //     for(auto && neighbor:faceToElemList[face])
  //     {
  //       if(neighbor == elem){
  //         continue;
  //       }
  //       if(!fracturedElements.contains(neighbor) && !m_baseCrackFront.contains(neighbor) && neighbor > -1)
  //       {
  //         GEOSX_LOG_LEVEL_RANK_0( 1, "base element " << neighbor << " being added to front "<<"\n" );
  //         m_baseCrackFront.insert(neighbor);
  //         frontIndicator[0][0][neighbor] = 1;
  //       }
  //     }
  //   }
  //   //remove from front
  //   GEOSX_LOG_LEVEL_RANK_0( 1, "base element " << elem << " removed from front "<<"\n" );
  //   m_baseCrackFront.remove(elem);
  //   frontIndicator[0][0][elem] = 0;
  // }

  efemGenerator.emptyCutList();

  MpiWrapper::barrier();

}      

real64 MultiResolutionHFSolver::splitOperatorStep( real64 const & time_n,
                                                   real64 const & dt,
                                                   integer const cycleNumber,
                                                   DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  //int const thisRank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  real64 dtReturn = dt;
  real64 dtReturnTemporary;

  SolidMechanicsEmbeddedFractures &
  baseSolver = this->getParent().getGroup< SolidMechanicsEmbeddedFractures >( m_baseSolverName );

  PhaseFieldFractureSolver &
  patchSolver = this->getParent().getGroup< PhaseFieldFractureSolver >( m_patchSolverName );

  // EmbeddedSurfaceGenerator &
  // efemGenerator = this->getParent().getGroup< EmbeddedSurfaceGenerator >( "SurfaceGenerator" ); //this is hard coded

  PhaseFieldDamageFEM &
  patchDamageSolver = *patchSolver.damageSolver();

  //SolidMechanicsLagrangianFEM &
  //patchSolidSolver = *patchSolver.solidMechanicsSolver();

  baseSolver.setupSystem( domain,
                          baseSolver.getDofManager(),
                          baseSolver.getLocalMatrix(),
                          baseSolver.getSystemRhs(),
                          baseSolver.getSystemSolution(),
                          true );

  baseSolver.implicitStepSetup( time_n, dt, domain );

  this->implicitStepSetup( time_n, dt, domain );

  map< std::pair< string, string >, array1d< string > > const & baseTargets = baseSolver.getReference< map< std::pair< string, string >, array1d< string > > >(
  SolverBase::viewKeyStruct::meshTargetsString());
  auto const baseTarget = baseTargets.begin()->first;
  map< std::pair< string, string >, array1d< string > > const & patchTargets = patchSolver.getReference< map< std::pair< string, string >, array1d< string > > >(
  SolverBase::viewKeyStruct::meshTargetsString());
  auto const patchTarget = patchTargets.begin()->first;
  MeshLevel & base = domain.getMeshBody( baseTarget.first ).getBaseDiscretization();
  MeshLevel & patch = domain.getMeshBody( patchTarget.first ).getBaseDiscretization();
  if(cycleNumber == 0){
      initializeCrackFront(base);
  }
  NonlinearSolverParameters & solverParams = getNonlinearSolverParameters();
  //although these iterations are not really Newton iterations, we will use this nomeclature to keep things consistent
  integer & iter = solverParams.m_numNewtonIterations;
  iter = 0;
  bool isConverged = false;
  bool isPatchConverged = false;
  while( iter < solverParams.m_maxIterNewton )
  {
    if( iter == 0 )
    {
      // reset the states of all slave solvers if any of them has been reset
      //TODO: this is potentially a code duplication since resetStateToBeginningOfStep(domain) already calls the slaves
      patchSolver.resetStateToBeginningOfStep( domain );
      baseSolver.resetStateToBeginningOfStep( domain );
      resetStateToBeginningOfStep( domain );
    }

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1 << ", BaseSolver: " );

    //we probably want to run a phase-field solve in the patch problem at timestep 0 to get a smooth initial crack. Also, re-run this
    // anytime the base crack changes


    // testElemMappingPatchToBase( domain.getMeshBody( baseTarget.first ).getBaseDiscretization(), 
    //                             domain.getMeshBody( patchTarget.first ).getBaseDiscretization() );

    // testElemMappingBaseToPatch( domain.getMeshBody( baseTarget.first ).getBaseDiscretization(), 
    //                             domain.getMeshBody( patchTarget.first ).getBaseDiscretization() );

    CRSMatrix< real64, globalIndex > & patchDamageLocalMatrix = patchDamageSolver.getLocalMatrix();
    this->setInitialCrackDamageBCs( patchDamageSolver.getDofManager(), patchDamageLocalMatrix.toViewConstSizes(), patch,
                                    base );
    patchDamageSolver.setInitialCrackNodes( m_nodeFixDamage );

    //now perform the subproblem run with no BCs on displacements, just to set the damage inital condition;
    if( iter == 0 )
    {
      //TODO: eventually, we will need to update the patch domain
      real64 dtUseless = patchSolver.solverStep( time_n,
                                                 dtReturn,
                                                 cycleNumber,
                                                 domain );
      GEOSX_UNUSED_VAR( dtUseless );
    }

    dtReturnTemporary = baseSolver.nonlinearImplicitStep( time_n,
                                                          dtReturn,
                                                          cycleNumber,
                                                          domain );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }
    //test for convergence of MR scheme, based on changes to the fracture topology
    int added = m_addedFractureElements;
    added = MpiWrapper::sum(added);
    GEOSX_LOG_LEVEL_RANK_0(1, "isPatchConverged "<<isPatchConverged);
    if ( added == 0 && iter > 0 && isPatchConverged)
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "***** The Global-Local iterative scheme has converged in " << iter << " iterations! *****\n" );
      isConverged = true;
      break;  
    }


    // if( baseSolver.getNonlinearSolverParameters().m_numNewtonIterations >= 1 && iter > 0 )
    // {
    //   GEOSX_LOG_LEVEL_RANK_0( 1, "***** The Global-Local iterative scheme has converged in " << iter << " iterations! *****\n" );
    //   isConverged = true;
    //   break;
    // }

    //here, before calling the nonlinarImplicitStep of the patch solver, we must prescribe the displacement boundary conditions
    //this->prepareSubProblemBCs( domain.getMeshBody( baseTarget.first ).getBaseDiscretization(), domain.getMeshBody( patchTarget.first ).getBaseDiscretization());

    //write disp BCs to local disp solver
    //TODO: m_nodeFixDisp and m_fixedDispList dont need to be members, they can be initialized at every time step, this is actually safer
    //since the size of the boundary can change
    //patchSolidSolver.setInternalBoundaryConditions( m_nodeFixDisp, m_fixedDispList );

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1 << ", PatchSolver: " );

    // dtReturnTemporary = patchSolver.solverStep( time_n,
    //                                             dtReturn,
    //                                             cycleNumber,
    //                                             domain );

    dtReturnTemporary = patchSolver.unitSequentiallyCoupledSolverStep( isPatchConverged,
                                                                       time_n,
                                                                       dtReturn,
                                                                       cycleNumber,
                                                                       domain );                                            
                             
    // this->findPhaseFieldTip( m_patchTip, domain.getMeshBody( patchTarget.first ).getBaseDiscretization());
    if( time_n > 0 )
    {
      //efemGenerator.propagationStep( domain, m_baseTip, m_patchTip, m_baseTipElementIndex );
      cutDamagedElements( base, patch );  
      baseSolver.setupSystem( domain,
                              baseSolver.getDofManager(),
                              baseSolver.getLocalMatrix(),
                              baseSolver.getSystemRhs(),
                              baseSolver.getSystemSolution(),
                              true );
      baseSolver.implicitStepSetup( time_n, dt, domain );
    }

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    ++iter;
  }
  
  GEOSX_UNUSED_VAR(isConverged);
  GEOSX_LOG_LEVEL_RANK_0(1, "MultiResolutionHFSolver::SplitOperatorStep() did not converge, accepting anyway.");
  //GEOSX_ERROR_IF( !isConverged, "MultiResolutionHFSolver::SplitOperatorStep() did not converge" );

  baseSolver.implicitStepComplete( time_n, dt, domain );
  patchSolver.implicitStepComplete( time_n, dt, domain );


  return dtReturn;
}

REGISTER_CATALOG_ENTRY( SolverBase, MultiResolutionHFSolver, string const &, Group * const )
} /* namespace geosx */
