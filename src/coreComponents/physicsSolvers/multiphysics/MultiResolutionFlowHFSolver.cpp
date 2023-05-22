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


#include "MultiResolutionFlowHFSolver.hpp"

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
#include "physicsSolvers/multiphysics/SinglePhasePoromechanicsEmbeddedFractures.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/surfaceGeneration/SurfaceGenerator.hpp"
#include "physicsSolvers/surfaceGeneration/EmbeddedSurfaceGenerator.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
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

MultiResolutionFlowHFSolver::MultiResolutionFlowHFSolver( const string & name,
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
    "Name of the EFEM/EDFM (SinglePhasePoromechanicsEmbeddedFractures) solver to be used as the base solver in the MultiResolution scheme" );

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

void MultiResolutionFlowHFSolver::buildBaseToPatchMaps(const MeshLevel & base, 
                                                       const MeshLevel & patch)
{
  //need to get base and patch elem managers
  ElementRegionManager const & baseElemManager  = base.getElemManager();
  FaceManager const & baseFaceManager = base.getFaceManager();
  ElementRegionManager const & patchElemManager = patch.getElemManager();
  arrayView2d<const real64> const & baseFaceCenters  = baseFaceManager.faceCenter();
  arrayView2d<const real64> const & baseFaceNormals  = baseFaceManager.faceNormal();
  //get elem to face map
  //elem -> elemToFaces -> loop over these faces only
  
  //get base element faces

  //here we will loop over all element of base and patch and write the relation between them in the maps
  //loop over all base subregions
  baseElemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & cellElementSubRegion )
  {
    array2d< localIndex > const & baseFaceList = cellElementSubRegion.faceList();
    arrayView2d< const real64 > const & baseElemCenters  = cellElementSubRegion.getElementCenter();
    arrayView1d< integer const > const baseGhostRank = cellElementSubRegion.ghostRank();
    //loop over all elements
    forAll< serialPolicy >( cellElementSubRegion.size(), [&] ( localIndex const k )
    {
      //skip ghosts from global
      if(baseGhostRank[k]<0)
      {
        arraySlice1d<const real64> baseCenter = baseElemCenters[k];
        patchElemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & patchCellElementSubRegion )
        {
          arrayView1d< integer const > const & patchGhostRank = patchCellElementSubRegion.ghostRank();
          arrayView2d< const real64> const & patchElemCenters  = patchCellElementSubRegion.getElementCenter();
          forAll< serialPolicy >( patchCellElementSubRegion.size(), [&] ( localIndex const n )
          {
            //skip ghosts
            if(patchGhostRank[n]<0)
            {              
                arraySlice1d<const real64> patchCenter = patchElemCenters[n];
                bool isInside = true;
                //loop over faces of base element
                for(int f=0; f<baseFaceList[k].size(); f++)
                {
                  arraySlice1d<const real64> faceCenter = baseFaceCenters[baseFaceList[k][f]];
                  arraySlice1d<const real64> faceNormal = baseFaceNormals[baseFaceList[k][f]];
                  //array1d<real64> baseCenterMod = {baseCenter[0], baseCenter[1], baseCenter[2]};
                  //array1d<real64> patchCenterMod = {patchCenter[0], patchCenter[1], patchCenter[2]};
                  //I need these copies because the tensorOps doesnt work on ArraySlices
                  R1Tensor baseCenterMod = {baseCenter[0], baseCenter[1], baseCenter[2]};
                  R1Tensor patchCenterMod = {patchCenter[0], patchCenter[1], patchCenter[2]};                  
                  //baseCenterMod[0] = baseCenter[0]; baseCenterMod[1] = baseCenter[1]; baseCenterMod[2] = baseCenter[2];
                  //patchCenterMod[0] = patchCenter[0]; patchCenterMod[1] = patchCenter[1]; patchCenterMod[2] = patchCenter[2];
                  LvArray::tensorOps::subtract< 3 >(baseCenterMod,faceCenter);
                  LvArray::tensorOps::subtract< 3 >(patchCenterMod,faceCenter);
                  real64 s1 = LvArray::tensorOps::AiBi< 3 >(baseCenterMod, faceNormal);
                  real64 s2 = LvArray::tensorOps::AiBi< 3 >(patchCenterMod, faceNormal);
                  //std::cout<<"Base Element: "<<K<<", Patch Element: "<<N<<std::endl;
                  bool debugMode = false;
                  if(k==0 && n<12 && debugMode)
                  {
                    //DEBUG INFO
                    std::cout<<"Base element "<<k<<"\n";
                    std::cout<<"Patch element "<<n<<"\n";
                    std::cout<<"Face "<<f<<", s1*s2 = "<<s1*s2<<std::endl;
                    std::cout<<"baseCenter: ["<<baseCenter[0]<<", "<<baseCenter[1]<<", "<<baseCenter[2]<<"]\n";
                    std::cout<<"patchCenter: ["<<patchCenter[0]<<", "<<patchCenter[1]<<", "<<patchCenter[2]<<"]\n";
                    std::cout<<"baseCenterMod: ["<<baseCenterMod[0]<<", "<<baseCenterMod[1]<<", "<<baseCenterMod[2]<<"]\n";
                    std::cout<<"patchCenterMod: ["<<patchCenterMod[0]<<", "<<patchCenterMod[1]<<", "<<patchCenterMod[2]<<"]\n";
                    std::cout<<"faceCenter: ["<<faceCenter[0]<<", "<<faceCenter[1]<<", "<<faceCenter[2]<<"]\n";
                    std::cout<<"faceNormal: ["<<faceNormal[0]<<", "<<faceNormal[1]<<", "<<faceNormal[2]<<"]\n";
                  }
                  if(s1*s2 <= 0)
                  {
                    isInside = false;
                  }
                }
                if(isInside)
                {
                  //get globalIndices
                  globalIndex N = patchCellElementSubRegion.localToGlobalMap()[n];
                  globalIndex K = cellElementSubRegion.localToGlobalMap()[k];
                  //std::cout<<"s1*s2>0 / Base Element: "<<K<<", Patch Element: "<<N<<std::endl;
                  m_patchToBaseElementRelation[N] = K;
                  //if K is already in the map, just add another patchElem to its set
                  if(m_baseToPatchElementRelation.find(K)!=m_baseToPatchElementRelation.end())
                  {
                    m_baseToPatchElementRelation[K].insert(N);
                  }
                  else
                  {
                    m_baseToPatchElementRelation[K] = {N};
                  }
                }
            }
          });

        });

        //}
      }

    });
  });

  //MPI communication to synchronize the maps accross ranks - straight from my buddy chatGPT
  // Gather the local maps from all ranks onto each rank
  // MPI_Comm comm = base.getComm();
  // int commSize;
  // MPI_Comm_size(MPI_COMM_GEOSX, &commSize);
  // std::vector<map<globalIndex, set<globalIndex>>> allBaseToPatchMaps(commSize);
  // std::vector<map<globalIndex, globalIndex>> allPatchToBaseMaps(commSize);

  // MPI_Allgather(&m_baseToPatchElementRelation, sizeof(m_baseToPatchElementRelation), MPI_BYTE,
  //               allBaseToPatchMaps.data(), sizeof(m_baseToPatchElementRelation), MPI_BYTE,
  //               MPI_COMM_GEOSX);
  // MPI_Allgather(&m_patchToBaseElementRelation, sizeof(m_patchToBaseElementRelation), MPI_BYTE,
  //               allPatchToBaseMaps.data(), sizeof(m_patchToBaseElementRelation), MPI_BYTE,
  //               MPI_COMM_GEOSX);

  // // Merge the gathered maps into a single map on each rank
  // map<globalIndex, set<globalIndex>> mergedBaseToPatchMap;
  // map<globalIndex, globalIndex> mergedPatchToBaseMap;

  // for (int rank = 0; rank < commSize; ++rank) {
  //   for (auto& entry : allBaseToPatchMaps[rank]) {
  //     const auto& globalIndex = entry.first;
  //     const auto& patchIndices = entry.second;
  //     if (mergedBaseToPatchMap.find(globalIndex) != mergedBaseToPatchMap.end()) {
  //       mergedBaseToPatchMap[globalIndex].insert(patchIndices.begin(), patchIndices.end());
  //     } else {
  //       mergedBaseToPatchMap[globalIndex] = patchIndices;
  //     }
  //   }

  //   for (auto& entry : allPatchToBaseMaps[rank]) {
  //     const auto& globalIndex = entry.first;
  //     const auto& baseIndex = entry.second;
  //     mergedPatchToBaseMap[globalIndex] = baseIndex;
  //   }
  // }

  // // Update the local maps to include the merged maps
  // m_baseToPatchElementRelation = mergedBaseToPatchMap;
  // m_patchToBaseElementRelation = mergedPatchToBaseMap;
  
}

void MultiResolutionFlowHFSolver::registerDataOnMesh( dataRepository::Group & meshBodies )
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

    const MeshLevel & base  = meshBodies.getGroup<MeshBody>(0).getBaseDiscretization();
    const MeshLevel & patch = meshBodies.getGroup<MeshBody>(1).getBaseDiscretization();

    buildBaseToPatchMaps(base, patch);

  // } );
}

MultiResolutionFlowHFSolver::~MultiResolutionFlowHFSolver()
{
  // TODO Auto-generated destructor stub
}

void MultiResolutionFlowHFSolver::implicitStepSetup( real64 const & time_n,
                                                     real64 const & dt,
                                                     DomainPartition & domain )
{
  m_baseSolver->implicitStepSetup( time_n, dt, domain );
  m_patchSolver->implicitStepSetup( time_n, dt, domain );
}

void MultiResolutionFlowHFSolver::postProcessInput()
{
  SolverBase::postProcessInput();
  m_baseSolver = &this->getParent().getGroup< SinglePhasePoromechanicsEmbeddedFractures >( m_baseSolverName );
  m_patchSolver = &this->getParent().getGroup< PhaseFieldFractureSolver >( m_patchSolverName );
}

void MultiResolutionFlowHFSolver::initializePostInitialConditionsPreSubGroups()
{}

void MultiResolutionFlowHFSolver::resetStateToBeginningOfStep( DomainPartition & domain )
{
  m_baseSolver->resetStateToBeginningOfStep( domain );
  m_patchSolver->resetStateToBeginningOfStep( domain );
}

real64 MultiResolutionFlowHFSolver::solverStep( real64 const & time_n,
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
void MultiResolutionFlowHFSolver::setInitialCrackDamageBCs( DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
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

void MultiResolutionFlowHFSolver::findNewlyDamagedElements(SortedArray<localIndex> GEOSX_UNUSED_PARAM(toCutElems), 
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
void MultiResolutionFlowHFSolver::findPhaseFieldTip( R1Tensor & tip,
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
void MultiResolutionFlowHFSolver::prepareSubProblemBCs( MeshLevel const & base,
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

void MultiResolutionFlowHFSolver::testElemMappingPatchToBase( MeshLevel & GEOSX_UNUSED_PARAM(base),
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
  patchElemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & cellElementSubRegion )
  {
    arrayView1d< integer const > const & patchGhostRank = cellElementSubRegion.ghostRank();
    forAll< serialPolicy >( cellElementSubRegion.size(), [&] ( localIndex const k )
    {
      if(patchGhostRank[k]<0)
      {                 
        globalIndex K = cellElementSubRegion.localToGlobalMap()[k];
        fineToCoarseMap[0][0][k] = m_patchToBaseElementRelation[K];
      }
    } );
  } );  
}

void MultiResolutionFlowHFSolver::testElemMappingBaseToPatch( MeshLevel & base,
                                                              MeshLevel & patch )
{
  //for every elem in patch
  //call fineToCoarseStructuredElemMap
  //write index of the associated base elem to all patch elems
  //this should be exactly the same as the function above
  ElementRegionManager & coarseElemManager = base.getElemManager();
  ElementRegionManager & patchElemManager = patch.getElemManager();
  CellElementSubRegion const & patchElementSubRegion = patchElemManager.getRegion(0).getSubRegion< CellElementSubRegion >( 0 );
  
  //get accessor to elemental field
  ElementRegionManager::ElementViewAccessor< arrayView1d< localIndex > > const coarseToFineMap =
  patchElemManager.constructViewAccessor< array1d< localIndex >, arrayView1d< localIndex > >( "coarseToFineMap" );

  coarseElemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & cellElementSubRegion )
  {
    arrayView1d< integer const > const & baseGhostRank = cellElementSubRegion.ghostRank();
    forAll< serialPolicy >( cellElementSubRegion.size(), [&] ( localIndex const k )
    {
      //loop over all fines in m_baseToPatchElementRelation[K]
      if(baseGhostRank[k]<0)
      {
        globalIndex K = cellElementSubRegion.localToGlobalMap()[k];
        for(auto fineElem:m_baseToPatchElementRelation[K])
        {
          // try{
          //   localIndex fineElemLocal = patchElementSubRegion.globalToLocalMap().at(fineElem);
          //   coarseToFineMap[0][0][fineElemLocal] = K;   
          //   throw;
          // }
          // catch (...){
          //   std::cout<<"rank "<<MpiWrapper::commRank( MPI_COMM_GEOSX )<<" failed in globalToLocalMap at patchElementSubRegion "<<fineElem<<std::endl;
          // }  
          localIndex fineElemLocal = patchElementSubRegion.globalToLocalMap().at(fineElem);
          coarseToFineMap[0][0][fineElemLocal] = K;          
        }
      }
    } );
  } );

}


real64 MultiResolutionFlowHFSolver::utilGetElemAverageDamage(globalIndex const patchElemNum,
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
      average_d = average_d + qp_damage( subRegion.globalToLocalMap().at(patchElemNum), q )/8.0;
    }
    //GEOSX_LOG_LEVEL( 2, "after damage call\n" );
    //std::cout<<"this rank is: "<<MpiWrapper::commRank( MPI_COMM_GEOSX )<<" is exiting the damageFunction"<<std::endl;
    return average_d;    
}


void MultiResolutionFlowHFSolver::initializeCrackFront( MeshLevel & base )
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

void MultiResolutionFlowHFSolver::cutDamagedElements( MeshLevel & base,
                                                      MeshLevel const & patch )
{

  EmbeddedSurfaceGenerator &
  efemGenerator = this->getParent().getGroup< EmbeddedSurfaceGenerator >( "SurfaceGenerator" ); //this is hard coded
  
  ElementRegionManager & baseElemManager = base.getElemManager();

  // Hard-coded region and subRegion
  ElementRegionBase const & elementRegion = baseElemManager.getRegion( 0 );
  CellElementSubRegion const & subRegion = elementRegion.getSubRegion< CellElementSubRegion >( 0 );

  //ghosts
  arrayView1d< integer const > const baseGhostRank = subRegion.ghostRank();
  
  //FaceManager const & baseFaceManager = base.getFaceManager();
  //array2d< localIndex > const & elemToFaceList = subRegion.faceList();
  //auto faceToElemList = baseFaceManager.elementList();
  //auto faceNormals = baseFaceManager.faceNormal();
  // Get domain
  MeshManager & meshManager = this->getGroupByPath< MeshManager >( "/Problem/Mesh");
  //integer Nx = meshManager.getGroup<InternalMeshGenerator>(0).getNx();
  integer Ny = meshManager.getGroup<InternalMeshGenerator>(0).getNy();
  integer Nz = meshManager.getGroup<InternalMeshGenerator>(0).getNz();
  //integer nx = meshManager.getGroup<InternalMeshGenerator>(1).getNx();
  integer ny = meshManager.getGroup<InternalMeshGenerator>(1).getNy();
  integer nz = meshManager.getGroup<InternalMeshGenerator>(1).getNz();

  //localIndex rx=nx/Nx;
  localIndex ry=ny/Ny;
  localIndex rz=nz/Nz;
  SortedArrayView< localIndex const > const fracturedElements = subRegion.fracturedElementsList();
  //get accessor to elemental field
  ElementRegionManager::ElementViewAccessor< arrayView1d< localIndex > > const frontIndicator =
  baseElemManager.constructViewAccessor< array1d< localIndex >, arrayView1d< localIndex > >( "frontIndicator" );

  SortedArray<localIndex> baseFrontCopy = m_baseCrackFront;

  for(auto && elem:baseFrontCopy)
  {
    //PARALLEL: split the work to the right rank
    //localIndex triplet[3];
    integer fracCount = 0;
    GEOSX_LOG_LEVEL( 3, "Entering coarseMap\n" );
    ////////version with new maps

    //loop over all fines in m_baseToPatchElementRelation[K]
    if(baseGhostRank[elem]<0)
    {
      for(auto fineElemGlobal:m_baseToPatchElementRelation[subRegion.localToGlobalMap()[elem]])
      { 
        //get averageDamage in fineElem
        real64 averageDamage = utilGetElemAverageDamage(fineElemGlobal, patch);
        if (averageDamage > 0.9)
        {
          fracCount++;  
        }        
      }
    }
    ///////////
    //////////other version
    // coarseToFineStructuredElemMap(subRegion.localToGlobalMap()[elem],Nx,Ny,Nz,rx,ry,rz,triplet);
    // GEOSX_LOG_LEVEL( 3, "Base elem index: " << elem << " triplet: ( "<<triplet[0]<<", "<<triplet[1]<<", "<<triplet[2]<<" )\n" );
    // for(int i=0; i<rx; i++){
    //   for(int j=0; j<ry; j++){
    //     for(int k=0; k<rz; k++){
    //       globalIndex fineElem = (triplet[0] + i)*ry*Ny*rz*Nz + (triplet[1] + j)*rz*Nz + (triplet[2] + k);
    //       if(fineElem > 0){
    //           GEOSX_LOG_LEVEL( 3, "fineElem: " << fineElem << " from base elem: "<<elem<<"\n" );
    //       }
    //       //get averageDamage in fineElem
    //       real64 averageDamage;
    //       //I dont think we need to cut the ghosts
    //       if(ghostRank[elem]>=0){
    //          averageDamage = 0.0;
    //       }
    //       else{
    //          averageDamage = utilGetElemAverageDamage(fineElem, patch);
    //       }
    //       if (averageDamage > 0.9)
    //       {
    //         fracCount++;  
    //       }
    //     }        
    //   }
    // }
    //////////
    if (fracCount >= ry*rz || fracturedElements.contains(elem))
    {
      efemGenerator.insertToCut(elem);
    }
  }

  m_addedFractureElements = efemGenerator.propagationStep3D();
  initializeCrackFront(base);

  efemGenerator.emptyCutList();

  MpiWrapper::barrier();

}      

void MultiResolutionFlowHFSolver::writeBasePressuresToPatch(MeshLevel & base,
                                                            MeshLevel & patch)
{
  using namespace fields::flow;

  ElementRegionManager & patchElemManager = patch.getElemManager();

  ElementRegionManager & baseElemManager = base.getElemManager();

  // Hard-coded region and subRegion
  ElementRegionBase const & baseMatrixElementRegion = baseElemManager.getRegion( 0 );
  ElementRegionBase const & baseFractureElementRegion = baseElemManager.getRegion( 1 );
  CellElementSubRegion const & patchSubRegion = patchElemManager.getRegion( 0 ).getSubRegion< CellElementSubRegion >( 0 );
  CellElementSubRegion const & baseMatrixSubRegion = baseMatrixElementRegion.getSubRegion< CellElementSubRegion >( 0 );
  EmbeddedSurfaceSubRegion const & baseFractureSubRegion = baseFractureElementRegion.getSubRegion< EmbeddedSurfaceSubRegion >( 0 );

  //get accessor to elemental field
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const patchMatrixPressure =
  patchElemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( "hardCodedPMatrixName" );

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const patchFracturePressure =
  patchElemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( "hardCodedPFractureName" );

  //get fields from SinglePhaseFlow problem (or maybe multiphase later)
  arrayView1d< real64 const > const baseMatrixPressure = baseMatrixSubRegion.getField< fields::flow::pressure >();

  //get fields from SinglePhaseFlow problem (or maybe multiphase later)
  arrayView1d< real64 const > const baseFracturePressure = baseFractureSubRegion.getField< fields::flow::pressure >();

  //loop over patch subRegion and fill values for matrixPressure and fracturePressure
  //TODO: make this a regular for so that all ranks will loop over all elements?
  //allgather fracture cell centers and pressure values - be careful with ghosts
  //get cell centers of all fracture elements
  arrayView2d< const real64 > allFracElemCenters = baseFractureSubRegion.getElementCenter();
  array2d< real64 > toSend(baseFractureSubRegion.size(),4);
  forAll< serialPolicy >(baseFractureSubRegion.size(), [&] (localIndex const k)
  {

    toSend[k][0] = allFracElemCenters[k][0];
    toSend[k][1] = allFracElemCenters[k][1];
    toSend[k][2] = allFracElemCenters[k][2];
    toSend[k][3] = baseFracturePressure[k];

  });
  // Exchange the sizes of the data across all ranks.
  array1d< int > dataSizes( MpiWrapper::commSize() );
  MpiWrapper::allGather( LvArray::integerConversion< int >( toSend.size(0)*toSend.size(1) ), dataSizes, MPI_COMM_GEOSX );

  int const totalDataSize = std::accumulate( dataSizes.begin(), dataSizes.end(), 0 );
  GEOSX_LOG_LEVEL_RANK_0( 1, "totalDataSize = "<<totalDataSize<<"\n");
  //array1d<T> allData( totalDataSize );
  array2d<real64> allData( totalDataSize/4, 4 );
  std::vector< int > mpiDisplacements( MpiWrapper::commSize(), 0 );
  std::partial_sum( dataSizes.begin(), dataSizes.end() - 1, mpiDisplacements.begin() + 1 );
  MpiWrapper::allgatherv( toSend.data(), 
                          toSend.size(0)*toSend.size(1), 
                          allData.data(),
                          dataSizes.data(), 
                          mpiDisplacements.data(), 
                          MPI_COMM_GEOSX );                         
  arrayView2d< const real64 > patchElemCenters = patchElemManager.getRegion(0).getSubRegion< CellElementSubRegion >(0).getElementCenter();
  arrayView1d< integer const > const & patchGhostRank = patchSubRegion.ghostRank();
  forAll< serialPolicy >( patchSubRegion.size(), [&] ( localIndex const k )
  {
    //convert patch element index to base - first convert patch to global index
    //this probably needs to be done at ghosts or require a sync operation at the end
    if(patchGhostRank[k]<0){
      globalIndex K = patchSubRegion.localToGlobalMap()[k];
      globalIndex baseK = m_patchToBaseElementRelation[K];
      //globalIndex baseK = fineToCoarseStructuredElemMap(K,Nx,Ny,Nz,rx,ry,rz);
      localIndex basek = 0;
      basek = baseMatrixSubRegion.globalToLocalMap().at(baseK);
      // try{
      //   basek = baseMatrixSubRegion.globalToLocalMap().at(baseK);
      //   throw;
      // }
      // catch (...){
      //   std::cout<<"rank "<<MpiWrapper::commRank( MPI_COMM_GEOSX )<<" failed in globalToLocalMap at baseGlobalElement "<<baseK<<std::endl;
      // }
      patchMatrixPressure[0][0][k] = baseMatrixPressure[basek];

      R1Tensor patchElemCenter;
      patchElemCenter[0] = patchElemCenters[k][0];
      patchElemCenter[1] = patchElemCenters[k][1];
      patchElemCenter[2] = patchElemCenters[k][2];
      localIndex minIndex = 0;
      real64 minDist = 1.0e20;
      for(localIndex j=0; j<allData.size(0); j++)
      {
        R1Tensor currFracElemCenter;
        currFracElemCenter[0] = allData[j][0];
        currFracElemCenter[1] = allData[j][1];
        currFracElemCenter[2] = allData[j][2];
        
        R1Tensor centerCopy = LVARRAY_TENSOROPS_INIT_LOCAL_3( patchElemCenter );
        LvArray::tensorOps::subtract< 3 >( centerCopy, currFracElemCenter );
        real64 dist = LvArray::tensorOps::l2Norm< 3 >( centerCopy );
        if(dist < minDist){
          minDist = dist;
          minIndex = j;
        }
      }

      patchFracturePressure[0][0][k] = allData[minIndex][3];
    } //TODO: needs to find closest fracture cell from base
  } );
      
}
                                         
real64 MultiResolutionFlowHFSolver::splitOperatorStep( real64 const & time_n,
                                                       real64 const & dt,
                                                       integer const cycleNumber,
                                                       DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  //int const thisRank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  real64 dtReturn = dt;
  real64 dtReturnTemporary;

  SinglePhasePoromechanicsEmbeddedFractures &
  baseSolver = this->getParent().getGroup< SinglePhasePoromechanicsEmbeddedFractures >( m_baseSolverName );

  PhaseFieldFractureSolver &
  patchSolver = this->getParent().getGroup< PhaseFieldFractureSolver >( m_patchSolverName );

  // EmbeddedSurfaceGenerator &
  // efemGenerator = this->getParent().getGroup< EmbeddedSurfaceGenerator >( "SurfaceGenerator" ); //this is hard coded

  PhaseFieldDamageFEM &
  patchDamageSolver = *patchSolver.damageSolver();

  SolidMechanicsLagrangianFEM &
  patchSolidSolver = *patchSolver.solidMechanicsSolver();

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
      // patchSolver.resetStateToBeginningOfStep( domain );
      // baseSolver.resetStateToBeginningOfStep( domain );
      // resetStateToBeginningOfStep( domain );
    }

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1 << ", BaseSolver: " );

    //we probably want to run a phase-field solve in the patch problem at timestep 0 to get a smooth initial crack. Also, re-run this
    // anytime the base crack changes


    testElemMappingPatchToBase( domain.getMeshBody( baseTarget.first ).getBaseDiscretization(), 
                                domain.getMeshBody( patchTarget.first ).getBaseDiscretization() );

    testElemMappingBaseToPatch( domain.getMeshBody( baseTarget.first ).getBaseDiscretization(), 
                                domain.getMeshBody( patchTarget.first ).getBaseDiscretization() );

    CRSMatrix< real64, globalIndex > & patchDamageLocalMatrix = patchDamageSolver.getLocalMatrix();
    this->setInitialCrackDamageBCs( patchDamageSolver.getDofManager(), patchDamageLocalMatrix.toViewConstSizes(), patch,
                                    base );
    patchDamageSolver.setInitialCrackNodes( m_nodeFixDamage );

    //now perform the subproblem run with no BCs on displacements, just to set the damage inital condition;
    // if( iter == 0 )
    // {
    //   //TODO: eventually, we will need to update the patch domain
    //   real64 dtUseless = patchSolver.solverStep( time_n,
    //                                              dtReturn,
    //                                              cycleNumber,
    //                                              domain );
    //   GEOSX_UNUSED_VAR( dtUseless );
    // }

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

    //if MR iteration = 0 or add = true
    if(iter==0 || added!=0){
      dtReturnTemporary = baseSolver.nonlinearImplicitStep( time_n,
                                                            dtReturn,
                                                            cycleNumber,
                                                            domain );  
    }                                            

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
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
    //THIS IS LIKELY TRANSFERING PRESSURE_N
    //patchSolver.implicitStepComplete( time_n, dt, domain );
    baseSolver.implicitStepComplete( time_n, dt, domain );
    patchSolidSolver.setInternalBoundaryConditions( m_nodeFixDisp, m_fixedDispList );
    writeBasePressuresToPatch(base, patch);

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
    isPatchConverged=true;                                                                                                               
                             
    // this->findPhaseFieldTip( m_patchTip, domain.getMeshBody( patchTarget.first ).getBaseDiscretization());
    if( time_n > 0 )
    {
      //efemGenerator.propagationStep( domain, m_baseTip, m_patchTip, m_baseTipElementIndex );
      //cutDamagedElements( base, patch );  
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
  GEOSX_LOG_LEVEL_RANK_0(1, "MultiResolutionFlowHFSolver::SplitOperatorStep() did not converge, accepting anyway.");
  //GEOSX_ERROR_IF( !isConverged, "MultiResolutionFlowHFSolver::SplitOperatorStep() did not converge" );

  baseSolver.implicitStepComplete( time_n, dt, domain );
  patchSolver.implicitStepComplete( time_n, dt, domain );


  return dtReturn;
}

REGISTER_CATALOG_ENTRY( SolverBase, MultiResolutionFlowHFSolver, string const &, Group * const )
} /* namespace geosx */
