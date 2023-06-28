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
 * @file NonplanarTracker.cpp
 *
 */

#include "NonplanarTracker.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

NonplanarTracker::NonplanarTracker( const string & name,
                                    Group * const parent ):
  Base( name, parent )
{}

NonplanarTracker::~NonplanarTracker()
{
  // TODO Auto-generated destructor stub
}

void NonplanarTracker::buildBaseToPatchMaps( const MeshLevel & base, 
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
                  //I need these copies because the tensorOps doesnt work on ArraySlices
                  R1Tensor baseCenterMod = {baseCenter[0], baseCenter[1], baseCenter[2]};
                  R1Tensor patchCenterMod = {patchCenter[0], patchCenter[1], patchCenter[2]};                  
                  LvArray::tensorOps::subtract< 3 >(baseCenterMod,faceCenter);
                  LvArray::tensorOps::subtract< 3 >(patchCenterMod,faceCenter);
                  real64 s1 = LvArray::tensorOps::AiBi< 3 >(baseCenterMod, faceNormal);
                  real64 s2 = LvArray::tensorOps::AiBi< 3 >(patchCenterMod, faceNormal);
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

      }

    });
  });

}

void NonplanarTracker::buildBaseToPatchEdgeMap( meshLevel const & base, meshLevel const & patch )
{
    // Get num cuts of each edge
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
    //loop over all base elems using map<globalIndex, set<globalIndex>> m_baseToPatchElementRelation;
    ElementRegionManager const & baseElemManager = base.getElemManager();
    //get elemsToNodes for patch
    NodeManager const & patchNodeManager = patch.getNodeManager();
    ElementRegionManager const & patchElemManager = patch.getElemManager();
    // Hard-coded patch region and subRegion
    ElementRegionBase const & patchElementRegion = patchElemManager.getRegion( 1 );
    CellElementSubRegion const & patchSubRegion = patchElementRegion.getSubRegion< CellElementSubRegion >( 1 );
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const & patchElemsToNodes = patchSubRegion.nodeList();   
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & patchNodePosition = patchNodeManager.referencePosition();
    baseElemManager.forElementSubRegions< CellElementSubRegion,
                                          FaceElementSubRegion >( [&] ( auto & baseElementSubRegion )
    {       
        arrayView1d< integer const > const baseGhostRank = baseElementSubRegion.ghostRank();  
        arrayView2d< const real64 > const & baseElemCenters  = baseElementSubRegion.getElementCenter();
        auto & cellToEdges = baseElementSubRegion.edgeList();
        forAll< serialPolicy >( baseElementSubRegion.size(), [&] ( localIndex const k ) 
        {
            if(baseGhostRank[k]<0)
            {
                //get element center
                real64 currBaseElemCenter[ 3 ];
                LvArray::tensorOps::copy< 3 >( currBaseElemCenter, baseElemCenters[k] );
                //build edgesToTriplet map
                auto elementEdges = cellToEdges[k];
                map<R1Tensor, localIndex> tripletToEdge;
                //iterate over edges and find their associated triplet - build edge to triplet map
                for(auto edge:elementEdges)
                {
                    //check triplet associated with edge
                    //get edgeToNode at edge
                    R1Tensor nodeA = edgeToNode[edge][0];
                    R1Tensor nodeB = edgeToNode[edge][1];
                    R1Tensor edgeMidPoint = (nodeA + nodeB)/2;
                    R1Tensor diff = edgeMidPoint - currBaseElemCenter;
                    //diff = (+- hx/2, +- hy/2, 0) or permutation
                    diff = diff/{hx/2, hy/2, hz/2};
                    tripletToEdge[diff] = edge;
                }
                for(auto patchElemGlobal:m_baseToPatchElementRelation[baseElementSubRegion.localToGlobalMap()[k]])
                { 
                    //get nodes
                    for(size_t i=0; i<NUM_NODES_PER_ELEM; ++i)
                    {
                        localIndex currNode = patchElemsToNodes[patchSubRegion.globalToLocalMap()[patchElemGlobal]][i];
                        real64 currNodePosition[ 3 ]; 
                        LvArray::tensorOps::copy< 3 >( currNodePosition, patchNodePosition[currNode] );
                        //compare currNodePosition to baseElemCenters[k];
                        R1Tensor nodalTriplet = {0,0,0};
                        nodalTiplet[0] = currNodePosition[0]-baseElemCenter[0]/hx;
                        nodalTiplet[1] = currNodePosition[1]-baseElemCenter[1]/hy;
                        nodalTiplet[2] = currNodePosition[2]-baseElemCenter[2]/hz;
                        if (tipletToEdge.contains(nodalTriplet))
                        {
                            //relate current node to edge tripletToEdge[nodalTriplet];
                            //find position using triplet
                            //position = if abs(nodalTriplet[0]) < 1 
                            //position = (1+nodalTriplet[0]/hx)*(Nx);
                            //same for nodalTriplet[1] por [2]
                            //insert to edge in the position order of SortedArray
                        }                                                                         
                    }
                }
            }      
        });
    });




}

void NonplanarTracker::registerDataOnMesh( dataRepository::Group & meshBodies )
{

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
    buildBaseToPatchEdgeMap(base, patch);

}

void NonplanarTracker::initializeFracturedElements( MeshLevel & base, MeshLevel & patch )
{}      

void NonplanarTracker::initializeCrackFront( MeshLevel & base )
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

void NonplanarTracker::cutDamagedElements( MeshLevel & base,
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

real64 NonplanarTracker::sequentiallyCoupledSolverStep( real64 const & time_n,
                                                        real64 const & dt,
                                                        int const cycleNumber,
                                                        DomainPartition & domain )
  {
    GEOS_MARK_FUNCTION;

    real64 dtReturn = dt;

    real64 dtReturnTemporary;

    Timestamp const meshModificationTimestamp = getMeshModificationTimestamp( domain );

    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {

      // Only build the sparsity pattern if the mesh has changed
      if( meshModificationTimestamp > solver->getSystemSetupTimestamp() )
      {
        solver->setupSystem( domain,
                             solver->getDofManager(),
                             solver->getLocalMatrix(),
                             solver->getSystemRhs(),
                             solver->getSystemSolution() );
        solver->setSystemSetupTimestamp( meshModificationTimestamp );
      }

      solver->implicitStepSetup( time_n, dt, domain );

    } );

    NonlinearSolverParameters & solverParams = getNonlinearSolverParameters();
    integer & iter = solverParams.m_numNewtonIterations;
    iter = 0;
    bool isConverged = false;
    /// Sequential coupling loop
    while( iter < solverParams.m_maxIterNewton )
    {
      if( iter == 0 )
      {
        // Reset the states of all solvers if any of them had to restart
        forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
        {
          solver->resetStateToBeginningOfStep( domain );
        } );
        resetStateToBeginningOfStep( domain );
      }

      // Increment the solver statistics for reporting purposes
      // Pass a "0" as argument (0 linear iteration) to skip the output of linear iteration stats at the end
      m_solverStatistics.logNonlinearIteration( 0 );

      // Solve the subproblems nonlinearly
      forEachArgInTuple( m_solvers, [&]( auto & solver, auto idx )
      {
        GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "  Iteration {:2}: {}", iter+1, solver->getName() ) );
        dtReturnTemporary = solver->nonlinearImplicitStep( time_n,
                                                           dtReturn,
                                                           cycleNumber,
                                                           domain );

        mapSolutionBetweenSolvers( domain, idx() );

        if( dtReturnTemporary < dtReturn )
        {
          iter = 0;
          dtReturn = dtReturnTemporary;
        }
      } );

      // Check convergence of the outer loop
      isConverged = checkSequentialConvergence( iter,
                                                time_n,
                                                dtReturn,
                                                domain );

      if( isConverged )
      {
        break;
      }
      // Add convergence check:
      ++iter;
    }

    GEOS_ERROR_IF( !isConverged, getName() << "::sequentiallyCoupledSolverStep did not converge!" );

    implicitStepComplete( time_n, dt, domain );

    return dtReturn;
  }

REGISTER_CATALOG_ENTRY( SolverBase, NonplanarTracker, string const &, Group * const )

} /* namespace geos */