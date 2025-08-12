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
{}

SolidMechanicsMortarContact::~SolidMechanicsMortarContact()
{}

void SolidMechanicsMortarContact::registerDataOnMesh( dataRepository::Group & meshBodies )
{
  GEOS_MARK_FUNCTION;
  ContactSolverBase::registerDataOnMesh( meshBodies );
}

void SolidMechanicsMortarContact::setupDofs( DomainPartition const & domain,
                                             DofManager & dofManager ) const
{
  GEOS_MARK_FUNCTION;
  SolidMechanicsLagrangianFEM::setupDofs( domain, dofManager );
}

void SolidMechanicsMortarContact::setupSystem( DomainPartition & domain,
                                                            DofManager & dofManager,
                                                            CRSMatrix< real64, globalIndex > & localMatrix,
                                                            ParallelVector & rhs,
                                                            ParallelVector & solution,
                                                            bool const setSparsity )
{

  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( domain, dofManager, localMatrix, rhs, solution, setSparsity );

}

void SolidMechanicsMortarContact::implicitStepSetup( real64 const & time_n,
                                                     real64 const & dt,
                                                     DomainPartition & domain )
{

  GEOS_MARK_FUNCTION;
  SolidMechanicsLagrangianFEM::implicitStepSetup( time_n, dt, domain );

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
    }
    
    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< FaceElementSubRegion >([&](const FaceElementSubRegion & subRegion )
    {
      // assign surface to master or slave depending on object path (again, naive)
      string surfacePath = subRegion.getPath();
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

  
  std::cout << "Processing master/slave connectivity" << std::endl;

  this->getConnectivityMap();
  
  //////////////////////////////////////////////////////////////////////////////////
  GEOS_ERROR("Mortar solver is not implemented yet. ");

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

void SolidMechanicsMortarContact::updateState( DomainPartition & domain )
{
  GEOS_UNUSED_VAR( domain );
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
  m_connectivityMap.resize(nSurfSlave);

  // perform contact search and populate connectivity map
  contactSearch(treeSlave,treeMaster);

    for (localIndex i=0; i<nSurfSlave; ++i)
  {
    localIndex N = m_connectivityMap.sizeOfArray(i);
    std::cout << "SLAVE ELEMENT: " << i << std::endl;
    std::cout << "MASTER ELEMENT:"; 
    for (localIndex j=0; j<N; ++j)
    {  
      std::cout << " " << m_connectivityMap[i][j];
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
    m_connectivityMap.emplaceBack(nodeSlave->leafId,nodeMaster->leafId);
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
} /* namespace geos */
