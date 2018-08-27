/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file TribolCoupling.cpp
 */

#include "TribolCoupling.hpp"
#include "dataRepository/KeyNames.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/MeshBody.hpp"

namespace geosx
{

vista::View *TribolCoupling::s_tribolProblem = nullptr ;
vista::View *TribolCoupling::s_tribolDomain = nullptr ;
vista::View *TribolCoupling::s_slideWorldSourceNodes = nullptr ;
vista::View *TribolCoupling::s_slideWorldSourceFaces = nullptr ;

void TribolCoupling::Initialize(dataRepository::ManagedGroup * eventManager, dataRepository::ManagedGroup * domain)
{
  real64& dt = *(eventManager->getData<real64>("dt"));
  integer& cycle = *(eventManager->getData<integer>("cycle"));

  // Currently we need the previous dt, but we assume fixed dt for now.
  // We probably want to pass in prevDt in the update anyway.
  s_tribolProblem = new vista::View("geosx", 0) ;
  s_tribolDomain = s_tribolProblem->viewCreate("DomainWorld", 0)->viewCreate("domain", 0, 0) ;

  MeshLevel * const meshLevel = domain->group_cast<DomainPartition*>()->getMeshBody(0)->getMeshLevel(0);
  CellBlockSubRegion * const subRegion = meshLevel->getElemManager()->GetRegion(0)->GetSubRegion(0);
  NodeManager * const nodeManager = meshLevel->getNodeManager();
  FaceManager * const faceManager = meshLevel->getFaceManager();

  int numFaces = faceManager->size() ;
  int numNodes = nodeManager->size() ;
  int numBricks = subRegion->size() ;
  const globalIndex_array faceMap = faceManager->m_localToGlobalMap ;
  const globalIndex_array nodeMap = nodeManager->m_localToGlobalMap ;
  const globalIndex_array brickMap = subRegion->m_localToGlobalMap ;

  const array2d<localIndex> & facesToElems = faceManager->elementSubRegionList();
  const OrderedVariableOneToManyRelation & facesToNodes = faceManager->nodeList();
  const integer_array isExternalFace = faceManager->m_isExternal ;

  globalID *extFaceMap ;
  int *extFacesToNodes ;
  int *extFacesToElement ;
  vista::MemAlloc(numFaces, &extFaceMap) ;
  vista::MemAlloc(numFaces*4, &extFacesToNodes) ;
  vista::MemAlloc(numFaces*2, &extFacesToElement) ;

  int numExtFaces = 0 ;

  for (int i = 0 ; i < numFaces ; ++i) {
     if (isExternalFace[i]) {
        extFaceMap[numExtFaces] = globalID((GIDTYPE)faceMap[i]) ;
        for (int j = 0 ; j < 4 ; ++j) {
           extFacesToNodes[numExtFaces*4+j] = facesToNodes[i][j] ;
        }
        for (int j = 0 ; j < 2 ; ++j) {
           extFacesToElement[numExtFaces*2+j] = facesToElems[i][j] ;
        }
        ++numExtFaces ;
     }
  }

  vista::MemRealloc(numExtFaces, &extFaceMap) ;
  vista::MemRealloc(numExtFaces*4, &extFacesToNodes) ;
  vista::MemRealloc(numExtFaces*2, &extFacesToElement) ;

  vista::View *srcFaces = s_tribolDomain->viewCreate("extFaces", 1, &numExtFaces, new vista::IndexSet(vista::VISTA_ACQUIRES, numExtFaces, extFaceMap)) ;
  vista::View *srcNodes = s_tribolDomain->viewCreate("nodes", 1, &numNodes, new vista::IndexSet(vista::VISTA_COPIES, numNodes, (globalID*) nodeMap.data())) ;
  vista::View *srcBricks = s_tribolDomain->viewCreate("bricks", 1, &numBricks, new vista::IndexSet(vista::VISTA_COPIES, numBricks, (globalID*) brickMap.data())) ;

  // include all nodes for now, this could be just the external nodes.
  int numSlideNodes = numNodes ;
  s_slideWorldSourceNodes = srcNodes->attach(new vista::View("LSslaveNodes0", 1, &numSlideNodes, new vista::IndexSet(numSlideNodes))) ;

  s_slideWorldSourceFaces = srcFaces->attach(new vista::View("LSslaveFaces0", 1, &numExtFaces, new vista::IndexSet(numExtFaces))) ;

  srcFaces->relationCreateFixed("facesToNodes", 4, extFacesToNodes, srcNodes, vista::VISTA_ACQUIRES) ;
  srcFaces->relationCreateFixed("facesToElement", 2, extFacesToElement, srcBricks, vista::VISTA_ACQUIRES) ;

  const array1d< R1Tensor > & X = nodeManager->referencePosition();
  
  real64 *x0 = srcNodes->fieldCreateReal("x0")->Real() ;
  real64 *y0 = srcNodes->fieldCreateReal("y0")->Real() ;
  real64 *z0 = srcNodes->fieldCreateReal("z0")->Real() ;

  const int *slideMap = s_slideWorldSourceNodes->map() ;

  for (int i = 0 ; i < numNodes ; ++i) {
     x0[i] = X[i][0] ;
     y0[i] = X[i][1] ;
     z0[i] = X[i][2] ;
  }

  s_slideWorldSourceNodes->fieldCreateReal("x") ;
  s_slideWorldSourceNodes->fieldCreateReal("y") ;
  s_slideWorldSourceNodes->fieldCreateReal("z") ;
  s_slideWorldSourceNodes->fieldCreateReal("fx") ;
  s_slideWorldSourceNodes->fieldCreateReal("fy") ;
  s_slideWorldSourceNodes->fieldCreateReal("fz") ;

  CopyPositionsToTribolSourceData(nodeManager) ;
  CopyForcesToTribolSourceData(nodeManager) ;

  srcBricks->relationCreateFixed("bricksToNodes", 8, (int*)subRegion->nodeList().data(), srcNodes, vista::VISTA_SHARES) ;

  SlideWorldAdapter::CreateWorld(MPI_COMM_GEOSX, MPI_COMM_WORLD,
                                 4000, // comm tag 
                                 3, // dimension
                                 0, // axisym
                                 1, // numSS
                                 4, // nodesPerFace
                                 8, // nodesPerElem
                                 &cycle,
                                 &dt,
                                 &dt, // prevDt
                                 1, // noParamInput
                                 s_tribolProblem) ;

  SlideWorldAdapter::SetSourceData(s_tribolProblem) ;
}

void TribolCoupling::CopyPositionsToTribolSourceData(NodeManager const * const nodeManager)
{
  const array1d< R1Tensor > & X = nodeManager->referencePosition();
  const array1d< R1Tensor > & u = nodeManager->getReference< array1d< R1Tensor > >("TotalDisplacement");
  const int numSlideNodes = s_slideWorldSourceNodes->length() ;
  const int *slideMap = s_slideWorldSourceNodes->map() ;
  real64 *x = s_slideWorldSourceNodes->fieldReal("x") ;
  real64 *y = s_slideWorldSourceNodes->fieldReal("y") ;
  real64 *z = s_slideWorldSourceNodes->fieldReal("z") ;

  for (int i = 0 ; i < numSlideNodes ; ++i) {
     int nodeIdx = slideMap[i] ;
     x[i] = X[nodeIdx][0] + u[nodeIdx][0] ;
     y[i] = X[nodeIdx][1] + u[nodeIdx][1] ;
     z[i] = X[nodeIdx][2] + u[nodeIdx][2] ;
  }
}

void TribolCoupling::CopyForcesToTribolSourceData(NodeManager const * const nodeManager)
{
  const array1d< R1Tensor > & a = nodeManager->getReference< array1d< R1Tensor > >( dataRepository::keys::Acceleration );
  const array1d< real64 > & nodalMass = nodeManager->getReference< array1d< real64 > >( dataRepository::keys::Mass );
  const int numSlideNodes = s_slideWorldSourceNodes->length() ;
  const int *slideMap = s_slideWorldSourceNodes->map() ;
  real64 *fx = s_slideWorldSourceNodes->fieldReal("fx") ;
  real64 *fy = s_slideWorldSourceNodes->fieldReal("fy") ;
  real64 *fz = s_slideWorldSourceNodes->fieldReal("fz") ;

  for (int i = 0 ; i < numSlideNodes ; ++i) {
     int nodeIdx = slideMap[i] ;
     real64 nmass = nodalMass[nodeIdx] ;
     fx[i] = a[nodeIdx][0]*nmass ;
     fy[i] = a[nodeIdx][1]*nmass ;
     fz[i] = a[nodeIdx][2]*nmass ;
  }
}

void TribolCoupling::CopyAccelerationsFromTribolSourceData(NodeManager * const nodeManager)
{
  array1d< R1Tensor > & a = nodeManager->getReference< array1d< R1Tensor > >( dataRepository::keys::Acceleration );
  const array1d< real64 > & nodalMass = nodeManager->getReference< array1d< real64 > >( dataRepository::keys::Mass );
  const int numSlideNodes = s_slideWorldSourceNodes->length() ;
  const int *slideMap = s_slideWorldSourceNodes->map() ;
  const real64 *fx = s_slideWorldSourceNodes->fieldReal("fx") ;
  const real64 *fy = s_slideWorldSourceNodes->fieldReal("fy") ;
  const real64 *fz = s_slideWorldSourceNodes->fieldReal("fz") ;

  for (int i = 0 ; i < numSlideNodes ; ++i) {
     int nodeIdx = slideMap[i] ;
     real64 nmass = nodalMass[nodeIdx] ;
     a[nodeIdx][0] = fx[i]/nmass ;
     a[nodeIdx][1] = fy[i]/nmass ;
     a[nodeIdx][2] = fz[i]/nmass ;
  }
}

void TribolCoupling::SyncTermination(int* terminate)
{
   SlideWorldAdapter::SyncTermination(terminate) ;
}

void TribolCoupling::SyncTimestep(real64* newDt)
{
    // Currently does not change GEOS timestep.
    // Will also be incorrect if timestep is adjusted by maxTime.
    SlideWorldAdapter::GetRequiredTimestep(newDt) ;

    // Ignore GEOS requests
    int plotThisCycle = 0 ;
    int dumpThisCycle = 0 ;
    int advectCountThisCycle = 0 ;
    int slideDeleteThisCycle = 0 ;

    SlideWorldAdapter::SyncControl(&slideDeleteThisCycle, &advectCountThisCycle, &plotThisCycle, &dumpThisCycle) ;

    if (slideDeleteThisCycle) {
        SlideWorldAdapter::Cleanup();
        SlideWorldAdapter::SetSourceData(s_tribolProblem, true);
    }
}

void TribolCoupling::Cleanup()
{
   SlideWorldAdapter::DestroyWorld() ;
   delete s_tribolProblem ;
   s_tribolProblem = nullptr ;
   s_tribolDomain = nullptr ;
   s_slideWorldSourceNodes = nullptr ;
   s_slideWorldSourceFaces = nullptr ;
}

} /* namespace geosx */
