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

/*
 * @file FiniteVolumeManager.cpp
 *
 */

#include "FiniteVolumeManager.hpp"

#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"

namespace geosx
{

using namespace dataRepository;


FiniteVolumeManager::FiniteVolumeManager(string const &name, ManagedGroup *const parent)
  : ManagedGroup(name, parent)
{

}

FiniteVolumeManager::~FiniteVolumeManager()
{

}

void FiniteVolumeManager::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName("FiniteVolume");
  docNode->setSchemaType("UniqueNode");
}

void FiniteVolumeManager::CreateChild(string const &childKey, string const &childName)
{
  std::unique_ptr<FluxApproximationBase> approx = FluxApproximationBase::CatalogInterface::Factory(childKey, childName, this);
  FluxApproximationBase * newApprox = this->RegisterGroup<FluxApproximationBase>(childName, std::move(approx));
}

FluxApproximationBase const * FiniteVolumeManager::getFluxApproximation(std::string const &name) const
{
  return this->GetGroup<FluxApproximationBase>(name);
}

void FiniteVolumeManager::FinalInitialization(ManagedGroup *const rootGroup)
{
  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);
  precomputeFiniteVolumeData(domain);
}

void FiniteVolumeManager::precomputeFiniteVolumeData(DomainPartition * const domain)
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();
  FaceManager * const faceManager = mesh->getFaceManager();

  r1_array const & X = nodeManager->referencePosition();

  auto elemCenter = elemManager->ConstructViewAccessor<r1_array>(CellBlock::
                                                                 viewKeyStruct::
                                                                 elementCenterString);

  auto elemVolume = elemManager->ConstructViewAccessor<real64_array>(CellBlock::
                                                                     viewKeyStruct::
                                                                     elementVolumeString);

  auto elemsToNodes = elemManager->ConstructViewAccessor<FixedOneToManyRelation>(CellBlockSubRegion::
                                                                                 viewKeyStruct::
                                                                                 nodeListString);


  // Loop over all the elements and calculate element centers, and element volumes
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k )->void
  {
    localIndex const * const nodeList = elemsToNodes[er][esr][k];
    localIndex const nodeListSize = elemsToNodes[er][esr].get().size(1);
    R1Tensor Xlocal[ElementRegionManager::maxNumNodesPerElem];

    R1Tensor & center = elemCenter[er][esr][k];
    center = 0.0;

    // TODO different center options
    for (localIndex a = 0; a < nodeListSize; ++a)
    {
      Xlocal[a] = X[nodeList[a]];
      center += Xlocal[a];
    }
    center /= nodeListSize;

    // TODO proper volumes for all shapes
    if( nodeListSize == 8 )
    {
        elemVolume[er][esr][k] = computationalGeometry::HexVolume(Xlocal);
    }
    else if( nodeListSize == 4)
    {
        elemVolume[er][esr][k] = computationalGeometry::TetVolume(Xlocal);
    }
    else if( nodeListSize == 6)
    {
        elemVolume[er][esr][k] = computationalGeometry::WedgeVolume(Xlocal);
    }
    else if ( nodeListSize == 5)
    {
        elemVolume[er][esr][k] = computationalGeometry::PyramidVolume(Xlocal);
    }
    else
    {
        GEOS_ERROR("GEOX does not support cells with " << nodeListSize << " nodes");
    }
  });

  r1_array & faceCenter = faceManager->getReference<r1_array>(FaceManager::viewKeyStruct::faceCenterString);
  array1d<array1d<localIndex>> const & faceToNodes = faceManager->nodeList();

  R1Tensor normal;
  for (localIndex kf = 0; kf < faceManager->size(); ++kf)
  {
    computationalGeometry::Centroid_3DPolygon(faceToNodes[kf], X, faceCenter[kf], normal);
  }
}

} // namespace geosx
