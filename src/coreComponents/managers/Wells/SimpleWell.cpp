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
 * @file SimpleWell.cpp
 *
 */

#include "SimpleWell.hpp"
#include "Perforation.hpp"

#include "managers/DomainPartition.hpp"
#include "math/TensorT/TensorBaseT.h"

namespace geosx
{

using namespace dataRepository;

SimpleWell::SimpleWell(string const & name, dataRepository::ManagedGroup * const parent)
  : WellBase(name, parent)
{
  RegisterViewWrapper( viewKeys.connectionElementRegion.Key(), &m_connectionElementRegion, false );
  RegisterViewWrapper( viewKeys.connectionElementSubregion.Key(), &m_connectionElementSubregion, false );
  RegisterViewWrapper( viewKeys.connectionElementIndex.Key(), &m_connectionElementIndex, false );
}

SimpleWell::~SimpleWell()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->AllocateChildNode( viewKeys.pressure.Key(),
                              viewKeys.pressure.Key(),
                              -1,
                              "real64_array",
                              "real64_array",
                              "Connection pressure",
                              "Connection pressure",
                              "",
                              getName(),
                              1,
                              0,
                              1 );

  docNode->AllocateChildNode( viewKeys.transmissibility.Key(),
                              viewKeys.transmissibility.Key(),
                              -1,
                              "real64_array",
                              "real64_array",
                              "Connection transmissibility",
                              "Connection transmissibility",
                              "",
                              getName(),
                              1,
                              0,
                              1 );

  docNode->AllocateChildNode( viewKeys.gravityDepth.Key(),
                              viewKeys.gravityDepth.Key(),
                              -1,
                              "real64_array",
                              "real64_array",
                              "Connection gravity-depth product",
                              "Connection gravity-depth product",
                              "",
                              getName(),
                              1,
                              0,
                              1 );
}

void SimpleWell::FillDocumentationNode()
{
  WellBase::FillDocumentationNode();
}

const string SimpleWell::getCatalogName() const
{
  return CatalogName();
}

void SimpleWell::InitializePostSubGroups(ManagedGroup * const problemManager)
{
  WellBase::InitializePostSubGroups(problemManager);

  // all array fields live on connections/perforations
  this->resize(numConnections());

  DomainPartition * domain = problemManager->GetGroup<DomainPartition>( keys::domain );
  ConnectToCells( domain );
}

void SimpleWell::ConnectToCells( DomainPartition const * domain )
{
  ElementRegionManager const * elemManager = domain->getMeshBody(0)->getMeshLevel(0)->getElemManager();

  for (localIndex iconn = 0; iconn < numConnections(); ++iconn)
  {
    Perforation const * perf = GetGroup<Perforation>(iconn);
    R1Tensor const & loc = perf->getLocation();

    elemManager->forCellBlocksComplete( [&] (localIndex er,
                                             localIndex esr,
                                             ElementRegion const * region,
                                             CellBlockSubRegion const * subRegion) -> void
    {
      // TODO locate cell
//    m_connectionElementRegion[iconn] = ...
//    m_connectionElementSubregion[iconn] = ...
//    m_connectionElementIndex[iconn] = ...
    });
  }
}

void SimpleWell::FinalInitialization(ManagedGroup * const problemManager)
{
  WellBase::FinalInitialization(problemManager);

  DomainPartition const * domain = problemManager->GetGroup<DomainPartition>( keys::domain );
  PrecomputeData( domain );
}

void SimpleWell::PrecomputeData(DomainPartition const * domain)
{
  R1Tensor const & gravity = getGravityVector();
  array1d<real64> gravDepth = getReference<array1d<real64>>(viewKeys.gravityDepth);

  for (localIndex iconn = 0; iconn < numConnections(); ++iconn)
  {
    Perforation const * perf = GetGroup<Perforation>(iconn);
    gravDepth[iconn] = Dot(perf->getLocation(), gravity);
  }
}

REGISTER_CATALOG_ENTRY( WellBase, SimpleWell, string const &, ManagedGroup * const )

} //namespace geosx