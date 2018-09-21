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

#include "constitutive/ConstitutiveManager.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "math/TensorT/TensorBaseT.h"
#include "mesh/MeshForLoopInterface.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

SimpleWell::SimpleWell(string const & name, dataRepository::ManagedGroup * const parent)
  : WellBase(name, parent)
{

}

SimpleWell::~SimpleWell()
{

}

void SimpleWell::FillDocumentationNode()
{
  WellBase::FillDocumentationNode();

  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName("SimpleWell");
  docNode->setSchemaType("Node");

  docNode->AllocateChildNode( viewKeysSimpleWell.pressure.Key(),
                              viewKeysSimpleWell.pressure.Key(),
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

  docNode->AllocateChildNode( viewKeysSimpleWell.gravityDepth.Key(),
                              viewKeysSimpleWell.gravityDepth.Key(),
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

const string SimpleWell::getCatalogName() const
{
  return CatalogName();
}

void SimpleWell::ConnectToCells( DomainPartition const * domain )
{
  MeshLevel const * mesh = domain->getMeshBody(0)->getMeshLevel(0);
  ElementRegionManager const * elemManager = mesh->getElemManager();

  auto elemCenter = elemManager->ConstructViewAccessor<array1d<R1Tensor>>( CellBlockSubRegion::
                                                                           viewKeyStruct::
                                                                           elementCenterString );

  PerforationManager const * perfManager = GetGroup<PerforationManager>( groupKeysSimpleWell.perforations );

  // TODO Until we can properly trace perforations to cells,
  // just connect to the nearest cell center (this is NOT correct in general)

  m_numConnections = 0;
  localIndex iconn_global = 0;
  perfManager->forSubGroups<Perforation>( [&] ( Perforation const * perf ) -> void
  {
    R1Tensor const & loc = perf->getLocation();

    auto ret = minLocOverElemsInMesh( mesh, [&] ( localIndex er,
                                                  localIndex esr,
                                                  localIndex ei ) -> real64
    {
      R1Tensor v = loc;
      v -= elemCenter[er][esr][ei];
      return v.L2_Norm();
    });

    m_connectionElementRegion[m_numConnections]    = std::get<0>(ret.second);
    m_connectionElementSubregion[m_numConnections] = std::get<1>(ret.second);
    m_connectionElementIndex[m_numConnections]     = std::get<2>(ret.second);
    m_connectionPerforationIndex[m_numConnections] = iconn_global++;

    // This will not be correct in parallel until we can actually check that
    // the perforation belongs to local mesh partition
    ++m_numConnections;
  });
}

void SimpleWell::FinalInitialization(ManagedGroup * const problemManager)
{
  WellBase::FinalInitialization(problemManager);

  DomainPartition const * domain = problemManager->GetGroup<DomainPartition>( keys::domain );

  // initially allocate enough memory for all (global) perforations
  resize( numConnectionsGlobal() );
  ConnectToCells( domain );
  resize( numConnectionsLocal() );

  // generate the "all" set to enable application of BC
  ManagedGroup * sets = GetGroup( keys::sets );
  set<localIndex> & setAll = sets->RegisterViewWrapper<set<localIndex>>("all")->reference();
  for (localIndex iconn = 0; iconn < numConnectionsLocal(); ++iconn)
  {
    setAll.insert( iconn );
  }

  PrecomputeData( domain );

  NumericalMethodsManager const * numericalMethodManager = domain->
    getParent()->GetGroup<NumericalMethodsManager>( keys::numericalMethodsManager );

  FiniteVolumeManager const * fvManager = numericalMethodManager->
    GetGroup<FiniteVolumeManager>( keys::finiteVolumeManager );

  FluxApproximationBase const * fluxApprox = fvManager->GetGroup<FluxApproximationBase>(0); // TODO HACK!

  auto vw = this->RegisterViewWrapper<FluxApproximationBase::WellStencil>( keys::FVstencil );
  vw->setRestartFlags( RestartFlags::NO_WRITE );
  FluxApproximationBase::WellStencil & stencil = vw->reference();
  stencil.reserve( numConnectionsLocal(), 2 );

  fluxApprox->computeWellStencil( domain, this, stencil );
}

void SimpleWell::PrecomputeData(DomainPartition const * domain)
{
  R1Tensor const & gravity = getGravityVector();
  array1d<real64> & gravDepth = getReference<array1d<real64>>(viewKeysSimpleWell.gravityDepth);

  PerforationManager const * perfManager = GetGroup<PerforationManager>( groupKeysSimpleWell.perforations );

  for (localIndex iconn = 0; iconn < numConnectionsLocal(); ++iconn)
  {
    Perforation const * perf = perfManager->GetGroup<Perforation>(iconn);
    gravDepth[iconn] = Dot( perf->getLocation(), gravity );
  }
}

void SimpleWell::UpdateConnectionPressure( DomainPartition const * domain, localIndex fluidIndex, bool gravityFlag )
{
  if ( !gravityFlag )
    return;

  MeshLevel const * mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager const * elemManager = mesh->getElemManager();

  ConstitutiveManager const * constitutiveManager =
    domain->GetGroup<ConstitutiveManager>( keys::ConstitutiveManager );

  auto & pres      = getReference<array1d<real64>>( viewKeysSimpleWell.pressure );
  auto & gravDepth = getReference<array1d<real64>>( viewKeysSimpleWell.gravityDepth );

  auto constitutiveRelations = elemManager->ConstructConstitutiveAccessor<ConstitutiveBase const>( constitutiveManager );

  R1Tensor const & gravity = getGravityVector();
  R1Tensor const & refDepth = { 0, 0, m_referenceDepth };
  real64 const refGravDepth = Dot( refDepth, gravity );

  // ECLIPSE 100/300 treatment: use density at BHP to compute hydrostatic head
  real64 dens, dummy;
  for (localIndex iconn = 0; iconn < numConnectionsLocal(); ++iconn)
  {
    constitutiveRelations[m_connectionElementRegion[iconn]]
                         [m_connectionElementSubregion[iconn]]
                         [fluidIndex]->FluidDensityCompute( pres[iconn], m_connectionElementIndex[iconn], dens, dummy );

    pres[iconn] += dens * (gravDepth[iconn] - refGravDepth);
  }
}

REGISTER_CATALOG_ENTRY( WellBase, SimpleWell, string const &, ManagedGroup * const )

} //namespace geosx
