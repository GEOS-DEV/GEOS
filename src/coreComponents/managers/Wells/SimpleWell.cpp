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

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

SimpleWell::SimpleWell(string const & name, dataRepository::ManagedGroup * const parent)
  : WellBase( name, parent ),
    m_bhp()
{
  RegisterViewWrapper( viewKeysSimpleWell.bhp.Key(), &m_bhp, false );
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

  cxx_utilities::DocumentationNode * const perfDocNode = m_perfManager.getDocumentationNode();

  perfDocNode->AllocateChildNode( viewKeysSimpleWell.pressure.Key(),
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

  perfDocNode->AllocateChildNode( viewKeysSimpleWell.flowRate.Key(),
                                  viewKeysSimpleWell.flowRate.Key(),
                                  -1,
                                  "real64_array",
                                  "real64_array",
                                  "Connection flow rate",
                                  "Connection flow rate",
                                  "",
                                  getName(),
                                  1,
                                  0,
                                  1 );

}

void SimpleWell::InitializationOrder(string_array & order)
{
  // Skip initializing PerforationManager, call it manually from FinalInitialization()
}

void SimpleWell::FinalInitializationPreSubGroups(ManagedGroup * const problemManager)
{
  WellBase::FinalInitializationPreSubGroups(problemManager);

  // fields owned by the well itself are scalar, but must be arrays for BC, so resize to 1
  resize(1);

  // generate the "all" set to enable application of BC
  ManagedGroup * sets = GetGroup( keys::sets );
  set<localIndex> & setAll = sets->RegisterViewWrapper<set<localIndex>>("all")->reference();
  setAll.insert(0);

  DomainPartition const * domain = problemManager->GetGroup<DomainPartition>( keys::domain );

  NumericalMethodsManager * numericalMethodManager =
    problemManager->GetGroup<NumericalMethodsManager>( keys::numericalMethodsManager );

  FiniteVolumeManager * fvManager =
    numericalMethodManager->GetGroup<FiniteVolumeManager>( keys::finiteVolumeManager );

  FluxApproximationBase * fluxApprox = fvManager->GetGroup<FluxApproximationBase>(0); // TODO HACK!

  auto vw = fluxApprox->GetGroup( fluxApprox->groupKeysFABase.wellStencils )
    ->RegisterViewWrapper<FluxApproximationBase::WellStencil>( getName() );

  vw->setRestartFlags( RestartFlags::NO_WRITE );
  FluxApproximationBase::WellStencil & stencil = vw->reference();
  stencil.reserve( numConnectionsLocal(), 2 );

  fluxApprox->computeWellStencil( domain, this, stencil );
}

void SimpleWell::StateUpdate( DomainPartition const * domain, localIndex fluidIndex )
{
  bool isGravityOn = getParent()->group_cast<WellManager *>()->getGravityFlag();

  MeshLevel const * mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager const * elemManager = mesh->getElemManager();

  ConstitutiveManager const * constitutiveManager =
    domain->GetGroup<ConstitutiveManager>( keys::ConstitutiveManager );

  auto & pres      = m_perfManager.getReference<array1d<real64>>( viewKeysSimpleWell.pressure );
  auto & gravDepth = m_perfManager.getReference<array1d<real64>>( m_perfManager.viewKeysPerfManager.gravityDepth );

  auto & elemRegion    = m_perfManager.getReference<array1d<localIndex>>( m_perfManager.viewKeysPerfManager.connectionElementRegion );
  auto & elemSubregion = m_perfManager.getReference<array1d<localIndex>>( m_perfManager.viewKeysPerfManager.connectionElementSubregion );
  auto & elemIndex     = m_perfManager.getReference<array1d<localIndex>>( m_perfManager.viewKeysPerfManager.connectionElementIndex );

  auto constitutiveRelations = elemManager->ConstructConstitutiveAccessor<ConstitutiveBase const>( constitutiveManager );

  R1Tensor const & gravity = getGravityVector();
  R1Tensor const & refDepth = { 0, 0, m_referenceDepth };
  real64 const refGravDepth = Dot( refDepth, gravity );

  // ECLIPSE 100/300 treatment: use density at BHP to compute hydrostatic head
  real64 dens, dummy;
  for (localIndex iconn = 0; iconn < numConnectionsLocal(); ++iconn)
  {
    constitutiveRelations[elemRegion[iconn]][elemSubregion[iconn]][fluidIndex]->FluidDensityCompute( pres[iconn],
                                                                                                     elemIndex[iconn],
                                                                                                     dens, dummy );

    pres[iconn] = pres[iconn] = m_bhp[0] + (isGravityOn ? dens * (gravDepth[iconn] - refGravDepth) : 0.0);
  }
}

real64 SimpleWell::GetTotalFlowRate()
{
  auto const & flowRate = m_perfManager.getReference<array1d<real64>>( viewKeysSimpleWell.flowRate );

  real64 totalRate = 0.0;
  for (localIndex iconn = 0; iconn < numConnectionsLocal(); ++iconn)
  {
    totalRate += flowRate[iconn];
  }

  return totalRate;
}

REGISTER_CATALOG_ENTRY( WellBase, SimpleWell, string const &, ManagedGroup * const )

} //namespace geosx
