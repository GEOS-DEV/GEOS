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
#include "constitutive/Fluid/SingleFluidBase.hpp"
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
  RegisterViewWrapper( viewKeyStruct::bhpString, &m_bhp, false );

  // Most perforations-based fields are registered by specific well models, rather than
  // PerforationManager itself, which is physics-agnostic
  PerforationManager * perfManager = GetGroup<PerforationManager>( groupKeyStruct::perforationsString );
  perfManager->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::pressureString );
  perfManager->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::flowRateString );
}

SimpleWell::~SimpleWell()
{

}

void SimpleWell::InitializePostSubGroups( ManagedGroup * const rootGroup )
{
  WellBase::InitializePostSubGroups( rootGroup );

  // vars owned by the well itself are scalar, but must be stored as arrays for BC, so resize to 1
  resize(1);

  // generate the "all" set to enable application of BC
  set<localIndex> & setAll = this->sets()->RegisterViewWrapper<set<localIndex>>("all")->reference();
  setAll.insert(0);
}

void SimpleWell::StateUpdate( DomainPartition const * domain, localIndex fluidIndex )
{
  bool isGravityOn = getParent()->group_cast<WellManager *>()->getGravityFlag();

  MeshLevel const * mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager const * elemManager = mesh->getElemManager();

  ConstitutiveManager const * constitutiveManager = domain->getConstitutiveManager();

  arrayView1d<real64> const & pres =
    m_perfManager.getReference<array1d<real64>>( viewKeyStruct::pressureString );

  arrayView1d<real64 const> const & gravDepth =
    m_perfManager.getReference<array1d<real64>>( PerforationManager::viewKeyStruct::gravityDepthString );

  arrayView1d<localIndex const> const & elemRegion =
    m_perfManager.getReference<array1d<localIndex>>( PerforationManager::viewKeyStruct::connectionElementRegionString );

  arrayView1d<localIndex const> const & elemSubregion =
    m_perfManager.getReference<array1d<localIndex>>( PerforationManager::viewKeyStruct::connectionElementSubregionString );

  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase const> constitutiveRelations =
    elemManager->ConstructConstitutiveAccessor<ConstitutiveBase const>( constitutiveManager );

  R1Tensor const & gravity = getGravityVector();
  R1Tensor const & refDepth = { 0, 0, m_referenceDepth };
  real64 const refGravDepth = Dot( refDepth, gravity );

  // ECLIPSE 100/300 treatment: use density at BHP to compute hydrostatic head
  real64 dens, dummy;
  for (localIndex iconn = 0; iconn < numConnectionsLocal(); ++iconn)
  {
    localIndex const er  = elemRegion[iconn];
    localIndex const esr = elemSubregion[iconn];

    SingleFluidBase const * fluid = constitutiveRelations[er][esr][fluidIndex]->group_cast<SingleFluidBase const *>();
    fluid->Compute( m_bhp[0], dens, dummy, dummy, dummy );

    pres[iconn] = m_bhp[0] + (isGravityOn ? dens * (gravDepth[iconn] - refGravDepth) : 0.0);
  }
}

real64 SimpleWell::GetTotalFlowRate()
{
  arrayView1d<real64 const> const & flowRate =
    m_perfManager.getReference<array1d<real64>>( viewKeyStruct::flowRateString );

  real64 totalRate = 0.0;
  for (localIndex iconn = 0; iconn < numConnectionsLocal(); ++iconn)
  {
    totalRate += flowRate[iconn];
  }

  return totalRate;
}

REGISTER_CATALOG_ENTRY( WellBase, SimpleWell, string const &, ManagedGroup * const )

} //namespace geosx
