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
 * @file FluxApproximationBase.cpp
 *
 */

#include "FluxApproximationBase.hpp"

#include "managers/BoundaryConditions/BoundaryConditionManager.hpp"

namespace geosx
{

using namespace dataRepository;

FluxApproximationBase::FluxApproximationBase(string const &name, ManagedGroup *const parent)
  : ManagedGroup(name, parent)
{
  this->RegisterViewWrapper(viewKeyStruct::fieldNameString, &m_fieldName, false);
  this->RegisterViewWrapper(viewKeyStruct::boundaryFieldNameString, &m_boundaryFieldName, false);
  this->RegisterViewWrapper(viewKeyStruct::coeffNameString, &m_coeffName, false);

  ViewWrapper<CellStencil> * stencil = this->RegisterViewWrapper<CellStencil>(viewKeyStruct::cellStencilString);
  stencil->setRestartFlags(RestartFlags::NO_WRITE);
}

FluxApproximationBase::CatalogInterface::CatalogType &
FluxApproximationBase::GetCatalog()
{
  static FluxApproximationBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void FluxApproximationBase::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->AllocateChildNode(viewKeyStruct::fieldNameString,
                             viewKeyStruct::fieldNameString,
                             -1,
                             "string",
                             "string",
                             "Name of primary solution field",
                             "Name of primary solution field",
                             "REQUIRED",
                             "",
                             0,
                             1,
                             0 );

  docNode->AllocateChildNode(viewKeyStruct::boundaryFieldNameString,
                             viewKeyStruct::boundaryFieldNameString,
                             -1,
                             "string",
                             "string",
                             "Name of boundary (face) field",
                             "Name of boundary (face) field",
                             "",
                             "",
                             0,
                             1,
                             0 );

  docNode->AllocateChildNode(viewKeyStruct::coeffNameString,
                             viewKeyStruct::coeffNameString,
                             -1,
                             "string",
                             "string",
                             "Name of coefficient field",
                             "Name of coefficient field",
                             "REQUIRED",
                             "",
                             0,
                             1,
                             0 );
}

void FluxApproximationBase::compute(DomainPartition * domain)
{
  // compute cell-cell stencil in the domain
  computeMainStencil(domain, getStencil());

  // compute face-cell stencils for boundary conditions
  BoundaryConditionManager * bcManager = BoundaryConditionManager::get();
  bcManager->ApplyBoundaryCondition( 0.0,
                                     domain,
                                     "faceManager",
                                     m_boundaryFieldName,
                                     [&] ( BoundaryConditionBase const * bc,
                                           string const & setName,
                                           set<localIndex> const & targetSet,
                                           ManagedGroup * targetGroup,
                                           string const & targetName) -> void
  {
    ViewWrapper<BoundaryStencil> * stencil = this->RegisterViewWrapper<BoundaryStencil>(setName);
    stencil->setRestartFlags(RestartFlags::NO_WRITE);
    computeBoundaryStencil(domain, targetSet, stencil->reference());
  });
}

FluxApproximationBase::CellStencil const &
FluxApproximationBase::getStencil() const
{
  return this->getReference<CellStencil>(viewKeyStruct::cellStencilString);
}

FluxApproximationBase::CellStencil &
FluxApproximationBase::getStencil()
{
  return this->getReference<CellStencil>(viewKeyStruct::cellStencilString);
}

FluxApproximationBase::BoundaryStencil const &
FluxApproximationBase::getBoundaryStencil(string const & setName) const
{
  return this->getReference<BoundaryStencil>(setName);
}

FluxApproximationBase::BoundaryStencil &
FluxApproximationBase::getBoundaryStencil(string const & setName)
{
  return this->getReference<BoundaryStencil>(setName);
}

bool FluxApproximationBase::hasBoundaryStencil(string const & setName) const
{
  return this->hasView(setName);
}

void FluxApproximationBase::FinalInitialization(ManagedGroup * const rootGroup)
{
  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);
  compute(domain);
}

} //namespace geosx
