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

#include "physicsSolvers/BoundaryConditions/BoundaryConditionManager.hpp"

namespace geosx
{

using namespace dataRepository;

FluxApproximationBase::FluxApproximationBase(string const &name, ManagedGroup *const parent)
  : ManagedGroup(name, parent)
{
  this->RegisterViewWrapper(viewKeyStruct::fieldNameString, &m_fieldName, false);
  this->RegisterViewWrapper(viewKeyStruct::coeffNameString, &m_coeffName, false);
  this->RegisterViewWrapper(viewKeyStruct::cellLocationString, &m_cellLocation, false);
  this->RegisterViewWrapper<CellStencil>(viewKeyStruct::cellStencilString);
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

  docNode->AllocateChildNode(viewKeyStruct::cellLocationString,
                             viewKeyStruct::cellLocationString,
                             -1,
                             "string",
                             "string",
                             "Option for cell collocation points",
                             "Option for cell collocation points",
                             CellBlock::viewKeyStruct::elementCenterString,
                             "",
                             0,
                             1,
                             0 );
}

void FluxApproximationBase::compute(DomainPartition * domain)
{
  computeMainStencil(domain, getStencil());

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  FaceManager * const faceManager = mesh->getFaceManager();
  BoundaryConditionManager * bcManager = BoundaryConditionManager::get();

  dataRepository::ManagedGroup const * sets = faceManager->GetGroup(dataRepository::keys::sets);

  bcManager->forSubGroups<BoundaryConditionBase>([&] (BoundaryConditionBase * bc) -> void
  {
    if (bc->initialCondition() || bc->GetFieldName() != m_fieldName)
      return;

    string_array setNames = bc->GetSetNames();
    for (auto & setName : setNames)
    {
      dataRepository::ViewWrapper<lSet> const * const setWrapper = sets->getWrapper<lSet>(setName);
      if (setWrapper != nullptr)
      {
        lSet const & set = setWrapper->reference();
        ViewWrapper<BoundaryStencil> * stencil = this->RegisterViewWrapper<BoundaryStencil>(setName);
        computeBoundaryStencil(domain, set, stencil->reference());
      }
    }
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

}