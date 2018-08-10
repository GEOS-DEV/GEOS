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
  m_boundarySetData = this->RegisterGroup(groupKeyStruct::boundarySetDataString);

  this->RegisterViewWrapper(viewKeyStruct::fieldNameString, &m_fieldName, false);
  this->RegisterViewWrapper(viewKeyStruct::coeffNameString, &m_coeffName, false);
  this->RegisterViewWrapper(viewKeyStruct::cellLocationString, &m_cellLocation, false);

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
    if (bc->initialCondition())
      return;

    string_array const objectPath = stringutilities::Tokenize(bc->GetObjectPath(), "/");

    if (objectPath.size() < 1 || objectPath[0] != MeshLevel::groupStructKeys::faceManagerString)
      return;

    // TODO this validation should really be in BoundaryConditionBase::ReadXML_Postprocess
    string fieldName = bc->GetFieldName();
    if (objectPath.size() > 1)
    {
      GEOS_ASSERT(fieldName.empty() || objectPath[1].empty(),
                  "field name specified in both fieldName (" << fieldName
                                                             << ") and objectPath (" << objectPath[1] << ")");

      GEOS_ASSERT(!(bc->GetFieldName().empty() && objectPath[3].empty()),
                  "field name not specified in either fieldName or objectPath");

      if (!objectPath[1].empty())
        fieldName = objectPath[1];
    }
    else
    {
      GEOS_ASSERT(!fieldName.empty(), "field name not specified in either fieldName or objectPath");
    }

    if (fieldName != m_fieldName)
      return;

    string_array setNames = bc->GetSetNames();
    for (auto & setName : setNames)
    {
      dataRepository::ViewWrapper<lSet> const * const setWrapper = sets->getWrapper<lSet>(setName);
      if (setWrapper != nullptr)
      {
        lSet const & set = setWrapper->reference();
        ViewWrapper<BoundaryStencil> * stencil = this->RegisterViewWrapper<BoundaryStencil>(setName);
        stencil->setRestartFlags(RestartFlags::NO_WRITE);
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

bool FluxApproximationBase::hasBoundaryStencil(string const & setName) const
{
  return this->hasView(setName);
}

}