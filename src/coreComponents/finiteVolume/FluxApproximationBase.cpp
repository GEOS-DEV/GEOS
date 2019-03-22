/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

#include "managers/FieldSpecification/FieldSpecificationManager.hpp"

namespace geosx
{

using namespace dataRepository;

FluxApproximationBase::FluxApproximationBase(string const &name, ManagedGroup *const parent)
  : ManagedGroup(name, parent),
    m_boundarySetData(nullptr),
    m_fieldName(),
    m_boundaryFieldName(),
    m_coeffName()
{
  setInputFlags(InputFlags::OPTIONAL_NONUNIQUE);

  m_boundarySetData = this->RegisterGroup(groupKeyStruct::boundarySetDataString);

  RegisterViewWrapper(viewKeyStruct::fieldNameString, &m_fieldName, false)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of primary solution field");

  RegisterViewWrapper(viewKeyStruct::boundaryFieldNameString, &m_boundaryFieldName, false)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Name of boundary (face) field");

  RegisterViewWrapper(viewKeyStruct::fratureRegionNameString, &m_fractureRegionName, false)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Names of the fracture region that will have a fracture stencil generated for them.");

  RegisterViewWrapper(viewKeyStruct::coeffNameString, &m_coeffName, false)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of coefficient field");

  RegisterViewWrapper(viewKeyStruct::areaRelativeToleranceString, &m_areaRelTol, false)->
    setInputFlag(InputFlags::OPTIONAL)->
    setApplyDefaultValue(1.0e-8)->
    setDescription("Relative tolerance for area calculations.");

  RegisterViewWrapper<CellStencil>(viewKeyStruct::cellStencilString)->
    setRestartFlags(RestartFlags::NO_WRITE);

  RegisterViewWrapper<CellStencil>(viewKeyStruct::fratureStencilString)->
    setRestartFlags(RestartFlags::NO_WRITE);

}

FluxApproximationBase::CatalogInterface::CatalogType &
FluxApproximationBase::GetCatalog()
{
  static FluxApproximationBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void FluxApproximationBase::compute(DomainPartition * domain)
{
  GEOSX_MARK_BEGIN("FluxApproximationBase::compute");

  computeMainStencil(domain, getStencil());

  if( !m_fractureRegionName.empty() )
  {
    computeFractureStencil( *domain,
                            this->getReference<CellStencil>(viewKeyStruct::fratureStencilString),
                            getStencil() );
  }

  FieldSpecificationManager * fsManager = FieldSpecificationManager::get();

  fsManager->Apply( 0.0,
                    domain,
                    "faceManager",
                    m_boundaryFieldName,
                    [&] ( FieldSpecificationBase const * bc,
                    string const & setName,
                    set<localIndex> const & targetSet,
                    ManagedGroup * targetGroup,
                    string const & targetName) -> void
  {
    ViewWrapper<BoundaryStencil> * stencil = this->RegisterViewWrapper<BoundaryStencil>(setName);
    stencil->setRestartFlags(RestartFlags::NO_WRITE);
    computeBoundaryStencil(domain, targetSet, stencil->reference());
  });

  GEOSX_MARK_END("FluxApproximationBase::compute");
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

void FluxApproximationBase::InitializePostInitialConditions_PreSubGroups(ManagedGroup * const rootGroup)
{
  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);
  compute(domain);
}

} //namespace geosx
