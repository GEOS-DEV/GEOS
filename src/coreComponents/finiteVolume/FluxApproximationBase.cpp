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
    m_fieldName(),
    m_boundaryFieldName(),
    m_coeffName()
{
  setInputFlags(InputFlags::OPTIONAL_NONUNIQUE);

  RegisterViewWrapper(viewKeyStruct::fieldNameString, &m_fieldName, false)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of primary solution field");

  RegisterViewWrapper(viewKeyStruct::boundaryFieldNameString, &m_boundaryFieldName, false)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Name of boundary (face) field");

  RegisterViewWrapper(viewKeyStruct::coeffNameString, &m_coeffName, false)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of coefficient field");

  RegisterViewWrapper(viewKeyStruct::targetRegionsString, &m_targetRegions, false)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("List of regions to build the stencil for");

  RegisterViewWrapper(viewKeyStruct::areaRelativeToleranceString, &m_areaRelTol, false)->
    setInputFlag(InputFlags::OPTIONAL)->
    setApplyDefaultValue(1.0e-8)->
    setDescription("Relative tolerance for area calculations.");

  RegisterViewWrapper<CellStencil>(viewKeyStruct::cellStencilString)->
    setRestartFlags(RestartFlags::NO_WRITE);

  RegisterViewWrapper<CellStencil>(viewKeyStruct::fractureStencilString)->
    setRestartFlags(RestartFlags::NO_WRITE);

}

FluxApproximationBase::CatalogInterface::CatalogType &
FluxApproximationBase::GetCatalog()
{
  static FluxApproximationBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void FluxApproximationBase::compute( DomainPartition const & domain )
{
  GEOSX_MARK_FUNCTION_SCOPED;

  computeCellStencil( domain, getStencil() );

//  computeFractureStencil( domain,
//                          this->getReference<CellStencil>(viewKeyStruct::fractureStencilString),
//                          getStencil() );

  FieldSpecificationManager * fsManager = FieldSpecificationManager::get();

  fsManager->Apply( 0.0,
                    const_cast<DomainPartition *>( &domain ), // hack, but guaranteed we won't modify it
                    "faceManager",
                    m_boundaryFieldName,
                    [&] ( FieldSpecificationBase const * bc,
                          string const & setName,
                          set<localIndex> const & targetSet,
                          ManagedGroup const * targetGroup,
                          string const & targetName)
  {
    ViewWrapper<BoundaryStencil> * stencil = this->RegisterViewWrapper<BoundaryStencil>( setName );
    stencil->setRestartFlags(RestartFlags::NO_WRITE);
    computeBoundaryStencil( domain, targetSet, stencil->reference() );
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

void FluxApproximationBase::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const rootGroup )
{
  DomainPartition const * domain = rootGroup->GetGroup<DomainPartition>( keys::domain );
  compute( *domain );
}

} //namespace geosx
