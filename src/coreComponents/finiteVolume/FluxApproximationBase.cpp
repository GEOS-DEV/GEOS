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
  this->RegisterGroup( groupKeysFABase.faceStencils );
  this->RegisterGroup( groupKeysFABase.wellStencils );

  RegisterViewWrapper(viewKeyStruct::fieldNameString, &m_fieldName, false)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of primary solution field");

  RegisterViewWrapper(viewKeyStruct::boundaryFieldNameString, &m_boundaryFieldName, false)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Name of boundary (face) field");

  RegisterViewWrapper(viewKeyStruct::coeffNameString, &m_coeffName, false)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of coefficient field");

  RegisterViewWrapper<CellStencil>(viewKeyStruct::cellStencilString)->
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
  // compute cell-cell stencil in the domain
  computeCellStencil( domain, getCellStencil() );

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
    ManagedGroup * faceStencils = this->GetGroup( groupKeysFABase.faceStencils );
    ViewWrapper<FaceStencil> * stencil = faceStencils->RegisterViewWrapper<FaceStencil>( setName );
    stencil->setRestartFlags( RestartFlags::NO_WRITE );
    computeFaceStencil( domain, targetSet, stencil->reference() );
  });
}

void FluxApproximationBase::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const rootGroup )
{
  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>( keys::domain );
  compute( domain );
}

FluxApproximationBase::CellStencil const &
FluxApproximationBase::getCellStencil() const
{
  return this->getReference<CellStencil>( viewKeyStruct::cellStencilString );
}

FluxApproximationBase::CellStencil &
FluxApproximationBase::getCellStencil()
{
  return this->getReference<CellStencil>( viewKeyStruct::cellStencilString );
}

FluxApproximationBase::FaceStencil const &
FluxApproximationBase::getFaceStencil( string const & setName ) const
{
  return this->GetGroup( groupKeysFABase.faceStencils )->getReference<FaceStencil>( setName );
}

FluxApproximationBase::FaceStencil &
FluxApproximationBase::getFaceStencil( string const & setName )
{
  return this->GetGroup( groupKeysFABase.faceStencils )->getReference<FaceStencil>( setName );
}

bool FluxApproximationBase::hasFaceStencil( string const & setName ) const
{
  return this->GetGroup( groupKeysFABase.faceStencils )->hasView( setName );
}

const FluxApproximationBase::WellStencil & FluxApproximationBase::getWellStencil( string const & wellName ) const
{
  return this->GetGroup( groupKeysFABase.wellStencils )->getReference<WellStencil>( wellName );
}

FluxApproximationBase::WellStencil & FluxApproximationBase::getWellStencil( string const & wellName )
{
  return this->GetGroup( groupKeysFABase.wellStencils )->getReference<WellStencil>( wellName );
}

bool FluxApproximationBase::hasWellStencil( string const & wellName ) const
{
  return this->GetGroup( groupKeysFABase.wellStencils )->hasView( wellName );
}

} //namespace geosx
