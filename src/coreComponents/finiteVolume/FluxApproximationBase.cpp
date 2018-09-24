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
  this->RegisterGroup( groupKeysFABase.faceStencils );
  this->RegisterGroup( groupKeysFABase.wellStencils );

  this->RegisterViewWrapper( viewKeysFABase.fieldName.Key(), &m_fieldName, false );
  this->RegisterViewWrapper( viewKeysFABase.boundaryFieldName.Key(), &m_boundaryFieldName, false );
  this->RegisterViewWrapper( viewKeysFABase.coeffName.Key(), &m_coeffName, false );

  this->RegisterViewWrapper<CellStencil>( viewKeysFABase.cellStencil )->setRestartFlags( RestartFlags::NO_WRITE );
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

void FluxApproximationBase::compute( DomainPartition * domain )
{
  // compute cell-cell stencil in the domain
  computeCellStencil( domain, getCellStencil() );

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
    ManagedGroup * faceStencils = this->GetGroup( groupKeysFABase.faceStencils );
    ViewWrapper<FaceStencil> * stencil = faceStencils->RegisterViewWrapper<FaceStencil>( setName );
    stencil->setRestartFlags( RestartFlags::NO_WRITE );
    computeFaceStencil( domain, targetSet, stencil->reference() );
  });
}

void FluxApproximationBase::FinalInitialization( ManagedGroup * const rootGroup )
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
