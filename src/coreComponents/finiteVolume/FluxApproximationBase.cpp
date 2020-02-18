/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FluxApproximationBase.cpp
 *
 */

#include "FluxApproximationBase.hpp"

#include "managers/FieldSpecification/FieldSpecificationManager.hpp"

namespace geosx
{

using namespace dataRepository;

FluxApproximationBase::FluxApproximationBase(string const &name, Group *const parent)
  : Group(name, parent),
    m_fieldName(),
    m_boundaryFieldName(),
    m_coeffName()
{
  setInputFlags(InputFlags::OPTIONAL_NONUNIQUE);

  registerWrapper(viewKeyStruct::fieldNameString, &m_fieldName, false)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of primary solution field");

  registerWrapper(viewKeyStruct::boundaryFieldNameString, &m_boundaryFieldName, false)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Name of boundary (face) field");

  registerWrapper(viewKeyStruct::coeffNameString, &m_coeffName, false)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of coefficient field");

  registerWrapper(viewKeyStruct::targetRegionsString, &m_targetRegions, false)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("List of regions to build the stencil for");

  registerWrapper(viewKeyStruct::areaRelativeToleranceString, &m_areaRelTol, false)->
    setInputFlag(InputFlags::OPTIONAL)->
    setApplyDefaultValue(1.0e-8)->
    setDescription("Relative tolerance for area calculations.");


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

  computeCellStencil( domain );

  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();

  fsManager.Apply( 0.0,
                   const_cast<DomainPartition *>( &domain ), // hack, but guaranteed we won't modify it
                   "faceManager",
                   m_boundaryFieldName,
                   [&] ( FieldSpecificationBase const * GEOSX_UNUSED_PARAM( bc ),
                         string const & setName,
                         SortedArray<localIndex> const & targetSet,
                         Group const * GEOSX_UNUSED_PARAM( targetGroup ),
                         string const & GEOSX_UNUSED_PARAM( targetName ))
  {
    Wrapper<BoundaryStencil> * stencil = this->registerWrapper<BoundaryStencil>( setName );
    stencil->setRestartFlags(RestartFlags::NO_WRITE);
    computeBoundaryStencil( domain, targetSet, stencil->reference() );
  });
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
  return this->hasWrapper( setName );
}

void FluxApproximationBase::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
{
  DomainPartition const * domain = rootGroup->GetGroup<DomainPartition>( keys::domain );
  compute( *domain );
}

} //namespace geosx
