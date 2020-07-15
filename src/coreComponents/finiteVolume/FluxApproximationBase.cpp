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
#include "mpiCommunications/CommunicationTools.hpp"

namespace geosx
{

using namespace dataRepository;

FluxApproximationBase::FluxApproximationBase( string const & name, Group * const parent )
  : Group( name, parent ),
  m_fieldName(),
  m_boundaryFieldName(),
  m_coeffName()
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::fieldNameString, &m_fieldName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of primary solution field" );

  registerWrapper( viewKeyStruct::boundaryFieldNameString, &m_boundaryFieldName )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Name of boundary (face) field" );

  registerWrapper( viewKeyStruct::coeffNameString, &m_coeffName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of coefficient field" );

  registerWrapper( viewKeyStruct::targetRegionsString, &m_targetRegions )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "List of regions to build the stencil for" );

  registerWrapper( viewKeyStruct::areaRelativeToleranceString, &m_areaRelTol )->
    setInputFlag( InputFlags::OPTIONAL )->
    setApplyDefaultValue( 1.0e-8 )->
    setDescription( "Relative tolerance for area calculations." );


}

FluxApproximationBase::CatalogInterface::CatalogType &
FluxApproximationBase::GetCatalog()
{
  static FluxApproximationBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void FluxApproximationBase::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
{
  GEOSX_MARK_FUNCTION;

  DomainPartition & domain = *rootGroup->GetGroup< DomainPartition >( keys::domain );

  // Compute the main cell-based stencil
  computeCellStencil( domain );

  // For each boundary condition on target field, create a boundary stencil
  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();
  fsManager.Apply( 0.0,
                   &domain,
                   "faceManager",
                   m_boundaryFieldName,
                   [&] ( FieldSpecificationBase const *,
                         string const & setName,
                         SortedArrayView< localIndex const > const & faceSet,
                         Group const *,
                         string const & )
  {
    computeBoundaryStencil( domain, setName, faceSet );
  } );
}

} //namespace geosx
