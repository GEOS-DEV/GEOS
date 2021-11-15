/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "FieldSpecificationBase.hpp"

#include "common/MpiWrapper.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "mesh/DomainPartition.hpp"


namespace geosx
{
using namespace dataRepository;

FieldSpecificationBase::FieldSpecificationBase( string const & name, Group * parent ):
  Group( name, parent ),
  m_normalizeBySetSize( false )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::setNamesString(), &m_setNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Name of sets that boundary condition is applied to." );

  registerWrapper( viewKeyStruct::objectPathString(), &m_objectPath ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Path to the target field" );

  registerWrapper( viewKeyStruct::fieldNameString(), &m_fieldName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of field that boundary condition is applied to." );

  registerWrapper( viewKeyStruct::componentString(), &m_component ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Component of field (if tensor) to apply boundary condition to." );

  registerWrapper( viewKeyStruct::directionString(), &m_direction ).
    setApplyDefaultValue( {0, 0, 0} ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Direction to apply boundary condition to." );

  registerWrapper( viewKeyStruct::functionNameString(), &m_functionName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of function that specifies variation of the boundary condition." );

  registerWrapper( viewKeyStruct::bcApplicationTableNameString(), &m_bcApplicationFunctionName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of table that specifies the on/off application of the boundary condition." );

  registerWrapper( viewKeyStruct::scaleString(), &m_scale ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Scale factor for value of the boundary condition." );

  registerWrapper( viewKeyStruct::initialConditionString(), &m_initialCondition ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Boundary condition is applied as an initial condition." );

  registerWrapper( viewKeyStruct::beginTimeString(), &m_beginTime ).
    setApplyDefaultValue( -1.0e99 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Time at which the boundary condition will start being applied." );

  registerWrapper( viewKeyStruct::endTimeString(), &m_endTime ).
    setApplyDefaultValue( 1.0e99 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Time at which the boundary condition will stop being applied." );

}


FieldSpecificationBase::~FieldSpecificationBase()
{}

FieldSpecificationBase::CatalogInterface::CatalogType &
FieldSpecificationBase::getCatalog()
{
  static FieldSpecificationBase::CatalogInterface::CatalogType catalog;
  return catalog;
}



void FieldSpecificationBase::postProcessInput()
{}



REGISTER_CATALOG_ENTRY( FieldSpecificationBase, FieldSpecificationBase, string const &, Group * const )

}
