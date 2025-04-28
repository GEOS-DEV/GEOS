/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "FieldSpecificationBase.hpp"

#include "common/MpiWrapper.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/LogLevelsInfo.hpp"

namespace geos
{
using namespace dataRepository;

FieldSpecificationBase::FieldSpecificationBase( string const & name, Group * parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::setNamesString(), &m_setNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription(
    "Name of sets  that boundary condition is applied to." \
    "Sets name can be." \
    "all, to select any mesh element," \
    "xpos, ypos, zpos to select all elements facing the positive x, y or z axis" \
    "xneg, yneg, zneg to select all elements facing the negative x, y or z axis" \
    "can be the name of an geometrical object in <Geometry> to select the mesh elements that are strictly contained in it." );

  registerWrapper( viewKeyStruct::objectPathString(), &m_objectPath ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Path to the target field" );

  registerWrapper( viewKeyStruct::fieldNameString(), &m_fieldName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of field that boundary condition is applied to." );

  registerWrapper( viewKeyStruct::componentString(), &m_component ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Component name of field (if tensor) to apply boundary condition to. Component name must use the order in which the phaseNames" );

  registerWrapper( viewKeyStruct::directionString(), &m_direction ).
    setApplyDefaultValue( {0, 0, 0} ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Direction to apply boundary condition to." );

  registerWrapper( viewKeyStruct::functionNameString(), &m_functionName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of function that specifies variation of the boundary condition." );

  registerWrapper( viewKeyStruct::bcApplicationTableNameString(), &m_bcApplicationFunctionName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of table that specifies the on/off application of the boundary condition." );

  registerWrapper( viewKeyStruct::scaleString(), &m_scale ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Scale factor by which the value of the initial condition is multiplied" );

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

  addLogLevel< logInfo::BoundaryCondition >();
  addLogLevel< logInfo::FaceBoundaryCondition >();
  addLogLevel< logInfo::SourceFluxFailure >();
}


FieldSpecificationBase::~FieldSpecificationBase()
{}

FieldSpecificationBase::CatalogInterface::CatalogType &
FieldSpecificationBase::getCatalog()
{
  static FieldSpecificationBase::CatalogInterface::CatalogType catalog;
  return catalog;
}



void FieldSpecificationBase::setMeshObjectPath( Group const & meshBodies )
{
  try
  {
    m_meshObjectPaths = std::make_unique< MeshObjectPath >( m_objectPath, meshBodies );
  }
  catch( std::exception const & e )
  {
    throw InputError( e, getWrapperDataContext( viewKeyStruct::objectPathString() ).toString() +
                      " is a wrong objectPath: " + m_objectPath + "\n" );
  }
}



REGISTER_CATALOG_ENTRY( FieldSpecificationBase, FieldSpecificationBase, string const &, Group * const )

}
