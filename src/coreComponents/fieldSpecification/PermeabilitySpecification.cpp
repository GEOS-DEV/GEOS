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

#include "PermeabilitySpecification.hpp"

namespace geos
{
using namespace dataRepository;

PermeabilitySpecification::PermeabilitySpecification( string const & name, Group * parent ):
  FieldSpecificationABC( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::setNamesString(), &m_setNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Names of sets that the boundary condition is applied to.\n"
                    "A set can contain heterogeneous elements in the mesh (volumes, nodes, faces, edges).\n"
                    "A set can be be defined by a 'Geometry' component, or correspond to imported sets in case of an external mesh" );

  registerWrapper( viewKeyStruct::regionNamesString(), &m_regionNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Names of the regions that boundary condition is applied to." );

  registerWrapper( viewKeyStruct::fieldNameString(), &m_fieldName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of field that boundary condition is applied to.\n"
                    "A field can represent a physical variable. (pressure, temperature, global composition fraction of the fluid, ...)" );

  registerWrapper( viewKeyStruct::functionNameString(), &m_functionName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of function that specifies variation of the boundary condition." );

  registerWrapper( viewKeyStruct::scalesString(), &m_scales ).
    // setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Apply a scaling factor for the value of the boundary condition." );
}


PermeabilitySpecification::~PermeabilitySpecification()
{}


void PermeabilitySpecification::postInputInitialization()
{
  GEOS_THROW_IF( m_scales.size() != 3,
                 getWrapperDataContext( viewKeyStruct::scalesString() ) << 
                 ": expected 3 components, got " << m_scales.size(),
                 InputError );
}


REGISTER_CATALOG_ENTRY( FieldSpecificationABC, PermeabilitySpecification, string const &, Group * const )

}