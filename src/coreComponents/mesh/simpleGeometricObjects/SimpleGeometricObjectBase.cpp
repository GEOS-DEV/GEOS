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

/**
 * @file SimpleGeometricObjectBase.cpp
 */

#include "SimpleGeometricObjectBase.hpp"
#include "dataRepository/InputFlags.hpp"
#include "dataRepository/ProblemManagerBase.hpp"
#include "mesh/DomainPartition.hpp"

namespace geos
{

SimpleGeometricObjectBase::SimpleGeometricObjectBase( string const & name,
                                                      Group * const parent ):
  Group( name, parent )
{
  setInputFlags( dataRepository::InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::epsilonString(), &m_epsilon ).
    setApplyDefaultValue( -1 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Absolute tolerance for coordinate checks. "
                    "If not specified, default value of 1e-6 * GlobalLengthScale will be used, where GlobalLengthScale is the length scale of the domain." );
}

SimpleGeometricObjectBase::CatalogInterface::CatalogType & SimpleGeometricObjectBase::getCatalog()
{
  static SimpleGeometricObjectBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void SimpleGeometricObjectBase::postInputInitialization()
{
  if( m_epsilon < 0.0 ) // if not specified in input, compute it
  {
    // determine m_epsilon
    m_epsilon = std::numeric_limits< real64 >::max();
    DomainPartition & domain = getProblemManagerBase( *this ).getDomainPartition();
    domain.forMeshBodies( [&]( MeshBody const & meshBody )
    {
      m_epsilon = std::min( m_epsilon, 1e-6 * meshBody.getGlobalLengthScale() );
    } );
  }
}

} /// namespace geos
