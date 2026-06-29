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
 * @file FieldSpecificationFactory.cpp
 */

#include "FieldSpecificationFactory.hpp"
#include "PermeabilitySpecification.hpp"

namespace geos
{

namespace
{

template< typename ... SPEC_TYPES, typename LAMBDA >
void forExpandableSpecifications( dataRepository::Group & manager,
                                  types::TypeList< SPEC_TYPES... >,
                                  LAMBDA && lambda )
{
  manager.template forSubGroups< SPEC_TYPES... >( std::forward< LAMBDA >( lambda ) );
}

}

void expandFieldSpecifications( dataRepository::Group & manager )
{
  forExpandableSpecifications( manager, ExpandableSpecTypes{}, [&]( auto const & spec )
  {
    generateFieldSpecifications( spec, manager );
  } );
}

} // namespace geos
