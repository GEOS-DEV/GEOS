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

#include "fieldSpecification/PermeabilitySpecification.hpp"
#include "fieldSpecification/FieldSpecification.hpp"

#include <gtest/gtest.h>
#include <conduit.hpp>

using namespace geos;
using namespace dataRepository;

namespace
{

void fillValidInput( PermeabilitySpecification & ps )
{
  ps.getReference< string_array >( PermeabilitySpecification::viewKeyStruct::setNamesString() ) = { "all" };
  ps.getReference< string_array >( PermeabilitySpecification::viewKeyStruct::regionNamesString() ) = { "region1", "region2" };
  ps.getReference< string >( PermeabilitySpecification::viewKeyStruct::fieldNameString() ) = "rockPerm_permeability";
  array1d< real64 > & scales = ps.getReference< array1d< real64 > >( PermeabilitySpecification::viewKeyStruct::scalesString() );
  scales.resize( 3 );
  scales[0] = 9.869233e-16;
  scales[1] = 9.869233e-16;
  scales[2] = 9.869233e-16;
}

}

TEST( PermeabilitySpecificationTest, ExpansionPropagatesAttributes )
{
  conduit::Node node;
  Group rootGroup( "root", node );
  PermeabilitySpecification ps( "perm", &rootGroup );

  fillValidInput( ps );

  EXPECT_NO_THROW( ps.postInputInitializationRecursive() );

  Group manager( "FieldSpecifications", &rootGroup );
  generateFieldSpecifications( ps, manager );

  FieldSpecification const & fs = manager.getGroup< FieldSpecification >( "perm_region2" );

  EXPECT_DOUBLE_EQ( fs.getScales()[0], 9.869233e-16 );
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
