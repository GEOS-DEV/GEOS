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
#include "fieldSpecification/FieldSpecificationManager.hpp"

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
  ps.getReference< string >( PermeabilitySpecification::viewKeyStruct::permeabilityModelNameString() ) = "rockPerm";
  array1d< real64 > & scales = ps.getReference< array1d< real64 > >( PermeabilitySpecification::viewKeyStruct::scalesString() );
  scales.resize( 3 );
  scales[0] = 9.869233e-16;
  scales[1] = 9.869233e-16;
  scales[2] = 9.869233e-16;
}

}

TEST( PermeabilitySpecificationTest, ExpansionPropagatesAttributes )
{
  using ProcessorRegistry = FieldSpecificationProcessorRegistry;

  conduit::Node node;
  Group root( "root", node );

  FieldSpecificationManager manager( "FieldSpecifications", &root );
  root.registerGroup( manager.getName(), &manager );

  PermeabilitySpecification permSpec( "permSpec", &manager );
  manager.registerGroup( permSpec.getName(), &permSpec );


  fillValidInput( permSpec );

  // verify that the permeability specification processor exists
  EXPECT_NE( ProcessorRegistry::getProcessors().find( PermeabilitySpecification::catalogName() ),
             ProcessorRegistry::getProcessors().end() ) << GEOS_FMT( "Processor of {} does not exist", PermeabilitySpecification::catalogName() );

  // indirectly call postInputInitialization() -> expandFieldSpecification() on permSpec
  EXPECT_NO_THROW( root.postInputInitializationRecursive() );

  FieldSpecification const * generatedFS = manager.getGroupPointer< FieldSpecification >( "permSpec_region2" );
  // verify that  the generated ("expanded") field specification exists
  EXPECT_NE( generatedFS, nullptr ) << "Field specification has not been generated.";
  // verify that the scale values have been well applied on the generated ("expanded") field specification
  EXPECT_DOUBLE_EQ( generatedFS->getScales()[0], 9.869233e-16 );
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
