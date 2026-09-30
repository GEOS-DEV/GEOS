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

/** @file testCellElementRegionSelector.cpp */
#include "common/MpiWrapper.hpp"
#include "mesh/CellElementRegionSelector.hpp"
#include "mesh/generators/CellBlockManager.hpp"

#include <gtest/gtest.h>
#include <conduit.hpp>

using namespace geos;

namespace
{
string const pyramids = "1_pyramids";
string const tetrahedra = "1_tetrahedra";
string const introduced = "1_pyramids__refined_tetrahedra";

class RefinedRegionSelection : public ::testing::Test
{
protected:
  conduit::Node repository;
  dataRepository::Group root{ "Problem", repository };
  CellBlockManager blocks{ "cellBlocks", &root };
  CellElementRegion pyramidMaterial{ "pyramidMaterial", &root };
  CellElementRegion tetrahedronMaterial{ "tetrahedronMaterial", &root };

  RefinedRegionSelection()
  {
    for( auto const & name : { pyramids, tetrahedra, introduced } )
      blocks.registerCellBlock( name, 1 ); // Empty local blocks still belong to the agreed schema.
    blocks.setSourceCellBlockDescendants( { { pyramids, { pyramids, introduced } }, { tetrahedra, { tetrahedra } } } );
  }
  CellElementRegionSelector selector()
  {
    return { blocks.getCellBlocks(), blocks.getRegionAttributesCellBlocks(), &blocks.getSourceCellBlockDescendants() };
  }
};
} // namespace

TEST_F( RefinedRegionSelection, ExactSourcesKeepDifferentMaterialsSeparate )
{
  auto selection = selector();
  pyramidMaterial.setCellBlockNames( std::set< string >{ pyramids } );
  tetrahedronMaterial.setCellBlockNames( std::set< string >{ tetrahedra } );
  EXPECT_EQ( selection.buildCellBlocksSelection( pyramidMaterial ), ( std::set< string >{ pyramids, introduced } ) );
  EXPECT_EQ( selection.buildCellBlocksSelection( tetrahedronMaterial ), ( std::set< string >{ tetrahedra } ) );
  EXPECT_NO_THROW( selection.checkSelectionConsistency() );
}

TEST_F( RefinedRegionSelection, PatternsMatchTheOriginalNamespace )
{
  auto selection = selector();
  pyramidMaterial.setCellBlockNames( std::set< string >{ "*_pyramids" } );
  tetrahedronMaterial.setCellBlockNames( std::set< string >{ "*_tetrahedra" } );
  EXPECT_EQ( selection.buildCellBlocksSelection( pyramidMaterial ), ( std::set< string >{ pyramids, introduced } ) );
  EXPECT_EQ( selection.buildCellBlocksSelection( tetrahedronMaterial ), ( std::set< string >{ tetrahedra } ) );
  EXPECT_NO_THROW( selection.checkSelectionConsistency() );
  pyramidMaterial.setCellBlockNames( std::set< string >{ introduced } );
  EXPECT_THROW( selection.buildCellBlocksSelection( pyramidMaterial ), InputError );
}

TEST_F( RefinedRegionSelection, AttributeAndWildcardSelectionsExpandAllSources )
{
  for( string const & qualifier : { string( "1" ), string( "*" ) } )
  {
    auto selection = selector();
    pyramidMaterial.setCellBlockNames( std::set< string >{ qualifier } );
    EXPECT_EQ( selection.buildCellBlocksSelection( pyramidMaterial ), ( std::set< string >{ pyramids, tetrahedra, introduced } ) );
    EXPECT_NO_THROW( selection.checkSelectionConsistency() );
  }
}

TEST_F( RefinedRegionSelection, CoverageAndConflictsAreCheckedAfterExpansion )
{
  auto selection = selector();
  pyramidMaterial.setCellBlockNames( std::set< string >{ pyramids } );
  selection.buildCellBlocksSelection( pyramidMaterial );
  EXPECT_THROW( selection.checkSelectionConsistency(), InputError );
  tetrahedronMaterial.setCellBlockNames( std::set< string >{ "*" } );
  selection.buildCellBlocksSelection( tetrahedronMaterial );
  EXPECT_THROW( selection.checkSelectionConsistency(), InputError );
}

TEST_F( RefinedRegionSelection, InvalidLineageCannotHideOrMergeFinalBlocks )
{
  for( SourceCellBlockDescendants const & bad :
       { SourceCellBlockDescendants{ { pyramids, { pyramids, introduced } } },
         SourceCellBlockDescendants{ { pyramids, { pyramids, introduced } }, { tetrahedra, { tetrahedra, introduced } } },
         SourceCellBlockDescendants{ { pyramids, { pyramids, introduced } }, { tetrahedra, { "missing" } } },
         SourceCellBlockDescendants{ { pyramids, {} }, { tetrahedra, { tetrahedra, introduced } } } } )
    EXPECT_THROW( CellElementRegionSelector( blocks.getCellBlocks(), blocks.getRegionAttributesCellBlocks(), &bad ), InputError );
}

TEST( CellElementRegionSelector, EmptyLineagePreservesLegacySelection )
{
  conduit::Node repository;
  dataRepository::Group root{ "Problem", repository };
  CellBlockManager blocks{ "cellBlocks", &root };
  blocks.registerCellBlock( pyramids, 1 );
  blocks.registerCellBlock( tetrahedra, 1 );
  EXPECT_TRUE( blocks.getSourceCellBlockDescendants().empty() );
  CellElementRegion region{ "material", &root };
  region.setCellBlockNames( std::set< string >{ "*" } );
  CellElementRegionSelector selector{ blocks.getCellBlocks(), blocks.getRegionAttributesCellBlocks(),
                                      &blocks.getSourceCellBlockDescendants() };
  EXPECT_EQ( selector.buildCellBlocksSelection( region ), ( std::set< string >{ pyramids, tetrahedra } ) );
  EXPECT_NO_THROW( selector.checkSelectionConsistency() );
}

int main( int argc, char ** argv )
{
  MpiWrapper::init( &argc, &argv );
  MPI_COMM_GEOS = MpiWrapper::commDup( MPI_COMM_WORLD );
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  MpiWrapper::finalize();
  return result;
}
