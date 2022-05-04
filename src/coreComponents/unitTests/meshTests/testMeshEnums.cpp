/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2020-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "mesh/ElementType.hpp"
#include "mesh/SurfaceElementRegion.hpp"

#include <gtest/gtest.h>

using namespace geosx;

// The `ENUM_STRING` implementation relies on consistency between the order of the `enum`,
// and the order of the `string` array provided. Since this consistency is not enforced, it can be corrupted anytime.
// This unit test aims at preventing from this implicit relationship to bring a bug.

TEST( MeshEnums, ElementType )
{
  using EnumType = ElementType;

  ASSERT_EQ( "Vertex", toString( EnumType::Vertex ) );
  ASSERT_EQ( "BEAM", toString( EnumType::Line ) );
  ASSERT_EQ( "C2D3", toString( EnumType::Triangle ) );
  ASSERT_EQ( "C2D4", toString( EnumType::Quadrilateral ) );
  ASSERT_EQ( "Polygon", toString( EnumType::Polygon ) );
  ASSERT_EQ( "C3D4", toString( EnumType::Tetrahedron ) );
  ASSERT_EQ( "C3D5", toString( EnumType::Pyramid ) );
  ASSERT_EQ( "C3D6", toString( EnumType::Wedge ) );
  ASSERT_EQ( "C3D8", toString( EnumType::Hexahedron ) );
  ASSERT_EQ( "C3D10", toString( EnumType::Prism5 ) );
  ASSERT_EQ( "C3D12", toString( EnumType::Prism6 ) );
  ASSERT_EQ( "Polyhedron", toString( EnumType::Polyhedron ) );

  ASSERT_EQ( numElementTypes(), 12 );
}


TEST( MeshEnums, SurfaceSubRegionType )
{
  using EnumType = SurfaceElementRegion::SurfaceSubRegionType;

  ASSERT_EQ( "faceElement", toString( EnumType::faceElement ) );
  ASSERT_EQ( "embeddedElement", toString( EnumType::embeddedElement ) );
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
