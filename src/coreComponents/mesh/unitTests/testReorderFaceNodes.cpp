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

#include "../generators/NewGlobalNumbering.hpp"

#include <gtest/gtest.h>

namespace geos
{

struct Data
{
  std::vector< int > const input;
  bool const isFlipped;
  std::uint8_t const start;
  std::vector< int > const expected;
};

class TestReorderFaceNodes : public ::testing::TestWithParam< Data >
{

};

template< class I >
std::vector< I > conv( std::vector< int > const & input )
{
  std::vector< I > res;
  res.reserve( std::size( input ) );
  for( int const & i: input )
  {
    res.emplace_back( I{ i } );
  }
  return res;
}

TEST_P( TestReorderFaceNodes, BB )
{
  Data const & d = GetParam();
  bool isFlipped;
  std::uint8_t start;
  std::vector< NodeGlbIdx > const res = ghosting::reorderFaceNodes( conv< NodeGlbIdx >( d.input ), isFlipped, start );
  EXPECT_EQ( std::make_tuple( conv< NodeGlbIdx >( d.expected ), d.isFlipped, d.start ), std::tie( res, isFlipped, start ) );
  EXPECT_EQ( conv< NodeLocIdx >( d.input ), ghosting::resetFaceNodes( conv< NodeLocIdx >( d.expected ), d.isFlipped, d.start ) );
}

INSTANTIATE_TEST_SUITE_P( instance, TestReorderFaceNodes, ::testing::Values(
  Data{ { 0, 1, 2 }, false, 0, { 0, 1, 2 } },
  Data{ { 2, 0, 1 }, false, 1, { 0, 1, 2 } },
  Data{ { 2, 1, 0 }, true, 2, { 0, 1, 2 } },
  Data{ { 0, 2, 1 }, true, 0, { 0, 1, 2 } },
  Data{ { 1, 0, 2, 3 }, true, 1, { 0, 1, 3, 2 } },
  Data{ { 2, 0, 1, 3 }, false, 1, { 0, 1, 3, 2 } }
) );

}  // end of namespace
