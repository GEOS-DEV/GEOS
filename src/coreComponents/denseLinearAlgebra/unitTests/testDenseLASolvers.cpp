/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "denseLinearAlgebra/denseLASolvers.hpp"

// TPL includes
#include <gtest/gtest.h>


using namespace geos;


constexpr real64 machinePrecision = 1.0e2 * LvArray::NumericLimits< real64 >::epsilon;

template< std::ptrdiff_t N >
struct Reference 
{};

template<>
struct Reference< 2 >
{

  static constexpr std::ptrdiff_t size() { return 2; }
  
  static constexpr real64 A[2][2] = { { 2.0, -1.0},
                                      { 0.0, 1.0 } };
  
  static constexpr real64 rhs[2] = {0.0, 1.0};

  static constexpr real64 solution[2] = {0.5, 1.0};
};

template<>
struct Reference< 3 >
{
  static constexpr std::ptrdiff_t size() { return 3; }

  static constexpr real64 A[3][3] = { { 4.0, 3.0, -2.0},
                                      { 2.0, -1.0, 5.0 },
                                      { -1.0, 3.0, -2.0 } };
  
  static constexpr real64 rhs[3] = {-8.0, 19.0, -13.0};

  static constexpr real64 solution[3] = {1.0, -2.0, 3.0};
};

template< typename MATRIX_TYPE, typename VECTOR_TYPE, typename REFERENCE >
void assemble( MATRIX_TYPE & A, VECTOR_TYPE & rhs )
{
  for( std::ptrdiff_t i=0; i < REFERENCE::size(); ++i )
  {
    rhs[i] = REFERENCE::rhs[i];
    for( std::ptrdiff_t j=0; j < REFERENCE::size(); ++j )
    {
      A[i][j] = REFERENCE::A[i][j];
    }
  }
}

template< typename REFERENCE >
void test_denseLASolve()
{
  constexpr std::ptrdiff_t N = REFERENCE::size();

  real64 A[N][N];
  real64 b[N];
  real64 sol[N];

  assemble< decltype(A), decltype(b), REFERENCE >( A, b );
     
  denseLinearAlgebra::solve< N >( A, b, sol );

  for( std::ptrdiff_t i = 0; i < N; ++i )
  {
    EXPECT_NEAR( sol[i],
                 REFERENCE::solution[i],
                 machinePrecision );
  }
}

TEST( denseLASolve, testTwoByTwo )
{
  test_denseLASolve< Reference< 2 > >();
}

TEST( denseLASolve, testThreeByThree )
{
  test_denseLASolve< Reference< 3 > >();
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}