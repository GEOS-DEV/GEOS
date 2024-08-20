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

#include <random>


using namespace geos;


constexpr real64 machinePrecision = 1.0e2 * LvArray::NumericLimits< real64 >::epsilon;

template< std::ptrdiff_t N >
struct LinearSystem
{
  static constexpr std::ptrdiff_t size() { return N; }

  LinearSystem( short int seed )
  {
    std::mt19937 generator( seed ); 
    std::uniform_real_distribution< real64 > distribution( -10.0, 10.0 );
    for( ptrdiff_t i=0; i<N; ++i )
    {
      solution[i] = distribution( generator );
      for( ptrdiff_t j=0; j<N; ++j )
      {
        matrix[i][j] =  distribution( generator );
      }
    }

    // Compute rhs as matrix * solution
    for( std::ptrdiff_t i = 0; i < N; ++i )
    {
      rhs[i] = 0.0;
      for( std::ptrdiff_t j = 0; j < N; ++j )
      {
        rhs[i] += matrix[i][j] * solution[j];
      }
    }
  }

  real64 matrix[N][N];
  real64 solution[N];
  real64 rhs[N];
};

template< typename LINEAR_SYSTEM >
void test_denseLASolve()
{
  constexpr std::ptrdiff_t N = LINEAR_SYSTEM::size();

  LINEAR_SYSTEM LS(2024);
  real64 sol[N];

  denseLinearAlgebra::solve< N >( LS.matrix, LS.rhs, sol );

  for( std::ptrdiff_t i = 0; i < N; ++i )
  {
    EXPECT_NEAR( sol[i],
                 LS.solution[i],
                 machinePrecision );
  }
}

TEST( denseLASolve, testTwoByTwo )
{
  test_denseLASolve< LinearSystem< 2 > >();
}

TEST( denseLASolve, testThreeByThree )
{
  test_denseLASolve< LinearSystem< 3 > >();
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
