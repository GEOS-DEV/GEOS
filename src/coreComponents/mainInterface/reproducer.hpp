#pragma once

#include "RAJA/RAJA.hpp"

#include <vector>
#include <iostream>

inline void testInline()
{
  std::vector< int > data { 0, 3, 1, 1, 0, 1, 0, 0, 0 };

  // Perform an inplace prefix-sum to get the unique edge offset.
  RAJA::inclusive_scan_inplace< RAJA::omp_parallel_for_exec >( RAJA::make_span( data.data(), data.size() ) );

  std::cout << " after scan = ";

  for( auto const & value : data )
  {
    std::cout << value << ", ";
  }

  std::cout << std::endl;
}

void test();
