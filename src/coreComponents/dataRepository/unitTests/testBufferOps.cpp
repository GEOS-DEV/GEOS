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

#include <gtest/gtest.h>

#include "dataRepository/BufferOps.hpp"

using namespace geos;
using namespace bufferOps;


TEST( testGeosxTraits, test_is_host_packable_object )
{
  static_assert( is_host_packable_object_v< int >, "Should be true." );
  static_assert( is_host_packable_object_v< double >, "Should be true." );
  static_assert( is_host_packable_object_v< R1Tensor >, "Should be true." );
  static_assert( is_host_packable_object_v< string >, "Should be true." );

  static_assert( !is_host_packable_object_v< void >, "Should be false." );
  static_assert( !is_host_packable_object_v< array1d< double > >, "Should be false." );
  static_assert( !is_host_packable_object_v< SortedArray< double > >, "Should be false." );
  static_assert( !is_host_packable_object_v< map< string, int > >, "Should be false." );
  static_assert( !is_host_packable_object_v< std::pair< string, int > >, "Should be false." );
}


TEST( testGeosxTraits, test_is_array_packable )
{
  static_assert( is_host_packable_array_v< array2d< real64, RAJA::PERM_IJ > >, "Should be true." );
  static_assert( is_host_packable_array_v< array2d< real64, RAJA::PERM_JI > >, "Should be true." );

  static_assert( !is_host_packable_array_v< int >, "Should be false." );
  static_assert( !is_host_packable_array_v< double >, "Should be false." );
  static_assert( !is_host_packable_array_v< void >, "Should be false." );
}


TEST( testGeosxTraits, test_is_host_packable_map )
{
  static_assert( is_host_packable_map_v< map< string, int > >, "Should be true." );
  static_assert( is_host_packable_map_v< map< string, array1d< int > > >, "Should be true." );
  static_assert( !is_host_packable_map_v< map< string, std::pair< int, int > > >, "Should be false" );
}

TEST( testGeosxTraits, test_is_host_packable_scalar )
{
  static_assert( is_host_packable_scalar_v< int >, "Should be true." );
  static_assert( is_host_packable_scalar_v< double >, "Should be true." );
  static_assert( !is_host_packable_scalar_v< R1Tensor >, "Should be false" );
  static_assert( !is_host_packable_scalar_v< string >, "Should be false" );
  static_assert( !is_host_packable_scalar_v< array1d< double > >, "Should be false." );
}
