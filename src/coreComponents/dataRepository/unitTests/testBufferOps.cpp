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
  static_assert( is_host_packable_object< int >, "Should be true." );
  static_assert( is_host_packable_object< double >, "Should be true." );
  static_assert( is_host_packable_object< R1Tensor >, "Should be true." );
  static_assert( is_host_packable_object< string >, "Should be true." );

  static_assert( !is_host_packable_object< void >, "Should be false." );
  static_assert( !is_host_packable_object< array1d< double > >, "Should be false." );
  static_assert( !is_host_packable_object< SortedArray< double > >, "Should be false." );
  static_assert( !is_host_packable_object< map< string, int > >, "Should be false." );
  static_assert( !is_host_packable_object< std::pair< string, int > >, "Should be false." );
}


TEST( testGeosxTraits, test_is_array_packable )
{
  static_assert( is_host_packable_array< array2d< real64, RAJA::PERM_IJ > >, "Should be true." );
  static_assert( is_host_packable_array< array2d< real64, RAJA::PERM_JI > >, "Should be true." );

  static_assert( !is_host_packable_array< int >, "Should be false." );
  static_assert( !is_host_packable_array< double >, "Should be false." );
  static_assert( !is_host_packable_array< void >, "Should be false." );
}


TEST( testGeosxTraits, test_is_host_packable_map )
{
  static_assert( is_host_packable_map< map< string, int > >, "Should be true." );
  static_assert( is_host_packable_map< map< string, array1d< int > > >, "Should be true." );
  static_assert( !is_host_packable_map< map< string, std::pair< int, int > > >, "Should be false" );
}

TEST( testGeosxTraits, test_is_host_packable_scalar )
{
  static_assert( is_host_packable_scalar< int >, "Should be true." );
  static_assert( is_host_packable_scalar< double >, "Should be true." );
  static_assert( !is_host_packable_scalar< R1Tensor >, "Should be false" );
  static_assert( !is_host_packable_scalar< string >, "Should be false" );
  static_assert( !is_host_packable_scalar< array1d< double > >, "Should be false." );
}
