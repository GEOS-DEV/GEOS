/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include <gtest/gtest.h>

#include "Array.hpp"
#include "dataRepository/DefaultValue.hpp"

#include <functional>
#include <string>
#include <typeindex>
#include <vector>

using namespace geosx;
using namespace dataRepository;


TEST( testDefaultValue, testScalar )
{
  EXPECT_TRUE( DefaultValue< int >::has_default_value == true );
  EXPECT_TRUE( DefaultValue< long int >::has_default_value == true );
  EXPECT_TRUE( DefaultValue< long long int >::has_default_value == true );
  EXPECT_TRUE( DefaultValue< double >::has_default_value == true );
}

TEST( testDefaultValue, testArray )
{
  using array1 = array_decl< double, 1 >;
  using array2 = array_decl< double, 2 >;
  using array3 = array_decl< double, 3 >;
  using array4 = array_decl< int, 1 >;
  using array5 = array_decl< long int, 1 >;
  using array6 = array_decl< long long int, 1 >;
  EXPECT_TRUE( DefaultValue< array1 >::has_default_value==true );
  EXPECT_TRUE( DefaultValue< array2 >::has_default_value==true );
  EXPECT_TRUE( DefaultValue< array3 >::has_default_value==true );
  EXPECT_TRUE( DefaultValue< array4 >::has_default_value==true );
  EXPECT_TRUE( DefaultValue< array5 >::has_default_value==true );
  EXPECT_TRUE( DefaultValue< array6 >::has_default_value==true );
}
