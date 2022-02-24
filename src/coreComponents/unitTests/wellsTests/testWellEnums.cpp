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

#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"

#include <gtest/gtest.h>

using namespace geosx;

// The `ENUM_STRING` implementation relies on consistency between the order of the `enum`,
// and the order of the `string` array provided. Since this consistency is not enforced, it can be corrupted anytime.
// This unit test aims at preventing from this implicit relationship to bring a bug.

TEST( WellControlsEnums, Type )
{
  using EnumType = WellControls::Type;

  ASSERT_EQ( "producer", toString( EnumType::PRODUCER ) );
  ASSERT_EQ( "injector", toString( EnumType::INJECTOR ) );
}


TEST( WellControlsEnums, Control )
{
  using EnumType = WellControls::Control;

  ASSERT_EQ( "BHP", toString( EnumType::BHP ) );
  ASSERT_EQ( "phaseVolRate", toString( EnumType::PHASEVOLRATE ) );
  ASSERT_EQ( "totalVolRate", toString( EnumType::TOTALVOLRATE ) );
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
