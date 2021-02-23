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

#include "gtest/gtest.h"

#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/solid/PropertyConversions.hpp"

using namespace geosx;
using namespace ::geosx::constitutive;

TEST( PropertyConversionTests, testElasticConversions )
{
  real64 bulk, poisson, youngs, shear;
  real64 bulk_out, poisson_out, youngs_out, shear_out;

  // bulk and poisson
  bulk        = 5.4321;
  poisson     = 0.2345;
  youngs      = conversions::BulkModAndPoissonRatio::toYoungsMod( bulk, poisson );
  shear       = conversions::BulkModAndPoissonRatio::toShearMod( bulk, poisson );

  EXPECT_DOUBLE_EQ( youngs, 3*bulk*(1-2*poisson) );
  EXPECT_DOUBLE_EQ( shear, 3*bulk*(1-2*poisson)/(2+2*poisson) );

  // bulk and shear
  poisson_out = conversions::BulkModAndShearMod::toPoissonRatio( bulk, shear );
  youngs_out = conversions::BulkModAndShearMod::toYoungsMod( bulk, shear );

  EXPECT_DOUBLE_EQ( poisson, poisson_out );
  EXPECT_DOUBLE_EQ( youngs, youngs_out );

  // bulk and youngs
  shear_out = conversions::BulkModAndYoungsMod::toShearMod( bulk, youngs );
  poisson_out = conversions::BulkModAndYoungsMod::toPoissonRatio( bulk, youngs );

  EXPECT_DOUBLE_EQ( poisson, poisson_out );
  EXPECT_DOUBLE_EQ( shear, shear_out );

  // youngs and shear
  bulk_out    = conversions::ShearModAndYoungsMod::toBulkMod( shear, youngs );
  poisson_out = conversions::ShearModAndYoungsMod::toPoissonRatio( shear, youngs );

  EXPECT_DOUBLE_EQ( bulk, bulk_out );
  EXPECT_DOUBLE_EQ( poisson, poisson_out );

  // youngs and poisson
  bulk_out = conversions::YoungsModAndPoissonRatio::toBulkMod( youngs, poisson );
  shear_out = conversions::YoungsModAndPoissonRatio::toShearMod( youngs, poisson );

  EXPECT_DOUBLE_EQ( bulk, bulk_out );
  EXPECT_DOUBLE_EQ( shear, shear_out );

  // shear and poisson
  bulk_out = conversions::ShearModAndPoissonRatio::toBulkMod( shear, poisson );
  youngs_out =conversions::ShearModAndPoissonRatio::toYoungsMod( shear, poisson );

  EXPECT_DOUBLE_EQ( bulk, bulk_out );
  EXPECT_DOUBLE_EQ( youngs, youngs_out );

  // first lame (bulk and shear input only)
  real64 lame = conversions::BulkModAndShearMod::toFirstLame( bulk, shear );

  EXPECT_DOUBLE_EQ( lame, bulk-2*shear/3 );
}
