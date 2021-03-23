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
  real64 bulk, poisson, young, shear;
  real64 bulk_out, poisson_out, young_out, shear_out;

  // bulk and poisson
  bulk        = 5.4321;
  poisson     = 0.2345;
  //young       = conversions::BulkModAndPoissonRatio::toYoungsMod( bulk, poisson );
  young = YoungModulus().
            setBulkModulus( bulk ).
            setPoissonRatio( poisson ).
            getValue();

  //shear       = conversions::BulkModAndPoissonRatio::toShearMod( bulk, poisson );
  shear = ShearModulus().
             setBulkModulus( bulk ).
             setPoissonRatio( poisson ).
             getValue();

  EXPECT_DOUBLE_EQ( young, 3*bulk*(1-2*poisson) );
  EXPECT_DOUBLE_EQ( shear, 3*bulk*(1-2*poisson)/(2+2*poisson) );

  // bulk and shear
  //poisson_out = conversions::BulkModAndShearMod::toPoissonRatio( bulk, shear );
  poisson_out = PoissonRatio().
                  setBulkModulus( bulk ).
                  setShearModulus( shear ).
                  getValue();

  //young_out = conversions::BulkModAndShearMod::toYoungsMod( bulk, shear );
  young_out = YoungModulus().
                 setBulkModulus( bulk ).
                 setShearModulus( shear ).
                 getValue();
  
  EXPECT_DOUBLE_EQ( poisson, poisson_out );
  EXPECT_DOUBLE_EQ( young, young_out );

  // bulk and youngs
  //shear_out = conversions::BulkModAndYoungsMod::toShearMod( bulk, young );
  shear_out = ShearModulus().
                setBulkModulus( bulk ).
                setYoungModulus( young ).
                getValue();

  //poisson_out = conversions::BulkModAndYoungsMod::toPoissonRatio( bulk, young );
  poisson_out = PoissonRatio().
                  setBulkModulus( bulk ).
                  setYoungModulus( young ).
                  getValue();

  EXPECT_DOUBLE_EQ( poisson, poisson_out );
  EXPECT_DOUBLE_EQ( shear, shear_out );

  // young and shear
  //bulk_out    = conversions::ShearModAndYoungsMod::toBulkMod( shear, young );
  //poisson_out = conversions::ShearModAndYoungsMod::toPoissonRatio( shear, young );
  bulk_out = BulkModulus().
               setYoungModulus( young ).
               setShearModulus( shear ).
               getValue();

  poisson_out = PoissonRatio().
                  setShearModulus( shear ).
                  setYoungModulus( young ).
                  getValue();

  EXPECT_DOUBLE_EQ( bulk, bulk_out );
  EXPECT_DOUBLE_EQ( poisson, poisson_out );

  // young and poisson
  //bulk_out = conversions::YoungsModAndPoissonRatio::toBulkMod( young, poisson );
  //shear_out = conversions::YoungsModAndPoissonRatio::toShearMod( young, poisson );
  bulk_out = BulkModulus().
               setYoungModulus( young ).
               setPoissonRatio( poisson ).
               getValue();

  shear_out = ShearModulus().
               setYoungModulus( young ).
               setPoissonRatio( poisson ).
               getValue();

  EXPECT_DOUBLE_EQ( bulk, bulk_out );
  EXPECT_DOUBLE_EQ( shear, shear_out );

  // shear and poisson
  //bulk_out = conversions::ShearModAndPoissonRatio::toBulkMod( shear, poisson );
  //young_out =conversions::ShearModAndPoissonRatio::toYoungsMod( shear, poisson );

  bulk_out = BulkModulus().
               setShearModulus( shear ).
               setPoissonRatio( poisson ).
               getValue();

  young_out = YoungModulus().
                 setShearModulus( shear ).
                 setPoissonRatio( poisson ).
                 getValue();

  EXPECT_DOUBLE_EQ( bulk, bulk_out );
  EXPECT_DOUBLE_EQ( young, young_out );

  // first lame (bulk and shear input only)
  //real64 lame = conversions::BulkModAndShearMod::toFirstLame( bulk, shear );
  real64 const lame = LameModulus().
                          setBulkModulus( bulk ).
                          setShearModulus( shear ).
                          getValue();

  EXPECT_DOUBLE_EQ( lame, bulk-2*shear/3 );
}
