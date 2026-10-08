/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "constitutive/cohesiveZone/SurfaceInformedPolymerCohesiveZone.hpp"

#include <gtest/gtest.h>
#include <cmath>
#include <initializer_list>

using namespace geos;
using namespace geos::constitutive;

namespace
{
using Keys = SurfaceInformedPolymerCohesiveZone::viewKeyStruct;
using Measure = PolymerCohesiveNormalStrainMeasure;

real64 constexpr thickness = 0.1;
real64 constexpr bulkModulus = 260.0;
real64 constexpr shearModulus = 5.0;
real64 constexpr constrainedModulus = bulkModulus + 4.0 * shearModulus / 3.0;

// Allocate a single interface directly: the MPM solver normally sizes these arrays.
struct Film
{
  conduit::Node node;
  dataRepository::Group root{ "root", node };
  SurfaceInformedPolymerCohesiveZone model{ "film", &root };

  Film()
  {
    model.getReference< real64 >( Keys::thicknessString() ) = thickness;
    model.getReference< real64 >( Keys::bulkModulusString() ) = bulkModulus;
    model.getReference< real64 >( Keys::shearModulusString() ) = shearModulus;
    model.getReference< real64 >( Keys::defaultYieldStrengthString() ) = 1.0e6;

    for( char const * key : { Keys::damageString(), Keys::temperatureString(),
                             Keys::previousLambdaString(), Keys::equivalentPlasticStrainString(),
                             Keys::plasticNormalStrainString(), Keys::plasticTangentialStrainString() } )
    {
      array1d< real64 > & state = model.getReference< array1d< real64 > >( key );
      state.resize( 1 );
      state[0] = 0.0;
    }
    state( Keys::temperatureString() ) = 300.0;
    state( Keys::previousLambdaString() ) = 1.0;
  }

  real64 & state( char const * key )
  {
    return model.getReference< array1d< real64 > >( key )[0];
  }

  void select( Measure const measure )
  {
    model.getReference< Measure >( Keys::normalStrainMeasureString() ) = measure;
  }

  void update( real64 const opening, real64 const sliding, real64 & normal, real64 & shear )
  {
    model.createKernelUpdates().jumpDisplacementUpdate( 0, opening, sliding, normal, shear );
  }
};
}

TEST( SurfaceInformedPolymerCohesiveZone, EngineeringDefault )
{
  Film film;
  EXPECT_EQ( film.model.getReference< Measure >( Keys::normalStrainMeasureString() ), Measure::Engineering );
  real64 normal, shear;
  film.update( thickness, 0.0, normal, shear );
  EXPECT_NEAR( normal, -constrainedModulus, 1.0e-10 );
}

TEST( SurfaceInformedPolymerCohesiveZone, ElasticNormalStrainMeasures )
{
  for( Measure const measure : { Measure::Engineering, Measure::Logarithmic } )
  {
    for( real64 const engineeringStrain : { -0.5, 0.1, 1.0 } )
    {
      Film film;
      film.select( measure );
      real64 normal, shear;
      film.update( thickness * engineeringStrain, 0.0, normal, shear );
      real64 const strain = measure == Measure::Logarithmic ? std::log1p( engineeringStrain ) : engineeringStrain;
      EXPECT_NEAR( normal, -constrainedModulus * strain, 1.0e-10 );
      EXPECT_NEAR( shear, 0.0, 1.0e-12 );
      EXPECT_NEAR( film.state( Keys::equivalentPlasticStrainString() ), 0.0, 1.0e-12 );
    }
  }
}

TEST( SurfaceInformedPolymerCohesiveZone, SameInitialTangentAndPureShear )
{
  for( Measure const measure : { Measure::Engineering, Measure::Logarithmic } )
  {
    Film film;
    film.select( measure );
    real64 positive, negative, shear;
    real64 constexpr opening = thickness * 1.0e-6;
    film.update( opening, 0.0, positive, shear );
    film.update( -opening, 0.0, negative, shear );
    EXPECT_NEAR( ( positive - negative ) / ( 2.0 * opening ), -constrainedModulus / thickness, 1.0e-5 );

    real64 normal;
    film.update( 0.0, thickness * 0.2, normal, shear );
    EXPECT_NEAR( normal, 0.0, 1.0e-12 );
    EXPECT_NEAR( shear, -shearModulus * 0.2, 1.0e-12 );
  }
}

TEST( SurfaceInformedPolymerCohesiveZone, CoupledReturnUsesSelectedPlasticNormalStrain )
{
  for( Measure const measure : { Measure::Engineering, Measure::Logarithmic } )
  {
    Film film;
    film.select( measure );
    film.model.getReference< real64 >( Keys::defaultYieldStrengthString() ) = 1.0;
    real64 normal, shear;
    film.update( thickness * 0.1, thickness * 0.2, normal, shear );
    real64 const strain = measure == Measure::Logarithmic ? std::log1p( 0.1 ) : 0.1;
    real64 const deviatoricNormal = -normal - bulkModulus * strain;
    EXPECT_NEAR( std::sqrt( 2.25 * deviatoricNormal * deviatoricNormal + 3.0 * shear * shear ), 1.0, 1.0e-10 );
    EXPECT_NEAR( film.state( Keys::plasticNormalStrainString() ), strain - deviatoricNormal / ( 4.0 * shearModulus / 3.0 ), 1.0e-12 );
    real64 const kappa = film.state( Keys::equivalentPlasticStrainString() );
    EXPECT_GT( kappa, 0.0 );

    real64 repeatedNormal, repeatedShear;
    film.update( thickness * 0.1, thickness * 0.2, repeatedNormal, repeatedShear );
    EXPECT_NEAR( repeatedNormal, normal, 1.0e-10 );
    EXPECT_NEAR( repeatedShear, shear, 1.0e-10 );
    EXPECT_NEAR( film.state( Keys::equivalentPlasticStrainString() ), kappa, 1.0e-12 );
  }
}

TEST( SurfaceInformedPolymerCohesiveZone, HardeningAndFailureUsePhysicalStretch )
{
  for( Measure const measure : { Measure::Engineering, Measure::Logarithmic } )
  {
    Film film;
    film.select( measure );
    film.model.getReference< real64 >( Keys::defaultYieldStrengthString() ) = 0.0;
    film.model.getReference< real64 >( Keys::strainHardeningSlopeString() ) = 1.0;
    real64 normal, shear;
    film.update( thickness, 0.0, normal, shear );
    real64 const strain = measure == Measure::Logarithmic ? std::log( 2.0 ) : 1.0;
    // Physical stretch is two, so the flow strength is 2^2 - 1/2 in both modes.
    EXPECT_NEAR( -normal - bulkModulus * strain, ( 2.0 / 3.0 ) * 3.5, 1.0e-10 );
    EXPECT_NEAR( film.state( Keys::previousLambdaString() ), 2.0, 1.0e-12 );

    film.model.getReference< real64 >( Keys::maximumStretchString() ) = 2.6;
    film.update( 2.0 * thickness, 0.0, normal, shear );
    EXPECT_NEAR( film.state( Keys::previousLambdaString() ), 3.0, 1.0e-12 );
    EXPECT_NEAR( film.state( Keys::damageString() ), 1.0, 1.0e-12 );
    EXPECT_NEAR( normal, 0.0, 1.0e-12 );
    EXPECT_NEAR( shear, 0.0, 1.0e-12 );
  }
}

TEST( SurfaceInformedPolymerCohesiveZone, LogarithmGuardIsFinite )
{
  Film film;
  film.select( Measure::Logarithmic );
  real64 normal, shear;
  film.update( -thickness, 0.0, normal, shear );
  EXPECT_TRUE( std::isfinite( normal ) );
  EXPECT_NEAR( normal, -constrainedModulus * std::log( 1.0e-16 ), 1.0e-8 );
}

int main( int argc, char ** argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
