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
#include "constitutive/solid/DruckerPrager.hpp"
#include "constitutive/solid/DruckerPragerExtended.hpp"
#include "constitutive/solid/InvariantDecompositions.hpp"

#include "dataRepository/xmlWrapper.hpp"

using namespace geosx;
using namespace ::geosx::constitutive;

TEST( DruckerPragerTests, testDruckerPrager )
{
  // create a Drucker-Prager model, and test xml input

  ConstitutiveManager constitutiveManager( "constitutive", nullptr );

  real64 const friction = 30.0; // will use later in checks

  string const inputStream =
    "<Constitutive>"
    "   <DruckerPrager"
    "      name=\"granite\" "
    "      defaultDensity=\"2700\" "
    "      defaultBulkModulus=\"1000.0\" "
    "      defaultShearModulus=\"1000.0\" "
    "      defaultFrictionAngle=\"" + std::to_string( friction )+ "\" "
                                                                  "      defaultDilationAngle=\"15.0\" "
                                                                  "      defaultHardeningRate=\"-4000.0\" "
                                                                  "      defaultCohesion=\"1\"/>"
                                                                  "</Constitutive>";

  xmlWrapper::xmlDocument xmlDocument;
  xmlWrapper::xmlResult xmlResult = xmlDocument.load_buffer( inputStream.c_str(),
                                                             inputStream.size() );
  if( !xmlResult )
  {
    GEOSX_LOG_RANK_0( "XML parsed with errors!" );
    GEOSX_LOG_RANK_0( "Error description: " << xmlResult.description());
    GEOSX_LOG_RANK_0( "Error offset: " << xmlResult.offset );
  }

  xmlWrapper::xmlNode xmlConstitutiveNode = xmlDocument.child( "Constitutive" );
  constitutiveManager.processInputFileRecursive( xmlConstitutiveNode );
  constitutiveManager.postProcessInputRecursive();

  localIndex constexpr numElem = 2;
  localIndex constexpr numQuad = 4;

  dataRepository::Group disc( "discretization", nullptr );
  disc.resize( numElem );

  DruckerPrager & cm = *(constitutiveManager.getConstitutiveRelation< DruckerPrager >( "granite" ));
  cm.allocateConstitutiveData( &disc, numQuad );

  // confirm allocation sizes

  EXPECT_EQ( cm.size(), numElem );
  EXPECT_EQ( cm.numQuadraturePoints(), numQuad );

  // use updates class to run a uniaxial compression test.  the material
  // begins with some cohesion, but it quickly degrades to zero under
  // plastic loading.  at the end of the loading, we confirm the
  // stress point lies on the correctly-positioned yield surface.

  DruckerPrager::KernelWrapper cmw = cm.createKernelUpdates();

  real64 strainIncrement[6] = {-1e-4, 0, 0, 0, 0, 0};
  real64 stress[6];
  real64 stiffness[6][6];

  for( localIndex loadstep=0; loadstep < 30; ++loadstep )
  {
    cmw.smallStrainUpdate( 0, 0, strainIncrement, stress, stiffness );
    cm.saveConvergedState();
  }

  // loading was set up to drive to total cohesion loss (c=0), at which
  // point the Q/P ratio should equal the slope of the DP yield surface:

  real64 invariantP, invariantQ;
  real64 deviator[6];

  twoInvariant::stressDecomposition( stress,
                                     invariantP,
                                     invariantQ,
                                     deviator );

  real64 phi = friction * M_PI / 180;
  real64 slope = -6 * sin( phi ) / ( 3 - sin( phi ) );
  EXPECT_TRUE( fabs( invariantQ / invariantP / slope - 1 ) < 1e-8 );

  // we now use a finite-difference check of tangent stiffness to confirm
  // the analytical form is working properly.

  EXPECT_TRUE( cmw.checkSmallStrainStiffness( 0, 0, strainIncrement ) );
}



TEST( DruckerPragerTests, testDruckerPragerExtended )
{
  // create an Extended-Drucker-Prager model, and test xml input

  ConstitutiveManager constitutiveManager( "constitutive", nullptr );

  real64 const initialFriction = 15.27;   // will use later in checks
  real64 const residualFriction = 23.05;
  real64 const cohesion = 1.0;

  string const inputStream =
    "<Constitutive>"
    "   <ExtendedDruckerPrager"
    "      name=\"granite\" "
    "      defaultDensity=\"2700\" "
    "      defaultBulkModulus=\"300\" "
    "      defaultShearModulus=\"300\" "
    "      defaultInitialFrictionAngle=\"" + std::to_string( initialFriction )+ "\" "
                                                                                "      defaultResidualFrictionAngle=\"" + std::to_string( residualFriction )+ "\" "
                                                                                                                                                              "      defaultDilationRatio=\"0.5\" "
                                                                                                                                                              "      defaultHardening=\"0.001\" "
                                                                                                                                                              "      defaultCohesion=\"" +
    std::to_string( cohesion )+ "\" "
                                "   />"
                                "</Constitutive>";

  xmlWrapper::xmlDocument xmlDocument;
  xmlWrapper::xmlResult xmlResult = xmlDocument.load_buffer( inputStream.c_str(),
                                                             inputStream.size() );
  if( !xmlResult )
  {
    GEOSX_LOG_RANK_0( "XML parsed with errors!" );
    GEOSX_LOG_RANK_0( "Error description: " << xmlResult.description());
    GEOSX_LOG_RANK_0( "Error offset: " << xmlResult.offset );
  }

  xmlWrapper::xmlNode xmlConstitutiveNode = xmlDocument.child( "Constitutive" );
  constitutiveManager.processInputFileRecursive( xmlConstitutiveNode );
  constitutiveManager.postProcessInputRecursive();

  localIndex constexpr numElem = 2;
  localIndex constexpr numQuad = 4;

  dataRepository::Group disc( "discretization", nullptr );
  disc.resize( numElem );

  DruckerPragerExtended & cm = *(constitutiveManager.getConstitutiveRelation< DruckerPragerExtended >( "granite" ));
  cm.allocateConstitutiveData( &disc, numQuad );

  // confirm allocation sizes

  EXPECT_EQ( cm.size(), numElem );
  EXPECT_EQ( cm.numQuadraturePoints(), numQuad );

  // use updates class to run a uniaxial compression test.  the material
  // begins with some cohesion, but it quickly degrades to zero under
  // plastic loading.  at the end of the loading, we confirm the
  // stress point lies on the correctly-positioned yield surface.

  DruckerPragerExtended::KernelWrapper cmw = cm.createKernelUpdates();

  real64 strainIncrement[6] = {-1e-3, 0, 0, 0, 0, 0};
  real64 stress[6];
  real64 stiffness[6][6];
  real64 invariantP, invariantQ;
  real64 deviator[6];
  real64 total = 0;

  //FILE* fp = fopen("pq.txt","w");
  for( localIndex loadstep=0; loadstep < 300; ++loadstep )
  {
    total += strainIncrement[0];
    cmw.smallStrainUpdate( 0, 0, strainIncrement, stress, stiffness );
    cm.saveConvergedState();

    twoInvariant::stressDecomposition( stress,
                                       invariantP,
                                       invariantQ,
                                       deviator );
    //fprintf(fp,"%.4e %.4e %.4e\n",invariantP,invariantQ,total);
  }
  //fclose(fp);

  // loading was set up to drive to residual strength, at which
  // point the Q/(P-P0) ratio should equal the slope of residual yield surface:

  real64 phi_r = residualFriction * M_PI / 180;
  real64 slope_r = 6 * sin( phi_r ) / ( 3 - sin( phi_r ) );
  real64 phi_i = initialFriction * M_PI / 180;
  real64 slope_i = 6 * sin( phi_i ) / ( 3 - sin( phi_i ) );
  real64 intercept = 6 * cohesion * cos( phi_i ) / ( 3 - sin( phi_i ) ) / slope_i;

  EXPECT_TRUE( fabs( invariantQ / (invariantP-intercept) / slope_r + 1 ) < 1e-2 );

  // we now use a finite-difference check of tangent stiffness to confirm
  // the analytical form is working properly.

  EXPECT_TRUE( cmw.checkSmallStrainStiffness( 0, 0, strainIncrement ) );
}
