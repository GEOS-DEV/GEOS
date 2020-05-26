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

#include "dataRepository/xmlWrapper.hpp"

using namespace geosx;
using namespace ::geosx::constitutive;

TEST( DruckerPragerTests, testAllocation )
{
  DruckerPrager cm( "model", nullptr );

  localIndex constexpr numElems = 2;
  localIndex constexpr numQuadraturePoints = 4;

  dataRepository::Group disc( "discretization", nullptr );
  disc.resize( numElems );
  cm.AllocateConstitutiveData( &disc, numQuadraturePoints );

  EXPECT_EQ( cm.size(), numElems );
  EXPECT_EQ( cm.numQuadraturePoints(), numQuadraturePoints );

  arrayView1d< real64 const > const & bulkModulus = cm.bulkModulus();
  arrayView1d< real64 const > const & shearModulus = cm.shearModulus();
  arrayView3d< real64 const, solid::STRESS_USD > const & stress = cm.getStress();

  EXPECT_EQ( bulkModulus.size(), numElems );
  EXPECT_EQ( shearModulus.size(), numElems );
  EXPECT_EQ( stress.size( 0 ), numElems );
  EXPECT_EQ( stress.size( 1 ), numQuadraturePoints );
  EXPECT_EQ( stress.size( 2 ), 6 );
}

TEST( DruckerPragerTests, testStateUpdatePoint )
{
/*
  DruckerPrager cm( "model", nullptr );
  real64 constexpr K = 2e10;
  real64 constexpr G = 1e10;
  cm.setDefaultBulkModulus( K );
  cm.setDefaultShearModulus( G );

  dataRepository::Group disc( "discretization", nullptr );
  disc.resize( 2 );
  cm.AllocateConstitutiveData( &disc, 2 );
  DruckerPrager::KernelWrapper cmw = cm.createKernelWrapper();

  arrayView3d< real64, solid::STRESS_USD > const & stress = cm.getStress();

  real64 const strain = 0.1;
  R2SymTensor Ddt;
  R2Tensor Rot;

  {
    Ddt( 0, 0 ) = strain;
    Rot( 0, 0 ) = 1;
    Rot( 1, 1 ) = 1;
    Rot( 2, 2 ) = 1;

    //cmw.SmallStrain(0, 0, Ddt);

    ASSERT_DOUBLE_EQ( stress( 0, 0, 0 ), (2.0/3.0*strain)*2*G + strain*K );
    ASSERT_DOUBLE_EQ( stress( 0, 0, 1 ), (-1.0/3.0*strain)*2*G + strain*K );
    ASSERT_DOUBLE_EQ( stress( 0, 0, 2 ), (-1.0/3.0*strain)*2*G + strain*K );
    ASSERT_DOUBLE_EQ( stress( 0, 0, 3 ), 0.0 );
    ASSERT_DOUBLE_EQ( stress( 0, 0, 4 ), 0.0 );
    ASSERT_DOUBLE_EQ( stress( 0, 0, 5 ), 0.0 );
  }
  */
}

TEST( DruckerPragerTests, testXML )
{
  ConstitutiveManager constitutiveManager( "constitutive", nullptr );
  DruckerPrager cm( "model", &constitutiveManager );

  string const inputStream =
    "<Constitutive>"
    "  <DruckerPrager name=\"granite\" "
    "  defaultDensity=\"2700\" "
    "  defaultBulkModulus=\"5.5556e9\" "
    "  defaultShearModulus=\"4.16667e9\"/>"
    "</Constitutive>";

  xmlWrapper::xmlDocument xmlDocument;
  xmlWrapper::xmlResult xmlResult = xmlDocument.load_buffer( inputStream.c_str(), inputStream.size() );
  if( !xmlResult )
  {
    GEOSX_LOG_RANK_0( "XML parsed with errors!" );
    GEOSX_LOG_RANK_0( "Error description: " << xmlResult.description());
    GEOSX_LOG_RANK_0( "Error offset: " << xmlResult.offset );
  }

  xmlWrapper::xmlNode xmlConstitutiveNode = xmlDocument.child( "Constitutive" );
  constitutiveManager.ProcessInputFileRecursive( xmlConstitutiveNode );
  constitutiveManager.PostProcessInputRecursive();
}
