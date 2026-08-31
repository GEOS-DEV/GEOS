/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include <gtest/gtest.h>

#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/solid/ElasticIsotropicFiniteStrain.hpp"
#include "constitutive/solid/SolidUtilities.hpp"
#include "dataRepository/xmlWrapper.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"

using namespace geos;
using namespace ::geos::constitutive;

TEST( ElasticFiniteStrainTests, testMaterialTangentFiniteDifference )
{
  // create a model and test xml input
  conduit::Node node;
  dataRepository::Group rootGroup( "root", node );
  ConstitutiveManager constitutiveManager( "constitutive", &rootGroup );

  string const inputStream =
    "<Constitutive>"
    "   <ElasticIsotropicFiniteStrain"
    "      name=\"lithium\" "
    "      defaultDensity=\"1.0\" "
    "      defaultBulkModulus=\"1.26666666667e5\" "
    "      defaultShearModulus=\"4.0e4\" "
    "   />"
    "</Constitutive>";

  xmlWrapper::xmlDocument xmlDocument;
  xmlWrapper::xmlResult xmlResult = xmlDocument.loadString( inputStream );
  if( !xmlResult )
  {
    GEOS_LOG_RANK_0( "XML parsed with errors!" );
    GEOS_LOG_RANK_0( "Error description: " << xmlResult.description());
    GEOS_LOG_RANK_0( "Error offset: " << xmlResult.offset );
  }

  xmlWrapper::xmlNode xmlConstitutiveNode = xmlDocument.getChild( "Constitutive" );
  constitutiveManager.processInputFileRecursive( xmlDocument, xmlConstitutiveNode );
  constitutiveManager.postInputInitializationRecursive();

  localIndex constexpr numElem = 1;
  localIndex constexpr numQuad = 1;

  dataRepository::Group disc( "discretization", &rootGroup );
  disc.resize( numElem );

  ElasticIsotropicFiniteStrain & constitutive_model = constitutiveManager.getConstitutiveRelation< ElasticIsotropicFiniteStrain >( "lithium" );
  constitutive_model.allocateConstitutiveData( disc, numQuad );

  // confirm allocation sizes
  EXPECT_EQ( constitutive_model.size(), numElem );
  EXPECT_EQ( constitutive_model.numQuadraturePoints(), numQuad );

  ElasticIsotropicFiniteStrain::KernelWrapper cmw = constitutive_model.createKernelUpdates();
  real64 deformationGradient[3][3] = {};
  deformationGradient[0][0] = 1.05;
  deformationGradient[0][1] = 0.02;
  deformationGradient[0][2] = 0.01;
  deformationGradient[1][0] = 0.0;
  deformationGradient[1][1] = 1.02;
  deformationGradient[1][2] = 0.0;
  deformationGradient[2][0] = 0.0;
  deformationGradient[2][1] = 0.0;
  deformationGradient[2][2] = 1.01;

  EXPECT_TRUE( SolidUtilities::checkFiniteStrainStiffness( cmw, 0, 0, deformationGradient, true ) );

  deformationGradient[0][0] = 1.5;
  deformationGradient[0][1] = 2.0;
  deformationGradient[0][2] = 3.0;
  deformationGradient[1][0] = 1.2;
  deformationGradient[1][1] = 1.3;
  deformationGradient[1][2] = 4.0;
  deformationGradient[2][0] = 2.5;
  deformationGradient[2][1] = 1.2;
  deformationGradient[2][2] = 1.4;

  EXPECT_TRUE( SolidUtilities::checkFiniteStrainStiffness( cmw, 0, 0, deformationGradient, true ) );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::GeosxState state( geos::basicSetup( argc, argv ) );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
