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

// Source includes
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/permeability/ParallelPlatesPermeability.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "dataRepository/xmlWrapper.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

using namespace geos;
using namespace ::geos::constitutive;


TEST( JointRoughnessCoefficientTests, testJRC5 )
{
  conduit::Node node;
  dataRepository::Group rootGroup( "root", node );
  ConstitutiveManager constitutiveManager( "constitutive", &rootGroup );

  real64 constexpr defaultAperture = 1.0e-4; // mechanical aperture defined in SurfaceElementRegion

  string const inputStream =
    "<Constitutive>"
    "   <ParallelPlatesPermeability"
    "      name=\"fracturePerm\" "
    "      jointRoughnessCoefficient=\"5.0\"/>"
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

  ParallelPlatesPermeability & cm = constitutiveManager.getConstitutiveRelation< ParallelPlatesPermeability >( "fracturePerm" );

  ParallelPlatesPermeability::KernelWrapper cmw = cm.createKernelWrapper();

  {    
    // When JRC = 5 in this case, the new aperture should be in the invalid region. 
    // So the new aperture should be equal to the default aperture.
    real64 newHydraulicAperture = cmw.computeApertureUsingJRC( defaultAperture );

    EXPECT_DOUBLE_EQ( newHydraulicAperture, defaultAperture ); 
  }
}


TEST( JointRoughnessCoefficientTests, testJRC10 )
{
  conduit::Node node;
  dataRepository::Group rootGroup( "root", node );
  ConstitutiveManager constitutiveManager( "constitutive", &rootGroup );

  real64 constexpr defaultAperture = 1.0e-4; // mechanical aperture defined in SurfaceElementRegion

  string const inputStream =
    "<Constitutive>"
    "   <ParallelPlatesPermeability"
    "      name=\"fracturePerm\" "
    "      jointRoughnessCoefficient=\"10.0\"/>"
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

  ParallelPlatesPermeability & cm = constitutiveManager.getConstitutiveRelation< ParallelPlatesPermeability >( "fracturePerm" );

  ParallelPlatesPermeability::KernelWrapper cmw = cm.createKernelWrapper();

  {    
    // When JRC = 5 in this case, the new aperture should be in the invalid region. 
    // So the new aperture should be equal to the default aperture.
    real64 newHydraulicAperture = cmw.computeApertureUsingJRC( defaultAperture );

    EXPECT_DOUBLE_EQ( newHydraulicAperture, 3.162277660168379e-05 ); 
  }
}
/*
TEST( BartonBandisStressPathDrivenTests, testPressure )
{
  conduit::Node node;
  dataRepository::Group rootGroup( "root", node );
  ConstitutiveManager constitutiveManager( "constitutive", &rootGroup );

  real64 constexpr referenceAperture = 1.0e-4; // in-situ
  std::stringstream ss;
  ss << std::scientific << std::setprecision(4) << referenceAperture;
  
  string const inputStream =
    "<Constitutive>"
    "   <BartonBandisStressPathDriven"
    "      name=\"hApertureModel\" "
    "      biot=\"1.0\" "
    "      poisson=\"0.3\" "
    "      normalStiffness=\"10.0e9\" "
    "      referenceAperture=\"" + ss.str() + "\" "
    "      referencePressure=\"1e5\" "
    "      referenceTotalStress=\"{ 40.0e6, 40.0e6, 20.0e6 }\"/>"
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

  BartonBandisStressPathDriven & cm = constitutiveManager.getConstitutiveRelation< BartonBandisStressPathDriven >( "hApertureModel" );

  BartonBandisStressPathDriven::KernelWrapper cmw = cm.createKernelWrapper();

  {    
    array1d < real64 > normal(3);
    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 1.0;
    
    real64 const newAperture = cmw.computeHydraulicAperture(13842265.4230671, normal);
    EXPECT_DOUBLE_EQ( newAperture, 0.00022328650392488669 );
  }

}
  */
