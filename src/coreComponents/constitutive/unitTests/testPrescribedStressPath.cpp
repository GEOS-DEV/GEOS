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
#include "constitutive/contact/BartonBandisStressPathDriven.hpp"
#include "dataRepository/xmlWrapper.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

using namespace geos;
using namespace ::geos::constitutive;


TEST( BartonBandisStressPathDrivenTests, testNewAperture )
{
  conduit::Node node;
  dataRepository::Group rootGroup( "root", node );
  ConstitutiveManager constitutiveManager( "constitutive", &rootGroup );

  real64 constexpr referenceAperture = 1.2421e-4; // in-situ
  std::stringstream ss;
  ss << std::scientific << std::setprecision(4) << referenceAperture;
  
  string const inputStream =
    /*"<ElementRegions>"
    "  <SurfaceElementRegion"
    "    name=\"Fracture\""
    "    faceBlock=\"embeddedSurfaceSubRegion\""
    "    defaultAperture=\"" + ss.str() + "\""
    "    materialList=\"{ hApertureModel }\""
    "    subRegionType=\"embeddedElement\"/>"
    "</ElementRegions>"*/
    "<Constitutive>"
    "   <BartonBandisStressPathDriven"
    "      name=\"hApertureModel\" "
    "      biot=\"1\" "
    "      poisson=\"0.3\" "
    "      normalStiffness=\"12.041e9\" "
    "      referenceAperture=\"" + ss.str() + "\" "
    "      referencePressure=\"1e5\" "
    "      referenceTotalStress=\"{ 85.0e6, 85.0e6, 105.0e6 }\"/>"
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
    
    real64 const newAperture = cmw.computeHydraulicAperture(1e5, normal);
    EXPECT_DOUBLE_EQ( newAperture, referenceAperture );
  }

}
