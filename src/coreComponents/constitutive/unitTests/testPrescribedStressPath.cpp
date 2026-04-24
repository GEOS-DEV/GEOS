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
#include "constitutive/permeability/BartonBandisPermeability.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "dataRepository/xmlWrapper.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

using namespace geos;
using namespace ::geos::constitutive;


TEST( BartonBandisPermeabilityTests, testNewAperture )
{
  conduit::Node node;
  dataRepository::Group rootGroup( "root", node );
  ConstitutiveManager constitutiveManager( "constitutive", &rootGroup );

  real64 constexpr referenceAperture = 1.2421e-4; // in-situ
  real64 constexpr referencePressure = 1.0e5; 
  real64 constexpr fractureStress = 1.049e+08; // only if Z component of normal is 1.0 and referenceTotalStress in Z is 105.0e6
  std::stringstream ss_aperture;
  ss_aperture << std::scientific << std::setprecision(4) << referenceAperture;
  std::stringstream ss_pressure;
  ss_pressure << std::scientific << std::setprecision(4) << referencePressure;
  
  string const inputStream =
    "<Constitutive>"
    "   <BartonBandisPermeability"
    "      name=\"fracturePerm\" "
    "      biot=\"1\" "
    "      poisson=\"0.3\" "
    "      normalStiffness=\"12.041e9\" "
    "      referenceAperture=\"" + ss_aperture.str() + "\" "
    "      referencePressure=\"" + ss_pressure.str() + "\" "
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

  BartonBandisPermeability & cm = constitutiveManager.getConstitutiveRelation< BartonBandisPermeability >( "fracturePerm" );

  BartonBandisPermeability::KernelWrapper cmw = cm.createKernelWrapper();

  {    
    array1d < real64 > normal(3);
    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 1.0;
    //real64 const normal[ 3 ] = {0.0, 0.0, 1.0};
    

    real64 dAperture_dStress = -1.0;
    real64 newHydraulicAperture = cmw.computeHydraulicAperture(referencePressure, fractureStress, normal, dAperture_dStress, 0);
    
    EXPECT_DOUBLE_EQ( newHydraulicAperture, referenceAperture );
  }
}

