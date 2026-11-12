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


class JointRoughnessCoefficientTest : public ::testing::Test
{
protected:

  JointRoughnessCoefficientTest():
    m_rootGroup( "root", m_node ),
    m_constitutiveManager( "constitutive", &m_rootGroup )
  {}

  ParallelPlatesPermeability & createModel( string const & jointRoughnessCoefficient )
  {
    string const inputStream =
      "<Constitutive>"
      "   <ParallelPlatesPermeability"
      "      name=\"fracturePerm\" "
      "      jointRoughnessCoefficient=\"" + jointRoughnessCoefficient + "\"/>"
      "</Constitutive>";

    xmlWrapper::xmlResult const xmlResult = m_xmlDocument.loadString( inputStream );
    EXPECT_TRUE( xmlResult ) << xmlResult.description();

    xmlWrapper::xmlNode xmlConstitutiveNode = m_xmlDocument.getChild( "Constitutive" );
    m_constitutiveManager.processInputFileRecursive( m_xmlDocument, xmlConstitutiveNode );
    m_constitutiveManager.postInputInitializationRecursive();

    return m_constitutiveManager.getConstitutiveRelation< ParallelPlatesPermeability >( "fracturePerm" );
  }

  conduit::Node m_node;
  dataRepository::Group m_rootGroup;
  ConstitutiveManager m_constitutiveManager;
  xmlWrapper::xmlDocument m_xmlDocument;
};

// mechanical aperture defined in SurfaceElementRegion
real64 constexpr defaultAperture = 1.0e-4;

TEST_F( JointRoughnessCoefficientTest, testJRC5 )
{
  ParallelPlatesPermeability::KernelWrapper const cmw = createModel( "5.0" ).createKernelWrapper();

  // With JRC = 5 the reduced aperture (1.79e-4) would exceed the input aperture,
  // so the input aperture is kept and the derivative is 1.
  real64 dAperture_dInput = 0.0;
  real64 const aperture = cmw.computeApertureUsingJRC( defaultAperture, dAperture_dInput );

  EXPECT_DOUBLE_EQ( aperture, defaultAperture );
  EXPECT_DOUBLE_EQ( dAperture_dInput, 1.0 );
}


TEST_F( JointRoughnessCoefficientTest, testJRC10 )
{
  ParallelPlatesPermeability::KernelWrapper const cmw = createModel( "10.0" ).createKernelWrapper();

  // With JRC = 10: e_h = E^2 / JRC^2.5 (in microns) = 100^2 / 10^2.5 um = 10^-4.5 m
  real64 const expectedAperture = 3.162277660168379e-05;
  real64 dAperture_dInput = 0.0;
  real64 const aperture = cmw.computeApertureUsingJRC( defaultAperture, dAperture_dInput );

  EXPECT_NEAR( aperture, expectedAperture, 1e-12 * expectedAperture );

  // analytical derivative: 2 e_h / E
  EXPECT_NEAR( dAperture_dInput, 2.0 * expectedAperture / defaultAperture, 1e-12 );

  // finite-difference check of the derivative
  real64 const eps = 1e-6 * defaultAperture;
  real64 const finiteDifference = ( cmw.computeApertureUsingJRC( defaultAperture + eps )
                                    - cmw.computeApertureUsingJRC( defaultAperture - eps ) ) / ( 2.0 * eps );
  EXPECT_NEAR( dAperture_dInput, finiteDifference, 1e-6 * dAperture_dInput );
}

TEST_F( JointRoughnessCoefficientTest, testInvalidJRC )
{
  EXPECT_THROW( createModel( "0.5" ), InputError );
}
