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

TEST( DruckerPragerTests, testModel )
{
  ConstitutiveManager constitutiveManager( "constitutive", nullptr );
  DruckerPrager cm( "model", &constitutiveManager );

  string const inputStream =
    "<Constitutive>"
    "   <DruckerPrager"
    "      name=\"granite\" "
    "      defaultDensity=\"2700\" "
    "      defaultBulkModulus=\"5e9\" "
    "      defaultShearModulus=\"5e9\" "
    "      defaultTanFrictionAngle=\"1.0\" "
    "      defaultHardeningRate=\"0.0\" "
    "      defaultCohesion=\"1.0e6\"/>"
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
  
  localIndex constexpr numElem = 2;
  localIndex constexpr numQuad = 4;
  
  dataRepository::Group disc( "discretization", nullptr );
  disc.resize( numElem );
  cm.AllocateConstitutiveData( &disc, numQuad );
  
  EXPECT_EQ( cm.size(), numElem );
  EXPECT_EQ( cm.numQuadraturePoints(), numQuad );
  
  DruckerPrager::KernelWrapper cmw = cm.createKernelWrapper();
  
  array2d< real64 > strainIncrement(1,6);
                    strainIncrement = 0;
                    strainIncrement[0][0] = 1e-4;
  
  array2d< real64 > stress(1,6);
  array3d< real64 > stiffness(1,6,6);
  
  for(localIndex loadstep=0; loadstep < 10; ++loadstep)
  {
    cmw.SmallStrainUpdate(0,0,strainIncrement[0],stress[0],stiffness[0]);
    cmw.SaveConvergedState(0,0);
  }

}



