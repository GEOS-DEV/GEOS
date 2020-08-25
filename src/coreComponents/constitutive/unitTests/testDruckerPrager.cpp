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
#include "constitutive/solid/InvariantDecompositions.hpp"

#include "dataRepository/xmlWrapper.hpp"

using namespace geosx;
using namespace ::geosx::constitutive;

TEST( DruckerPragerTests, testModel )
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
    "      defaultFrictionAngle=\"" + std::to_string(friction)+ "\" "
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
  constitutiveManager.ProcessInputFileRecursive( xmlConstitutiveNode );
  constitutiveManager.PostProcessInputRecursive();
  
  localIndex constexpr numElem = 2;
  localIndex constexpr numQuad = 4;
  
  dataRepository::Group disc( "discretization", nullptr );
  disc.resize( numElem );
  
  DruckerPrager & cm = *(constitutiveManager.GetConstitutiveRelation<DruckerPrager>("granite"));
  cm.allocateConstitutiveData( &disc, numQuad );
  
  // confirm allocation sizes
  
  EXPECT_EQ( cm.size(), numElem );
  EXPECT_EQ( cm.numQuadraturePoints(), numQuad );
  
  // use updates class to run a uniaxial compression test.  the material
  // begins with some cohesion, but it quickly degrades to zero under
  // plastic loading.  at the end of the loading, we confirm the
  // stress point lies on the correctly-positioned yield surface.
  
  DruckerPrager::KernelWrapper cmw = cm.createKernelWrapper();
  
  real64 strainIncrement[6] = {-1e-4,0,0,0,0,0};
  real64 stress[6];
  real64 stiffness[6][6];
  
  for(localIndex loadstep=0; loadstep < 30; ++loadstep)
  {
    cmw.smallStrainUpdate(0,0,strainIncrement,stress,stiffness);
    cmw.saveConvergedState();
  }
  
  // loading was set up to drive to total cohesion loss (c=0), at which
  // point the Q/P ratio should equal the slope of the DP yield surface:
  
  real64 invariantP, invariantQ;
  real64 deviator[6];
  
  twoInvariant::stressDecomposition(stress,
                                    invariantP,
                                    invariantQ,
                                    deviator);

  real64 phi = friction * M_PI / 180;
  real64 slope = -6 * cos(phi) / ( 3 + sin(phi) );
  
  EXPECT_DOUBLE_EQ( invariantQ / invariantP , slope );

    
  // we now use a finite-difference check of tangent stiffness to confirm
  // the analytical form is working properly.
  
  real64 fd_stiffness[6][6]; // finite difference approximation
  real64 pstress[6];         // perturbed stress
  real64 pstiffness[6][6];   // perturbed stiffness
  real64 eps = 1e-10;        // finite difference perturbation
  
  cmw.smallStrainUpdate(0,0,strainIncrement,stress,stiffness);
  
  for(localIndex i=0; i<6; ++i)
  {
    strainIncrement[i] += eps;
    
    if(i>0)
    {
      strainIncrement[i-1] -= eps;
    }
    
    cmw.smallStrainUpdate(0,0,strainIncrement,pstress,pstiffness);
    
    for(localIndex j=0; j<6; ++j)
    {
      fd_stiffness[j][i] = (pstress[j]-stress[j])/eps;
    }
  }
  
  // confirm entrywise differences are less than a desired tolerance
  
  for(localIndex i=0; i<6; ++i)
  {
    for(localIndex j=0; j<6; ++j)
    {
      EXPECT_TRUE( fabs(fd_stiffness[i][j]-stiffness[i][j]) < 1e-3 );
    }
  }
}



