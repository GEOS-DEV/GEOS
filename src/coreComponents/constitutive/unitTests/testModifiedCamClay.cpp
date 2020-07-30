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
#include "constitutive/solid/ModifiedCamClay.hpp"

#include "dataRepository/xmlWrapper.hpp"

using namespace geosx;
using namespace ::geosx::constitutive;


// NOTE: using this for debugging, will set up proper unit tests later

TEST( ModifiedCamClayTests, testModel1 )
{
  ConstitutiveManager constitutiveManager( "constitutive", nullptr );

  string const inputStream =
    "<Constitutive>"
    "   <ModifiedCamClay"
    "      name=\"granite\" "
    "      defaultDensity=\"2700.0\" "
    "      defaultRefPInvariant=\"-90000.0\" "
    "      defaultRefElasticStrainVolumetric=\"0.0\" "
    "      defaultRefShearModulus=\"0.0\" "
    "      defaultShearModulusEvolution=\"120.0\" "
    "      defaultVirginCompressionIndex=\"0.13\" "
    "      defaultRecompressionIndex=\"0.018\" "
    "      defaultCriticalStateSlope=\"1.05\" "
    "      defaultAssociativity=\"1.0\" "
    "      defaultPreconsolidationPressure=\"-90000.0\"/>"
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
  
  ModifiedCamClay & cm = *(constitutiveManager.GetConstitutiveRelation<ModifiedCamClay>("granite"));
  cm.AllocateConstitutiveData( &disc, numQuad );
  
  EXPECT_EQ( cm.size(), numElem );
  EXPECT_EQ( cm.numQuadraturePoints(), numQuad );
  
  ModifiedCamClay::KernelWrapper cmw = cm.createKernelWrapper();
  
  real64 totalIncrement = -5.0e-2;
  real64 stepsIncrement = 10;
  localIndex finalLoadStep = stepsIncrement;
  real64 inc = totalIncrement / stepsIncrement;
  real64 total = 0;
    
  array2d< real64 > strainIncrement(1,6);
                    strainIncrement = 0;
                    strainIncrement[0][3] = inc;
                    strainIncrement[0][4] = inc;
                    strainIncrement[0][5] = inc;
  
  array2d< real64 > stress(1,6);
  array3d< real64 > stiffness(1,6,6);
  array2d< real64 > deviator(1,6);

  array1d< real64 > refSolutionP(10);
  array1d< real64 > refSolutionQ(10);
    
  refSolutionP[0] = -76467.89973;
  refSolutionP[1] = -68266.62166;
  refSolutionP[2] = -63159.61333;
  refSolutionP[3] = -59738.81023;
  refSolutionP[4] = -57341.81925;
  refSolutionP[5] = -55611.39243;
  refSolutionP[6] = -54335.95526;
  refSolutionP[7] = -53381.69161;
  refSolutionP[8] = -52659.78902;
  refSolutionP[9] = -52109.11973;

  refSolutionQ[0] = 36950.97166;
  refSolutionQ[1] = 44529.43167;
  refSolutionQ[2] = 47821.15333;
  refSolutionQ[3] = 49556.30214;
  refSolutionQ[4] = 50573.13629;
  refSolutionQ[5] = 51211.87426;
  refSolutionQ[6] = 51633.55619;
  refSolutionQ[7] = 51922.52678;
  refSolutionQ[8] = 52126.32914;
  refSolutionQ[9] = 52273.33077;
  
  real64 relError = 1e-04;
  real64 absError;
  
  for(localIndex loadstep=0; loadstep < finalLoadStep; ++loadstep)
  {
    cmw.SmallStrainUpdate(0,0,strainIncrement[0],stress[0],stiffness[0]);
    cmw.SaveConvergedState(0,0);
    total += inc;
      
    real64 mean = (stress[0][0]+stress[0][1]+stress[0][2])/3.;
    
    for(localIndex i=0; i<3; ++i)
    {
      deviator[0][i] = stress[0][i] - mean;
      deviator[0][i+3] = stress[0][i+3];
    }
  
    real64 invariantQ = 0;
    for(localIndex i=0; i<3; ++i)
    {
      invariantQ += deviator[0][i]*deviator[0][i];
      invariantQ += 2 * deviator[0][i+3]*deviator[0][i+3];
    }
    invariantQ = std::sqrt(invariantQ) + 1e-15; // perturbed to avoid divide by zero when Q=0;
  
    invariantQ *= std::sqrt(3./2.);
    
    absError = std::fabs(relError * refSolutionP[loadstep]);
    EXPECT_NEAR( refSolutionP[loadstep], mean, absError );
    absError = std::fabs(relError * refSolutionQ[loadstep]);
    EXPECT_NEAR( refSolutionQ[loadstep], invariantQ, absError );
  }
}

//Test 2
TEST( ModifiedCamClayTests, testModel2 )
{
  ConstitutiveManager constitutiveManager( "constitutive", nullptr );

  string const inputStream =
    "<Constitutive>"
    "   <ModifiedCamClay"
    "      name=\"granite\" "
    "      defaultDensity=\"2700.0\" "
    "      defaultRefPInvariant=\"-90000.0\" "
    "      defaultRefElasticStrainVolumetric=\"0.0\" "
    "      defaultRefShearModulus=\"0.0\" "
    "      defaultShearModulusEvolution=\"120.0\" "
    "      defaultVirginCompressionIndex=\"0.13\" "
    "      defaultRecompressionIndex=\"0.018\" "
    "      defaultCriticalStateSlope=\"1.05\" "
    "      defaultAssociativity=\"1.0\" "
    "      defaultPreconsolidationPressure=\"-90000.0\"/>"
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
  
  ModifiedCamClay & cm = *(constitutiveManager.GetConstitutiveRelation<ModifiedCamClay>("granite"));
  cm.AllocateConstitutiveData( &disc, numQuad );
  
  EXPECT_EQ( cm.size(), numElem );
  EXPECT_EQ( cm.numQuadraturePoints(), numQuad );
  
  ModifiedCamClay::KernelWrapper cmw = cm.createKernelWrapper();
  
  real64 totalIncrement = -5.0e-2;
  real64 stepsIncrement = 100;
  localIndex finalLoadStep = stepsIncrement;
  real64 inc = totalIncrement / stepsIncrement;
  real64 total = 0;

  array2d< real64 > strainIncrement(1,6);
  strainIncrement = 0;
  strainIncrement[0][0] = inc/3.;
  strainIncrement[0][1] = inc/3.;
  strainIncrement[0][2] = inc/3.;
  strainIncrement[0][3] = inc;
  strainIncrement[0][4] = inc;
  strainIncrement[0][5] = inc;
  
  array2d< real64 > stress(1,6);
  array3d< real64 > stiffness(1,6,6);
  array2d< real64 > deviator(1,6);

  array1d< real64 > refSolutionP(10);
  array1d< real64 > refSolutionQ(10);
  
  refSolutionP[0] = -83528.97804;
  refSolutionP[1] = -83764.55661;
  refSolutionP[2] = -85898.4324;
  refSolutionP[3] = -88788.3358;
  refSolutionP[4] = -92064.26773;
  refSolutionP[5] = -95584.67539;
  refSolutionP[6] = -99293.43385;
  refSolutionP[7] = -103169.6063;
  refSolutionP[8] = -107207.4099;
  refSolutionP[9] = -111407.7736;

  refSolutionQ[0] = 33110.31071;
  refSolutionQ[1] = 38500.48404;
  refSolutionQ[2] = 41380.23559;
  refSolutionQ[3] = 43548.40176;
  refSolutionQ[4] = 45486.51095;
  refSolutionQ[5] = 47369.76996;
  refSolutionQ[6] = 49270.70093;
  refSolutionQ[7] = 51221.72081;
  refSolutionQ[8] = 53238.53964;
  refSolutionQ[9] = 55329.74335;
  
  real64 relError = 1e-04;
  real64 absError;
  localIndex count = 0;
  
  for(localIndex loadstep=0; loadstep < finalLoadStep; ++loadstep)
  {
    cmw.SmallStrainUpdate(0,0,strainIncrement[0],stress[0],stiffness[0]);
    cmw.SaveConvergedState(0,0);
    total += inc;
      
    real64 mean = (stress[0][0]+stress[0][1]+stress[0][2])/3.;
    
    for(localIndex i=0; i<3; ++i)
    {
      deviator[0][i] = stress[0][i] - mean;
      deviator[0][i+3] = stress[0][i+3];
    }
  
    real64 invariantQ = 0;
    for(localIndex i=0; i<3; ++i)
    {
      invariantQ += deviator[0][i]*deviator[0][i];
      invariantQ += 2 * deviator[0][i+3]*deviator[0][i+3];
    }
    invariantQ = std::sqrt(invariantQ) + 1e-15;
    invariantQ *= std::sqrt(3./2.);
    
    int64_t reminder = (loadstep + 1) % 10;
    if (reminder == 0)
    {
      absError = std::fabs(relError * refSolutionP[count]);
      EXPECT_NEAR( refSolutionP[count], mean, absError );
      absError = std::fabs(relError * refSolutionQ[count]);
      EXPECT_NEAR( refSolutionQ[count], invariantQ, absError );
      count++;
    }
  }
  
  // finite-difference check of tangent stiffness
  
  array2d< real64 > fd_stiffness(6,6);
  array2d< real64 > pstress(1,6);
  array3d< real64 > pstiffness(1,6,6);
  
  real64 eps = -1e-12;
  
  cmw.SmallStrainUpdate(0,0,strainIncrement[0],stress[0],stiffness[0]);
    
  for(localIndex i=0; i<6; ++i)
  {
    strainIncrement[0][i] += eps;
    
    if(i>0)
    {
      strainIncrement[0][i-1] -= eps;
    }
    
    cmw.SmallStrainUpdate(0,0,strainIncrement[0],pstress[0],pstiffness[0]);
      
    for(localIndex j=0; j<6; ++j)
    {
      fd_stiffness[j][i] = (pstress[0][j]-stress[0][j])/eps;
    }
  }
    
  for(localIndex i=0; i<6; ++i)
  {
    for(localIndex j=0; j<6; ++j)
    {
      absError = std::fabs(relError * stiffness[0][i][j]);
      EXPECT_NEAR( fd_stiffness[i][j], stiffness[0][i][j], absError );
    }
  }
  
}



