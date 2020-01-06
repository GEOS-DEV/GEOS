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
#include "constitutive/solid/LinearElasticAnisotropic.hpp"

#include "dataRepository/xmlWrapper.hpp"
using namespace geosx;
using namespace ::geosx::constitutive;

TEST( LinearElasticAnisotropicTests, testAllocation )
{
  LinearElasticAnisotropic cm( "model", nullptr );

  localIndex constexpr numElems = 2;
  localIndex constexpr numQuadraturePoints = 3;

  dataRepository::Group disc( "discretization", nullptr );
  disc.resize(numElems);
  cm.AllocateConstitutiveData( &disc, numQuadraturePoints );

  EXPECT_EQ( cm.size(), numElems );
  EXPECT_EQ( cm.numQuadraturePoints(), numQuadraturePoints );

//  arrayView1d<LinearElasticAnisotropic::StiffnessTensor const> const &
//  stiffness = cm.stiffness() ;

  arrayView2d<R2SymTensor const> const & stress = cm.getStress();

//  EXPECT_EQ( stiffness.size(), numElems );
  EXPECT_EQ( stress.size(0), numElems );
  EXPECT_EQ( stress.size(1), numQuadraturePoints );

}

void stressCalc( real64 const c[6][6], R2SymTensor const Ddt, real64 stressVoigt[6] )
{
  real64 const DdtVoigt[6] = { Ddt(0,0), Ddt(1,1), Ddt(2,2), 2*Ddt(1,2), 2*Ddt(0,2), 2*Ddt(0,1) };

  for( int i=0 ; i<6 ; ++i )
  {
    stressVoigt[i] = 0.0;
    for( int j=0 ; j<6 ; ++j )
    {
      stressVoigt[i] += c[i][j] * DdtVoigt[j];
    }
  }
}

void stressCheck( R2SymTensor const & stress, real64 const stressV[6] )
{
  ASSERT_DOUBLE_EQ( stress(0,0) , stressV[0] );
  ASSERT_DOUBLE_EQ( stress(1,1) , stressV[1] );
  ASSERT_DOUBLE_EQ( stress(2,2) , stressV[2] );
  ASSERT_DOUBLE_EQ( stress(1,2) , stressV[3] );
  ASSERT_DOUBLE_EQ( stress(0,2) , stressV[4] );
  ASSERT_DOUBLE_EQ( stress(0,1) , stressV[5] );
}

TEST( LinearElasticAnisotropicTests, testStateUpdatePoint )
{
  LinearElasticAnisotropic cm( "model", nullptr );

  LinearElasticAnisotropic::StiffnessTensor c { { { 1.0e11, 0.1e10, 0.2e10, 0.3e10, 0.4e10, 0.5e10 },
                                                  { 1.0e10, 1.1e11, 1.2e10, 1.3e10, 1.4e10, 1.5e10 },
                                                  { 2.0e10, 2.1e10, 2.2e11, 2.3e10, 2.4e10, 2.5e10 },
                                                  { 3.0e10, 3.1e10, 3.2e10, 3.3e10, 3.4e10, 3.5e10 },
                                                  { 4.0e10, 4.1e10, 4.2e10, 4.3e10, 4.4e10, 4.5e10 },
                                                  { 5.0e10, 5.1e10, 5.2e10, 5.3e10, 5.4e10, 5.5e10 }
                                              } };

  cm.setDefaultStiffness( c );

  dataRepository::Group disc( "discretization", nullptr );
  disc.resize(2);
  cm.AllocateConstitutiveData( &disc, 2 );

//  arrayView1d<LinearElasticAnisotropic::StiffnessTensor> const &
//  stiffness = cm.stiffness();


  arrayView2d<R2SymTensor> const & stress = cm.getStress();

  real64 const strain = 0.1;
  R2SymTensor Ddt;
  real64 stressV[6] = {0.0};
  R2Tensor Rot;
  R2SymTensor zero;

  {
    Ddt(0,0) = strain;
    Rot(0,0) = 1;
    Rot(1,1) = 1;
    Rot(2,2) = 1;

    cm.StateUpdatePoint( 0, 0, Ddt, Rot, 0 );
    stressCalc( c.m_data, Ddt, stressV );
    stressCheck( stress[0][0], stressV );
  }

  {
    stress = zero;
    Ddt = 0;

    Ddt(1,1) = strain;
    Rot(0,0) = 1;
    Rot(1,1) = 1;
    Rot(2,2) = 1;

    cm.StateUpdatePoint( 0, 0, Ddt, Rot, 0 );
    stressCalc( c.m_data, Ddt, stressV );
    stressCheck( stress[0][0], stressV );

  }

  {
    stress = zero;
    Ddt = 0;

    Ddt(2,2) = strain;
    Rot(0,0) = 1;
    Rot(1,1) = 1;
    Rot(2,2) = 1;

    cm.StateUpdatePoint( 0, 0, Ddt, Rot, 0 );
    stressCalc( c.m_data, Ddt, stressV );
    stressCheck( stress[0][0], stressV );

  }

  {
    stress = zero;
    Ddt = 0;

    Ddt(0,1) = strain;
    Rot(0,0) = 1;
    Rot(1,1) = 1;
    Rot(2,2) = 1;

    cm.StateUpdatePoint( 0, 0, Ddt, Rot, 0 );
    stressCalc( c.m_data, Ddt, stressV );
    stressCheck( stress[0][0], stressV );

  }

  {
    stress = zero;
    Ddt = 0;

    Ddt(0,2) = strain;
    Rot(0,0) = 1;
    Rot(1,1) = 1;
    Rot(2,2) = 1;

    cm.StateUpdatePoint( 0, 0, Ddt, Rot, 0 );
    stressCalc( c.m_data, Ddt, stressV );
    stressCheck( stress[0][0], stressV );

  }

  {
    stress = zero;
    Ddt = 0;

    Ddt(1,2) = strain;
    Rot(0,0) = 1;
    Rot(1,1) = 1;
    Rot(2,2) = 1;

    cm.StateUpdatePoint( 0, 0, Ddt, Rot, 0 );
    stressCalc( c.m_data, Ddt, stressV );
    stressCheck( stress[0][0], stressV );

  }
}



TEST( LinearElasticAnisotropicTests, testXML )
{
  ConstitutiveManager constitutiveManager("constitutive",nullptr);

  string const inputStream =
  "<Constitutive>"
  "  <LinearElasticAnisotropic name=\"granite\""
  "                            defaultDensity=\"2700\""
  "                            defaultC11=\"1.0e10\" defaultC12=\"1.1e9\"  defaultC13=\"1.2e9\"  defaultC14=\"1.3e9\" defaultC15=\"1.4e9\" defaultC16=\"1.5e9\""
  "                            defaultC21=\"2.0e9\"  defaultC22=\"2.1e10\" defaultC23=\"2.2e9\"  defaultC24=\"2.3e9\" defaultC25=\"2.4e9\" defaultC26=\"2.5e9\""
  "                            defaultC31=\"3.0e9\"  defaultC32=\"3.1e9\"  defaultC33=\"3.2e10\" defaultC34=\"3.3e9\" defaultC35=\"3.4e9\" defaultC36=\"3.5e9\""
  "                            defaultC41=\"4.0e9\"  defaultC42=\"4.1e9\"  defaultC43=\"4.2e9\"  defaultC44=\"4.3e9\" defaultC45=\"4.4e9\" defaultC46=\"4.5e9\""
  "                            defaultC51=\"5.0e9\"  defaultC52=\"5.1e9\"  defaultC53=\"5.2e9\"  defaultC54=\"5.3e9\" defaultC55=\"5.4e9\" defaultC56=\"5.5e9\""
  "                            defaultC61=\"6.0e9\"  defaultC62=\"6.1e9\"  defaultC63=\"6.2e9\"  defaultC64=\"6.3e9\" defaultC65=\"6.4e9\" defaultC66=\"6.5e9\" />"
  "</Constitutive>";

  xmlWrapper::xmlDocument xmlDocument;
  xmlWrapper::xmlResult xmlResult = xmlDocument.load_buffer( inputStream.c_str(), inputStream.size() );
  if (!xmlResult)
  {
    GEOSX_LOG_RANK_0("XML parsed with errors!");
    GEOSX_LOG_RANK_0("Error description: " << xmlResult.description());
    GEOSX_LOG_RANK_0("Error offset: " << xmlResult.offset);
  }

  xmlWrapper::xmlNode xmlConstitutiveNode = xmlDocument.child("Constitutive");
  constitutiveManager.ProcessInputFileRecursive( xmlConstitutiveNode );
  constitutiveManager.PostProcessInputRecursive();

  LinearElasticAnisotropic * const model = constitutiveManager.GetConstitutiveRelation<LinearElasticAnisotropic>("granite");
  dataRepository::Group disc( "discretization", nullptr );
  disc.resize(1);
  model->AllocateConstitutiveData( &disc, 1 );


  LinearElasticAnisotropic::KernelWrapper kernelWrapper = model->createKernelWrapper();

  real64 stiffness[6][6];
  kernelWrapper.GetStiffness(0,stiffness);

  real64 c[6][6] = { { 1.0e10, 1.1e9,  1.2e9,  1.3e9,  1.4e9,  1.5e9 },
                     { 2.0e9,  2.1e10, 2.2e9,  2.3e9,  2.4e9,  2.5e9 },
                     { 3.0e9,  3.1e9,  3.2e10, 3.3e9,  3.4e9,  3.5e9 },
                     { 4.0e9,  4.1e9,  4.2e9,  4.3e9,  4.4e9,  4.5e9 },
                     { 5.0e9,  5.1e9,  5.2e9,  5.3e9,  5.4e9,  5.5e9 },
                     { 6.0e9,  6.1e9,  6.2e9,  6.3e9,  6.4e9,  6.5e9 }
                   };

  for( int i=0 ; i<6 ; ++i )
  {
    for( int j=0 ; j<6 ; ++j )
    {
      ASSERT_DOUBLE_EQ( stiffness[i][j] , c[i][j] );
    }
  }

}
