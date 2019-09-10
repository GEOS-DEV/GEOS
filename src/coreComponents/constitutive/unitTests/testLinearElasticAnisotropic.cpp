/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
#include "gtest/gtest.h"

#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/Solid/LinearElasticAnisotropic.hpp"

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

  arrayView1d<LinearElasticAnisotropic::StiffnessTensor const> const &
  stiffness = cm.stiffness() ;

  arrayView2d<real64 const>      const & meanStress = cm.meanStress();

  arrayView2d<R2SymTensor const> const & deviatorStress = cm.deviatorStress();

  EXPECT_EQ( stiffness.size(), numElems );
  EXPECT_EQ( meanStress.size(0), numElems );
  EXPECT_EQ( meanStress.size(1), numQuadraturePoints );
  EXPECT_EQ( deviatorStress.size(0), numElems );
  EXPECT_EQ( deviatorStress.size(1), numQuadraturePoints );

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

void stressCheck( real64 const meanStress, R2SymTensor const & deviatorStress, real64 const stressV[6] )
{
  real64 p = ( stressV[0] + stressV[1] + stressV[2] )/3.0;
  ASSERT_DOUBLE_EQ( meanStress , p );
  ASSERT_DOUBLE_EQ( deviatorStress(0,0) , stressV[0] - p );
  ASSERT_DOUBLE_EQ( deviatorStress(1,1) , stressV[1] - p );
  ASSERT_DOUBLE_EQ( deviatorStress(2,2) , stressV[2] - p );
  ASSERT_DOUBLE_EQ( deviatorStress(1,2) , stressV[3] );
  ASSERT_DOUBLE_EQ( deviatorStress(0,2) , stressV[4] );
  ASSERT_DOUBLE_EQ( deviatorStress(0,1) , stressV[5] );
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

  dataRepository::Group disc( "discretization", nullptr );
  disc.resize(2);
  cm.AllocateConstitutiveData( &disc, 2 );

  arrayView1d<LinearElasticAnisotropic::StiffnessTensor> const &
  stiffness = cm.stiffness();

  stiffness[0] = c;

  arrayView2d<real64>      const & meanStress = cm.meanStress();
  arrayView2d<R2SymTensor> const & deviatorStress = cm.deviatorStress();

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
    stressCheck( meanStress[0][0], deviatorStress[0][0], stressV );
  }

  {
    meanStress = 0.0;
    deviatorStress = zero;
    Ddt = 0;

    Ddt(1,1) = strain;
    Rot(0,0) = 1;
    Rot(1,1) = 1;
    Rot(2,2) = 1;

    cm.StateUpdatePoint( 0, 0, Ddt, Rot, 0 );
    stressCalc( c.m_data, Ddt, stressV );
    stressCheck( meanStress[0][0], deviatorStress[0][0], stressV );

  }

  {
    meanStress = 0.0;
    deviatorStress = zero;
    Ddt = 0;

    Ddt(2,2) = strain;
    Rot(0,0) = 1;
    Rot(1,1) = 1;
    Rot(2,2) = 1;

    cm.StateUpdatePoint( 0, 0, Ddt, Rot, 0 );
    stressCalc( c.m_data, Ddt, stressV );
    stressCheck( meanStress[0][0], deviatorStress[0][0], stressV );

  }

  {
    meanStress = 0.0;
    deviatorStress = zero;
    Ddt = 0;

    Ddt(0,1) = strain;
    Rot(0,0) = 1;
    Rot(1,1) = 1;
    Rot(2,2) = 1;

    cm.StateUpdatePoint( 0, 0, Ddt, Rot, 0 );
    stressCalc( c.m_data, Ddt, stressV );
    stressCheck( meanStress[0][0], deviatorStress[0][0], stressV );

  }

  {
    meanStress = 0.0;
    deviatorStress = zero;
    Ddt = 0;

    Ddt(0,2) = strain;
    Rot(0,0) = 1;
    Rot(1,1) = 1;
    Rot(2,2) = 1;

    cm.StateUpdatePoint( 0, 0, Ddt, Rot, 0 );
    stressCalc( c.m_data, Ddt, stressV );
    stressCheck( meanStress[0][0], deviatorStress[0][0], stressV );

  }

  {
    meanStress = 0.0;
    deviatorStress = zero;
    Ddt = 0;

    Ddt(1,2) = strain;
    Rot(0,0) = 1;
    Rot(1,1) = 1;
    Rot(2,2) = 1;

    cm.StateUpdatePoint( 0, 0, Ddt, Rot, 0 );
    stressCalc( c.m_data, Ddt, stressV );
    stressCheck( meanStress[0][0], deviatorStress[0][0], stressV );

  }
}



TEST( LinearElasticAnisotropicTests, testXML )
{
  ConstitutiveManager constitutiveManager("constitutive",nullptr);

  string const inputStream =
  "<?xml version=\"1.0\" ?>"
  "<Constitutive xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"geos_v0.0.xsd\">"
  "  <LinearElasticAnisotropic name=\"granite\""
  "                            defaultDensity=\"2700\""
  "                            c11=\"1.0e10\" c12=\"1.1e9\"  c13=\"1.2e9\"  c14=\"1.3e9\" c15=\"1.4e9\" c16=\"1.5e9\""
  "                            c21=\"2.0e9\"  c22=\"2.1e10\" c23=\"2.2e9\"  c24=\"2.3e9\" c25=\"2.4e9\" c26=\"2.5e9\""
  "                            c31=\"3.0e9\"  c32=\"3.1e9\"  c33=\"3.2e10\" c34=\"3.3e9\" c35=\"3.4e9\" c36=\"3.5e9\""
  "                            c41=\"4.0e9\"  c42=\"4.1e9\"  c43=\"4.2e9\"  c44=\"4.3e9\" c45=\"4.4e9\" c46=\"4.5e9\""
  "                            c51=\"5.0e9\"  c52=\"5.1e9\"  c53=\"5.2e9\"  c54=\"5.3e9\" c55=\"5.4e9\" c56=\"5.5e9\""
  "                            c61=\"6.0e9\"  c62=\"6.1e9\"  c63=\"6.2e9\"  c64=\"6.3e9\" c65=\"6.4e9\" c66=\"6.5e9\" />"
  "</Constitutive>";

  xmlWrapper::xmlDocument xmlDocument;
  xmlWrapper::xmlResult xmlResult = xmlDocument.load_buffer( inputStream.c_str(), inputStream.size() );
  if (!xmlResult)
  {
    GEOS_LOG_RANK_0("XML parsed with errors!");
    GEOS_LOG_RANK_0("Error description: " << xmlResult.description());
    GEOS_LOG_RANK_0("Error offset: " << xmlResult.offset);
  }

  xmlWrapper::xmlNode xmlConstitutiveNode = xmlDocument.child("Constitutive");
  constitutiveManager.ProcessInputFileRecursive( xmlConstitutiveNode );
  constitutiveManager.PostProcessInputRecursive();

  LinearElasticAnisotropic * const model = constitutiveManager.GetConstitutiveRelation<LinearElasticAnisotropic>("granite");
  dataRepository::Group disc( "discretization", nullptr );
  disc.resize(1);
  model->AllocateConstitutiveData( &disc, 1 );

  arrayView1d<LinearElasticAnisotropic::StiffnessTensor const> const &
  stiffness = model->stiffness() ;

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
      ASSERT_DOUBLE_EQ( stiffness[0].m_data[i][j] , c[i][j] );
    }
  }

}
