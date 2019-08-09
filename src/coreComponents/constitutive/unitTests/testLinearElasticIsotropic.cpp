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
#include "constitutive/Solid/LinearElasticIsotropic.hpp"

#include "dataRepository/xmlWrapper.hpp"
using namespace geosx;
using namespace ::geosx::constitutive;

TEST( LinearElasticIsotropicTests, testAllocation )
{
  LinearElasticIsotropic cm( "model", nullptr );

  localIndex constexpr numElems = 2;
  localIndex constexpr numQuadraturePoints = 3;

  dataRepository::Group disc( "discretization", nullptr );
  disc.resize(numElems);
  cm.AllocateConstitutiveData( &disc, numQuadraturePoints );

  EXPECT_EQ( cm.size(), numElems );
  EXPECT_EQ( cm.numQuadraturePoints(), numQuadraturePoints );

  arrayView1d<real64 const> const & bulkModulus = cm.bulkModulus() ;
  arrayView1d<real64 const> const & shearModulus = cm.shearModulus() ;
  arrayView2d<real64 const>      const & meanStress = cm.meanStress();
  arrayView2d<R2SymTensor const> const & deviatorStress = cm.deviatorStress();

  EXPECT_EQ( bulkModulus.size(), numElems );
  EXPECT_EQ( shearModulus.size(), numElems );
  EXPECT_EQ( meanStress.size(0), numElems );
  EXPECT_EQ( meanStress.size(1), numQuadraturePoints );
  EXPECT_EQ( deviatorStress.size(0), numElems );
  EXPECT_EQ( deviatorStress.size(1), numQuadraturePoints );

}

TEST( LinearElasticIsotropicTests, testStateUpdatePoint )
{
  LinearElasticIsotropic cm( "model", nullptr );

  real64 constexpr K = 2e10;
  real64 constexpr G = 1e10;
  cm.setDefaultBulkModulus(K);
  cm.setDefaultShearModulus(G);

  dataRepository::Group disc( "discretization", nullptr );
  disc.resize(2);
  cm.AllocateConstitutiveData( &disc, 2 );

//  cm.bulkModulus() = cm.setDefaultBulkModulus();
//  cm.shearModulus() = cm.setDefaultShearModulus();

  arrayView2d<real64>      const & meanStress = cm.meanStress();
  arrayView2d<R2SymTensor> const & deviatorStress = cm.deviatorStress();

  real64 const strain = 0.1;
  R2SymTensor Ddt;
  R2Tensor Rot;
  R2SymTensor zero;

  {
    Ddt(0,0) = strain;
    Rot(0,0) = 1;
    Rot(1,1) = 1;
    Rot(2,2) = 1;

    cm.StateUpdatePoint( 0, 0, Ddt, Rot, 0 );

    ASSERT_DOUBLE_EQ( meanStress[0][0] , strain*K );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](0,0) , (2.0/3.0*strain)*2*G );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](1,1) , (-1.0/3.0*strain)*2*G );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](2,2) , (-1.0/3.0*strain)*2*G );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](1,2) , 0.0 );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](0,2) , 0.0 );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](0,1) , 0.0 );
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

    ASSERT_DOUBLE_EQ( meanStress[0][0] , strain*K );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](0,0) , (-1.0/3.0*strain)*2*G );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](1,1) ,  (2.0/3.0*strain)*2*G );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](2,2) , (-1.0/3.0*strain)*2*G );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](1,2) , 0.0 );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](0,2) , 0.0 );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](0,1) , 0.0 );
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

    ASSERT_DOUBLE_EQ( meanStress[0][0] , strain*K );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](0,0) , (-1.0/3.0*strain)*2*G );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](1,1) , (-1.0/3.0*strain)*2*G );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](2,2) ,  (2.0/3.0*strain)*2*G );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](1,2) , 0.0 );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](0,2) , 0.0 );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](0,1) , 0.0 );
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

    ASSERT_DOUBLE_EQ( meanStress[0][0] , 0 );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](0,0) , 0.0 );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](1,1) , 0.0 );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](2,2) , 0.0 );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](1,2) , 0.0 );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](0,2) , 0.0 );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](0,1) , strain*2*G );
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

    ASSERT_DOUBLE_EQ( meanStress[0][0] , 0 );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](0,0) , 0.0 );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](1,1) , 0.0 );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](2,2) , 0.0 );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](1,2) , 0.0 );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](0,2) , strain*2*G );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](0,1) , 0.0 );
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

    ASSERT_DOUBLE_EQ( meanStress[0][0] , 0 );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](0,0) , 0.0 );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](1,1) , 0.0 );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](2,2) , 0.0 );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](1,2) , strain*2*G );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](0,2) , 0.0 );
    ASSERT_DOUBLE_EQ( deviatorStress[0][0](0,1) , 0.0 );
  }
}



TEST( LinearElasticIsotropicTests, testXML )
{
  ConstitutiveManager constitutiveManager("constitutive",nullptr);
  LinearElasticIsotropic cm( "model", &constitutiveManager );

  string const inputStream =
  "<?xml version=\"1.0\" ?>"
  "  <Constitutive xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"geos_v0.0.xsd\">"
  "  <LinearElasticIsotropic name=\"granite\" "
  "  defaultDensity=\"2700\" "
  "  defaultBulkModulus=\"5.5556e9\" "
  "  defaultShearModulus=\"4.16667e9\"/>"
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

}
