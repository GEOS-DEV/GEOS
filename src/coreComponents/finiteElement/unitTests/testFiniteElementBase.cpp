/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file testFiniteElementBase.cpp
 */

#include "finiteElement/elementFormulations/FiniteElementBase.hpp"
#include "managers/initialization.hpp"
#include "gtest/gtest.h"
#include "testFiniteElementHelpers.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include <chrono>


using namespace geosx;
using namespace finiteElement;

static real64 inverse( real64 (& J)[3][3] )
{
  real64 scratch[3][3];
  scratch[0][0] = J[1][1]*J[2][2] - J[1][2]*J[2][1];
  scratch[1][0] = J[0][2]*J[2][1] - J[0][1]*J[2][2];
  scratch[2][0] = J[0][1]*J[1][2] - J[0][2]*J[1][1];
  scratch[0][1] = J[1][2]*J[2][0] - J[1][0]*J[2][2];
  scratch[1][1] = J[0][0]*J[2][2] - J[0][2]*J[2][0];
  scratch[2][1] = J[0][2]*J[1][0] - J[0][0]*J[1][2];
  scratch[0][2] = J[1][0]*J[2][1] - J[1][1]*J[2][0];
  scratch[1][2] = J[0][1]*J[2][0] - J[0][0]*J[2][1];
  scratch[2][2] = J[0][0]*J[1][1] - J[0][1]*J[1][0];

  real64 const detJ = J[0][0] * scratch[0][0] + J[1][0] * scratch[1][0] + J[2][0] * scratch[2][0];
  real64 const invDet = 1 / detJ;

  for( int i=0; i<3; ++i )
  {
    for( int j=0; j<3; ++j )
    {
      J[i][j] = scratch[j][i] * invDet;
    }
  }

  return detJ;
}



template< int NUM_SUPPORT_POINTS >
static void value( real64 const (&N)[NUM_SUPPORT_POINTS],
                               real64 const (&var)[NUM_SUPPORT_POINTS],
                               real64 & value )
{
  value = 0;
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    value += N[a] * var[a];
  }
}

template< int NUM_SUPPORT_POINTS,
          int NUM_COMPONENTS >
static void value( real64 const (&N)[NUM_SUPPORT_POINTS],
                               real64 const (&var)[NUM_SUPPORT_POINTS][NUM_COMPONENTS],
                               real64 (& value)[NUM_COMPONENTS] )
{
  for( int i=0; i<NUM_COMPONENTS; ++i )
  {
    value[i] = 0 ;
  }
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    for( int i=0; i<NUM_COMPONENTS; ++i )
    {
      value[i] += N[a] * var[a][i];
    }
  }
}


template< int NUM_SUPPORT_POINTS >
static void symmetricGradient( real64 const (& gradN)[NUM_SUPPORT_POINTS][3],
                                           real64 const (&var)[NUM_SUPPORT_POINTS][3],
                                           real64 (& gradVar)[6] )
{
  for( int i=0; i<6; ++i )
  {
    gradVar[i] = 0.0;
  }
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    gradVar[0] = gradVar[0] + gradN[a][0] * var[ a ][0];
    gradVar[1] = gradVar[1] + gradN[a][1] * var[ a ][1];
    gradVar[2] = gradVar[2] + gradN[a][2] * var[ a ][2];
    gradVar[3] = gradVar[3] + gradN[a][2] * var[ a ][1] + gradN[a][1] * var[ a ][2];
    gradVar[4] = gradVar[4] + gradN[a][2] * var[ a ][0] + gradN[a][0] * var[ a ][2];
    gradVar[5] = gradVar[5] + gradN[a][1] * var[ a ][0] + gradN[a][0] * var[ a ][1];
  }
}

template< int NUM_SUPPORT_POINTS >
static void gradient( real64 const (& gradN)[NUM_SUPPORT_POINTS][3],
                                  real64 const (&var)[NUM_SUPPORT_POINTS],
                                  real64 (& gradVar)[3] )
{
  for( int i=0; i<3; ++i )
  {
    gradVar[i] = 0.0;
  }
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    for( int i = 0; i < 3; ++i )
    {
      gradVar[i] = gradVar[i] + var[ a ] * gradN[a][i];
    }
  }
}

template< int NUM_SUPPORT_POINTS >
static void gradient( real64 const (& gradN)[NUM_SUPPORT_POINTS][3],
                                  real64 const (&var)[NUM_SUPPORT_POINTS][3],
                                  real64 (& gradVar)[3][3] )
{
  for( int i = 0; i < 3; ++i )
  {
    for( int j = 0; j < 3; ++j )
    {
      gradVar[i][j] = 0.0;
    }
  }
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        gradVar[i][j] = gradVar[i][j] + var[ a ][i] * gradN[a][j];
      }
    }
  }
}



template< int NUM_SUPPORT_POINTS >
static void valueAndGradient( real64 const (&N)[NUM_SUPPORT_POINTS],
                              real64 const (& gradN)[NUM_SUPPORT_POINTS][3],
                                          real64 const (&var)[NUM_SUPPORT_POINTS],
                                          real64 & value,
                                          real64 (& gradVar)[3] )
{
  value = 0;
  for( int i = 0; i < 3; ++i )
  {
    gradVar[i] = 0;
  }

  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    value += N[a] * var[a];
    for( int i = 0; i < 3; ++i )
    {
      gradVar[i] = gradVar[i] + var[ a ] * gradN[a][i];
    }
  }
}



template< int NUM_SUPPORT_POINTS >
static void gradNajAij( real64 const (& gradN)[NUM_SUPPORT_POINTS][3],
                                    real64 const (&var_detJxW)[6],
                                    real64 (& R)[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    R[a][0] += var_detJxW[0] * gradN[a][0] + var_detJxW[5] * gradN[a][1] + var_detJxW[4] * gradN[a][2];
    R[a][1] += var_detJxW[5] * gradN[a][0] + var_detJxW[1] * gradN[a][1] + var_detJxW[3] * gradN[a][2];
    R[a][2] += var_detJxW[4] * gradN[a][0] + var_detJxW[3] * gradN[a][1] + var_detJxW[2] * gradN[a][2];
  }
}


template< int NUM_SUPPORT_POINTS >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void gradNajAij( real64 const (& gradN)[NUM_SUPPORT_POINTS][3],
                                    real64 const (&var_detJxW)[3][3],
                                    real64 (& R)[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    R[a][0] += var_detJxW[0][0] * gradN[a][0] + var_detJxW[0][1] * gradN[a][1] + var_detJxW[0][2] * gradN[a][2];
    R[a][1] += var_detJxW[1][0] * gradN[a][0] + var_detJxW[1][1] * gradN[a][1] + var_detJxW[1][2] * gradN[a][2];
    R[a][2] += var_detJxW[2][0] * gradN[a][0] + var_detJxW[2][1] * gradN[a][1] + var_detJxW[2][2] * gradN[a][2];
  }
}

template< int NUM_SUPPORT_POINTS >
static void NaFi( real64 const (&N)[NUM_SUPPORT_POINTS],
                              real64 const (&var_detJxW)[3],
                              real64 ( & R )[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    R[a][0] += var_detJxW[0] * N[a];
    R[a][1] += var_detJxW[1] * N[a];
    R[a][2] += var_detJxW[2] * N[a];
  }
}


template< int NUM_SUPPORT_POINTS >
static void gradNajAij_plus_NaFi( real64 const (& gradN)[NUM_SUPPORT_POINTS][3],
                                              real64 const (&var_detJxW)[6],
                                              real64 const (&N)[NUM_SUPPORT_POINTS],
                                              real64 const (&forcingTerm_detJxW)[3],
                                              real64 (& R)[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    R[a][0] += var_detJxW[0] * gradN[a][0] + var_detJxW[5] * gradN[a][1] + var_detJxW[4] * gradN[a][2] + forcingTerm_detJxW[0] * N[a];
    R[a][1] += var_detJxW[5] * gradN[a][0] + var_detJxW[1] * gradN[a][1] + var_detJxW[3] * gradN[a][2] + forcingTerm_detJxW[1] * N[a];
    R[a][2] += var_detJxW[4] * gradN[a][0] + var_detJxW[3] * gradN[a][1] + var_detJxW[2] * gradN[a][2] + forcingTerm_detJxW[2] * N[a];
  }
}

template< int NUM_SUPPORT_POINTS >
static void gradNajAij_plus_NaFi( real64 const (& gradN)[NUM_SUPPORT_POINTS][3],
                                              real64 const (&var_detJxW)[3][3],
                                              real64 const (&N)[NUM_SUPPORT_POINTS],
                                              real64 const (&forcingTerm_detJxW)[3],
                                              real64 (& R)[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    R[a][0] += var_detJxW[0][0] * gradN[a][0] + var_detJxW[0][1] * gradN[a][1] + var_detJxW[0][2] * gradN[a][2] + forcingTerm_detJxW[0] * N[a];
    R[a][1] += var_detJxW[1][0] * gradN[a][0] + var_detJxW[1][1] * gradN[a][1] + var_detJxW[1][2] * gradN[a][2] + forcingTerm_detJxW[1] * N[a];
    R[a][2] += var_detJxW[2][0] * gradN[a][0] + var_detJxW[2][1] * gradN[a][1] + var_detJxW[2][2] * gradN[a][2] + forcingTerm_detJxW[2] * N[a];
  }
}

struct Jacobian
{
  real64 data[3][3] = { {  1.19167, -0.0372008, -0.0766346, },
    {  0.0599679, 1.19167, 0.0205342, },
    { -0.0438996, 0.00610042, 1.1378, } };
};

template< typename POLICY >
void testInverseDriver()
{
  Jacobian J;
  Jacobian invJ;
  real64 detJ;

  forAll< POLICY >( 1,
                    [&] ( localIndex const )
  {
    detJ = FiniteElementBase::inverse( invJ.data );
  } );


  real64 const detJ_Ref = inverse( J.data );

  EXPECT_FLOAT_EQ( detJ, detJ_Ref );

  forAll< serialPolicy >( 1,
                          [=] ( localIndex const )
  {
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        EXPECT_FLOAT_EQ( J.data[i][j], invJ.data[i][j] );
      }
    }
  } );
}


template< typename POLICY >
void testDetJDriver()
{
  Jacobian J;
  real64 detJ;

  forAll< POLICY >( 1,
                    [&] ( localIndex const )
  {
    detJ = FiniteElementBase::detJ( J.data );
  } );

  real64 const detJ_Ref = inverse( J.data );

  EXPECT_FLOAT_EQ( detJ, detJ_Ref );

}


TEST( FiniteElementBase, testJacobianInverseHost )
{
  testInverseDriver< serialPolicy >();
}

TEST( FiniteElementBase, testDetJHost )
{
  testDetJDriver< serialPolicy >();
}


TEST( FiniteElementBase, setGradNView )
{
}

TEST( FiniteElementBase, setDetJView )
{

}

TEST( FiniteElementBase, testGetGradN )
{

}


TEST( FiniteElementBase, testValue )
{
  srand( 1234 );
  constexpr int NUM_SUPPORT_POINTS = 9;
  real64 N[NUM_SUPPORT_POINTS];
  real64 scalar[NUM_SUPPORT_POINTS];
  real64 vector[NUM_SUPPORT_POINTS][3];

  for( int q=0; q<20; ++q )
  {
    randomShape(N);
    randomSupportVar(scalar);
    randomSupportVar(vector);

    real64 referenceScalarValue = -1.0;
    real64 feBaseScalarValue = -1.0;
    value( N, scalar, referenceScalarValue );
    FiniteElementBase::value( N, scalar, feBaseScalarValue );

    EXPECT_FLOAT_EQ( feBaseScalarValue, referenceScalarValue );


    real64 referenceVectorValue[3] = {-1,-1,-1};
    real64 feBaseVectorValue[3] = {-1,-1,-1};
    value( N, vector, referenceVectorValue );
    FiniteElementBase::value( N, vector, feBaseVectorValue );

    EXPECT_FLOAT_EQ( feBaseVectorValue[0], referenceVectorValue[0] );
    EXPECT_FLOAT_EQ( feBaseVectorValue[1], referenceVectorValue[1] );
    EXPECT_FLOAT_EQ( feBaseVectorValue[2], referenceVectorValue[2] );
  }
}


TEST( FiniteElementBase, testSymmetricGradient )
{
  srand( 1234 );
  constexpr int NUM_SUPPORT_POINTS = 9;
  real64 gradN[NUM_SUPPORT_POINTS][3];
  real64 vector[NUM_SUPPORT_POINTS][3];

  for( int q=0; q<20; ++q )
  {
    randomShapeGradient(gradN);
    randomSupportVar(vector);

    real64 referenceVectorGradient[6] = {-1,-1,-1,-1,-1,-1};
    real64 feBaseVectorGradient[6] = {-1,-1,-1,-1,-1,-1};
    symmetricGradient( gradN, vector, referenceVectorGradient );
    FiniteElementBase::symmetricGradient( gradN, vector, feBaseVectorGradient );

    for( int i=0; i<6; ++i )
    {
      EXPECT_FLOAT_EQ( feBaseVectorGradient[i], referenceVectorGradient[i] );
    }
  }
}

TEST( FiniteElementBase, testGradient )
{
  srand( 1234 );
  constexpr int NUM_SUPPORT_POINTS = 9;
  real64 gradN[NUM_SUPPORT_POINTS][3];
  real64 scalar[NUM_SUPPORT_POINTS];
  real64 vector[NUM_SUPPORT_POINTS][3];

  for( int q=0; q<20; ++q )
  {
    randomShapeGradient(gradN);
    randomSupportVar(scalar);
    randomSupportVar(vector);

    real64 referenceScalarGradient[3] = {-1,-1,-1};
    real64 feBaseScalarGradient[3] = {-1,-1,-1};
    gradient( gradN, scalar, referenceScalarGradient );
    FiniteElementBase::gradient( gradN, scalar, feBaseScalarGradient );


    real64 referenceVectorGradient[3][3] = {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};
    real64 feBaseVectorGradient[3][3] = {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};
    gradient( gradN, vector, referenceVectorGradient );
    FiniteElementBase::gradient( gradN, vector, feBaseVectorGradient );

    for( int i=0; i<3; ++i )
    {
      for( int j=0; j<3; ++j )
      {
        EXPECT_FLOAT_EQ( feBaseVectorGradient[i][j], referenceVectorGradient[i][j] );
      }
    }
  }
}


TEST( FiniteElementBase, testValueAndGradient )
{
  srand( 1234 );
  constexpr int NUM_SUPPORT_POINTS = 9;
  real64 N[NUM_SUPPORT_POINTS];
  real64 gradN[NUM_SUPPORT_POINTS][3];
  real64 scalar[NUM_SUPPORT_POINTS];
  real64 vector[NUM_SUPPORT_POINTS][3];

  for( int q=0; q<20; ++q )
  {
    randomShape(N);
    randomShapeGradient(gradN);
    randomSupportVar(scalar);
    randomSupportVar(vector);



    real64 referenceScalarValue = -1.0;
    real64 feBaseScalarValue = -1.0;
    real64 referenceScalarGradient[3] = {-1,-1,-1};
    real64 feBaseScalarGradient[3] = {-1,-1,-1};
    valueAndGradient( N, gradN, scalar, referenceScalarValue, referenceScalarGradient );
    FiniteElementBase::valueAndGradient( N, gradN, scalar, feBaseScalarValue, feBaseScalarGradient );

    EXPECT_FLOAT_EQ( feBaseScalarValue, referenceScalarValue );
    for( int i=0; i<3; ++i )
    {
      EXPECT_FLOAT_EQ( feBaseScalarGradient[i], referenceScalarGradient[i] );
    }

    // NOT IMPLEMENTED YET
//    real64 referenceVectorValue[3] = {-1,-1,-1};
//    real64 feBaseVectorValue[3] = {-1,-1,-1};
//    real64 referenceVectorGradient[3][3] = {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};
//    real64 feBaseVectorGradient[3][3] = {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};
//    valueAndGradient( N, gradN, vector, referenceVectorValue, referenceVectorGradient );
//    FiniteElementBase::valueAndGradient( N, gradN, vector, feBaseVectorValue, feBaseVectorGradient );
//
//
//    for( int i=0; i<3; ++i )
//    {
//      EXPECT_FLOAT_EQ( feBaseVectorValue[i], referenceVectorValue[i] );
//      for( int j=0; j<3; ++j )
//      {
//        EXPECT_FLOAT_EQ( feBaseVectorGradient[i][j], referenceVectorGradient[i][j] );
//      }
//    }


  }
}

TEST( FiniteElementBase, testGradNajAij )
{
  srand( 1234 );
  constexpr int NUM_SUPPORT_POINTS = 9;
  real64 N[NUM_SUPPORT_POINTS];
  real64 gradN[NUM_SUPPORT_POINTS][3];

  real64 r2Tensor[3][3];
  real64 r2SymmTensor[6];
  real64 f[3];

  real64 baselineResult[NUM_SUPPORT_POINTS][3] = {0};
  real64 feResult[NUM_SUPPORT_POINTS][3] = {0};
  real64 baselineResultSym[NUM_SUPPORT_POINTS][3] = {0};
  real64 feResultSym[NUM_SUPPORT_POINTS][3] = {0};


  real64 baselineResult2[NUM_SUPPORT_POINTS][3] = {0};
  real64 feResult2[NUM_SUPPORT_POINTS][3] = {0};
  real64 baselineResultSym2[NUM_SUPPORT_POINTS][3] = {0};
  real64 feResultSym2[NUM_SUPPORT_POINTS][3] = {0};

  for( int q=0; q<20; ++q )
  {
    randomShape(N);
    randomShapeGradient(gradN);
    randomVar(r2Tensor);
    randomVar(r2SymmTensor);
    randomVar(f);

    gradNajAij( gradN, r2Tensor, baselineResult );
    FiniteElementBase::gradNajAij( gradN, r2Tensor, feResult );

    gradNajAij( gradN, r2SymmTensor, baselineResultSym );
    FiniteElementBase::gradNajAij( gradN, r2SymmTensor, feResultSym );

    gradNajAij_plus_NaFi( gradN, r2Tensor, N, f, baselineResult2 );
    FiniteElementBase::gradNajAij_plus_NaFi( gradN, r2Tensor, N, f, feResult2 );

    gradNajAij_plus_NaFi( gradN, r2SymmTensor, N, f, baselineResultSym2 );
    FiniteElementBase::gradNajAij_plus_NaFi( gradN, r2SymmTensor, N, f, feResultSym2 );
  }

  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    for( int i=0; i<3; ++i )
    {
      EXPECT_FLOAT_EQ( feResult[a][i], baselineResult[a][i] );
      EXPECT_FLOAT_EQ( feResultSym[a][i], baselineResultSym[a][i] );
      EXPECT_FLOAT_EQ( feResult2[a][i], baselineResult2[a][i] );
      EXPECT_FLOAT_EQ( feResultSym2[a][i], baselineResultSym2[a][i] );
    }
  }

}


using namespace geosx;
int main( int argc, char * argv[] )
{
  testing::InitGoogleTest();

  basicSetup( argc, argv, false );

  int const result = RUN_ALL_TESTS();

  basicCleanup();

  return result;
}
