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


#include "finiteElement/elementFormulations/FiniteElementBase.hpp"
#include "finiteElement/elementFormulations/FiniteElementOperators.hpp"
#include "gtest/gtest.h"
#include "testFiniteElementHelpers.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include <chrono>


using namespace geos;
using namespace finiteElement;

//***** TEST VIEW SETTERS/GETTERS *****************************************************************

class TestFiniteElementBase final : public FiniteElementBase,
  public FiniteElementBase_impl< 8, 8, 8 >
{
  TestFiniteElementBase():
    FiniteElementBase( 8, 8, 8 ),
    FiniteElementBase_impl< 8, 8, 8 >()
  {}

  template< typename SUBREGION_TYPE >
  static void fillMeshData( NodeManager const & GEOS_UNUSED_PARAM( nodeManager ),
                            EdgeManager const & GEOS_UNUSED_PARAM( edgeManager ),
                            FaceManager const & GEOS_UNUSED_PARAM( faceManager ),
                            CellElementSubRegion const & GEOS_UNUSED_PARAM( cellSubRegion ),
                            MeshData< SUBREGION_TYPE > & GEOS_UNUSED_PARAM( meshData )
                            )
  {}

  template< typename SUBREGION_TYPE >
  // GEOS_HOST_DEVICE
  static void setupStack( localIndex const & GEOS_UNUSED_PARAM( cellIndex ),
                          MeshData< SUBREGION_TYPE > const & GEOS_UNUSED_PARAM( meshData ),
                          StackVariables & GEOS_UNUSED_PARAM( stack ) )
  {}
};



//***** TEST value() ******************************************************************************

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
    value[i] = 0;
  }
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    for( int i=0; i<NUM_COMPONENTS; ++i )
    {
      value[i] += N[a] * var[a][i];
    }
  }
}

TEST( FiniteElementBase, test_value )
{
  srand( 1234 );
  constexpr int NUM_SUPPORT_POINTS = 9;
  real64 N[NUM_SUPPORT_POINTS];
  real64 scalar[NUM_SUPPORT_POINTS];
  real64 vector[NUM_SUPPORT_POINTS][3];

  for( int q=0; q<20; ++q )
  {
    randomShape( N );
    randomSupportVar( scalar );
    randomSupportVar( vector );

    real64 referenceScalarValue = -1.0;
    real64 feBaseScalarValue = -1.0;
    value( N, scalar, referenceScalarValue );
    finiteElement::feOps::value( N, scalar, feBaseScalarValue );

    EXPECT_FLOAT_EQ( feBaseScalarValue, referenceScalarValue );


    real64 referenceVectorValue[3] = {-1, -1, -1};
    real64 feBaseVectorValue[3] = {-1, -1, -1};
    value( N, vector, referenceVectorValue );
    finiteElement::feOps::value( N, vector, feBaseVectorValue );

    EXPECT_FLOAT_EQ( feBaseVectorValue[0], referenceVectorValue[0] );
    EXPECT_FLOAT_EQ( feBaseVectorValue[1], referenceVectorValue[1] );
    EXPECT_FLOAT_EQ( feBaseVectorValue[2], referenceVectorValue[2] );
  }
}

//***** TEST symmetricGradient() ******************************************************************


template< int NUM_SUPPORT_POINTS >
static void symmetricGradient( real64 const (&gradN)[NUM_SUPPORT_POINTS][3],
                               real64 const (&var)[NUM_SUPPORT_POINTS][3],
                               real64 (& gradVar)[6] )
{
  for( int i=0; i<6; ++i )
  {
    gradVar[i] = 0.0;
  }
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    gradVar[0] += gradN[a][0] * var[ a ][0];
    gradVar[1] += gradN[a][1] * var[ a ][1];
    gradVar[2] += gradN[a][2] * var[ a ][2];
    gradVar[3] += gradN[a][2] * var[ a ][1] + gradN[a][1] * var[ a ][2];
    gradVar[4] += gradN[a][2] * var[ a ][0] + gradN[a][0] * var[ a ][2];
    gradVar[5] += gradN[a][1] * var[ a ][0] + gradN[a][0] * var[ a ][1];
  }
}


TEST( FiniteElementBase, test_symmetricGradient )
{
  srand( 1234 );
  constexpr int NUM_SUPPORT_POINTS = 9;
  real64 gradN[NUM_SUPPORT_POINTS][3];
  real64 vector[NUM_SUPPORT_POINTS][3];

  for( int q=0; q<20; ++q )
  {
    randomShapeGradient( gradN );
    randomSupportVar( vector );

    real64 referenceVectorGradient[6] = {-1, -1, -1, -1, -1, -1};
    real64 feBaseVectorGradient[6] = {-1, -1, -1, -1, -1, -1};
    symmetricGradient( gradN, vector, referenceVectorGradient );
    finiteElement::feOps::symmetricGradient( gradN, vector, feBaseVectorGradient );

    for( int i=0; i<6; ++i )
    {
      EXPECT_FLOAT_EQ( feBaseVectorGradient[i], referenceVectorGradient[i] );
    }
  }
}

//***** TEST gradient() ******************************************************************

template< int NUM_SUPPORT_POINTS >
static void gradient( real64 const (&gradN)[NUM_SUPPORT_POINTS][3],
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
      gradVar[i] += var[ a ] * gradN[a][i];
    }
  }
}

template< int NUM_SUPPORT_POINTS >
static void gradient( real64 const (&gradN)[NUM_SUPPORT_POINTS][3],
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
        gradVar[i][j] += var[ a ][i] * gradN[a][j];
      }
    }
  }
}

TEST( FiniteElementBase, test_gradient )
{
  srand( 1234 );
  constexpr int NUM_SUPPORT_POINTS = 9;
  real64 gradN[NUM_SUPPORT_POINTS][3];
  real64 scalar[NUM_SUPPORT_POINTS];
  real64 vector[NUM_SUPPORT_POINTS][3];

  for( int q=0; q<20; ++q )
  {
    randomShapeGradient( gradN );
    randomSupportVar( scalar );
    randomSupportVar( vector );

    real64 referenceScalarGradient[3] = {-1, -1, -1};
    real64 feBaseScalarGradient[3] = {-1, -1, -1};
    gradient( gradN, scalar, referenceScalarGradient );
    finiteElement::feOps::gradient( gradN, scalar, feBaseScalarGradient );


    real64 referenceVectorGradient[3][3] = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}};
    real64 feBaseVectorGradient[3][3] = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}};
    gradient( gradN, vector, referenceVectorGradient );
    finiteElement::feOps::gradient( gradN, vector, feBaseVectorGradient );

    for( int i=0; i<3; ++i )
    {
      for( int j=0; j<3; ++j )
      {
        EXPECT_FLOAT_EQ( feBaseVectorGradient[i][j], referenceVectorGradient[i][j] );
      }
    }
  }
}


//***** TEST valueAndGradient() *******************************************************************

template< int NUM_SUPPORT_POINTS >
static void valueAndGradient( real64 const (&N)[NUM_SUPPORT_POINTS],
                              real64 const (&gradN)[NUM_SUPPORT_POINTS][3],
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
      gradVar[i] +=  var[ a ] * gradN[a][i];
    }
  }
}

TEST( FiniteElementBase, test_valueAndGradient )
{
  srand( 1234 );
  constexpr int NUM_SUPPORT_POINTS = 9;
  real64 N[NUM_SUPPORT_POINTS];
  real64 gradN[NUM_SUPPORT_POINTS][3];
  real64 scalar[NUM_SUPPORT_POINTS];
  real64 vector[NUM_SUPPORT_POINTS][3];

  for( int q=0; q<20; ++q )
  {
    randomShape( N );
    randomShapeGradient( gradN );
    randomSupportVar( scalar );
    randomSupportVar( vector );



    real64 referenceScalarValue = -1.0;
    real64 feBaseScalarValue = -1.0;
    real64 referenceScalarGradient[3] = {-1, -1, -1};
    real64 feBaseScalarGradient[3] = {-1, -1, -1};
    valueAndGradient( N, gradN, scalar, referenceScalarValue, referenceScalarGradient );
    finiteElement::feOps::valueAndGradient( N, gradN, scalar, feBaseScalarValue, feBaseScalarGradient );

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
//    finiteElement::feOps::valueAndGradient( N, gradN, vector, feBaseVectorValue, feBaseVectorGradient );
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

//***** TEST plusGradNajAij() ********************************************************************

template< int NUM_SUPPORT_POINTS >
static void plusGradNajAij( real64 const (&gradN)[NUM_SUPPORT_POINTS][3],
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
// GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void plusGradNajAij( real64 const (&gradN)[NUM_SUPPORT_POINTS][3],
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



TEST( FiniteElementBase, test_plusGradNajAij )
{
  srand( 1234 );
  constexpr int NUM_SUPPORT_POINTS = 9;
  real64 gradN[NUM_SUPPORT_POINTS][3];

  real64 r2Tensor[3][3];
  real64 r2SymmTensor[6];

  real64 baselineResult[NUM_SUPPORT_POINTS][3] = {{0}};
  real64 feResult[NUM_SUPPORT_POINTS][3] = {{0}};
  real64 baselineResultSym[NUM_SUPPORT_POINTS][3] = {{0}};
  real64 feResultSym[NUM_SUPPORT_POINTS][3] = {{0}};

  for( int q=0; q<20; ++q )
  {
    randomShapeGradient( gradN );
    randomVar( r2Tensor );
    randomVar( r2SymmTensor );

    plusGradNajAij( gradN, r2Tensor, baselineResult );
    finiteElement::feOps::plusGradNajAij( gradN, r2Tensor, feResult );

    plusGradNajAij( gradN, r2SymmTensor, baselineResultSym );
    finiteElement::feOps::plusGradNajAij( gradN, r2SymmTensor, feResultSym );
  }

  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    for( int i=0; i<3; ++i )
    {
      EXPECT_FLOAT_EQ( feResult[a][i], baselineResult[a][i] );
      EXPECT_FLOAT_EQ( feResultSym[a][i], baselineResultSym[a][i] );
    }
  }
}

//***** TEST plusNaFi() ********************************************************************

template< int NUM_SUPPORT_POINTS >
static void plusNaFi( real64 const (&N)[NUM_SUPPORT_POINTS],
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

TEST( FiniteElementBase, test_plusNaFi )
{
  srand( 1234 );
  constexpr int NUM_SUPPORT_POINTS = 9;
  real64 N[NUM_SUPPORT_POINTS];
  real64 f[3];

  real64 baselineResult[NUM_SUPPORT_POINTS][3] = {{0}};
  real64 feResult[NUM_SUPPORT_POINTS][3] = {{0}};

  for( int q=0; q<20; ++q )
  {
    randomShape( N );
    randomVar( f );

    plusNaFi( N, f, baselineResult );
    finiteElement::feOps::plusNaFi( N, f, feResult );
  }

  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    for( int i=0; i<3; ++i )
    {
      EXPECT_FLOAT_EQ( feResult[a][i], baselineResult[a][i] );
    }
  }
}



//***** TEST plusGradNajAijPlusNaFi() **********************************************************


template< int NUM_SUPPORT_POINTS >
static void plusGradNajAijPlusNaFi( real64 const (&gradN)[NUM_SUPPORT_POINTS][3],
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
static void plusGradNajAijPlusNaFi( real64 const (&gradN)[NUM_SUPPORT_POINTS][3],
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

TEST( FiniteElementBase, test_plusGradNajAijPlusNaFi )
{
  srand( 1234 );
  constexpr int NUM_SUPPORT_POINTS = 9;
  real64 N[NUM_SUPPORT_POINTS];
  real64 gradN[NUM_SUPPORT_POINTS][3];

  real64 r2Tensor[3][3];
  real64 r2SymmTensor[6];
  real64 f[3];

  real64 baselineResult[NUM_SUPPORT_POINTS][3] = {{0}};
  real64 feResult[NUM_SUPPORT_POINTS][3] = {{0}};
  real64 baselineResultSym[NUM_SUPPORT_POINTS][3] = {{0}};
  real64 feResultSym[NUM_SUPPORT_POINTS][3] = {{0}};

  for( int q=0; q<20; ++q )
  {
    randomShape( N );
    randomShapeGradient( gradN );
    randomVar( r2Tensor );
    randomVar( r2SymmTensor );
    randomVar( f );

    plusGradNajAijPlusNaFi( gradN, r2Tensor, N, f, baselineResult );
    finiteElement::feOps::plusGradNajAijPlusNaFi( gradN, r2Tensor, N, f, feResult );

    plusGradNajAijPlusNaFi( gradN, r2SymmTensor, N, f, baselineResultSym );
    finiteElement::feOps::plusGradNajAijPlusNaFi( gradN, r2SymmTensor, N, f, feResultSym );
  }

  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    for( int i=0; i<3; ++i )
    {
      EXPECT_FLOAT_EQ( feResult[a][i], baselineResult[a][i] );
      EXPECT_FLOAT_EQ( feResultSym[a][i], baselineResultSym[a][i] );
    }
  }
}
