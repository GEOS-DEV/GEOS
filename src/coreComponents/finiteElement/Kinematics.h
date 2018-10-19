/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

#ifndef _KINEMATICS_H_
#define _KINEMATICS_H_

#include "common/DataTypes.hpp"
#include "assert.h"

//*****************************************************************************
//***** DECLARATIONS **********************************************************
//*****************************************************************************

namespace geosx
{
void IncrementalKinematics( const R2TensorT<3>& A,
                            R2SymTensorT<3>& Dadt,
                            R2TensorT<3>& Rhat );

void IncrementalRotation( const R2TensorT<3>& A,
                          R2TensorT<3>& Rot );

inline void CalculateGradient( R2TensorT<3>& Gradient,
                               const int* bConnectivity,
                               const array1d<R1TensorT<3> >& disp,
                               const array1d<R1TensorT<3> >& dNdX )

{
  Gradient = 0.0;
  for( int a=0 ; a<8 ; ++a )
    Gradient.plus_dyadic_ab( disp(bConnectivity[a]), dNdX(a));

}

inline void CalculateGradient( R2Tensor& Gradient,
                               const array1d<R1Tensor >& disp,
                               const array1d<R1Tensor >& dNdX )

{
  //Gradient = 0.0;
  //for( int a=1 ; a<=8 ; ++a )
  //Gradient.plus_dyadic_ab( disp(bConnectivity[a-1]) , dNdX(a));

  assert( disp.size() == dNdX.size() );

  Gradient.dyadic_ab( disp(0), dNdX(0) );
  for( auto a=1 ; a<disp.size() ; ++a )
  {
    Gradient.plus_dyadic_ab( disp(a), dNdX(a) );
  }
}

inline void CalculateGradient( R2Tensor& Gradient,
                               const array1d<R1Tensor >& disp,
                               const R1Tensor* const dNdX )

{

  Gradient.dyadic_ab( disp(0), dNdX[0] );
  for( auto a=1 ; a<disp.size() ; ++a )
  {
    Gradient.plus_dyadic_ab( disp(a), dNdX[a] );
  }
}


inline void CalculateGradient(R2Tensor& Gradient, const R1Tensor * disp,
                              const arraySlice1d<R1Tensor> dNdX,
                              const localIndex numNodes)
{
  Gradient.dyadic_ab( disp[0], dNdX[0] );
  for( auto a=1 ; a<numNodes ; ++a )
  {
    Gradient.plus_dyadic_ab( disp[a], dNdX[a] );
  }
}

inline void CalculateGradient(R2Tensor& Gradient,
                              const R1Tensor * disp,
                              const R1TensorT<3> * dNdX,
                              const localIndex numNodes)
{
  Gradient.dyadic_ab( disp[0], dNdX[0] );
  for( auto a=1 ; a<numNodes ; ++a )
  {
    Gradient.plus_dyadic_ab( disp[a], dNdX[a] );
  }
}

template< int N >
inline void CalculateGradient(R2Tensor& Gradient,
                              const R1Tensor * disp,
                              const R1TensorT<3> * dNdX )
{
  Gradient.dyadic_ab( disp[0], dNdX[0] );
  for( auto a=1 ; a<N ; ++a )
  {
    Gradient.plus_dyadic_ab( disp[a], dNdX[a] );
  }
}

template< int N >
inline void CalculateGradients( R2Tensor& Gradient0,
                                R2Tensor& Gradient1,
                                R1Tensor const * restrict const var0,
                                R1Tensor const * restrict const var1,
                                R1TensorT<3> const * restrict const dNdX )
{
  Gradient0.dyadic_ab( var0[0], dNdX[0] );
  Gradient1.dyadic_ab( var1[0], dNdX[0] );
  for( localIndex a=1 ; a<N ; ++a )
  {
    Gradient0.plus_dyadic_ab( var0[a], dNdX[a] );
    Gradient1.plus_dyadic_ab( var1[a], dNdX[a] );
  }
}
void CalculatePhantomGradient( R2TensorT<3>& Gradient,
                               const int* bConnectivity,
                               const array1d<R1TensorT<3> >& disp,
                               const array2d<R1TensorT<3> >& dNdX );


//*****************************************************************************


}

#endif
