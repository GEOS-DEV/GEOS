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

#ifndef _KINEMATICS_H_
#define _KINEMATICS_H_

#include "common/DataTypes.hpp"
#include "assert.h"

//*****************************************************************************
//***** DECLARATIONS **********************************************************
//*****************************************************************************

namespace geosx
{
void IncrementalKinematics( const R2Tensor& A,
                            R2SymTensor& Dadt,
                            R2Tensor& Rhat );

void IncrementalRotation( const R2Tensor& A,
                          R2TensorT<3>& Rot );

inline void CalculateGradient( R2Tensor& Gradient,
                               const int* bConnectivity,
                               arraySlice1d<R1Tensor> const & disp,
                               arraySlice1d<R1Tensor> const & dNdX )
{
  Gradient = 0.0;
  for( localIndex a=0 ; a<8 ; ++a )
    Gradient.plus_dyadic_ab( disp[bConnectivity[a]], dNdX[a]);
}

inline void CalculateGradient(R2Tensor& Gradient,
                              arraySlice1d<R1Tensor const> const & disp,
                              arraySlice1d<R1Tensor const> const & dNdX,
                              localIndex numNodes)
{
  Gradient.dyadic_ab( disp[0], dNdX[0] );
  for( localIndex a=1 ; a<numNodes ; ++a )
  {
    Gradient.plus_dyadic_ab( disp[a], dNdX[a] );
  }
}

template< int N >
inline void CalculateGradient(R2Tensor& Gradient,
                              arraySlice1d<R1Tensor const> const & disp,
                              arraySlice1d<R1Tensor const> const & dNdX )
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
                                arraySlice1d<R1Tensor const> const & var0,
                                arraySlice1d<R1Tensor const> const & var1,
                                arraySlice1d<R1Tensor const> const & dNdX )
{
  Gradient0.dyadic_ab( var0[0], dNdX[0] );
  Gradient1.dyadic_ab( var1[0], dNdX[0] );
  for( localIndex a=1 ; a<N ; ++a )
  {
    Gradient0.plus_dyadic_ab( var0[a], dNdX[a] );
    Gradient1.plus_dyadic_ab( var1[a], dNdX[a] );
  }
}

}

#endif
