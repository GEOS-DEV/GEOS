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

inline void HughesWinget( R2Tensor &Rot, R2SymTensor & Dadt, R2Tensor const & G)
{

  real64 * restrict const Dadt_data = Dadt.Data();
  real64 * restrict const Rot_data = Rot.Data();
  real64 const * restrict const G_data = G.Data();


  //Dadt = 0.5*(G + GT);
  Dadt_data[0] = G_data[0];

  Dadt_data[1] = 0.5*(G_data[1] + G_data[3]);
  Dadt_data[2] = G_data[4];

  Dadt_data[3] = 0.5*(G_data[6] + G_data[2]);
  Dadt_data[4] = 0.5*(G_data[7] + G_data[5]);
  Dadt_data[5] = G_data[8];


  //Omega = 0.5*(G - GT);
  real64 const w12 = 0.5*(G_data[1] - G_data[3]);
  real64 const w13 = 0.5*(G_data[2] - G_data[6]);
  real64 const w23 = 0.5*(G_data[5] - G_data[7]);

  real64 const w12w12div4 = 0.25*w12*w12;
  real64 const w13w13div4 = 0.25*w13*w13;
  real64 const w23w23div4 = 0.25*w23*w23;
  real64 const w12w13div2 = 0.5*(w12*w13);
  real64 const w12w23div2 = 0.5*(w12*w23);
  real64 const w13w23div2 = 0.5*(w13*w23);
  real64 const invDetIplusOmega = 1.0 / ( 1 + ( w12w12div4 + w13w13div4 + w23w23div4 ) );

  Rot_data[0] = ( 1.0 + (-w12w12div4 - w13w13div4 + w23w23div4) ) * invDetIplusOmega;
  Rot_data[1] = ( w12 - w13w23div2 ) * invDetIplusOmega;
  Rot_data[2] = ( w13 + w12w23div2 ) * invDetIplusOmega;

  Rot_data[3] = (-w12 - w13w23div2 ) * invDetIplusOmega;
  Rot_data[4] = ( 1.0 + (-w12w12div4 + w13w13div4 - w23w23div4) ) * invDetIplusOmega;
  Rot_data[5] = ( w23 - w12w13div2 ) * invDetIplusOmega;

  Rot_data[6] = (-w13 + w12w23div2 ) * invDetIplusOmega;
  Rot_data[7] = (-w23 - w12w13div2 ) * invDetIplusOmega;
  Rot_data[8] = ( 1.0 + ( w12w12div4 - w13w13div4 - w23w23div4) ) * invDetIplusOmega;


}

}

#endif
