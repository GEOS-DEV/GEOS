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

#ifndef _KINEMATICS_H_
#define _KINEMATICS_H_

#include "common/DataTypes.hpp"
#include "common/GeosxMacros.hpp"
#include "assert.h"

//*****************************************************************************
//***** DECLARATIONS **********************************************************
//*****************************************************************************

namespace geosx
{
void IncrementalKinematics( const R2Tensor & A,
                            R2SymTensor & Dadt,
                            R2Tensor & Rhat );

void IncrementalRotation( const R2Tensor & A,
                          R2TensorT< 3 > & Rot );

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void CalculateGradient( R2Tensor & Gradient,
                        const int * bConnectivity,
                        arraySlice1d< R1Tensor > const & disp,
                        arraySlice1d< R1Tensor > const & dNdX )
{
  Gradient = 0.0;
  for( localIndex a=0; a<8; ++a )
    Gradient.plus_dyadic_ab( disp[bConnectivity[a]], dNdX[a] );
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void CalculateGradient( R2Tensor & Gradient,
                        arraySlice1d< R1Tensor const > const & disp,
                        arraySlice1d< R1Tensor const > const & dNdX,
                        localIndex numNodes )
{
  Gradient.dyadic_ab( disp[0], dNdX[0] );
  for( localIndex a=1; a<numNodes; ++a )
  {
    Gradient.plus_dyadic_ab( disp[a], dNdX[a] );
  }
}

template< int N >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void CalculateGradient( R2Tensor & Gradient,
                        arraySlice1d< R1Tensor const > const & disp,
                        arraySlice1d< R1Tensor const > const & dNdX )
{
  Gradient.dyadic_ab( disp[0], dNdX[0] );
  for( auto a=1; a<N; ++a )
  {
    Gradient.plus_dyadic_ab( disp[a], dNdX[a] );
  }
}

template< int N >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void CalculateGradients( R2Tensor & Gradient0,
                         R2Tensor & Gradient1,
                         R1Tensor const * GEOSX_RESTRICT const var0,
                         R1Tensor const * GEOSX_RESTRICT const var1,
                         arraySlice1d< R1Tensor const > const & dNdX )
{
  Gradient0.dyadic_ab( var0[0], dNdX[0] );
  Gradient1.dyadic_ab( var1[0], dNdX[0] );
  for( localIndex a=1; a<N; ++a )
  {
    Gradient0.plus_dyadic_ab( var0[a], dNdX[a] );
    Gradient1.plus_dyadic_ab( var1[a], dNdX[a] );
  }
}

template< int N >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void CalculateGradients( R2Tensor & Gradient0,
                         R2Tensor & Gradient1,
                         R1Tensor const * GEOSX_RESTRICT const var0,
                         R1Tensor const * GEOSX_RESTRICT const var1,
                         real64 const (&dNdX)[8][3] )
{
  Gradient0 = 0;
  Gradient1 = 0;
  real64 * const GEOSX_RESTRICT g0 = Gradient0.Data();
  real64 * const GEOSX_RESTRICT g1 = Gradient1.Data();

  for( int a=0; a<N; ++a )
  {
    real64 const * const GEOSX_RESTRICT v0 = var0[a].Data();
    real64 const * const GEOSX_RESTRICT v1 = var1[a].Data();
    g0[0] += v0[0]*dNdX[a][0];
    g0[1] += v0[0]*dNdX[a][1];
    g0[2] += v0[0]*dNdX[a][2];

    g0[3] += v0[1]*dNdX[a][0];
    g0[4] += v0[1]*dNdX[a][1];
    g0[5] += v0[1]*dNdX[a][2];

    g0[6] += v0[2]*dNdX[a][0];
    g0[7] += v0[2]*dNdX[a][1];
    g0[8] += v0[2]*dNdX[a][2];


    g1[0] += v1[0]*dNdX[a][0];
    g1[1] += v1[0]*dNdX[a][1];
    g1[2] += v1[0]*dNdX[a][2];

    g1[3] += v1[1]*dNdX[a][0];
    g1[4] += v1[1]*dNdX[a][1];
    g1[5] += v1[1]*dNdX[a][2];

    g1[6] += v1[2]*dNdX[a][0];
    g1[7] += v1[2]*dNdX[a][1];
    g1[8] += v1[2]*dNdX[a][2];
  }
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void HughesWinget( R2Tensor & Rot, R2SymTensor & Dadt, R2Tensor const & G )
{

  real64 * GEOSX_RESTRICT const Dadt_data = Dadt.Data();
  real64 * GEOSX_RESTRICT const Rot_data = Rot.Data();
  real64 const * GEOSX_RESTRICT const G_data = G.Data();


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
