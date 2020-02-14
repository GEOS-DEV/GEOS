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

GEOSX_HOST_DEVICE
inline void CalculateGradient(R2Tensor& Gradient,
                              R1Tensor const * const disp,
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
inline void CalculateGradients( R2Tensor & Gradient0,
                                R2Tensor & Gradient1,
                                R1Tensor const * GEOSX_RESTRICT const var0,
                                R1Tensor const * GEOSX_RESTRICT const var1,
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
