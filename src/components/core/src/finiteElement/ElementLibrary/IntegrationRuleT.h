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

/**
 * @file IntegrationRuleT.h
 * @author settgast1
 * @date Dec 6, 2010
 */

#ifndef INTEGRATIONRULET_H_
#define INTEGRATIONRULET_H_

#include "Common/Common.h"

class IntegrationRuleT
{
public:
  IntegrationRuleT();
  virtual ~IntegrationRuleT();


  void CalculateShapeFunctionDerivatives( const array<R1Tensor>& X_local,
                                          array<R1Tensor>& dNdX,
                                          realT& detJ );


};


void CalculateShapeFunctionDerivative( const realT y[8],
                                       const realT z[8],
                                       realT b[8] );



inline void CalculateShapeFunctionDerivative( const realT y[8],
                                              const realT z[8],
                                              realT b[8] )
{
  const realT y0 = y[0];
  const realT y1 = y[1];
  const realT y2 = y[2];
  const realT y3 = y[3];
  const realT y4 = y[4];
  const realT y5 = y[5];
  const realT y6 = y[6];
  const realT y7 = y[7];

  const realT z0 = z[0];
  const realT z1 = z[1];
  const realT z2 = z[2];
  const realT z3 = z[3];
  const realT z4 = z[4];
  const realT z5 = z[5];
  const realT z6 = z[6];
  const realT z7 = z[7];
  const realT twelfth = 1.0/12.0;

  b[0] = ( y1*((z5-z2)-(z3-z4))
           +y2*(z1-z3)
           +y3*((z2-z7)-(z4-z1))
           +y4*((z7-z5)-(z1-z3))
           +y5*(z4-z1)
           +y7*(z3-z4) )*twelfth;

  b[1] = ( y2*((z6-z3)-(z0-z5))
           +y3*(z2-z0)
           +y0*((z3-z4)-(z5-z2))
           +y5*((z4-z6)-(z2-z0))
           +y6*(z5-z2)
           +y4*(z0-z5) )*twelfth;

  b[2] = ( y3*((z7-z0)-(z1-z6))
           +y0*(z3-z1)
           +y1*((z0-z5)-(z6-z3))
           +y6*((z5-z7)-(z3-z1))
           +y7*(z6-z3)
           +y5*(z1-z6))*twelfth;

  b[3] = ( y0*((z4-z1)-(z2-z7))
           +y1*(z0-z2)
           +y2*((z1-z6)-(z7-z0))
           +y7*((z6-z4)-(z0-z2))
           +y4*(z7-z0)
           +y6*(z2-z7))*twelfth;

  b[4] = ( y7*((z3-z6)-(z5-z0))
           +y6*(z7-z5)
           +y5*((z6-z1)-(z0-z7))
           +y0*((z1-z3)-(z7-z5))
           +y3*(z0-z7)
           +y1*(z5-z0))*twelfth;

  b[5] = ( y4*((z0-z7)-(z6-z1))
           +y7*(z4-z6)
           +y6*((z7-z2)-(z1-z4))
           +y1*((z2-z0)-(z4-z6))
           +y0*(z1-z4)
           +y2*(z6-z1))*twelfth;

  b[6] = ( y5*((z1-z4)-(z7-z2))
           +y4*(z5-z7)
           +y7*((z4-z3)-(z2-z5))
           +y2*((z3-z1)-(z5-z7))
           +y1*(z2-z5)
           +y3*(z7-z2))*twelfth;

  b[7] = ( y6*((z2-z5)-(z4-z3))
           +y5*(z6-z4)
           +y4*((z5-z0)-(z3-z6))
           +y3*((z0-z2)-(z6-z4))
           +y2*(z3-z6)
           +y0*(z4-z3))*twelfth;
}


/*
   static void CalculateFBHourGlassModes( const realT* const x,
                                       const realT* const y,
                                       const realT* const z,
                                       const realT b[][8],
                                       const realT vol,
                                       realT hourgamma[][8] )
   {
   const realT gamma[4][8] =  { { 1,  1, -1, -1, -1, -1, 1,  1},
                                { 1, -1, -1,  1, -1,  1, 1, -1},
                                { 1, -1,  1, -1,  1, -1, 1, -1},
                                {-1,  1, -1,  1,  1, -1, 1, -1} };

   const realT volinv = 1.0 / vol;;


   // compute the hourglass modes
   for( int mode=0 ; mode<4 ; ++mode )
   {
    const realT hourmodx = x[0] * gamma[mode][0]
 + x[1] * gamma[mode][1]
 + x[2] * gamma[mode][2]
 + x[3] * gamma[mode][3]
 + x[4] * gamma[mode][4]
 + x[5] * gamma[mode][5]
 + x[6] * gamma[mode][6]
 + x[7] * gamma[mode][7] ;

    const realT hourmody = y[0] * gamma[mode][0]
 + y[1] * gamma[mode][1]
 + y[2] * gamma[mode][2]
 + y[3] * gamma[mode][3]
 + y[4] * gamma[mode][4]
 + y[5] * gamma[mode][5]
 + y[6] * gamma[mode][6]
 + y[7] * gamma[mode][7] ;

    const realT hourmodz = z[0] * gamma[mode][0]
 + z[1] * gamma[mode][1]
 + z[2] * gamma[mode][2]
 + z[3] * gamma[mode][3]
 + z[4] * gamma[mode][4]
 + z[5] * gamma[mode][5]
 + z[6] * gamma[mode][6]
 + z[7] * gamma[mode][7] ;

    for( int a=0 ; a<8 ; ++a )
    {
      hourgamma[mode][a] = gamma[mode][a] - volinv * ( b[0][a] * hourmodx
 + b[1][a] * hourmody
 + b[2][a] * hourmodz );
    }
   }
   }


   static void
   CalcFBHourForce( const realT* const xd,
                 const realT* const yd,
                 const realT* const zd,
                 const realT hgamma[][8],
                 const realT coefficient,
                 realT* const hgfx,
                 realT* const hgfy,
                 realT* const hgfz )
   {

   realT hx[4] = { 0.0, 0.0, 0.0, 0.0 } ;
   realT hy[4] = { 0.0, 0.0, 0.0, 0.0 } ;
   realT hz[4] = { 0.0, 0.0, 0.0, 0.0 } ;

   for( int mode=0 ; mode<4 ; ++mode )
   {

    for( int a=0 ; a<8 ; ++a )
    {
      hx[mode] += hgamma[mode][a] * xd[a];
      hy[mode] += hgamma[mode][a] * yd[a];
      hz[mode] += hgamma[mode][a] * zd[a];
    }

   }

   for( int a=0 ; a<8 ; ++a )
   {

    hgfx[a] = coefficient * ( hgamma[0][a] * hx[0]
 + hgamma[1][a] * hx[1]
 + hgamma[2][a] * hx[2]
 + hgamma[3][a] * hx[3] );

    hgfy[a] = coefficient * ( hgamma[0][a] * hy[0]
 + hgamma[1][a] * hy[1]
 + hgamma[2][a] * hy[2]
 + hgamma[3][a] * hy[3] );

    hgfz[a] = coefficient * ( hgamma[0][a] * hz[0]
 + hgamma[1][a] * hz[1]
 + hgamma[2][a] * hz[2]
 + hgamma[3][a] * hz[3] );


   }
   }*/

#endif /* INTEGRATIONRULET_H_ */
