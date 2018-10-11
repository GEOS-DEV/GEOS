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
 * @file Kinematics.cpp
 * @author settgast1
 * @date Dec 20, 2010
 */


#include "Kinematics.h"
namespace geosx
{
void IncrementalKinematics( const R2TensorT<3>& A,
                            R2SymTensorT<3>& Dadt,
                            R2TensorT<3>& Rhat )
{
  R2SymTensorT<3> C_I;
  R2SymTensorT<3> C_I2;
  //static R2SymTensorT<3> C_I3 ;
//  static R2TensorT<3> Amod ;
//  static R2TensorT<3> UhatA ;

  C_I.AijAkj_m_Aik_m_Aki(A);
  C_I2.AijAjk(C_I);
  //C_I3.AijBjk(C_I2,C_I);

  C_I *= -0.5;
  C_I2 *= 0.25;
  //C_I3 *= -1.0/6.0;

  Dadt = C_I;
  Dadt += C_I2;
  //Dadt += C_I3;


/*
   // ***** Precondition A for production of Rhat *****

   // Estimate Uhat-I store in C_I
   C_I2 *= 1.5;
   C_I3 *= 15.0/8.0;

   C_I += C_I2;
   C_I += C_I3;

   // Add I-Uhat to Amod
   Amod = C_I;
   Amod *= -1.0;

   // Make Uhat
   C_I.PlusIdentity(1.0);

   //store Uhat*A
   UhatA.AijBjk(C_I,A);

   // Construct Modified A
   Amod += UhatA;

   IncrementalRotation( Amod , Rhat );
 */
  IncrementalRotation( A, Rhat );


}

void IncrementalRotation( const R2TensorT<3>& A,
                          R2TensorT<3>& Rot )
{
//  realT alpha[3];
  R1TensorT<3> alpha;

  alpha.eijkAjk(A);
//  alpha *=-1.0;

  realT Q = 0.25 * Dot(alpha,alpha);


  realT trA = A.Trace();
  realT trFhatinv_1 = ( 2.0 - trA );
  realT P = 0.25 * pow( trFhatinv_1, 2 );

  realT P2 = pow(P,2);
  realT P3 = P2*P;

  realT Q2 = pow(Q,2);

  realT PpQ = P + Q;
  realT PpQ2 = pow(PpQ,2);
  realT PpQ3 = PpQ2*PpQ;

  realT term1 = sqrt( std::max( P + 3*P2*(1-PpQ)/PpQ2  - 2*P3*(1-PpQ)/PpQ3, 0.0 ) );

  realT term2 = 0;

  if( fabs(Q) > 0.01 )
  {
    if( trFhatinv_1 > 0.0 )
      term2 = (1.0 - trFhatinv_1 / fabs(trFhatinv_1) * term1) / ( 4*Q );
    else
      term2 = (1.0) / ( 4*Q );

  }
  else
    term2 = 0.125 + Q * ( P2-12.0*(P-1.0) ) / (32 * P2)
            + Q2*( (P-2.0)*(P2-10*P+32) ) / (64*P3);

  realT term3 = 0.5 * sqrt( ( P*Q*(3.0-Q) + P3 + Q2 ) / PpQ3 );

  Rot.dyadic_aa(alpha);
  Rot *= term2;
  Rot.PlusIdentity(term1);

  alpha *= term3;

  Rot(0,1) += alpha(2);
  Rot(0,2) -= alpha(1);
  Rot(1,2) += alpha(0);
  Rot(1,0) -= alpha(2);
  Rot(2,0) += alpha(1);
  Rot(2,1) -= alpha(0);

/*
    Rot.t_data[1] += alpha.t_data[2];
    Rot.t_data[2] -= alpha.t_data[1];
    Rot.t_data[5] += alpha.t_data[0];
    Rot.t_data[3] -= alpha.t_data[2];
    Rot.t_data[6] += alpha.t_data[1];
    Rot.t_data[7] -= alpha.t_data[0];
 */

}



void CalculatePhantomGradient( R2TensorT<3>& Gradient,
                               const int* bConnectivity,
                               const array1d<R1TensorT<3> >& disp,
                               const array2d<R1TensorT<3> >& dNdX )

{
  Gradient = 0.0;
  for( int side=1 ; side<=2 ; ++side )
    for( int a=1 ; a<=4 ; ++a )
      Gradient.plus_dyadic_ab( disp(bConnectivity[(side-1)*4+a-1]), dNdX(side,a));

}

}
