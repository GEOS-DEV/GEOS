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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/* File: Constants.h */

/*
 * created      : RRS
 */


#ifndef TENSOR_T_H_
#define TENSOR_T_H_
#include "R1TensorT.h"
#include "R2TensorT.h"
#include "R2SymTensorT.h"


using R1Tensor    = R1TensorT<3>;
using R2Tensor    = R2TensorT<3>;
using R2SymTensor = R2SymTensorT<3>;



template< int T_dim >
inline void Decompose( const R2TensorT<T_dim>& F, R2TensorT<T_dim>& R, R2SymTensorT<T_dim>& U );


//template< int T_dim >
//void Construct_Q( const R1TensorT<T_dim>& a , const R1TensorT<T_dim>& b );


template< int T_dim >
void Decompose( const R2TensorT<T_dim>& F, R2TensorT<T_dim>& R, R2SymTensorT<T_dim>& U )
{
  R2SymTensorT<T_dim> Uinv;
  U.AjiAjk(F);
  U.Sqrt();
  Uinv.Inverse(U);
  R.AijBjk(F,Uinv);

}



//  R1TensorT<3> R1ZeroTensorRank3;
//  R2TensorT<3> R2ZeroTensorRank3;
//  R2SymTensorT<3> R2SZeroTensorRank3;

// Construct a Proper Orthogonal Tensor which "rotates" a vector from the
// direction a to b.
/*template< int T_dim >
   void Construct_Q( const R1TensorT<T_dim>& a , const R1TensorT<T_dim>& b )
   {
   static R1TensorT<T_dim> p;
   static R1TensorT<T_dim> n_a;
   static R1TensorT<T_dim> n_b;

   static R2TensorT<T_dim> pp;
   static R2TensorT<T_dim> qq;
   static R2TensorT<T_dim> rr;
   static R2TensorT<T_dim> qr;

   n_a = a;
   n_a /= a.L2_Norm();
   n_b = b;
   n_b /= b.L2_Norm();

   p.Cross(n_a,n_b);

   qq.dyadic_aa(n_a);
   rr.dyadic_aa(n_b);

   }*/


#endif
