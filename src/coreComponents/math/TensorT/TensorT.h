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

/* File: Constants.h */

/*
 * created      : RRS
 */


#ifndef TENSOR_T_H_
#define TENSOR_T_H_
#include "R1TensorT.h"
#include "R2TensorT.h"
#include "R2SymTensorT.h"


using R1Tensor    = R1TensorT< 3 >;
using R2Tensor    = R2TensorT< 3 >;
using R2SymTensor = R2SymTensorT< 3 >;



template< int T_dim >
inline void Decompose( const R2TensorT< T_dim > & F, R2TensorT< T_dim > & R, R2SymTensorT< T_dim > & U );


//template< int T_dim >
//void Construct_Q( const R1TensorT<T_dim>& a , const R1TensorT<T_dim>& b );


template< int T_dim >
void Decompose( const R2TensorT< T_dim > & F, R2TensorT< T_dim > & R, R2SymTensorT< T_dim > & U )
{
  R2SymTensorT< T_dim > Uinv;
  U.AjiAjk( F );
  U.Sqrt();
  Uinv.Inverse( U );
  R.AijBjk( F, Uinv );

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
