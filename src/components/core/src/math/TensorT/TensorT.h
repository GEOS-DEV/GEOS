//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security,
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

// Construct a Proper Orthogonal Tensor which "rotates" a vector from the direction a to b.
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
