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

#ifndef __MATRIX_MATH_HPP_
#define __MATRIX_MATH_HPP_

#include "RAJA/RAJA.hpp"
#include "Layout.hpp"

struct P_Wrapper{
  real64 array[8][3][8];
  RAJA_HOST_DEVICE RAJA_INLINE real64 &operator()(localIndex e, localIndex a, localIndex c){
    return array[e][a][c];
  }  
};

struct Sym_Mat
{
  real64 array[6];

  //Accessor for non diagonal entries
  RAJA_HOST_DEVICE real64 &operator()(localIndex r, localIndex c){
    return array[c + 2*r];
  }

  //Accesor for diagonal entries;
  RAJA_HOST_DEVICE real64 &operator()(localIndex i){
    return array[2*i];
  }
};

template<typename T>
RAJA_HOST_DEVICE
RAJA_INLINE real64 det(T F[local_dim][local_dim]){

  real64 detF =  F[0][0]*(F[1][1]*F[2][2]-F[2][1]*F[1][2])
    - F[0][1]*(F[1][0]*F[2][2]-F[1][2]*F[2][0])
    + F[0][2]*(F[1][0]*F[2][1]-F[1][1]*F[2][0]);

  return detF;
}


template<typename T>
RAJA_HOST_DEVICE
RAJA_INLINE void Finverse(T F[local_dim][local_dim], T Finv[local_dim][local_dim] ){

  real64 detF = det<real64>(F);
  real64 idetF = (1.0/detF);

  Finv[0][0] =  idetF*(F[1][1]*F[2][2]-F[2][1]*F[1][2]);
  Finv[0][1] = -idetF*(F[0][1]*F[2][2]-F[0][2]*F[2][1]);
  Finv[0][2] =  idetF*(F[0][1]*F[1][2]-F[0][2]*F[1][1]);

  Finv[1][0] = -idetF*(F[1][0]*F[2][2]-F[1][2]*F[2][0]);
  Finv[1][1] =  idetF*(F[0][0]*F[2][2]-F[0][2]*F[2][0]);
  Finv[1][2] = -idetF*(F[0][0]*F[1][2]-F[1][0]*F[0][2]);

  Finv[2][0] =  idetF*(F[1][0]*F[2][1]-F[2][0]*F[1][1]);
  Finv[2][1] = -idetF*(F[0][0]*F[2][1]-F[2][0]*F[0][1]);
  Finv[2][2] =  idetF*(F[0][0]*F[1][1]-F[1][0]*F[0][1]);
}

template <typename T>
RAJA_HOST_DEVICE
RAJA_INLINE void AijBjk (T A[local_dim][local_dim], T B[local_dim][local_dim], T C [local_dim][local_dim] ) {

  for(localIndex row=0; row<local_dim; ++row) {
    for(localIndex col=0; col<local_dim; ++col) {
      T dot=0.;
      for(localIndex c=0; c<local_dim; ++c){
        dot += A[row][c]*B[c][col];
      }

      C[row][col] = dot;

    }
  }
}

#if 0
//symmetric matrix multiply
template <typename T>
RAJA_HOST_DEVICE
RAJA_INLINE void SAijBjk (Sym_Mat A, T B[local_dim][local_dim], T C [local_dim][local_dim] ) {

  C[0][0] = A(0) * B[0][0]
    + A(0,1) * B[1][0]
    + A(0,2) * B[2][0];
  
  C[0][1] = A(0) * B[0][1]
    + A(0,1) * B[1][1]
    + A(0,2) * B[2][1];

  C[0][2] = A(0) * B[0][2]
    + A(0,1) * B[1][2]
    + A(0,2) * B[2][2];
  
  C[1][0] = A(0,1) * B[0][0]
    + A(1)   * B[1][0]
    + A(1,2) * B[2][0];
  
  C[1][1] = A(0,1) * B[0][1]
    + A(1)   * B[1][1]
    + A(1,2) * B[2][1];
  
  C[1][2] = A(0,1) * B[0][2]
    + A(1)   * B[1][2]
    + A(1,2) * B[2][2];

  
  C[2][0] = A(0,2) * B[0][0]
    + A(1,2) * B[1][0]
    + A(2)   * B[2][0];

  C(2,1) = A(0,2) * B[0][1]
    + A(1,2) * B[1][1]
    + A(2) * B[2][1];

  C(2,2) = A(0,2) * B[0][2]
    + A(1,2) * B[1][2]
    + A(2)   * B[2][2];  
}
#endif


//matrix multiply with tranposed B                                                                                               
template<typename T>
RAJA_HOST_DEVICE
RAJA_INLINE void AijBkj (T A[local_dim][local_dim], T B[local_dim][local_dim], T C[local_dim][local_dim] ) {

  for(localIndex row=0; row<local_dim; ++row){
    for(localIndex col=0; col<local_dim; ++col){

      T dot=0.0;
      for(localIndex c=0; c<local_dim; ++c){
        dot += A[row][c]*B[col][c];
      }

      C[row][col] = dot;
    }
  }
}
                             

//Copy Global to Local;
RAJA_HOST_DEVICE
RAJA_INLINE void GlobalToLocal(const localIndex nodeList[8], localIndex k, 
                               real64 u_local[24], real64 uhat_local[24],
                               geosxData iu, geosxData iuhat)
{  
  for(localIndex a=0; a<inumNodesPerElement; ++a)
    {
      localIndex id = nodeList[a];
      for(localIndex i=0; i<local_dim; ++i)
        {
          u_local[i + local_dim*a] = iu[i+local_dim*id];
          uhat_local[i + local_dim*a] = iuhat[i + local_dim*id];
        }
    }    
}

//Copy Global to Local;
RAJA_HOST_DEVICE
RAJA_INLINE void GlobalToLocal(const localIndex nodeList[8], localIndex k, 
                               real64 u_local_x[8], real64 u_local_y[8],real64 u_local_z[8],
                               real64 uhat_local_x[8], real64 uhat_local_y[8], real64 uhat_local_z[8],
                               geosxData iu_x, geosxData iu_y, geosxData iu_z,
                               geosxData iuhat_x, geosxData iuhat_y, geosxData iuhat_z)
{  
  for(localIndex a=0; a<inumNodesPerElement; ++a)
    {
      localIndex id = nodeList[a];

      u_local_x[a] = iu_x[id];
      u_local_y[a] = iu_y[id];
      u_local_z[a] = iu_z[id];
      
      uhat_local_x[a] = iuhat_x[id];
      uhat_local_y[a] = iuhat_y[id];
      uhat_local_z[a] = iuhat_z[id];    
    }    
}


//Integrate Function
RAJA_HOST_DEVICE
RAJA_INLINE void Integrate(real64 f_local[local_dim*inumNodesPerElement],real64 idetJ, real64 detF, real64 Finv[local_dim][local_dim],
                           real64 TotalStress[local_dim][local_dim], geosxData idNdX, localIndex k, localIndex q, localIndex noElem)
  
{
  //---------[Integrate - Function]---------------------
  real64 const integrationFactor = idetJ*detF;
  real64 P[local_dim][local_dim];
  
  AijBkj(TotalStress,Finv,P);
  for(localIndex ty=0; ty<local_dim; ++ty){
    for(localIndex tx=0; tx<local_dim; ++tx){
      P[ty][tx] *= integrationFactor;
    }
  } 
  
  //--------------------------------------------------
  for(int a=0; a<inumNodesPerElement; ++a){
    
    f_local[0 + local_dim*a] -= P[0][0]*idNdX(k,q,a,0) + P[0][1]*idNdX(k,q,a,1) + P[0][2]*idNdX(k,q,a,2);
    f_local[1 + local_dim*a] -= P[1][0]*idNdX(k,q,a,0) + P[1][1]*idNdX(k,q,a,1) + P[1][2]*idNdX(k,q,a,2);
    f_local[2 + local_dim*a] -= P[2][0]*idNdX(k,q,a,0) + P[2][1]*idNdX(k,q,a,1) + P[2][2]*idNdX(k,q,a,2);
  }
}

//Integrate Function
RAJA_HOST_DEVICE
RAJA_INLINE void Integrate(real64 f_local_x[inumNodesPerElement], real64 f_local_y[inumNodesPerElement], real64 f_local_z[inumNodesPerElement],
                           real64 idetJ, real64 detF, real64 Finv[local_dim][local_dim],
                           real64 TotalStress[local_dim][local_dim], geosxData idNdX_x, geosxData idNdX_y, geosxData idNdX_z,
                           localIndex k, localIndex q, localIndex noElem)
  
{
  //---------[Integrate - Function]---------------------
  real64 const integrationFactor = idetJ*detF;
  real64 P[local_dim][local_dim];
  
  AijBkj(TotalStress,Finv,P);
  for(localIndex ty=0; ty<local_dim; ++ty){
    for(localIndex tx=0; tx<local_dim; ++tx){
      P[ty][tx] *= integrationFactor;
    }
  } 
  
  //--------------------------------------------------
  for(int a=0; a<inumNodesPerElement; ++a){
    
    f_local_x[a] -= P[0][0]*idNdX_x(k,q,a) + P[0][1]*idNdX_y(k,q,a) + P[0][2]*idNdX_z(k,q,a);
    f_local_y[a] -= P[1][0]*idNdX_x(k,q,a) + P[1][1]*idNdX_y(k,q,a) + P[1][2]*idNdX_z(k,q,a);
    f_local_z[a] -= P[2][0]*idNdX_x(k,q,a) + P[2][1]*idNdX_y(k,q,a) + P[2][2]*idNdX_z(k,q,a);
  }
}

//Integrate Function
RAJA_HOST_DEVICE
RAJA_INLINE void Integrate(real64 f_local_x[inumNodesPerElement], real64 f_local_y[inumNodesPerElement], real64 f_local_z[inumNodesPerElement],
                           real64 idetJ, real64 detF, geosxData Finv_ptr,
                           real64 TotalStress[local_dim][local_dim], geosxData idNdX_x, geosxData idNdX_y, geosxData idNdX_z,
                           localIndex k, localIndex q, localIndex noElem)
  
{
  real64 const integrationFactor = idetJ*detF;
  real64 P[local_dim][local_dim];
  real64 Finv[local_dim][local_dim];
  
  for(localIndex r=0; r<local_dim; ++r){
    for(localIndex c=0; c<local_dim; ++c){
      Finv[r][c] = Finv_ptr(k,q,r,c);
    }
  }
  
  AijBkj(TotalStress,Finv,P);
  for(localIndex ty=0; ty<local_dim; ++ty){
    for(localIndex tx=0; tx<local_dim; ++tx){
      P[ty][tx] *= integrationFactor;
    }
  } 
  
  //--------------------------------------------------
  for(int a=0; a<inumNodesPerElement; ++a){
    f_local_x[a] -= P[0][0]*idNdX_x(k,q,a) + P[0][1]*idNdX_y(k,q,a) + P[0][2]*idNdX_z(k,q,a);
    f_local_y[a] -= P[1][0]*idNdX_x(k,q,a) + P[1][1]*idNdX_y(k,q,a) + P[1][2]*idNdX_z(k,q,a);
    f_local_z[a] -= P[2][0]*idNdX_x(k,q,a) + P[2][1]*idNdX_y(k,q,a) + P[2][2]*idNdX_z(k,q,a);
  }
  
}


//Integrate Function
RAJA_HOST_DEVICE
RAJA_INLINE void Integrate(real64 f_local[local_dim * inumNodesPerElement],
                           real64 idetJ, real64 detF, geosxData Finv_ptr,
                           real64 TotalStress[local_dim][local_dim], geosxData idNdX,
                           localIndex k, localIndex q, localIndex noElem)
  
{
  real64 const integrationFactor = idetJ*detF;
  real64 P[local_dim][local_dim];
  real64 Finv[local_dim][local_dim];
  
  for(localIndex r=0; r<local_dim; ++r){
    for(localIndex c=0; c<local_dim; ++c){
      Finv[r][c] = Finv_ptr(k,q,r,c);
    }
  }
  
  AijBkj(TotalStress,Finv,P);
  for(localIndex ty=0; ty<local_dim; ++ty){
    for(localIndex tx=0; tx<local_dim; ++tx){
      P[ty][tx] *= integrationFactor;
    }
  } 
  
  //--------------------------------------------------
  for(int a=0; a<inumNodesPerElement; ++a){
    f_local[0 + local_dim*a] -= P[0][0]*idNdX(k,q,a,0) + P[0][1]*idNdX(k,q,a,1) + P[0][2]*idNdX(k,q,a,2);
    f_local[1 + local_dim*a] -= P[1][0]*idNdX(k,q,a,0) + P[1][1]*idNdX(k,q,a,1) + P[1][2]*idNdX(k,q,a,2);
    f_local[2 + local_dim*a] -= P[2][0]*idNdX(k,q,a,0) + P[2][1]*idNdX(k,q,a,1) + P[2][2]*idNdX(k,q,a,2);
  }
  
}


//Integrate Functions
RAJA_HOST_DEVICE
RAJA_INLINE void Integrate(real64 f_local[local_dim*inumNodesPerElement],real64 idetJ, real64 detF, real64 Finv[local_dim][local_dim],
                           real64 TotalStress[local_dim][local_dim], real64 dNdX[inumQuadraturePoints][inumNodesPerElement][local_dim],
                           localIndex q, localIndex noElem)
  
{
  //---------[Integrate - Function]---------------------
  real64 const integrationFactor = idetJ*detF;
  real64 P[local_dim][local_dim];
  
  AijBkj(TotalStress,Finv,P);
  for(localIndex ty=0; ty<local_dim; ++ty){
    for(localIndex tx=0; tx<local_dim; ++tx){
      P[ty][tx] *= integrationFactor;
    }
  } 
  
  //--------------------------------------------------
  for(int a=0; a<inumNodesPerElement; ++a){
    
    f_local[0 + local_dim*a] -= P[0][0]*dNdX[q][a][0] + P[0][1]*dNdX[q][a][1] + P[0][2]*dNdX[q][a][2];
    f_local[1 + local_dim*a] -= P[1][0]*dNdX[q][a][0] + P[1][1]*dNdX[q][a][1] + P[1][2]*dNdX[q][a][2];
    f_local[2 + local_dim*a] -= P[2][0]*dNdX[q][a][0] + P[2][1]*dNdX[q][a][1] + P[2][2]*dNdX[q][a][2];
  }
}

//Integrate Functions
RAJA_HOST_DEVICE
RAJA_INLINE void Integrate(real64 f_local_x[inumNodesPerElement],real64 f_local_y[inumNodesPerElement],real64 f_local_z[inumNodesPerElement],
                           real64 idetJ, real64 detF, real64 Finv[local_dim][local_dim],
                           real64 TotalStress[local_dim][local_dim],real64 dNdX_x[inumQuadraturePoints][inumNodesPerElement],
                           real64 dNdX_y[inumQuadraturePoints][inumNodesPerElement],real64 dNdX_z[inumQuadraturePoints][inumNodesPerElement],
                           localIndex q) 
{
  //---------[Integrate - Function]---------------------
  real64 const integrationFactor = idetJ*detF;
  real64 P[local_dim][local_dim];
  
  AijBkj(TotalStress,Finv,P);
  for(localIndex ty=0; ty<local_dim; ++ty){
    for(localIndex tx=0; tx<local_dim; ++tx){
      P[ty][tx] *= integrationFactor;
    }
  } 
  
  //--------------------------------------------------
  for(int a=0; a<inumNodesPerElement; ++a)
    {    
      f_local_x[a] -= P[0][0]*dNdX_x[q][a] + P[0][1]*dNdX_y[q][a] + P[0][2]*dNdX_z[q][a];
      f_local_y[a] -= P[1][0]*dNdX_x[q][a] + P[1][1]*dNdX_y[q][a] + P[1][2]*dNdX_z[q][a];
      f_local_z[a] -= P[2][0]*dNdX_x[q][a] + P[2][1]*dNdX_y[q][a] + P[2][2]*dNdX_z[q][a];
    }
}


template<typename atomicPol, typename T>
RAJA_HOST_DEVICE
RAJA_INLINE void AddLocalToGlobal(T nodeList, real64 f_local[local_dim*inumNodesPerElement], geosxData iacc)
{
  for(localIndex a=0; a<inumNodesPerElement; ++a)
    {          
      localIndex id = nodeList[a];
      RAJA::atomic::atomicAdd<atomicPol>(&iacc[0 + local_dim*id],f_local[0 + local_dim*a]);
      RAJA::atomic::atomicAdd<atomicPol>(&iacc[1 + local_dim*id],f_local[1 + local_dim*a]);
      RAJA::atomic::atomicAdd<atomicPol>(&iacc[2 + local_dim*id],f_local[2 + local_dim*a]);
    }
}


template<typename atomicPol, typename T>
RAJA_HOST_DEVICE
RAJA_INLINE void AddLocalToGlobal(T nodeList,
                                  real64 f_local_x[inumNodesPerElement],real64 f_local_y[inumNodesPerElement],real64 f_local_z[inumNodesPerElement],
                                  geosxData iacc_x, geosxData iacc_y, geosxData iacc_z)
{
  for(localIndex a=0; a<inumNodesPerElement; ++a)
    {          
      localIndex id = nodeList[a];
      RAJA::atomic::atomicAdd<atomicPol>(&iacc_x[id],f_local_x[a]);
      RAJA::atomic::atomicAdd<atomicPol>(&iacc_y[id],f_local_y[a]);
      RAJA::atomic::atomicAdd<atomicPol>(&iacc_z[id],f_local_z[a]);
    }
}


RAJA_HOST_DEVICE
RAJA_INLINE void CalculateGradient(real64 dUdX[3][3], real64 u_local[24],
                                   geosxData idNdX, localIndex k, localIndex q, localIndex noElem){
  for(localIndex a=0; a<inumNodesPerElement; ++a)
    {              
      dUdX[0][0] += u_local[0 + local_dim*a]*idNdX(k,q,a,0);
      dUdX[0][1] += u_local[0 + local_dim*a]*idNdX(k,q,a,1);
      dUdX[0][2] += u_local[0 + local_dim*a]*idNdX(k,q,a,2);
      
      dUdX[1][0] += u_local[1 + local_dim*a]*idNdX(k,q,a,0);
      dUdX[1][1] += u_local[1 + local_dim*a]*idNdX(k,q,a,1);
      dUdX[1][2] += u_local[1 + local_dim*a]*idNdX(k,q,a,2);
      
      dUdX[2][0] += u_local[2 + local_dim*a]*idNdX(k,q,a,0);
      dUdX[2][1] += u_local[2 + local_dim*a]*idNdX(k,q,a,1);
      dUdX[2][2] += u_local[2 + local_dim*a]*idNdX(k,q,a,2);
    }
    
}

RAJA_HOST_DEVICE
RAJA_INLINE void CalculateGradient(real64 dUdX[3][3],
                                   real64 u_local_x[8], real64 u_local_y[8], real64 u_local_z[8],
                                   geosxData idNdX_x,geosxData idNdX_y,geosxData idNdX_z,
                                   localIndex k, localIndex q, localIndex noElem){
  
  for(localIndex a=0; a<inumNodesPerElement; ++a)
    {              
      dUdX[0][0] += u_local_x[a]*idNdX_x(k,q,a);
      dUdX[0][1] += u_local_x[a]*idNdX_y(k,q,a);
      dUdX[0][2] += u_local_x[a]*idNdX_z(k,q,a);
      
      dUdX[1][0] += u_local_y[a]*idNdX_x(k,q,a);
      dUdX[1][1] += u_local_y[a]*idNdX_y(k,q,a);
      dUdX[1][2] += u_local_y[a]*idNdX_z(k,q,a);
      
      dUdX[2][0] += u_local_z[a]*idNdX_x(k,q,a);
      dUdX[2][1] += u_local_z[a]*idNdX_y(k,q,a);
      dUdX[2][2] += u_local_z[a]*idNdX_z(k,q,a);
    }
    
}

RAJA_HOST_DEVICE
RAJA_INLINE void CalculateGradient(real64 dUdX[3][3],
                                   real64 u_local_x[8], real64 u_local_y[8], real64 u_local_z[8],
                                   real64 idNdX_x[inumQuadraturePoints][inumNodesPerElement],
                                   real64 idNdX_y[inumQuadraturePoints][inumNodesPerElement],
                                   real64 idNdX_z[inumQuadraturePoints][inumNodesPerElement],
                                   localIndex k, localIndex q){
  
  for(localIndex a=0; a<inumNodesPerElement; ++a)
    {              
      dUdX[0][0] += u_local_x[a]*idNdX_x[q][a];
      dUdX[0][1] += u_local_x[a]*idNdX_y[q][a];
      dUdX[0][2] += u_local_x[a]*idNdX_z[q][a];
      
      dUdX[1][0] += u_local_y[a]*idNdX_x[q][a];
      dUdX[1][1] += u_local_y[a]*idNdX_y[q][a];
      dUdX[1][2] += u_local_y[a]*idNdX_z[q][a];
      
      dUdX[2][0] += u_local_z[a]*idNdX_x[q][a];
      dUdX[2][1] += u_local_z[a]*idNdX_y[q][a];
      dUdX[2][2] += u_local_z[a]*idNdX_z[q][a];
    }
    
}


RAJA_HOST_DEVICE
RAJA_INLINE void CalculateGradient(real64 dUdX[3][3], real64 u_local[24],
                                   real64 dNdX[inumQuadraturePoints][inumNodesPerElement][local_dim], localIndex q){

  for(localIndex a=0; a<inumNodesPerElement; ++a)
    {              
      dUdX[0][0] += u_local[0 + local_dim*a]*dNdX[q][a][0];
      dUdX[0][1] += u_local[0 + local_dim*a]*dNdX[q][a][1];
      dUdX[0][2] += u_local[0 + local_dim*a]*dNdX[q][a][2];
      
      dUdX[1][0] += u_local[1 + local_dim*a]*dNdX[q][a][0];
      dUdX[1][1] += u_local[1 + local_dim*a]*dNdX[q][a][1];
      dUdX[1][2] += u_local[1 + local_dim*a]*dNdX[q][a][2];
      
      dUdX[2][0] += u_local[2 + local_dim*a]*dNdX[q][a][0];
      dUdX[2][1] += u_local[2 + local_dim*a]*dNdX[q][a][1];
      dUdX[2][2] += u_local[2 + local_dim*a]*dNdX[q][a][2];
    }
}


#endif
