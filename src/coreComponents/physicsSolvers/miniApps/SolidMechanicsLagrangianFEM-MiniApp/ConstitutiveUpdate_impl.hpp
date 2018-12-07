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

#ifndef __CONSTITUTIVE_UPDATE_HPP_
#define __CONSTITUTIVE_UPDATE_HPP_

#include "Layout.hpp"
#include "MatrixMath_impl.hpp"
#include "RAJA/RAJA.hpp"

template<typename T>
RAJA_HOST_DEVICE
RAJA_INLINE void HughesWinget(T Rot[LOCAL_DIM][LOCAL_DIM], T Dadt[LOCAL_DIM][LOCAL_DIM], T L[LOCAL_DIM][LOCAL_DIM], real64 dt){

  double Omega[LOCAL_DIM][LOCAL_DIM];
  //Omega = 0.5*(L-LT)
  //Dadt  = 0.5*(L+LT)

  for(localIndex ty=0; ty<LOCAL_DIM; ++ty)
    {
      for(localIndex tx=0; tx<LOCAL_DIM; ++tx)
        {
          Dadt[ty][tx] = 0.5*(L[ty][tx] + L[tx][ty]);
          Omega[ty][tx] = 0.5*(L[ty][tx] - L[tx][ty]);
        }
    }

  double IpR[LOCAL_DIM][LOCAL_DIM];
  double ImR[LOCAL_DIM][LOCAL_DIM];
  double ImRinv[LOCAL_DIM][LOCAL_DIM];
  
  for(localIndex ty=0; ty<LOCAL_DIM; ++ty)
    {
      for(localIndex tx=0; tx<LOCAL_DIM; ++tx)
        {
          IpR[ty][tx] =  0.5*Omega[ty][tx];
          ImR[ty][tx] = -0.5*Omega[ty][tx];
        }
    }
  
  for(localIndex tx=0; tx<LOCAL_DIM; ++tx)
    {
      IpR[tx][tx] = 1.0 + IpR[tx][tx];
      ImR[tx][tx] = 1.0 + ImR[tx][tx];
    }


  real64 detImR  = det(ImR);
  Finverse<real64> (ImR, ImRinv);
  
  AijBjk(ImRinv, ImR, Rot);

}


///Const update
#if defined(USE_GEOSX_ARRAY)
RAJA_HOST_DEVICE
RAJA_INLINE void UpdateStatePoint(real64 D[LOCAL_DIM][LOCAL_DIM], real64 Rot[LOCAL_DIM][LOCAL_DIM],
                                  localIndex m, localIndex q, globalIndex k,
                                  LvArray::ArrayView<real64,3,localIndex> idevStressData,
                                  LvArray::ArrayView<real64,1,localIndex> imeanStress,
                                  real64 shearModulus, real64 bulkModulus, localIndex noElem)
//                                  real64 shearModulus, real64 bulkModulus, localIndex noElem)
#elif defined(USE_RAJA_VIEW)
RAJA_HOST_DEVICE
RAJA_INLINE void UpdateStatePoint(real64 D[LOCAL_DIM][LOCAL_DIM], real64 Rot[LOCAL_DIM][LOCAL_DIM],
                                  localIndex m, localIndex q, globalIndex k,
                                  RAJA::View<real64,RAJA::Layout<3,localIndex,2>> idevStressData,
                                  RAJA::View<real64,RAJA::Layout<1,localIndex,0>> imeanStress,
                                  real64 shearModulus, real64 bulkModulus, localIndex noElem)
#else
RAJA_HOST_DEVICE
RAJA_INLINE void UpdateStatePoint(real64 D[LOCAL_DIM][LOCAL_DIM], real64 Rot[LOCAL_DIM][LOCAL_DIM],
                                  localIndex m, localIndex q, globalIndex k, geosxData idevStressData,
                                  geosxData imeanStress,
                                  real64 shearModulus, real64 bulkModulus, localIndex noElem)
#endif
{

  real64 volumeStrain = D[0][0] + D[1][1] + D[2][2];
  imeanStress(m) += volumeStrain * bulkModulus;
  //imeanStress(k,q) += volumeStrain * bulkModulus;  

  real64 temp[LOCAL_DIM][LOCAL_DIM];
  for(localIndex i=0; i<3; ++i)
    {
      for(localIndex j=0; j<3; ++j)
        {
          temp[i][j] = D[i][j];
        }
      temp[i][i] -= volumeStrain / 3.0;
    }

  for(localIndex ty=0; ty<3; ++ty)
    {
      for(localIndex tx=0; tx<3; ++tx)
        {
          temp[ty][tx] *= 2.0 * shearModulus;
        }
    }

  real64 localDevStress[LOCAL_DIM][LOCAL_DIM];

  //
  idevStressData(k,q,0) += temp[0][0];
  idevStressData(k,q,1) += temp[1][0];
  idevStressData(k,q,2) += temp[1][1];
  idevStressData(k,q,3) += temp[2][0];
  idevStressData(k,q,4) += temp[2][1];
  idevStressData(k,q,5) += temp[2][2];

  localDevStress[0][0] = idevStressData(k,q,0);
  localDevStress[1][0] = idevStressData(k,q,1);
  localDevStress[1][1] = idevStressData(k,q,2);
  localDevStress[2][0] = idevStressData(k,q,3);
  localDevStress[2][1] = idevStressData(k,q,4);
  localDevStress[2][2] = idevStressData(k,q,5);

  localDevStress[0][1] = localDevStress[1][0];
  localDevStress[0][2] = localDevStress[2][0];
  localDevStress[1][2] = localDevStress[2][1];


 //QijAjkQlk
  AijBjk(Rot,localDevStress,temp);
  AijBkj(temp,Rot,localDevStress);


  idevStressData(k,q,0) = localDevStress[0][0];
  idevStressData(k,q,1) = localDevStress[1][0];
  idevStressData(k,q,2) = localDevStress[1][1];
  idevStressData(k,q,3) = localDevStress[2][0];
  idevStressData(k,q,4) = localDevStress[2][1];
  idevStressData(k,q,5) = localDevStress[2][2];

}

//----------
RAJA_HOST_DEVICE
RAJA_INLINE void structuredElemToNodes(localIndex nodeList[8], localIndex k, localIndex nx, localIndex ny, localIndex nz)
{
  
  //Convert 1D Index to 3D                                                                                    
  localIndex x = k % nx;
  localIndex y = (k / nx) % ny;
  localIndex z = k / (nx * ny);
  
  //Compute node ids
  nodeList[0] = x + nx*(y + z*nx);
  nodeList[1] = (x+1) + nx*(y + z*nx);
  nodeList[2] = (x+1) + nx*((y+1) + z*nx);
  nodeList[3] = x + nx*((y+1) + z*nx);
  
  nodeList[4] = x + nx*(y + (z+1)*nx);
  nodeList[5] = (x+1) + nx*(y + (z+1)*nx);
  nodeList[6] = (x+1) + nx*((y+1) + (z+1)*nx);
  nodeList[7] = x + nx*((y+1) + (z+1)*nx);  
}


//----------
//Setup function pointers for stateUpdate

#if defined(USE_GEOSX_ARRAY)

typedef void (*constUpdate)(real64 D[LOCAL_DIM][LOCAL_DIM], real64 Rot[LOCAL_DIM][LOCAL_DIM],
                            localIndex m, localIndex q, globalIndex k, 
                            LvArray::ArrayView<real64,3,localIndex> devStressData,
                            LvArray::ArrayView<real64,1,localIndex> meanStress,
                            //real64 * meanStress,
                            real64 shearModulus, real64 bulkModulus, localIndex NoElem);

#if defined (USE_GPU)
__device__ constUpdate deviceUpdate = UpdateStatePoint;
#endif

#elif defined(USE_RAJA_VIEW)

typedef void (*constUpdate)(real64 D[LOCAL_DIM][LOCAL_DIM], real64 Rot[LOCAL_DIM][LOCAL_DIM],
                            localIndex m, localIndex q, globalIndex k, 
                            RAJA::View<real64,RAJA::Layout<3,localIndex,2>> devStressData,
                            RAJA::View<real64,RAJA::Layout<1,localIndex,0>> meanStress,
                            //real64 * meanStress,
                            real64 shearModulus, real64 bulkModulus, localIndex NoElem);

#if defined (USE_GPU)
__device__ constUpdate deviceUpdate = UpdateStatePoint;
#endif


#else
//Created a type 
typedef void (*constUpdate)(real64 D[LOCAL_DIM][LOCAL_DIM], real64 Rot[LOCAL_DIM][LOCAL_DIM],
                            localIndex m, localIndex q, globalIndex k,
                            geosxData devStressData, geosxData meanStress,
                            real64 shearModulus, real64 bulkModulus, localIndex NoElem);
#if defined (USE_GPU)
__device__ constUpdate deviceUpdate = UpdateStatePoint;
#endif

#endif



#endif
