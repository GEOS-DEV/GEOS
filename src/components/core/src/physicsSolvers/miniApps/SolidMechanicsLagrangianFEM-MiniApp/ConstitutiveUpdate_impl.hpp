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
RAJA_INLINE void HughesWinget(T Rot[local_dim][local_dim], T Dadt[local_dim][local_dim], T L[local_dim][local_dim], real64 dt){

  double Omega[local_dim][local_dim];
  //Omega = 0.5*(L-LT)
  //Dadt  = 0.5*(L+LT)*dt

  for(localIndex ty=0; ty<local_dim; ++ty)
    {
      for(localIndex tx=0; tx<local_dim; ++tx)
        {
          Dadt[ty][tx] = 0.5*(L[ty][tx] + L[tx][ty])*dt;
          Omega[ty][tx] = 0.5*(L[ty][tx] - L[tx][ty]);
        }
    }

  double IpR[local_dim][local_dim];
  double ImR[local_dim][local_dim];
  double ImRinv[local_dim][local_dim];
  
  for(localIndex ty=0; ty<local_dim; ++ty)
    {
      for(localIndex tx=0; tx<local_dim; ++tx)
        {
          IpR[ty][tx] =  0.5*Omega[ty][tx];
          ImR[ty][tx] = -0.5*Omega[ty][tx];
        }
    }
  
  for(localIndex tx=0; tx<local_dim; ++tx)
    {
      IpR[tx][tx] = 1.0 + IpR[tx][tx];
      ImR[tx][tx] = 1.0 + ImR[tx][tx];
    }


  real64 detImR  = det(ImR);
  Finverse<real64> (ImR, ImRinv);
  
  AijBjk(ImRinv, ImR, Rot);
    
}


///Const update
RAJA_HOST_DEVICE
RAJA_INLINE void UpdateStatePoint(real64 D[local_dim][local_dim], real64 Rot[local_dim][local_dim],
                                  localIndex m, localIndex q, globalIndex k, geosxData idevStressData,
                                  geosxData imeanStress,
                                  real64 shearModulus, real64 bulkModulus, localIndex noElem)
{

  real64 volumeStrain = D[0][0] + D[1][1] + D[2][2];
  //imeanStress(k,q) += volumeStrain * bulkModulus;
  imeanStress[m] += volumeStrain * bulkModulus;
  
  real64 temp[local_dim][local_dim];
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

  real64 localDevStress[local_dim][local_dim];
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
  globalIndex x = k % nx;
  globalIndex y = (k / nx) % ny;
  globalIndex z = k / (nx * ny);
  
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

//Created a type 
typedef void (*constUpdate)(real64 D[local_dim][local_dim], real64 Rot[local_dim][local_dim],
                            localIndex m, localIndex q, globalIndex k, geosxData devStressData,
                            geosxData meanStress,
                            real64 shearModulus, real64 bulkModulus, localIndex NoElem);
#if defined (USE_CUDA)
__device__ constUpdate deviceUpdate = UpdateStatePoint;
#endif


#endif
