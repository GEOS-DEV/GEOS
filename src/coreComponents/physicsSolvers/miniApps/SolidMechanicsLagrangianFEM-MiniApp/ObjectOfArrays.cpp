/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

#include "RAJA/RAJA.hpp"
#include <iostream>
#include <cstring>

#include "MemoryManager_impl.hpp"
#include "RAJA/RAJA.hpp"
#include "MatrixMath_impl.hpp"
#include "Layout.hpp"
#include "MeshGen_impl.hpp"
#include "ConstitutiveUpdate_impl.hpp"
#include "ShapeFun_impl.hpp"
#include "omp.h"

#include "../../src/SolidMechanicsLagrangianFEMKernels_impl.hpp"

//
//Driver for the GEOSX proxy app
//Assumes an Objects of Array data layout
//

int main(int argc, char* const argv[])
{

  if(argc != 3)
    {
      std::cout<<"usage ./main NoElem Iter"<<std::endl;
      exit(-1);
    }

  size_t dataAllocated = 0.0;
  Index_type Kx      = atoi(argv[1]); 
  Index_type Niter   =  atoi(argv[2]);

  std::cout<<"GEOSX mini-app: Object of Arrays data structures"<<std::endl;
  printParameters();

  real64 dt = 0.125; //step size
  const Index_type nx = Kx+1; //Number of nodes in a cartesian dimension
  const Index_type NoElem = Kx*Kx*Kx; //Total number of elements
  const Index_type numNodes = nx*nx*nx; //Total number of nodes

  //
  //Generate an element list
  //
  Index_type * const elementList = memoryManager::allocate<Index_type>(NoElem, dataAllocated);
  for(Index_type i=0; i<NoElem; ++i) elementList[i] = i;
  
  //
  //Allocate space for constitutive map
  //
  Index_type *  constitutiveMap = memoryManager::allocate<Index_type>(inumQuadraturePoints*NoElem, dataAllocated);

  //
  //Generate space for an element to node list
  //
  Index_type  *  elemsToNodes = memoryManager::allocate<Index_type>(inumNodesPerElement*NoElem, dataAllocated);

  //
  //Allocate space for a list of vertices, generate a mesh, and populate the constitutive map
  //
  geosxData VX = memoryManager::allocate<real64>(numNodes*local_dim, dataAllocated);
  meshGen(VX, elemsToNodes,constitutiveMap,Kx);

  //
  //Precompute evaluate the quadrature points at every basis function
  //with respect to the parent basis function
  //
  P_Wrapper P; 
  generateP(P, inumNodesPerElement, inumQuadraturePoints);
  
  
  //
  //Allocate space for shape function derivatives and compute their values
  //
  geosxData dNdX_x = memoryManager::allocate<real64>(inumNodesPerElement*inumQuadraturePoints*NoElem, dataAllocated);
  geosxData dNdX_y = memoryManager::allocate<real64>(inumNodesPerElement*inumQuadraturePoints*NoElem, dataAllocated);
  geosxData dNdX_z = memoryManager::allocate<real64>(inumNodesPerElement*inumQuadraturePoints*NoElem, dataAllocated);
  make_dNdX(dNdX_x, dNdX_y, dNdX_z, VX, elemsToNodes, NoElem, inumQuadraturePoints, inumNodesPerElement);  

  //
  //Allocate space for nodal degrees of freedom as Object of Arrays
  //
  geosxData u_x = memoryManager::allocate<real64>(numNodes, dataAllocated);
  geosxData u_y = memoryManager::allocate<real64>(numNodes, dataAllocated);
  geosxData u_z = memoryManager::allocate<real64>(numNodes, dataAllocated);
  
  geosxData uhat_x = memoryManager::allocate<real64>(numNodes, dataAllocated);
  geosxData uhat_y = memoryManager::allocate<real64>(numNodes, dataAllocated);
  geosxData uhat_z = memoryManager::allocate<real64>(numNodes, dataAllocated);

  geosxData acc_x = memoryManager::allocate<real64>(numNodes, dataAllocated);
  geosxData acc_y = memoryManager::allocate<real64>(numNodes, dataAllocated);
  geosxData acc_z = memoryManager::allocate<real64>(numNodes, dataAllocated);

  std::memset(u_x,0,numNodes*sizeof(real64));
  std::memset(u_y,0,numNodes*sizeof(real64));
  std::memset(u_z,0,numNodes*sizeof(real64));

  std::memset(uhat_x,0,numNodes*sizeof(real64));
  std::memset(uhat_y,0,numNodes*sizeof(real64));
  std::memset(uhat_z,0,numNodes*sizeof(real64));

  std::memset(acc_x,0,numNodes*sizeof(real64));
  std::memset(acc_y,0,numNodes*sizeof(real64));
  std::memset(acc_z,0,numNodes*sizeof(real64));

  
  //
  //In the case of three kernel launches, we need to allocate
  //extra memory for intermediate results
  //
#if defined(THREE_KERNEL_UPDATE)
  geosxData Dadt = memoryManager::allocate<real64>(localMatSz*inumQuadraturePoints*NoElem, dataAllocated);
  geosxData Rot  = memoryManager::allocate<real64>(localMatSz*inumQuadraturePoints*NoElem, dataAllocated);
  geosxData detF  = memoryManager::allocate<real64>(inumQuadraturePoints*NoElem, dataAllocated);
  geosxData inverseF  = memoryManager::allocate<real64>(localMatSz*inumQuadraturePoints*NoElem, dataAllocated);
  
  std::memset(Dadt, 0, localMatSz*inumQuadraturePoints*NoElem*sizeof(real64));
  std::memset(Rot, 0, localMatSz*inumQuadraturePoints*NoElem*sizeof(real64));
  std::memset(detF, 0, inumQuadraturePoints*NoElem*sizeof(real64));
  std::memset(inverseF, 0, localMatSz*inumQuadraturePoints*NoElem*sizeof(real64));
#endif


  //
  //Allocate and set physical parameters
  //
  Index_type noSymEnt = 6;
  real64 bulkModulus = 10, shearModulus=20;
  geosxData detJ            = memoryManager::allocate<real64>(inumQuadraturePoints*NoElem, dataAllocated);
  geosxData meanStress      = memoryManager::allocate<real64>(inumQuadraturePoints*NoElem, dataAllocated);
  geosxData devStressData   = memoryManager::allocate<real64>(noSymEnt*inumQuadraturePoints*NoElem, dataAllocated);
  
  std::memset(detJ, 1 ,inumQuadraturePoints * NoElem);
  std::memset(meanStress, 1 ,inumQuadraturePoints * NoElem);
  std::memset(devStressData, 1 , noSymEnt * inumQuadraturePoints * NoElem);
  
  //
  //Set up function pointer to constitutive relationship
  //
#if defined(USE_CUDA)
  constUpdate myUpdate;
  cudaMemcpyFromSymbol(&myUpdate,deviceUpdate,sizeof(deviceUpdate));
#else
  constUpdate myUpdate = UpdateStatePoint;
#endif  

  //-----set up timer-----
  double start, end, diff;
  double myMin = 100000000;

  std::cout<<"Allocated "<<dataAllocated*1e-9<<" GBs of data"<<std::endl;  
  std::cout<<"Launching kernel. . . "<<std::endl;  
  for(Index_type it=0; it<Niter; ++it)
    {
      start = omp_get_wtime();


#if !defined(COMPUTE_SHAPE_FUN) && !defined(THREE_KERNEL_UPDATE)
      //Monolithic Kernel: Assumes precomputed shape function derivatives
      SolidMechanicsLagrangianFEMKernels::ObjectOfArraysKernel<kernelPol>(NoElem, elementList, dt, elemsToNodes,
                                                                          u_x, u_y, u_z, uhat_x, uhat_y, uhat_z,
                                                                          dNdX_x, dNdX_y, dNdX_z, constitutiveMap, devStressData,
                                                                          meanStress,shearModulus, bulkModulus, detJ, acc_x, acc_y, acc_z, myUpdate, nx, nx, nx);      
#elif defined(COMPUTE_SHAPE_FUN) && !defined(THREE_KERNEL_UPDATE)

      //Monolithic Kernel which computes shape function derivatives 
      SolidMechanicsLagrangianFEMKernels::ObjectOfArraysKernel_Shape<kernelPol>(NoElem,elementList, dt, elemsToNodes,
                                                                                u_x, u_y, u_z, uhat_x, uhat_y, uhat_z, VX, P,
                                                                                constitutiveMap, devStressData, meanStress,shearModulus,
                                                                                bulkModulus, detJ, acc_x, acc_y, acc_z, myUpdate, nx, nx, nx);

#elif defined(THREE_KERNEL_UPDATE)
      //Kinematic step
      SolidMechanicsLagrangianFEMKernels::ObjectOfArrays_KinematicKernel<kernelPol>(NoElem,elementList, dt, elemsToNodes,
                                                                                    u_x, u_y, u_z, uhat_x, uhat_y, uhat_z,
                                                                                    dNdX_x, dNdX_y, dNdX_z, constitutiveMap, devStressData, meanStress,shearModulus,
                                                                                    bulkModulus, detJ, acc_x, acc_y, acc_z, Dadt, Rot, detF, inverseF);
      //Constitutive update
      SolidMechanicsLagrangianFEMKernels::ConstitutiveUpdateKernel<kernelPol>(NoElem, elementList, Dadt, Rot, constitutiveMap, devStressData,
                                                                              meanStress, shearModulus, bulkModulus);

      //Integration step
      SolidMechanicsLagrangianFEMKernels::ObjectOfArrays_IntegrationKernel<kernelPol>(NoElem,elementList, dt, elemsToNodes,
                                                                                      u_x, u_y, u_z, uhat_x, uhat_y, uhat_z, 
                                                                                      dNdX_x, dNdX_y, dNdX_z, constitutiveMap, devStressData, meanStress,shearModulus,
                                                                                      bulkModulus, detJ, acc_x, acc_y, acc_z,
                                                                                      Dadt, Rot, detF, inverseF);
#endif

#if defined (USE_CUDA)
      cudaDeviceSynchronize();
#endif      
      end = omp_get_wtime();
      diff = end-start; 
      if( diff < myMin) myMin = diff;
    }

  
#if defined(USE_CUDA)
    std::cout<<"Computing on the GPU"<<std::endl;
#else
    std::cout<<"Computing on the CPU"<<std::endl;
#endif 

  std::cout<<"Run time = "<<myMin<<" sec"<<std::endl;
  std::cout<<"No of Elements:= "<<NoElem<<std::endl;

  //  
    //Free the data
    //
    memoryManager::deallocate(elementList);
    memoryManager::deallocate(constitutiveMap);
    memoryManager::deallocate(elemsToNodes);
    memoryManager::deallocate(VX);

#if !defined(COMPUTE_SHAPE_FUN)    
    memoryManager::deallocate(dNdX_x);
    memoryManager::deallocate(dNdX_y);
    memoryManager::deallocate(dNdX_z);
#endif    

    memoryManager::deallocate(u_x);
    memoryManager::deallocate(u_y);
    memoryManager::deallocate(u_z);
    memoryManager::deallocate(uhat_x);
    memoryManager::deallocate(uhat_y);
    memoryManager::deallocate(uhat_z);
    memoryManager::deallocate(acc_x);
    memoryManager::deallocate(acc_y);
    memoryManager::deallocate(acc_z);
    
    memoryManager::deallocate(detJ);
    memoryManager::deallocate(meanStress);
    memoryManager::deallocate(devStressData);
        
#if defined(THREE_KERNEL_UPDATE)
    memoryManager::deallocate(Dadt);
    memoryManager::deallocate(Rot);
    memoryManager::deallocate(detF);
    memoryManager::deallocate(inverseF);
#endif  
   
  return 0;
}
