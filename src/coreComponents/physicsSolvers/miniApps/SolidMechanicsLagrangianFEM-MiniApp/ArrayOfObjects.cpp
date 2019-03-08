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

//Stored the kernels
#include "../../src/SolidMechanicsLagrangianFEMKernels_brief_impl.hpp"

//
//Driver for the GEOSX proxy app
//Assumes an Array of Objects data layout
//

int main(int argc, char* const argv[])
{

#if defined(USE_GEOSX_ARRAY)
  std::cout<<"Using GEOSX ARRAY"<<std::endl;
#elif defined(USE_RAJA_VIEW)
  std::cout<<"Using RAJA VIEW"<<std::endl;
#else
  std::cout<<"Using RAW PTRS"<<std::endl;
#endif

  if(argc != 3)
    {
      std::cout<<"usage ./main NoElem Iter"<<std::endl;
      exit(-1);
    }

  size_t dataAllocated   = 0.0; 
  localIndex Kx      = atoi(argv[1]); 
  localIndex Niter   =  atoi(argv[2]);

  std::cout<<"GEOSX mini-app: Array of Objects data structures"<<std::endl;

  real64 dt = 0.125; //step size
  const localIndex nx = Kx+1; //Number of nodes in a cartesian dimension
  const localIndex NoElem = Kx*Kx*Kx; //Total number of elements
  const localIndex numNodes = nx*nx*nx; //Total number of nodes    

  //
  //Generate an element list
  localIndex * const elementList = memoryManager::allocate<localIndex>(NoElem, dataAllocated);

  for(localIndex i=0; i<NoElem; ++i) elementList[i] = i;

  //Allocate space for constitutive map
  //
#if defined(USE_GEOSX_ARRAY)
  LvArray::Array<localIndex, 2, localIndex> _constitutiveMap(NoElem, NUMQUADPTS);
  LvArray::ArrayView<localIndex,2,localIndex> & constitutiveMap = _constitutiveMap;

#elif defined(USE_RAJA_VIEW)

  localIndex *  _constitutiveMap = memoryManager::allocate<localIndex>(NUMQUADPTS*NoElem, dataAllocated); 
  RAJA::View<localIndex, RAJA::Layout<2,localIndex,1>> constitutiveMap(_constitutiveMap, NoElem, NUMQUADPTS);
#else

  localIndex *  constitutiveMap = memoryManager::allocate<localIndex>(NUMQUADPTS*NoElem, dataAllocated); 
#endif

  //
  //Generate space for an element to node list
  //
#if defined(USE_GEOSX_ARRAY)
  LvArray::Array<localIndex, 2, localIndex> _elemsToNodes(NoElem, NODESPERELEM);
  LvArray::ArrayView<localIndex,2,localIndex> & elemsToNodes = _elemsToNodes;

#elif defined(USE_RAJA_VIEW)
  localIndex  *  elemsToNodes = memoryManager::allocate<localIndex>(NODESPERELEM*NoElem, dataAllocated);
  //RAJA::View<localIndex, RAJA::Layout<2, localIndex, 1>> elemsToNodes(_elemsToNodes, NoElem, NODESPERELEM);
#else
  localIndex  *  elemsToNodes = memoryManager::allocate<localIndex>(NODESPERELEM*NoElem, dataAllocated);
#endif


  //
  //Allocate space for a list of vertices, generate a mesh, and populate the constitutive map
  //
#if defined(USE_GEOSX_ARRAY)
  LvArray::Array<real64, 2, localIndex> _VX(numNodes, LOCAL_DIM);
  LvArray::ArrayView<real64,2, localIndex> & VX = _VX;
#elif defined(USE_RAJA_VIEW)

  geosxData _VX = memoryManager::allocate<real64>(numNodes*LOCAL_DIM, dataAllocated);
  RAJA::View<real64, RAJA::Layout<2, localIndex, 1> > VX(_VX, numNodes, LOCAL_DIM);
#else
  geosxData VX = memoryManager::allocate<real64>(numNodes*LOCAL_DIM, dataAllocated);
#endif


  meshGen(VX, elemsToNodes,constitutiveMap,Kx);

  //
  //Precompute evaluate the quadrature points at every basis function
  //with respect to the parent basis function
  //

  P_Wrapper P; 
  generateP(P, NODESPERELEM, NUMQUADPTS);

  //
  //Allocate space for shape function derivatives and compute their values
  //
#if defined(USE_GEOSX_ARRAY)
  LvArray::Array<real64, 4, localIndex> _dNdX(NoElem, NUMQUADPTS, NODESPERELEM, LOCAL_DIM);
  LvArray::ArrayView<real64,4, localIndex> & dNdX = _dNdX;
#elif defined(USE_RAJA_VIEW)
  geosxData _dNdX = memoryManager::allocate<real64>(NODESPERELEM*NUMQUADPTS*NoElem*LOCAL_DIM, dataAllocated);
  RAJA::View<real64, RAJA::Layout<4, localIndex, 3> > dNdX(_dNdX, NoElem, NUMQUADPTS, NODESPERELEM, LOCAL_DIM);
#else
  geosxData dNdX = memoryManager::allocate<real64>(NODESPERELEM*NUMQUADPTS*NoElem*LOCAL_DIM, dataAllocated);
#endif
  
  make_dNdX(dNdX, VX, elemsToNodes, NoElem, NUMQUADPTS, NODESPERELEM);

  ///
  //Allocate space for nodal degrees of freedom as an Array of Objects
  ///
#if defined(USE_GEOSX_ARRAY)
  LvArray::Array<real64, 2, localIndex> _u(numNodes, LOCAL_DIM);
  LvArray::Array<real64, 2, localIndex> _uhat(numNodes, LOCAL_DIM);
  LvArray::Array<real64, 2, localIndex> _acc(numNodes, LOCAL_DIM);

  LvArray::ArrayView<real64, 2, localIndex> & u = _u;
  LvArray::ArrayView<real64, 2, localIndex> & uhat = _uhat;
  LvArray::ArrayView<real64, 2, localIndex> & acc = _acc;
  std::memset(u.data(),0, numNodes*LOCAL_DIM*sizeof(real64));
  std::memset(uhat.data(),0, numNodes*LOCAL_DIM*sizeof(real64));
  std::memset(acc.data(),0, numNodes*LOCAL_DIM*sizeof(real64));

#elif defined(USE_RAJA_VIEW)

  geosxData _u = memoryManager::allocate<real64>(numNodes*LOCAL_DIM, dataAllocated);
  geosxData _uhat = memoryManager::allocate<real64>(numNodes*LOCAL_DIM, dataAllocated);
  geosxData _acc = memoryManager::allocate<real64>(numNodes*LOCAL_DIM, dataAllocated);
  
  RAJA::View<real64, RAJA::Layout<2, localIndex, 1>> u(_u, numNodes, LOCAL_DIM);
  RAJA::View<real64, RAJA::Layout<2, localIndex, 1>> uhat(_uhat, numNodes, LOCAL_DIM);
  RAJA::View<real64, RAJA::Layout<2, localIndex, 1>> acc(_acc, numNodes, LOCAL_DIM);  

  std::memset(_u,0, numNodes*LOCAL_DIM*sizeof(real64));
  std::memset(_uhat,0, numNodes*LOCAL_DIM*sizeof(real64));
  std::memset(_acc,0, numNodes*LOCAL_DIM*sizeof(real64));
 
#else
  geosxData u = memoryManager::allocate<real64>(numNodes*LOCAL_DIM, dataAllocated);
  geosxData uhat = memoryManager::allocate<real64>(numNodes*LOCAL_DIM, dataAllocated);
  geosxData acc = memoryManager::allocate<real64>(numNodes*LOCAL_DIM, dataAllocated);

  std::memset(u,0, numNodes*LOCAL_DIM*sizeof(real64));
  std::memset(uhat,0, numNodes*LOCAL_DIM*sizeof(real64));
  std::memset(acc,0, numNodes*LOCAL_DIM*sizeof(real64));
#endif

  //
  //Allocate and set physical parameters
  //
  localIndex noSymEnt = 6;
  real64 bulkModulus = 10, shearModulus=20;

#if defined(USE_GEOSX_ARRAY)
  
  LvArray::Array<real64,2,localIndex> _detJ(NoElem, NUMQUADPTS);
  LvArray::Array<real64,3,localIndex> _devStressData(NoElem, NUMQUADPTS, noSymEnt);
  LvArray::Array<real64,1,localIndex> _meanStress(NoElem*NUMQUADPTS); //Reformulate to 1D

  LvArray::ArrayView<real64,2,localIndex> & detJ          = _detJ;
  LvArray::ArrayView<real64,3,localIndex> & devStressData = _devStressData;  
  LvArray::ArrayView<real64,1,localIndex> & meanStress    = _meanStress;

  std::memset(detJ.data(), 1 ,NUMQUADPTS * NoElem * sizeof(real64));
  std::memset(meanStress.data(), 1 ,NUMQUADPTS * NoElem * sizeof(real64));
  std::memset(devStressData.data(), 1 , noSymEnt * NUMQUADPTS * NoElem * sizeof(real64));

#elif defined(USE_RAJA_VIEW)
  
  geosxData _detJ            = memoryManager::allocate<real64>(NUMQUADPTS*NoElem, dataAllocated);
  geosxData _meanStress      = memoryManager::allocate<real64>(NUMQUADPTS*NoElem, dataAllocated);
  geosxData _devStressData   = memoryManager::allocate<real64>(noSymEnt*NUMQUADPTS*NoElem, dataAllocated);

  RAJA::View<real64, RAJA::Layout<2,localIndex,1>> detJ(_detJ, NoElem, NUMQUADPTS);
  RAJA::View<real64, RAJA::Layout<3,localIndex,2>> devStressData(_devStressData, NoElem, NUMQUADPTS, noSymEnt);
  RAJA::View<real64, RAJA::Layout<1,localIndex,0>> meanStress(_meanStress, NoElem*NUMQUADPTS);

  std::memset(_detJ, 1 ,NUMQUADPTS * NoElem * sizeof(real64));
  std::memset(_meanStress, 1 ,NUMQUADPTS * NoElem * sizeof(real64));
  std::memset(_devStressData, 1 , noSymEnt * NUMQUADPTS * NoElem * sizeof(real64));  

#else
  geosxData detJ            = memoryManager::allocate<real64>(NUMQUADPTS*NoElem, dataAllocated);
  geosxData meanStress      = memoryManager::allocate<real64>(NUMQUADPTS*NoElem, dataAllocated);
  geosxData devStressData   = memoryManager::allocate<real64>(noSymEnt*NUMQUADPTS*NoElem, dataAllocated);

  std::memset(detJ, 1 ,NUMQUADPTS * NoElem * sizeof(real64));
  std::memset(meanStress, 1 ,NUMQUADPTS * NoElem * sizeof(real64));
  std::memset(devStressData, 1 , noSymEnt * NUMQUADPTS * NoElem * sizeof(real64));
#endif



  //
  //Set up function pointer for constitutive relationship
  //
#if defined(USE_GPU)
  std::cout<<"Using CUDA"<<std::endl;
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
  //-----------------[Kernel Launch]---------------------------------    
  for(localIndex it=0; it<Niter; ++it)
    {
      start = omp_get_wtime();


      //Standard Monolithic Kernel
      SolidMechanicsLagrangianFEMKernels::ArrayOfObjectsKernel<kernelPol>(NoElem,elementList, dt, elemsToNodes,
                                                                          u, uhat, dNdX, constitutiveMap, devStressData, meanStress,shearModulus,
                                                                          bulkModulus, detJ, acc, myUpdate, nx, nx, nx);

#if defined (USE_GPU)
      cudaDeviceSynchronize();
#endif      
      end = omp_get_wtime();
      diff = end-start; 
      if( diff < myMin) myMin = diff;
    }

#if defined(USE_GPU)
    std::cout<<"Computing on the GPU"<<std::endl;
#else
    std::cout<<"Computing on the CPU"<<std::endl;
#endif
  
  std::cout<<"Run time = "<<myMin<<" sec"<<std::endl;
  std::cout<<"No of Elements:= "<<NoElem<<std::endl;


  //
  //Free data
  //
  /*
#if !defined(USE_GEOSX_ARRAY)
  memoryManager::deallocate(elementList);
  memoryManager::deallocate(constitutiveMap);
  memoryManager::deallocate(elemsToNodes);
  memoryManager::deallocate(VX);

#if !defined(COMPUTE_SHAPE_FUN)  
  memoryManager::deallocate(dNdX);
#endif

  memoryManager::deallocate(u);
  memoryManager::deallocate(uhat);
  memoryManager::deallocate(acc);
  
  memoryManager::deallocate(detJ);
  memoryManager::deallocate(meanStress);
  memoryManager::deallocate(devStressData);
  
#if defined(THREE_KERNEL_UPDATE)
  memoryManager::deallocate(Dadt);
  memoryManager::deallocate(Rot);
  memoryManager::deallocate(detF);
  memoryManager::deallocate(inverseF);
#endif  
  */


  return 0;
}
