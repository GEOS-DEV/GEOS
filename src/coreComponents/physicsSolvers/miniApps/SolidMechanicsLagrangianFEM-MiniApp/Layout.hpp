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

#ifndef __DATA_LAYOUT_HPP_
#define __DATA_LAYOUT_HPP_

#include "RAJA/RAJA.hpp"
#include "Array.hpp"
#include "ChaiVector.hpp"
#include <iostream>
#include "MemoryManager_impl.hpp"

//#include "common/DataTypes.hpp"
#include "miniDataTypes.hpp"

//compile time constants
#define local_dim 3
#define inumQuadraturePoints 8
#define inumNodesPerElement 8
#define localMatSz 9

//define constants
//CPU friendly
using localIndex = geosx::localIndex;
using globalIndex = geosx::globalIndex;
using real64 = geosx::real64;

//using localIndex = RAJA::Index_type;
//using globalIndex = RAJA::Index_type;
//using real64 = double;

using geosxData   = real64 * const RAJA_RESTRICT;
using geosxData_const   = real64 const * const RAJA_RESTRICT;
using geosxIndex = const localIndex * const RAJA_RESTRICT;

//#define STRUCTURED_GRID

//
//#define EXTERNAL_KERNELS

//#define OBJECT_OF_ARRAYS_LAYOUT
#define ARRAY_OF_OBJECTS_LAYOUT

//#define COMPUTE_SHAPE_FUN
//#define PRE_COMPUTE_P

//#define THREE_KERNEL_UPDATE

#define USE_CUDA

//#define SHAPE_FUN_FAST_INDEX_ELEM
//#define STRESS_FUN_FAST_INDEX_ELEM
//#define INVERT_REST_FAST_INDEX_ELEM
//#define THREE_KERNEL_DATA_FAST_INDEX_ELEM

//Select between raw ptr vs array
#define USE_GEOSX_ARRAY



void printParameters()
{  
#if defined(STRUCTURED_GRID)
  std::cout<<"Assuming a structured grid"<<std::endl;
#endif  
#if defined(COMPUTE_SHAPE_FUN)
  std::cout<<"Computing shape function derivatives at each time step"<<std::endl;
#endif  
#if defined(THREE_KERNEL_UPDATE)
  std::cout<<"Using a three kernel update"<<std::endl;
#endif  
#if defined(USE_CUDA)
  std::cout<<"computing on the GPU"<<std::endl;
#endif
#if defined(ARRAY_OF_OBJECTS_LAYOUT)
  std::cout<<"Data layout: Array of Objects"<<std::endl;  
#endif
#if defined(OBJECT_OF_ARRAYS_LAYOUT)
  std::cout<<"Data layout: Object of Arrays"<<std::endl;
#endif
#if defined(PRE_COMPUTE_P)
  std::cout<<"Pre-computing P (partial assembly of shape fun derivatives)"<<std::endl;
#endif  

#if defined(SHAPE_FUN_FAST_INDEX_ELEM)
  std::cout<<"SHAPE_FUN_FAST_INDEX_ELEM"<<std::endl;
#endif

#if defined(STRESS_FUN_FAST_INDEX_ELEM)
  std::cout<<"STRESS_FUN_FAST_INDEX_ELEM"<<std::endl;
#endif  

}


#if defined(USE_CUDA)
using kernelPol = RAJA::cuda_exec<256>;
using atomicPol = RAJA::atomic::cuda_atomic;
#elif defined(RAJA_ENABLE_OPENMP) 
using kernelPol = RAJA::omp_parallel_for_exec;
using atomicPol = RAJA::atomic::omp_atomic;
#else
using kernelPol = RAJA::loop_exec;
using atomicPol = RAJA::atomic::loop_atomic;
#endif


#if !defined(USE_GEOSX_ARRAY)

#define iu(k,i) iu[i + local_dim*id]
#define iuhat(k,i) iuhat[i + local_dim*id]



//constitutive update
#if defined(THREE_KERNEL_DATA_FASTEST_INDEX_ELEM)
#define Dadt_ptr(k, q, r, c) Dadt_ptr[k + noElem*(c + 3*(r + 3*q))]
#define Rot_ptr(k, q, r, c) Rot_ptr[k + noElem*(c + 3*(r + 3*q))]

#define detF_ptr(k, q) detF_ptr[k + noElem*q]
#define Finv_ptr(k, q, r, c) Finv_ptr[k + noElem*(c + 3*(r + 3*q))]
#else
#define Dadt_ptr(k, q, r, c) Dadt_ptr[c + 3*(r + 3*(q + inumQuadraturePoints*k))]
#define Rot_ptr(k, q, r, c) Rot_ptr[c + 3*(r + 3*(q + inumQuadraturePoints*k))]

#define detF_ptr(k, q) detF_ptr[q + inumQuadraturePoints*k]
#define Finv_ptr(k, q, r, c) Finv_ptr[c + 3*(r + 3*(q + inumQuadraturePoints*k))]
#endif


#if defined(SHAPE_FUN_FAST_INDEX_ELEM)
#define idNdX_x(k,q,a) idNdX_x[k + noElem*(a + inumNodesPerElement*q)]
#define idNdX_y(k,q,a) idNdX_y[k + noElem*(a + inumNodesPerElement*q)]
#define idNdX_z(k,q,a) idNdX_z[k + noElem*(a + inumNodesPerElement*q)]
#define idNdX(k,q,a,i) idNdX[k + noElem*(i + local_dim*(a + inumNodesPerElement*q))]
#else
#define idNdX_x(k,q,a) idNdX_x[a + inumNodesPerElement*(q + inumQuadraturePoints*k)]
#define idNdX_y(k,q,a) idNdX_y[a + inumNodesPerElement*(q + inumQuadraturePoints*k)]
#define idNdX_z(k,q,a) idNdX_z[a + inumNodesPerElement*(q + inumQuadraturePoints*k)]
#define idNdX(k,q,a,i) idNdX[i + local_dim*(a + inumNodesPerElement*(q + inumQuadraturePoints*k))]
#endif


#if defined(STRESS_FUN_FAST_INDEX_ELEM)
#define idevStressData(k,q,i) idevStressData[k + noElem*(i + 6*q)]
#else
#define idevStressData(k,q,i) idevStressData[i + 6*(q + inumQuadraturePoints*k)]
#endif

#if defined(INVERT_REST_FAST_INDEX_ELEM)
#define imeanStress(k,q) imeanStress[k + q*noElem]
#define iconstitutiveMap(k,q) iconstitutiveMap[k + q*noElem]
#define idetJ(k,q) idetJ[k + q*noElem]
#else
#define imeanStress(k,q) imeanStress[q + inumQuadraturePoints*k]
#define iconstitutiveMap(k,q) iconstitutiveMap[q + inumQuadraturePoints*k]
#define idetJ(k,q) idetJ[q + inumQuadraturePoints*k]
#endif

#endif

#endif
