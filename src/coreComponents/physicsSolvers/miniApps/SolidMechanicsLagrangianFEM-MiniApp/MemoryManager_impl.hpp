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

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Copyright (c) 2016-18, Lawrence Livermore National Security, LLC.
//
// Produced at the Lawrence Livermore National Laboratory
//
// LLNL-CODE-689114
//
// All rights reserved.
//
// This file is part of RAJA.
//
// For details about use and distribution, please read RAJA/LICENSE.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#ifndef EXAMPLES_MEMORYMANAGER_HPP
#define EXAMPLES_MEMORYMANAGER_HPP

#include "RAJA/RAJA.hpp"
#include "Layout.hpp"
#if defined(RAJA_ENABLE_CHAI)
#include "chai/ManagedArray.hpp"
#include "chai/util/forall.hpp"
#endif

/*
  As RAJA does not manage memory we include a general purpose memory
  manager which may be used to perform c++ style allocation/deallocation 
  or allocate/deallocate CUDA unified memory. The type of memory allocated 
  is dependent on how RAJA was configured.  
*/
namespace memoryManager{

  template <typename T>
  T *allocate(RAJA::Index_type size, size_t &dataAlloc)
  {
    T *ptr;
    dataAlloc += size * sizeof(T);
#if defined(RAJA_ENABLE_CUDA)
    cudaErrchk(cudaMallocManaged((void **)&ptr, sizeof(T) * size, cudaMemAttachGlobal));               
#else
    ptr = new T[size];
#endif
    return ptr;
  }
  
  template <typename T>
  void deallocate(T ptr)
  {
    if (ptr) {
#if defined(RAJA_ENABLE_CUDA)
      cudaFree(ptr);
#else
      delete[] ptr;
#endif
      ptr = nullptr;
    }    
  }
  
}
#endif
