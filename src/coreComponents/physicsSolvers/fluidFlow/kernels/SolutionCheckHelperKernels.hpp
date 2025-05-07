/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SolutionCheckKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SOLUTIONCHECKHELPERKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SOLUTIONCHECKHELPERKERNELS_HPP

#include "physicsSolvers/fluidFlow/kernels/compositional/SolutionScalingAndCheckingKernelBase.hpp"

namespace geos
{

namespace flowSolverBaseKernels
{

//fonction inutile
template< typename KernelStackArray, typename IdType, typename IdCountType >
GEOS_HOST_DEVICE
void collectKernelId( KernelStackArray & idList,
                      IdCountType & idsCounter,
                      IdType id )
{
  if( idsCounter < idList.capacity() )
  {
    idList[idsCounter] = id; // ne peut avoir qu'une occurence
  }

  idsCounter += 1;
}

/**
 * @brief
 * @param stackBuffer
 * @param kernelIds
 * @param kernelIdsCount Strictly positive.
 */
template< typename OutputStackArray, typename KernelStackArray, typename IdCountType, typename ReducePolicy >
GEOS_HOST_DEVICE
void aggregateKernelsIds( OutputStackArray & outputBuffer,
                          T & outputIdsCounter,
                          KernelStackArray const & kernelIds,
                          IdCountType const kernelIdsCount )
{
  IdCountType const outputBufferStart = RAJA::atomicAdd< RAJA::auto_atomic >( &outputIdsCounter,
                                                                              kernelIdsCount );
  IdCountType const maxNbIdToAdd = IdCountType( outputBuffer.capacity() - outputBufferStart );
  IdCountType const nbIdToAdd = LvArray::math::min( kernelIdsCount, maxNbIdToAdd );
  for( IdCountType i = 0; i < nbIdToAdd; ++i )
  {
    outputBuffer[outputBufferStart + i] = kernelIds[i];
  }
}

template< typename OutputDynamicArray, typename InputArray, typename IdCountType >
void aggregateKernelsIds( OutputDynamicArray & outputBuffer,
                          IdCountType & outputIdsCounter,
                          InputArray const & ids,
                          IdCountType const idsCount )
{
  outputIdsCounter += idsCount;
  IdCountType const numIdsToAdd = std::min( idsCount, IdCountType( ids.capacity() ) );
  for( int i = 0; i < numIdsToAdd; ++i )
  {
    outputBuffer.emplace_back( ids[i] );   // todo local -> global
  }
}

} // namespace flowSolverBaseKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_SOLUTIONCHECKKERNEL_HPP
