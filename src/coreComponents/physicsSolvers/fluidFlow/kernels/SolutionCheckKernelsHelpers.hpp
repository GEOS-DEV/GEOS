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
 * @file SolutionCheckKernelsHelpers.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SOLUTIONCHECKKERNELSHELPER_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SOLUTIONCHECKKERNELSHELPER_HPP

#include "common/DataTypes.hpp"

namespace geos
{

class IdReporterBuffer;

/**
 * @brief TODO
 * @tparam IdCountType TODO
 * @tparam IdType TODO
 */
class IdReporterCollector
{
  friend class IdReporterBuffer;
public:

  using IdCountType = int32_t;
  using IdType = globalIndex;

  // TODO : proper docs. can be copied & moved as this class only has views to the internal chai memory buffers
  IdReporterCollector( IdReporterCollector const & other ) = default;
  IdReporterCollector( IdReporterCollector && other ) = default;
  IdReporterCollector & operator=( IdReporterCollector const & other ) = default;
  IdReporterCollector & operator=( IdReporterCollector && other ) = default;

  static IdReporterCollector disabled() // -> cpp
  {
    return IdReporterCollector( arrayView1d< IdCountType >(),
                                arrayView1d< IdType >(),
                                arrayView1d< globalIndex >() );
  }

  /**
   * @brief TODO
   * @tparam CollectorAtomicPolicy The policy of the atomic increment on the ids counter.
   * @param m_idsCounter The ids counter to increment with an atomic operation.
   * @param m_idsBuffer The output id buffer, in the same memory space as m_idsCounter.
   *                  If its size is 0 (= disabled output) or not not large enought, the buffer is not filled.
   * @param id The Id to add to the buffer.
   */
  template< typename CollectorAtomicPolicy >
  GEOS_HOST_DEVICE
  void collectId( CollectorAtomicPolicy, IdType id ) const
  {
    if( !m_idsCounter.empty() )
    {
      IdCountType const outputStart = RAJA::atomicInc< CollectorAtomicPolicy >( &m_idsCounter[0] );

      if( outputStart < m_idsBuffer.size() )
      {
        m_idsBuffer[outputStart] = m_localToGlobalId[id];
      }
    }
  }

  // // currently unused version for adding multiple ids from a given kernel
  // template< typename AddedArray, typename IdCountType >
  // GEOS_HOST_DEVICE
  // void collectIds( AddedArray const & newIds, IdCountType newIdsCount )
  // {
  //   IdCountType const outputStart = RAJA::atomicAdd< CollectorAtomicPolicy >( &m_idsCounter[0], newIdsCount );
  //   IdCountType const maxNbIdToAdd = IdCountType( m_idsBuffer.size() - outputStart );
  //   newIdsCount = LvArray::math::min( newIdsCount, maxNbIdToAdd );
  //   for( IdCountType i = 0; i < newIdsCount; ++i )
  //   {
  //     m_idsBuffer[outputStart + i] = newIds[i];
  //   }
  // }

private:

  // array of one element to get benefit of chai managed memory.
  arrayView1d< IdCountType > m_idsCounter;

  // ids of detected elements, quantity limited to 'maxIdsCount'
  arrayView1d< IdType > m_idsBuffer;

  arrayView1d< globalIndex > m_localToGlobalId;

  IdReporterCollector( arrayView1d< IdCountType > const & idsCounter,
                       arrayView1d< IdType > const & idsArray,
                       arrayView1d< globalIndex > const & localToGlobalId ):
    m_idsCounter( idsCounter ),
    m_idsBuffer( idsArray ),
    m_localToGlobalId( localToGlobalId )
  {}

};

// /**
//  * @brief TODO
//  * @tparam AtomicPolicy The policy of the atomic increment on idsCounter.
//  * @tparam KernelStackArray TODO
//  * @tparam IdType TODO
//  * @tparam IdCountType TODO
//  * @param idsCounter The ids counter to increment with an atomic operation.
//  * @param idsBuffer The output id buffer, in the same memory space as idsCounter.
//  *                  If its size is 0 (= disabled output) or not not large enought, the buffer is not filled.
//  * @param id The Id to add to the buffer.
//  * @return TODO
//  */
// template< typename AtomicPolicy, typename OutputStackArray, typename IdType, typename IdCountType >
// GEOS_HOST_DEVICE
// void collectKernelId( IdCountType & idsCounter, OutputStackArray & idsBuffer, IdType id )
// {
//   static constexpr IdCountType addingRequest = 1;
//   IdCountType const idsBufferStart = RAJA::atomicAdd< AtomicPolicy >( &idsCounter, addingRequest );
//   if( idsBufferStart < idsBuffer.size() )
//   {
//     idsBuffer[idsBufferStart] = kernelIds[i];
//   }
// }

// // // currently unused version for adding multiple ids from a given kernel
// // template< typename ReducePolicy, typename OutputStackArray, typename KernelStackArray, typename IdCountType >
// // GEOS_HOST_DEVICE
// // void collectKernelIds( OutputStackArray & outputBuffer,
// //                        IdCountType & outputIdsCounter,
// //                        KernelStackArray const & kernelIds,
// //                        IdCountType const kernelIdsCount )
// // {
// //   IdCountType const outputBufferStart = RAJA::atomicAdd< ReducePolicy >( &outputIdsCounter, kernelIdsCount );
// //   IdCountType const maxNbIdToAdd = IdCountType( outputBuffer.capacity() - outputBufferStart );
// //   IdCountType const nbIdToAdd = LvArray::math::min( kernelIdsCount, maxNbIdToAdd );
// //   for( IdCountType i = 0; i < nbIdToAdd; ++i )
// //   {
// //     outputBuffer[outputBufferStart + i] = kernelIds[i];
// //   }
// // }

// template< typename OutputDynamicArray, typename InputArray, typename IdCountType >
// void aggregateIdsBuffers( OutputDynamicArray & outputBuffer,
//                           IdCountType & outputIdsCounter,
//                           InputArray const & ids,
//                           IdCountType const idsCount )
// {
//   outputIdsCounter += idsCount;
//   IdCountType const numIdsToAdd = std::min( idsCount, IdCountType( ids.capacity() ) );
//   for( int i = 0; i < numIdsToAdd; ++i )
//   {
//     outputBuffer.emplace_back( ids[i] );   // todo local -> global
//   }
// }

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SOLUTIONCHECKKERNELSHELPER_HPP
