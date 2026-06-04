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
 * @file ElementReporterKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_ELEMENTREPORTERKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_ELEMENTREPORTERKERNELS_HPP

#include "common/DataTypes.hpp"

namespace geos
{

class ElementReporterBuffer;

struct ElementReport
{
  /// the global id of the reported element
  globalIndex m_id;
  /// a value to report for the given element (i.e. a negative pressure, a density... Or if needed
  /// could be of a templated type for composite values)
  real64 m_value;
};

/**
 * @brief Collects and reports elements ids and data using an atomic counter.
 *        This class provides functionality to collect data from multiple threads safely by incrementing
 *        through an atomic counter for each reported element's ID. The collected IDs are stored in a
 *        size limited buffer, which can be used later for reporting or analysis purposes.
 */
class ElementsReporterCollector
{
  friend class ElementReporterBuffer;
public:

  using ElementCount = int32_t;

  /**
   * @name Constructors
   * @brief This object can be copied and moved as it only provides views to internal memory buffers.
   */
  ///@{
  /** @cond DO_NOT_DOCUMENT */

  ElementsReporterCollector( ElementsReporterCollector const & other ) = default;

  ElementsReporterCollector( ElementsReporterCollector && other ) = default;

  ElementsReporterCollector & operator=( ElementsReporterCollector const & other ) = default;

  ElementsReporterCollector & operator=( ElementsReporterCollector && other ) = default;

  /** @endcond */
  ///@}

  static ElementsReporterCollector disabled()
  {
    return ElementsReporterCollector( arrayView1d< ElementCount >(),
                                      arrayView1d< ElementReport >(),
                                      arrayView1d< globalIndex const >() );
  }

  /**
   * @brief Collects a single element report and adds its ID to the output buffer if not disabled and
   *        there are available slots in the buffer.
   * @tparam CollectorAtomicPolicy The atomic increment operation to use for thread-safe counter increments.
   * @param report A constant reference to an `ElementReport` object containing data from a single element
   */
  template< typename CollectorAtomicPolicy >
  GEOS_HOST_DEVICE
  void collectElement( CollectorAtomicPolicy, ElementReport const & report ) const
  {
    if( !m_elementsCounter.empty() )
    {
      ElementCount const outputStart = RAJA::atomicInc< CollectorAtomicPolicy >( &m_elementsCounter[0] );

      if( outputStart < m_elementsBuffer.size() )
      {
        m_elementsBuffer[outputStart].m_id = m_localToGlobalId[report.m_id];
        m_elementsBuffer[outputStart].m_value = report.m_value;
      }
    }
  }

  // // currently unused version for adding multiple ids from a given kernel
  // template< typename AddedArray, typename ElementCount >
  // GEOS_HOST_DEVICE
  // void collectIds( AddedArray const & newIds, ElementCount newIdsCount )
  // {
  //   ElementCount const outputStart = RAJA::atomicAdd< CollectorAtomicPolicy >( &m_elementsCounter[0], newIdsCount );
  //   ElementCount const maxNbIdToAdd = ElementCount( m_elementsBuffer.size() - outputStart );
  //   newIdsCount = LvArray::math::min( newIdsCount, maxNbIdToAdd );
  //   for( ElementCount i = 0; i < newIdsCount; ++i )
  //   {
  //     m_elementsBuffer[outputStart + i] = newIds[i];
  //   }
  // }

private:

  /// array of one element to get benefit of chai managed memory.
  arrayView1d< ElementCount > m_elementsCounter;

  /// ids of detected elements, quantity limited to 'maxIdsCount'
  arrayView1d< ElementReport > m_elementsBuffer;

  /// Maps local element IDs to their respective global indices.
  arrayView1d< globalIndex const > m_localToGlobalId;

  ElementsReporterCollector( arrayView1d< ElementCount > const & elementsCounter,
                             arrayView1d< ElementReport > const & elementsBuffer,
                             arrayView1d< globalIndex const > const & localToGlobalId ):
    m_elementsCounter( elementsCounter ),
    m_elementsBuffer( elementsBuffer ),
    m_localToGlobalId( localToGlobalId )
  {}

};

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_ELEMENTREPORTERKERNELS_HPP
