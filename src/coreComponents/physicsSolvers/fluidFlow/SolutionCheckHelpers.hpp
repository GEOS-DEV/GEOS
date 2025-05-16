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
 * @file SolutionCheckHelpers.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SOLUTIONCHECKHELPERS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SOLUTIONCHECKHELPERS_HPP

#include "physicsSolvers/fluidFlow/kernels/SolutionCheckKernelsHelpers.hpp"
#include "common/Units.hpp"

namespace geos
{

class IdReporterOutput
{
public:

  IdReporterOutput( IdReporterBuffer const & buffer );

  void outputWrongValues( string_view linesPrefix,
                          string_view valueNaming,
                          real64 minValue,
                          units::Unit valueUnit ) const;

private:

  IdReporterBuffer const & m_buffer;

};

class IdReporterBuffer
{
public:

  using IdCountType = int32_t;
  using IdType = globalIndex;

  /**
   * @brief Construct a preallocated buffer to collect a limited quantity of ids in kernels.
   * @param maxCollectionSize Limit of the buffer.
   *                          If 0, the buffering functionnality is disabled and only the counting is enabled.
   */
  IdReporterBuffer( bool enabled, IdCountType maxCollectionSize );

  // TODO: Proper docs. can be moved without any issue.
  IdReporterBuffer( IdReporterBuffer && other ) = default;
  IdReporterBuffer & operator=( IdReporterBuffer && other ) = default;

  // TODO: Proper docs. copying prevented has it doesn't seem useful / relevant.
  IdReporterBuffer( IdReporterBuffer const & other ) = delete;
  IdReporterBuffer & operator=( IdReporterBuffer const & other ) = delete;

  IdCountType getSignaledIdsCount() const
  { return m_idsCounter.empty() ? 0 : m_idsCounter[0]; }

  IdCountType getCollectedIdsCount() const
  { return LvArray::math::min( getSignaledIdsCount(), m_idsBuffer.size() ); }

  auto begin() const
  { return m_idsBuffer.begin(); }

  auto end() const
  { return m_idsBuffer.begin() + getCollectedIdsCount(); }

  bool enabled() const
  { return !m_idsCounter.empty(); }

  bool empty() const
  { return getCollectedIdsCount() == 0; }

  bool isComplete() const
  { return getCollectedIdsCount() < getSignaledIdsCount(); }

  /**
   * @return A view on the ids array owned by the instance. -> change comment to explain the interest for kernels
   */
  IdReporterCollector createCollector( arrayView1d< globalIndex > const & localToGlobalId ) const;

  IdReporterOutput createOutput() const;

private:

  // array of one element to get benefit of managed host-device memory.
  array1d< IdCountType > m_idsCounter;

  // ids of detected elements
  array1d< IdType > m_idsBuffer;

};

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SOLUTIONCHECKHELPER_HPP
