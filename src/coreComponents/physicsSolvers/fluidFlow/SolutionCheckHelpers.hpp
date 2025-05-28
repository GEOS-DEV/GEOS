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

class ElementsReporterOutput
{
public:

  using ElementCount = ElementsReporterCollector::ElementCount;

  ElementsReporterOutput( ElementsReporterBuffer const & buffer );

  ElementCount getRanksSignaledIdsCount() const
  { return m_ranksSignaledElementsCount; }

  ElementCount getRanksCollectedIdsCount() const
  { return m_ranksCollectedElementsCount; }

  void outputTooLowValues( string_view linesPrefix,
                           string_view valueNaming,
                           real64 minValue,
                           units::Unit valueUnit ) const;

private:

  ElementsReporterBuffer const & m_buffer;

  ElementCount m_ranksSignaledElementsCount;

  ElementCount m_ranksCollectedElementsCount;

};

class ElementsReporterBuffer
{
public:

  using ElementCount = ElementsReporterCollector::ElementCount;

  /**
   * @brief Construct a preallocated buffer to collect a limited quantity of ids in kernels.
   * @param maxCollectionSize Limit of the buffer.
   *                          If 0, the buffering functionnality is disabled and only the counting is enabled.
   */
  ElementsReporterBuffer( bool enabled, ElementCount maxCollectionSize );

  // TODO: Proper docs. can be moved without any issue.
  ElementsReporterBuffer( ElementsReporterBuffer && other ) = default;
  ElementsReporterBuffer & operator=( ElementsReporterBuffer && other ) = default;

  // TODO: Proper docs. copying prevented has it doesn't seem useful / relevant.
  ElementsReporterBuffer( ElementsReporterBuffer const & other ) = delete;
  ElementsReporterBuffer & operator=( ElementsReporterBuffer const & other ) = delete;

  ElementCount getSignaledElementsCount() const
  { return m_elementsCounter.empty() ? 0 : m_elementsCounter[0]; }

  ElementCount getCollectedElementsCount() const
  { return LvArray::math::min( getSignaledElementsCount(), m_elementsBuffer.size() ); }

  ElementReport operator[]( ElementCount id ) const
  { return m_elementsBuffer[id]; }

  auto begin() const
  { return m_elementsBuffer.begin(); }

  auto end() const
  { return m_elementsBuffer.begin() + getCollectedElementsCount(); }

  bool enabled() const
  { return !m_elementsCounter.empty(); }

  bool empty() const
  { return getCollectedElementsCount() == 0; }

  bool isComplete() const
  { return getCollectedElementsCount() < getSignaledElementsCount(); }

  /**
   * @return A view on the ids array owned by the instance. -> change comment to explain the interest for kernels
   */
  ElementsReporterCollector createCollector( arrayView1d< globalIndex const > const & localToGlobalId ) const;

  ElementsReporterOutput createOutput() const;

private:

  // array of one element to get benefit of managed host-device memory.
  array1d< ElementCount > m_elementsCounter;

  // ids of detected elements
  array1d< ElementReport > m_elementsBuffer;

};

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SOLUTIONCHECKHELPER_HPP
