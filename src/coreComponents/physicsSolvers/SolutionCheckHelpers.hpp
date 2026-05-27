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

#ifndef GEOS_PHYSICSSOLVERS_SOLUTIONCHECKHELPERS_HPP
#define GEOS_PHYSICSSOLVERS_SOLUTIONCHECKHELPERS_HPP

#include "physicsSolvers/SolutionCheckKernelsHelpers.hpp"
#include "common/Units.hpp"

namespace geos
{

/**
 * @brief A class to report elements collected by the solver.
 */
class ElementsReporterOutput
{
public:

  /// Type alias for elements count (e.g., localIndex, globalIndex).
  using ElementCount = ElementsReporterCollector::ElementCount;

  /**
   * @brief Construct a preallocated buffer for collecting element ids in kernels.
   * @param buffer The buffer that will be utilized for counting & collecting elements IDs during kernel execution.
   */
  ElementsReporterOutput( ElementsReporterBuffer const & buffer );

  /**
   * @return The number of ranks that have signaled an id.
   */
  ElementCount getRanksSignaledIdsCount() const
  { return m_ranksSignaledElementsCount; }

  /**
   * @return The total count of collected elements across all ranks for signaling ids.
   */
  ElementCount getRanksCollectedIdsCount() const
  { return m_ranksCollectedElementsCount; }

  /**
   * @brief Report elements with values below a specified threshold in the log:
   *        Outputs lines indicating which variables have collected element ids whose corresponding
   *        solution components are too low, potentially signaling underflow or numerical instability.
   * @param linesPrefix Prefix for the line of text to be printed
   * @param valueNaming The name used when referring to variables within this context (e.g., "pressure", "density").
   * @param minValue Minimum acceptable solution component values. Values below this threshold are reported.
   * @param valueUnit Unit in which `minValue` is expressed.
   */
  void outputTooLowValues( string_view linesPrefix,
                           string_view valueNaming,
                           real64 minValue,
                           units::Unit valueUnit ) const;

private:

  /// Preallocated buffer for collecting ids.
  ElementsReporterBuffer const & m_buffer;

  /// Count of signaled elements per rank.
  ElementCount m_ranksSignaledElementsCount;

  /// Total collected signaling id count across ranks.
  ElementCount m_ranksCollectedElementsCount;

};

/**
 * @brief A buffer to count and store element ids during kernel execution.
 *        This facilitates the reporting mechanism by allowing a preallocated space for storing & counting elements.
 */
class ElementsReporterBuffer
{
public:

  /// Type alias for elements count (e.g., localIndex, globalIndex).
  using ElementCount = ElementsReporterCollector::ElementCount;

  /**
   * @brief Construct a preallocated buffer to collect a limited quantity of ids in kernels.
   * @param maxCollectionSize Limit of the buffer.
   *                          If 0, the buffering functionnality is disabled and only the counting is enabled.
   */
  ElementsReporterBuffer( bool enabled, ElementCount maxCollectionSize );

  /**
   * @brief Transfers ownership of an ElementsReporterBuffer to another instance (move semantics).
   */
  ElementsReporterBuffer( ElementsReporterBuffer && other ) = default;

  /**
   * @brief Transfers ownership of an ElementsReporterBuffer to another instance (move semantics).
   */
  ElementsReporterBuffer & operator=( ElementsReporterBuffer && other ) = default;

  /**
   * @brief Copying prevented as it doesn't seem relevant / useful.
   */
  ElementsReporterBuffer( ElementsReporterBuffer const & ) = delete;

  /**
   * @brief Copying prevented as it doesn't seem relevant / useful.
   */
  ElementsReporterBuffer & operator=( ElementsReporterBuffer const & ) = delete;

  /**
   * @return the count of signaled elements.
   */
  ElementCount getSignaledElementsCount() const
  { return m_elementsCounter.empty() ? 0 : m_elementsCounter[0]; }

  /**
   * @return the collected elements that could effectivly be stored (zero if no collection is enabled).
   */
  ElementCount getCollectedElementsCount() const
  { return LvArray::math::min( getSignaledElementsCount(), m_elementsBuffer.size() ); }

  /**
   * @return a reference to an element report by its ID within the buffer (0 -> collected count-1).
   */
  ElementReport const & operator[]( ElementCount id ) const
  { return m_elementsBuffer[id]; }

  /**
   * @return iterator pointing at beginning of collected elements in the buffer.
   */
  auto begin() const
  { return m_elementsBuffer.begin(); }

  /**
   * @return iterator pointing after the last collected element in the buffer.
   */
  auto end() const
  { return m_elementsBuffer.begin() + getCollectedElementsCount(); }

  /**
   * @return true when the collection of elements is enabled.
   */
  bool enabled() const
  { return !m_elementsCounter.empty(); }

  /**
   * @return true when there are no elements collected (always false when enabled() is false).
   */
  bool empty() const
  { return getCollectedElementsCount() == 0; }

  /**
   * @return true if the collection of elements completely fills the buffer.
   */
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

  // Preallocated array of ids of detected elements
  array1d< ElementReport > m_elementsBuffer;

};

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_SOLUTIONCHECKHELPERS_HPP
