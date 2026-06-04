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
 * @file ElementReporter.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_ELEMENTREPORTER_HPP
#define GEOS_PHYSICSSOLVERS_ELEMENTREPORTER_HPP

#include "physicsSolvers/ElementReporterKernels.hpp"
#include "common/Units.hpp"

namespace geos
{

/**
 * @brief A buffer to count and store element ids during kernel execution.
 *        This facilitates the reporting mechanism by allowing a preallocated space for storing & counting elements.
 */
class ElementReporterBuffer
{
public:

  /// Type alias for elements count (e.g., localIndex, globalIndex).
  using ElementCount = ElementsReporterCollector::ElementCount;

  /**
   * @brief Construct a preallocated buffer to collect a limited quantity of ids in kernels.
   * @param maxCollectionSize Limit of the buffer.
   *                          If 0, the buffering functionnality is disabled and only the counting is enabled.
   */
  ElementReporterBuffer( bool enabled, ElementCount maxCollectionSize );

  /**
   * @brief Transfers ownership of an ElementReporterBuffer to another instance (move semantics).
   */
  ElementReporterBuffer( ElementReporterBuffer && other ) = default;

  /**
   * @brief Transfers ownership of an ElementReporterBuffer to another instance (move semantics).
   */
  ElementReporterBuffer & operator=( ElementReporterBuffer && other ) = default;

  /**
   * @brief Copying prevented as it doesn't seem relevant / useful.
   */
  ElementReporterBuffer( ElementReporterBuffer const & ) = delete;

  /**
   * @brief Copying prevented as it doesn't seem relevant / useful.
   */
  ElementReporterBuffer & operator=( ElementReporterBuffer const & ) = delete;

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
   * @return the maximum elements that can effectivly be stored (zero if no collection is enabled).
   */
  ElementCount getMaxCollectionSize() const
  { return m_elementsBuffer.size(); }

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

  /**
   * @brief Signal that writings are completed for this collection. Must be called before reading / copying
   *        on the host. Called automatically when doing ElementReporterOutput
   */
  void finalizeCollection();

  /**
   * @brief Empty all collection.
   */
  void clear();

private:

  // array of one element to get benefit of managed host-device memory.
  array1d< ElementCount > m_elementsCounter;

  // Preallocated array of ids of detected elements
  array1d< ElementReport > m_elementsBuffer;

};

/**
 * @brief A class to report elements collected by the solver.
 */
class ElementReporterOutput
{
public:

  /// Type alias for elements count (e.g., localIndex, globalIndex).
  using ElementCount = ElementsReporterCollector::ElementCount;

  /// Type to store reports in maps
  using ReportsMap = stdMap< string, ElementReporterOutput >;

  /**
   * @brief Construct a preallocated buffer for collecting element ids in kernels.
   */
  ElementReporterOutput();

  /**
   * @brief Set the values metadata of this reporter
   * @param valueNaming The name used when referring to variables within this context (e.g., "pressure", "density").
   * @param valueUnit Unit in which values of this report are expressed.
   */
  ElementReporterOutput & setValueMetadata( string_view valuesNaming, units::Unit valuesUnit );

  /**
   * @brief Optionally set the prefix to apply to the first line of the log table output
   * @param linesPrefix Prefix for the line of text to be printed
   */
  ElementReporterOutput & setLogPrefix( string_view logPrefix );

  /**
   * @brief Set the acceptable ranges values of the report.
   * @param minValue Allow to optionally show the minimum acceptable values in the report.
   * @param maxValue Allow to optionally show the maximum acceptable values in the report.
   */
  ElementReporterOutput & setRanges( std::optional< real64 > minValue,
                                     std::optional< real64 > maxValue );

  /**
   * @brief Set the values of the report for the current rank.
   * @param buffer The buffer which contains the number and values to show from this rank.
   * @param rankMinValue Optionally store the minimum encountered value of this rank.
   * @param rankMaxValue Optionally store the maximum encountered value of this rank.
   */
  ElementReporterOutput & setValues( ElementReporterBuffer && buffer,
                                     std::optional< real64 > rankMinValue,
                                     std::optional< real64 > rankMaxValue );

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
   * @return The total count of collected elements across all ranks for signaling ids.
   */
  ElementReporterBuffer const & getRankBuffer() const
  { return m_buffer; }

  /**
   * @return The naming of the stored values, .
   */
  string_view getValuesNaming() const
  { return m_valuesNaming; }

  /**
   * @brief Output the elements which have been reported, potentially signaling underflow/overflow or
   *        numerical instability.
   * @param eventTitle Title of the report output, to identify when it occurs during the problem solving.
   */
  void outputReport( string_view eventTitle ) const;

  /**
   * @brief Empty all collection results.
   */
  void clear();

private:

  /// Preallocated buffer for collecting ids.
  ElementReporterBuffer m_buffer;

  /// Count of signaled elements per rank.
  ElementCount m_ranksSignaledElementsCount = 0;

  /// Total collected signaling id count across ranks.
  ElementCount m_ranksCollectedElementsCount = 0;

  string m_valuesNaming = "";

  units::Unit m_valuesUnit = units::Unknown;

  std::optional< real64 > m_minEncounteredValue;

  std::optional< real64 > m_maxEncounteredValue;

  std::optional< real64 > m_minAcceptableValue;

  std::optional< real64 > m_maxAcceptableValue;

  string m_logPrefix = "";
};

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_ELEMENTREPORTER_HPP
