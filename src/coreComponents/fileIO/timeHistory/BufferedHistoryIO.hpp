/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_FILEIO_TIMEHISTORY_BUFFEREDHISTORYIO_HPP_
#define GEOS_FILEIO_TIMEHISTORY_BUFFEREDHISTORYIO_HPP_

#include "common/DataTypes.hpp"

namespace geos
{

/**
 * @class BufferedHistoryIO
 * @brief An abstract class for performing buffered history output.
 */
class BufferedHistoryIO
{
public:
  /// Destructor
  virtual ~BufferedHistoryIO() {}

  /**
   * @brief Get the head of the internal history storage buffer.
   * @return The head of the internal history storage buffer to be written to.
   */
  buffer_unit_type * getBufferHead();

  /**
   * @brief Update the number of items being stored for IO in this object.
   * @param count [in] The new number of items being collected
   */
  void updateCollectingCount( localIndex count );

  /**
   * @brief Query the number of history states currently stored in the internal buffer.
   * @return The number of discrete time history records buffered to be written.
   * @note Since the size of each discrete time history can change, this should not be used
   *       to calculate size, but is useful to check for consistency between collectors
   *       that should be operating at the same cadence (ie the time collector and the data collector).
   */
  localIndex getBufferedCount()
  { return m_bufferedCount; }

  /**
   * @brief Perform and intialization needed for time-history output.
   * @param[in] existsOkay Whether it is acceptable for the intended output target to already exist
   *                       ( false on start from scratch, true on restart ).
   */
  virtual void init( bool existsOkay ) = 0;

  /**
   * @brief Write the buffered history data to the output target.
   */
  virtual void write() = 0;

  /**
   * @brief Do any final cleanup necessary to complete the TimeHistory output.
   */
  virtual void finalize() = 0;

  /**
   * @brief Whether the underlying file being written to can support multiple
   *       TimeHistories being written to it, or only accomodate a single history output.
   */
  virtual bool shareable() const = 0;

protected:
  BufferedHistoryIO( std::type_index typeIdx,
                     localIndex rank,
                     std::vector< localIndex > const & dims );

  /**
   * @brief Get the size in bytes the buffer is currently set to hold per collection operation.
   * @return The size in bytes.
   */
  size_t getRowBytes();

  /// @brief Empty the history collection buffer
  void emptyBuffer();

  /// @brief Resize the buffer to accomodate additional history collection.
  void resizeBuffer();

  /// How much to scale the internal and file allocations by when room runs out
  static inline constexpr localIndex m_overallocMultiple = 2;

  /// The current number of records in the buffer
  localIndex m_bufferedCount;
  /// The write head of the buffer
  buffer_unit_type * m_bufferHead;
  /// The data buffer containing the history info
  buffer_type m_dataBuffer;
  ///
  std::vector< size_t > m_bufferedLocalIdxCounts;


  // history metadata
  /// The underlying data type for this history data set
  std::type_index m_typeIdx;
  /// The size in byte of the data type
  size_t m_typeSize;
  /// The number of variables / atoms of the underlying data type in each element to be collected
  size_t m_typeCount;   // prod(dims[0:n])
  /// Whether the size of the collected data has changed between writes to file
  int m_sizeChanged;

    /// The rank of the data set
  size_t m_rank;
  /// The dimensions of the data set
  std::vector< size_t > m_dims;
};

}
#endif
