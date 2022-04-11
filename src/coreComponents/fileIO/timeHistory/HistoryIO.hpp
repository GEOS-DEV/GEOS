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

/**
 * @file HistoryIO.hpp
 */

#ifndef GEOSX_FILEIO_TIMEHISTORY_HISTORYIO_HPP_
#define GEOSX_FILEIO_TIMEHISTORY_HISTORYIO_HPP_

#include "dataRepository/Group.hpp"

namespace geosx
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
   * @note Depends on the virtual function resizeBuffer() being implemented correctly in
   *        an inheriting class.
   */
  virtual buffer_unit_type * getBufferHead() = 0;

  /**
   * @brief Perform and intialization needed for time-history output.
   * @param existsOkay Whether it is acceptable for the intended output target to already exist ( false on start from scratch, true on
   *                    restart ).
   */
  virtual void init( bool existsOkay ) = 0;

  /**
   * @brief Write the buffered history data to the output target.
   */
  virtual void write() = 0;

  /**
   * @brief Ensure the repressentation of the data in the output target is dense and terse.
   * @note Typically the file will be oversized to receive data writes without constant resizing.
   */
  virtual void compressInFile() = 0;

  /**
   * @brief Update the number of items being stored for IO in this object.
   * @param count [in] The new number of items being collected
   */
  virtual void updateCollectingCount( localIndex count ) = 0;

  /**
   * @brief Query the number of history states currently stored in the internal buffer.
   * @return The number of discrete time history records buffered to be written.
   * @note Since the size of each discrete time history can change, this should not be used
   *        to calculate size, but is useful to check for consistency between collectors
   *        that should be operating at the same cadence (ie the time collector and the data
   *        collector).
   */
  virtual localIndex getBufferedCount() = 0;
};

}
#endif
