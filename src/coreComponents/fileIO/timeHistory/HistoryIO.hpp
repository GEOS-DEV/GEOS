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
};

}
#endif
