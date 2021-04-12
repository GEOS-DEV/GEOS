/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
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
class DataSpec;

/**
 * @class BufferedHistoryIO
 * @brief An abstract class for performing buffered history output.
 */
class BufferedHistoryIO
{
public:

  /// Constructor
  BufferedHistoryIO():
    m_bufferedCount( 0 ),
    m_bufferHead( nullptr ),
    m_dataBuffer( 0 )
  {}

  /// Destructor
  virtual ~BufferedHistoryIO() {}

  /**
   * @brief Get the head of the internal history storage buffer.
   * @return The head of the internal history storage buffer to be written to.
   * @note Depends on the virtual function resizeBuffer() being implemented correctly in
   *        an inheriting class.
   */
  buffer_unit_type * getBufferHead( )
  {
    resizeBuffer();
    m_bufferedCount++;
    return m_bufferHead;
  }

  /**
   * @brief Perform and intialization needed for time-history output.
   * @param existsOkay Whether it is acceptable for the intended output target to already exist ( false on start from scratch, true on
   * restart ).
   * data).
   */
  virtual void init( bool existsOkay ) = 0;

  /**
   * @brief Write the buffered history data to the output target.
   */
  virtual void write( ) = 0;

  /**
   * @brief Ensure the repressentation of the data in the output target is dense and terse.
   * @note Typically the file will be oversized to receive data writes without constant resizing.
   */
  virtual void compressInFile( ) = 0;

  /**
   * @brief Query the number of history states currently stored in the internal buffer.
   * @return The number of discrete time history records buffered to be written.
   */
  localIndex getBufferedCount( ) { return m_bufferedCount; }
protected:
  /// @brief Resize the buffer to accomodate additional history collection.
  virtual void resizeBuffer( ) = 0;

  /// @brief Empty the history collection buffer
  void emptyBuffer( )
  {
    m_bufferedCount = 0;
    m_bufferHead = &m_dataBuffer[0];
  }

  /// The current number of records in the buffer
  localIndex m_bufferedCount;
  /// The write head of the buffer
  buffer_unit_type * m_bufferHead;
  /// The data buffer containing the history info
  buffer_type m_dataBuffer;
};

}
#endif
