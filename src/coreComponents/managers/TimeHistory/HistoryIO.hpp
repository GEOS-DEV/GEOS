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

#ifndef GEOSX_HISTORY_IO_HPP_
#define GEOSX_HISTORY_IO_HPP_

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
    m_buffered_count( 0 ),
    m_buffer_head( nullptr ),
    m_data_buffer( 0 )
  {}

  /// Destructor
  virtual ~BufferedHistoryIO() {}

  /**
   * @brief Get the head of the internal history storage buffer.
   * @return The head of the internal history storage buffer to be written to.
   * @note Depends on the virtual function resizeBuffer() being implemented correctly in
   *        an inheriting class.
   */
  buffer_unit_type * GetBufferHead( )
  {
    resizeBuffer();
    m_buffered_count++;
    return m_buffer_head;
  }

  /**
   * @brief Perform and intialization needed for time-history output.
   * @param exists_okay Whether it is acceptable for the intended output target to already exist ( false on start from scratch, true on
   * restart ).
   * data).
   */
  virtual void Init( bool exists_okay ) = 0;

  /**
   * @brief Write the buffered history data to the output target.
   */
  virtual void Write( ) = 0;

  /**
   * @brief Ensure the repressentation of the data in the output target is dense and terse.
   * @note Typically the file will be oversized to receive data writes without constant resizing.
   */
  virtual void CompressInFile( ) = 0;

  /**
   * @brief Get the offset for this history data across the communicator.
   * @return The first index in a theoretical contiguous parallel array being collected on this process.
   */
  virtual globalIndex GetRankOffset( ) = 0;

  /**
   * @brief Query the number of history states currently stored in the internal buffer.
   * @return The number of discrete time history records buffered to be written.
   */
  localIndex GetBufferedCount( ) { return m_buffered_count; }
protected:
  /// @brief Resize the buffer to accomodate additional history collection.
  virtual void resizeBuffer( ) = 0;

  /// @brief Empty the history collection buffer
  void EmptyBuffer( )
  {
    m_buffered_count = 0;
    m_buffer_head = &m_data_buffer[0];
  }

  /// The current number of records in the buffer
  localIndex m_buffered_count;
  /// The write head of the buffer
  buffer_unit_type * m_buffer_head;
  /// The data buffer containing the history info
  buffer_type m_data_buffer;
};

}
#endif
