/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron 
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file LogHistory.hpp
 */

#ifndef GEOS_COMMON_LOGHISTORY_HPP
#define GEOS_COMMON_LOGHISTORY_HPP

#include "common/LogMsg.hpp"


namespace geos
{ // TODO document


class LogHistory {
public:

  LogHistory() = default;
  LogHistory( LogHistory && other ) = default;
  LogHistory( LogHistory const & other ) = delete;
  LogHistory & operator=( LogHistory const & other ) = delete;
  LogHistory & operator=( LogHistory && other ) = delete;


  void setTimeStepsToKeep( integer numberOfTimeSteps );

  void storeMessage( LogMsg msg );

private:

  // settings //
  integer m_numberOfTimeSteps = 2;


  // storing //
  using MsgBuffer = std::vector< LogMsg >;
  using MsgBufferList = std::vector< MsgBuffer >;
  using MsgBufferByTimeSteps = std::map< real64, MsgBuffer * >;

  /// @brief List of buffer of messages from the same time-step.
  /// When needed, a buffer can be recycled to store the new time-step messages.
  MsgBufferList m_msgBuffers;

  /// @brief A map allowing to get the right message buffer for a given time-step.
  MsgBufferByTimeSteps m_buffersByTimeStep;


  void addNewTimeStep( real64 newTimeStep ); // can keep that private ? can do that on store ?

  // void flushOldMessages(); // useless ?

}


} // namespace geos

#endif /* GEOS_COMMON_LOGHISTORY_HPP */