/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file LoggingObject.hpp
 */

#ifndef GEOS_COMMON_LOGGINGOBJECT_HPP
#define GEOS_COMMON_LOGGINGOBJECT_HPP

#include "common/LogMsg.hpp"


namespace geos
{


/**
 * @brief Interface (implementation contract) of any class that can send log messages (GeneralLogger, Group)
 */
// TODO : adapt the new call convention everywhere
// TODO : Group class should implement it to log contextualized messages (implement logMsg() as private)
// TODO : check / add documentation
template< typename CONTEXT_T >
class LoggingObject
{

  void setLogLevel( LogLevel value );


  virtual void logRank0( LogMsg::CallParams loc, string_view msg, LogRouter router = LogRouter::m_main ) const = 0;

  // TODO not exactly sure how to implement that (avoid formating if cond is false, keeping condition as metadata).
  // virtual void logRank0If( bool cond, LogMsg::CallParams loc, string_view msg, LogRouter router = LogRouter::m_main ) const = 0;

  virtual void log( LogMsg::CallParams loc, string_view msg, LogRouter router = LogRouter::m_main ) const = 0;

  // TODO not exactly sure how to implement that (avoid formating if cond is false, keeping condition as metadata).
  // virtual void logIf( bool cond, LogMsg::CallParams loc, string_view msg, LogRouter router = LogRouter::m_main ) const = 0;

  //TODO : logDevice has been excluded as ... document why commented and should not exist (bad perfs, workaround)
  // logDevice( SrcCodeLoc loc, INPUTS ... inputs )

protected:

  LogMsg::Context createMessageContext( LogMsg::CallParams callParams ) const =0;

  CONTEXT_T getLoggingContext() const =0;

};
LogMsg::

class GeneralLogger final : public LoggingObject<LogMsg::Context::General>
{
public:

  void setLogLevel( LogLevel value ) override;

  void setRank( integer value );

  void setRankCount( integer value );

  void setLogSectionTitle( string value );

  void setTimeStepStart( real64 value );


  void logRank0( LogMsg:CallParams loc, string_view msg, LogRouter router = LogRouter::m_main ) const override;

  void log( LogMsg:CallParams loc, string_view msg, LogRouter router = LogRouter::m_main ) const override;

private:

  LogLevel m_generalLogLevel;

  int m_rank;

  int m_ranksCount;

  string m_logSectionTitle;

  real64 m_timeStepStart;
  
  LogMsg::Context::General m_currentContext;


  LogMsg::Context createMessageContext( LogMsg:CallParams callParams ) const override;

  LogMsg::Context::General getLoggingContext() const override;

};
extern GeneralLogger g_logger;


class TargetedLogger final : public LoggingObject<LogMsg::Context::Target>
{
public:

  void setLogLevel( LogLevel value ) override;

  void setTargetName( string value );

  void setTargetDataContext( string value );


  void logRank0( LogMsg:CallParams loc, string_view msg, LogRouter router = LogRouter::m_main ) const override;

  void log( LogMsg:CallParams loc, string_view msg, LogRouter router = LogRouter::m_main ) const override;

private:
  
  LogLevel m_logLevel;
  
  LogMsg::Context::Target m_currentContext;


  LogMsg::Context createMessageContext( LogMsg:CallParams callParams ) const override;

  LogMsg::Context::Target getLoggingContext() const override;

};


} // namespace geos

#endif /* GEOS_COMMON_LOGGINGOBJECT_HPP */