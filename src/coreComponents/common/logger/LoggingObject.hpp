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
 * @file LoggingObject.hpp
 */

#ifndef GEOS_COMMON_LOGGINGOBJECT_HPP
#define GEOS_COMMON_LOGGINGOBJECT_HPP

#include "common/LogMsg.hpp"

namespace geos
{

namespace logging
{

/**
 * @brief Interface (implementation contract) of any class that can send log messages (GeneralLogger, Group).
 * @tparam CONTEXT_T Data structure of the context of this class (LogMsgGeneralContext, LogMsgTargetContext).
 * @todo adapt the new call convention everywhere
 * @todo Group class should implement it to log contextualized messages (implement logMsg() as private)
 * @todo check / add documentation
 */
template< typename CONTEXT_T >
class LoggingObject
{
protected:

  /**
   * @brief setting / logging methods. Meant to be public in inheriting classes
   */
  ///@{

  void setLogLevel( LogLevel value );

  virtual void logRank0( LogMsg::CallParams callParams, INPUTS ... msgs ) const = 0;
  virtual void logRank0( LogMsg::CallParams callParams, LogRouter router, INPUTS ... msg ) const = 0;
  virtual void logRank0Fmt( LogMsg::CallParams callParams, string_view msg, INPUTS ... fmtParams ) const = 0;
  virtual void logRank0Fmt( LogMsg::CallParams callParams, LogRouter router, string_view msg, INPUTS ... fmtParams ) const = 0;

  virtual void log( LogMsg::CallParams callParams, INPUTS ... msgs ) const = 0;
  virtual void log( LogMsg::CallParams callParams, LogRouter router, INPUTS ... msg ) const = 0;
  virtual void logFmt( LogMsg::CallParams callParams, string_view msg, INPUTS ... fmtParams ) const = 0;
  virtual void logFmt( LogMsg::CallParams callParams, LogRouter router, string_view msg, INPUTS ... fmtParams ) const = 0;

  // TODO not exactly sure how to implement that (avoid formating if cond is false, keeping condition as metadata). Should we drop this in favor of classic compact if() expressions?
  // virtual void logRank0If( bool cond, LogMsg::CallParams callParams, string_view msg, LogRouter router = LogRouter::m_main ) const = 0;

  // TODO not exactly sure how to implement that (avoid formating if cond is false, keeping condition as metadata). Should we drop this in favor of classic compact if() expressions?
  // virtual void logIf( bool cond, LogMsg::CallParams callParams, string_view msg, LogRouter router = LogRouter::m_main ) const = 0;

  //TODO : analyze how we could support this in an non-intrusive way (about performances / CUDA cache usage) or as a developer limited feature
  // logDevice( SrcCodeLoc callParams, INPUTS ... inputs )

  ///@}

  /**
   * @brief Internal functionnal methods.
   */
  ///@{

  LogMsg::Context createMessageContext( LogMsg::CallParams callParams ) const = 0;

  CONTEXT_T getLoggingContext() const = 0;

  ///@}

};

/**
 * @brief LoggingObject implementations
 */
///@{

/**
 * @todo docs
 */
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

/**
 * @todo docs
 */
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

///@}

} /* namespace logging */

/**
 * @brief Global logger for general messages
 */
extern GeneralLogger g_logger;

} /* namespace geos */

#endif /* GEOS_COMMON_LOGGINGOBJECT_HPP */