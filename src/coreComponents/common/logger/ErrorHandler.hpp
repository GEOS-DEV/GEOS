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
 * @file ErrorHandler.hpp
 */

#ifndef INITIALIZATION_ERROR_HANDLER_HPP
#define INITIALIZATION_ERROR_HANDLER_HPP

#include "ErrorLogger.hpp"

namespace geos
{

class ErrorHandler
{
public:

  /**
   * @brief Singleton pattern (as we can only have one strategy of handling at the same time).
   *        Update the current error handling strategy. Once setup, the instance is not meant to be mutable.
   * @return ErrorHandler const&
   */
  static void setupErrorHandlingStrategy( ErrorHandler && instance );

  /**
   * @brief Singleton pattern (as we can only have one strategy of handling at the same time).
   * @return Handler instance currently setup.
   */
  static ErrorHandler const & getInstance();

  /**
   * @brief Construct a new Error Handler object. To use it as the current error handling strategy,
   *        setupErrorHandlingStrategy() must be called on this instance after setting it up.
   */
  ErrorHandler();

  /**
   * @brief Enable pipe deviation of external stderr messages,
   * @see ExternalErrorHandler
   */
  void enableExternalErrorPipeDeviation( bool enable )
  { m_isPipeDeviationEnabled = enable; }

  /**
   * @brief Set the callback to call when post-error program abortion is requested.
   * @param abortingFunctor a void() functor that terminates the program in a safe way.
   */
  void setProgramAborter( std::function< void() > const & abortingFunctor )
  { m_abortingFunctor = abortingFunctor; }

  /**
   * @brief the ErrorLogger where error messages are rooted
   * @param logger the instance of the logger, ErrorLogger::instance() or a local one.
   */
  void setLogger( ErrorLogger * logger )
  { m_logger = logger; }

  /**
   * @return the ErrorLogger where error messages are rooted
   */
  ErrorLogger & getLogger() const
  { return *m_logger; }

  /**
   * @brief Post error-handling function that terminates the program.
   */
  void abortProgram() const
  { m_abortingFunctor(); }

private:

  static ErrorHandler s_errorHandlerInstance;

  ErrorLogger * m_logger;

  bool m_isPipeDeviationEnabled = false;

  std::function< void() > m_abortingFunctor;

  /**
   * @brief set external error handling behaviour
   */
  void setupExternalErrorManagment();

  /**
   * @brief set signal handling behaviour
   */
  void setupSignalHandling();

  /**
   * @brief Set LvArray to use the GEOS error behaviour in a way that can be traced
   */
  void setupLvArrayErrorHandling();

};

} /* namespace geos */

#endif
