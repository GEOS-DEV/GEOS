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
 * @file ErrorHandler.cpp
 */

#include "ErrorHandler.hpp"

#include "ExternalErrorHandler.hpp"
#include "Logger.hpp"

namespace geos
{

ErrorHandler ErrorHandler::s_errorHandlerInstance;

ErrorHandler::ErrorHandler()
  : m_logger( &ErrorLogger::global() )
  , m_abortingFunctor( [this]() { std::abort(); } )
{}

ErrorHandler const & ErrorHandler::getInstance()
{
  return s_errorHandlerInstance;
}

void ErrorHandler::setupErrorHandlingStrategy( ErrorHandler && instance )
{
  s_errorHandlerInstance = instance;
  s_errorHandlerInstance.setupExternalErrorManagment();
  s_errorHandlerInstance.setupSignalHandling();
  s_errorHandlerInstance.setupLvArrayErrorHandling();
}

void ErrorHandler::setupExternalErrorManagment()
{
  ExternalErrorHandler::instance().enableStderrPipeDeviation( true );

  ExternalErrorHandler::instance().setErrorHandling( [&]( string_view errorMsg,
                                                          string_view detectionLocation )
  {
    // Filter out INFO level messages from external libraries (e.g., VTK)
    // ( error / signal lambda would calls either an error function or an info function, depending on a filtering function )
    if( ExternalErrorHandler::isNotAnErrorMsg( errorMsg ) )
    {
      // Just print the message without error formatting
      GEOS_LOG( errorMsg );
      return;
    }
    else
    {
      std::string const stackHistory = LvArray::system::stackTrace( true );
      DiagnosticMsg diagnosticMsg;
      m_logger->flushErrorMsg( DiagnosticMsgBuilder::init( diagnosticMsg,
                                                           MsgType::Error, errorMsg,
                                                           ::geos::logger::internal::g_rank )
                                 .addCallStackInfo( stackHistory )
                                 .addDetectionLocation( detectionLocation )
                                 .setCause( "Error pipe output from a dependency" )
                                 .getDiagnosticMsg() );

      // we do not terminate the program as 1. the error could be non-fatal, 2. there may be more messages to output.
    }
  } );
}

void ErrorHandler::setupSignalHandling()
{
  LvArray::system::setSignalHandling( []( int const signal )
  {
    // setSignalHandling() cannot take capturing lambda, so we get the instance with the static function.
    ErrorHandler const & instance = ErrorHandler::getInstance();

    // Disable signal handling to prevent catching exit signal (infinite loop)
    LvArray::system::setSignalHandling( nullptr );

    // first of all, there can be external error that await to be output
    ExternalErrorHandler::instance().flush( "before signal error output" );

    // error message output
    std::string const stackHistory = LvArray::system::stackTrace( true );
    DiagnosticMsg diagnosticMsg;
    DiagnosticMsgBuilder::init( diagnosticMsg,
                                MsgType::ExternalError, "",
                                ::geos::logger::internal::g_rank )
      .addSignal( signal )
      .addCallStackInfo( stackHistory );
    instance.m_logger->flushErrorMsg( diagnosticMsg );

    // call program termination
    instance.abortProgram();
  } );
}

void ErrorHandler::setupLvArrayErrorHandling()
{
  LvArray::system::setErrorHandler( [&]()
  {
    // first of all, there can be external error that await to be output
    ExternalErrorHandler::instance().flush( "before LvArray error handling" );

    // default error message output
    DiagnosticMsg diagnosticMsg;
    DiagnosticMsgBuilder::init( diagnosticMsg,
                                MsgType::ExternalError,
                                "LvArray Runtime Error",
                                ::geos::logger::internal::g_rank )
      .addCallStackInfo( LvArray::system::stackTrace( true ) )
      .addDetectionLocation( "LvArray Error Handler" );
    m_logger->flushErrorMsg( diagnosticMsg );

    // call program termination
    abortProgram();
  } );
}


} /* namespace geos */
