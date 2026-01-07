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
 * @file Logger.hpp
 */


#include "common/logger/ErrorHandling.hpp"

namespace geos
{

/**
 * @brief Geos Exception used in GEOS_THROW
 */
struct Exception : public std::exception
{
public:
  Exception() = default;

  /**
   * @brief Constructor
   * @param what The error message describing the exception
   */
  Exception( std::string const & what ):
    std::exception( )
  {
    m_errorMsg.m_msg = what;
  }

  /**
   * @return Reference to the ErrorMsg object
   */
  DiagnosticMsg & getErrorMsg()
  { return m_errorMsg; }

  /**
   * @brief System fallback to get description content if error system does not achieve to output the ErrorMsg
     can also provide exception information to debuggers
   * @return Returns The exception message
   * @note We does not allow to override what(), it's the GEOS_THROW responsability to write-in the exception message
   */
  virtual char const * what() const noexcept override final
  {
    return m_cachedWhat.c_str();
  }

  /**
   * @brief Prepare and cache the formatted exception message
   * @param msg The main error message
   * @param cause The cause of the exception
   * @param file The source file where the exception occurred
   * @param line The line number where the exception occurred
   * @param rank The MPI rank
   * @param stackTrace The stack trace at the point of exception
   */
  void prepareWhat( std::string const & msg,
                    std::string const & cause,
                    char const * file,
                    int line,
                    int rank,
                    string_view stackTrace ) noexcept
  {
    try
    {
      std::ostringstream oss;
      oss << "***** GEOS Exception\n";
      oss << "***** LOCATION: " << file << " l." << line << "\n";
      oss << "***** " << cause << "\n";
      oss << "***** Rank  " << rank << ": "<< msg <<"\n\n";
      oss << stackTrace;
      m_cachedWhat = oss.str();
    } catch( ... )
    {
      m_cachedWhat = "GEOS Exception (formatting failed)";
    }

  }

private:
  /// Formatted exception message for what() method
  string m_cachedWhat;
  /// Error message logger for structured error reporting
  DiagnosticMsg m_errorMsg;
};

/**
 * @brief Exception class used to report errors in user input.
 */
struct RuntimeError : public geos::Exception
{
  /**
   * @brief Constructor
   * @param what the error message
   */
  RuntimeError( std::string const & what ):
    geos::Exception( what )
  {}

  RuntimeError(): geos::Exception(){}
};

/**
 * @brief Exception class used to report bad GEOS state
 */
struct LogicError : public geos::Exception
{
  /**
   * @brief Constructor
   * @param what the error message
   */
  LogicError( std::string const & what ):
    geos::Exception( what )
  {}

  /**
   * @brief Default constructor
   */
  LogicError(): geos::Exception(){}
};

/**
 * @brief Exception class used to report errors in user input.
 */
struct DomainError : public geos::Exception
{
  /**
   * @brief Constructor
   * @param what the error message
   */
  DomainError( std::string const & what ):
    geos::Exception( what )
  {}

  /**
   * @brief Default constructor
   */
  DomainError(): geos::Exception(){}
};
/**
 * @brief Exception class used to report errors in user input.
 */
struct InputError : public geos::Exception
{
  /**
   * @brief Constructor
   * @param what the error message
   */
  InputError( std::string const & what ):
    geos::Exception( what )
  {}

  /**
   * @brief Default constructor
   */
  InputError(): geos::Exception(){}

  /**
   * @brief Constructor
   * @param what the error message
   */
  InputError( char const * const what ):
    geos::Exception( what )
  {}

  /**
   * @brief Constructs an InputError from an underlying exception.
   * @param subException The exception on which the created one is based.
   * @param msgToInsert The error message that will be inserted in the subException error message.
   */
  InputError( std::exception const & subException, std::string const & msgToInsert );
};

/**
 * @brief Exception class used to report errors in user input.
 */
struct SimulationError : public geos::Exception
{
  /**
   * @brief Constructor
   * @param what the error message
   */
  SimulationError( std::string const & what ):
    geos::Exception( what )
  {}

  /**
   * @brief Default constructor
   */
  SimulationError(): geos::Exception(){}

  /**
   * @brief Constructor
   * @param what the error message
   */
  SimulationError( char const * const what ):
    geos::Exception( what )
  {}

  /**
   * @brief Construct a SimulationError from an underlying exception.
   * @param subException An exception to base this new one on.
   * @param msgToInsert The error message.
   * It will be inserted before the error message inside of subException.
   */
  SimulationError( std::exception const & subException, std::string const & msgToInsert );
};

/**
 * @brief Exception class used to report errors from type conversion
 * @todo (ErrorManager EPIC #2940) Consider adding a way to precise custom exception parameters, to add
 * expected & encountered typeid for this one (in order to manage the exception output more precisely).
 * We could also manage this by having:    BadTypeErrorABC <|--- BadTypeError< T >      /!\ compilation time
 */
struct BadTypeError : public geos::Exception
{
  /**
   * @brief Constructor
   * @param what the error message
   */
  BadTypeError( std::string const & what ):
    geos::Exception( what )
  {}

  /**
   * @brief Default constructor
   */
  BadTypeError(): geos::Exception(){}
};

/**
 * @brief Exception class used for special control flow.
 */
class NotAnError : public geos::Exception
{};

}
