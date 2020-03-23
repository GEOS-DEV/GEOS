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
 * @file CommandLineParser.cpp
 */

#include <stdlib.h>
#include <stdio.h>
#include "CommandLineParser.hpp"

namespace geosx
{


struct Arg : public option::Arg
{
  static void printError( const char * msg1, const option::Option & opt, const char * msg2 )
  {
    fprintf( stderr, "%s", msg1 );
    fwrite( opt.name, static_cast< size_t >(opt.namelen), 1, stderr );
    fprintf( stderr, "%s", msg2 );
  }

  static option::ArgStatus Unknown( const option::Option & option, bool msg )
  {
    if( msg ) printError( "Unknown option '", option, "'\n" );
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus Required( const option::Option & option, bool msg )
  {
    if( option.arg != nullptr )
      return option::ARG_OK;

    if( msg ) printError( "Option '", option, "' requires an argument\n" );
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus NonEmpty( const option::Option & option, bool msg )
  {
    if( option.arg != nullptr && option.arg[0] != 0 )
      return option::ARG_OK;

    if( msg ) printError( "Option '", option, "' requires a non-empty argument\n" );
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus Numeric( const option::Option & option, bool msg )
  {
    char * endptr = nullptr;
    if( option.arg != nullptr && strtol( option.arg, &endptr, 10 )) {};
    if( endptr != option.arg && *endptr == 0 )
      return option::ARG_OK;

    if( msg ) printError( "Option '", option, "' requires a numeric argument\n" );
    return option::ARG_ILLEGAL;
  }
};


CommandLineParser::CommandLineParser()
{
  // TODO Auto-generated constructor stub

}

CommandLineParser::~CommandLineParser()
{
  // TODO Auto-generated destructor stub
}

} /* namespace geosx */
