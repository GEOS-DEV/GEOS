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

#include "LogLevelsRegistry.hpp"

void LogLevelsRegistry::addEntry( string string, string description )
{
  m_logLevelsDescriptions[ string( levelCondition ) ].push_back( string( logDescription ) );
}

/**
 * @brief Construct the log level string description for a wrapper
 */
void LogLevelsRegistry::buildLogLevelDescription( std::map< std::string, std::vector< std::string > > const & logLevelsDescriptions )
{
  std::ostringstream description;
  description << "Sets the level of information to write in the standard output (the console typically).\n"
                 "A level of 0 outputs minimal information, higher levels require more.";
  for( auto const & [logLevel, logDescriptions] : logLevelsDescriptions )
  {
    description << GEOS_FMT( "\n{}\n", logLevel );
    for( size_t idxDesc = 0; idxDesc< logDescriptions.size(); idxDesc++ )
    {
      description << " - " << logDescriptions[idxDesc];
      if( idxDesc != logDescriptions.size() - 1 )
        description << '\n';
    }
  }
  appendDescription( description.str() );
}