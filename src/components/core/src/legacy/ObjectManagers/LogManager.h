/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


/**
 * @file ChemistryManager.h
 * @author walsh24
 * @date April 4, 2013
 */

#ifndef LOGMANAGER_H_
#define LOGMANAGER_H_

#include "Common/typedefs.h"
#include "Utilities/Utilities.h"
#include "DataStructures/Tables/Table.h"

#include <map>
#include <utility>

#include "../IO/ticpp/HierarchicalDataNode.h"

#if GPAC_MPI
#include <mpi.h>
#endif

/// Log level typedef
enum logLevelT
{
  LOG_NOTHING,
  LOG_CRITICAL,
  LOG_ERROR,
  LOG_WARNING,
  LOG_INFO,
  LOG_DEBUG
};

class LogManager
{
public:

  static LogManager& Instance()
  {
    static LogManager theLogManager;

    return theLogManager;
  }
  void SetLocalLevel(logLevelT level){m_local_level = level;};
  void RestoreGlobalLevel(void){m_local_level = LOG_NOTHING;}; // default to
                                                               // global level
  void SetGlobalLevel(logLevelT level){m_global_level = level;};
  logLevelT GetLevel(void){return std::max(m_local_level,m_global_level);};

private:

  LogManager():
    m_local_level(LOG_NOTHING),
    m_global_level(LOG_NOTHING)
  {
    /*empty*/
  };
  ~LogManager() {}
  LogManager( const LogManager& );
  LogManager& operator=( const LogManager& );

  logLevelT m_local_level;
  logLevelT m_global_level;

};

namespace LOG_UTILS
{
inline bool StrMessageLogger(logLevelT level, const std::string& msg){
  LogManager& logManager = LogManager::Instance();
  bool wasCalled = false;
  if ( level <= logManager.GetLevel() )
  {
    std::cout << msg << std::endl;
    wasCalled = true;
  }
  return wasCalled;
}
}

// Wrapper function. StrMessageLogger is not used directly.
// Log<LOG_WARNING>("Abandon hope.");
template <logLevelT level> bool Log(const std::string& msg) {
  return LOG_UTILS::StrMessageLogger( level, msg );
}


inline
std::istream& operator >>(std::istream &is, logLevelT& level)
{
  int val;
  is >> val;

  level = (logLevelT) val;

  return is;
}



#endif /* LOGMANAGER_H_ */
