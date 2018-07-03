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



#ifndef TIMEHISTORYFILEMANAGER_H_
#define TIMEHISTORYFILEMANAGER_H_

#include "Common/Common.h"
#include <fstream>

class TimeHistoryFileManager
{
public:
  TimeHistoryFileManager();
  ~TimeHistoryFileManager();

  void AddField( const int numFields );

private:

//  std::map<std::string,std::ofstream> m_timeHistoryFiles;



};

#endif /* TIMEHISTORYFILEMANAGER_H_ */
