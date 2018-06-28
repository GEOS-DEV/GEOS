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

/*
 * LogStream.h
 *
 *  Created on: Feb 24, 2011
 *      Author: johnson346
 */

#ifndef LOGSTREAM_H_
#define LOGSTREAM_H_

// ***** Included Headers *****************************************************
#include <fstream>
#include <iostream>
#include <string>
#include "Common/intrinsic_typedefs.h"

namespace LogInfo
{
enum LogEnum
{
  noInfo                  =  0,
  debugInfo               =  1,
  inputInfo               =  2,
  numLogEnums
};
}

// ****************************************************************************
// ***** LOGSTREAM CLASS DECLARATION ******************************************
// ****************************************************************************
class LogStream
{
public:
  LogStream(void);
  ~LogStream(void);

  static inline std::ostream& Write(const std::string& msg, const LogInfo::LogEnum type)
  {
    Initialize();
    (*m_logEnumToStream[(INDEX)type]) << msg << std::endl;
    return std::cout;
  }

private:
  static std::ostream* m_logEnumToStream[LogInfo::numLogEnums];

  static inline void Initialize()
  {
    static bool initialized = false;
    if(!initialized)
    {
      for(INDEX i = 0 ; i < LogInfo::numLogEnums ; i++)
        m_logEnumToStream[i] = &std::cout;
      initialized = true;
    }
  }
};

LogStream::LogStream(){}
#endif /* LOGSTREAM_H_ */
