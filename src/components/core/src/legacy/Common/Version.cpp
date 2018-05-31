// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
/*
 * Version.cpp
 *
 *  created on Apr 25, 2014
 *      Author: Stuart Walsh
 */

#include "Version.h"

#include <string>


void DisplayVersion(int rank)
{
  if( rank==0 )
  {
    std::cout << "GEOS Version 3.0.1" << std::endl;
    std::cout << "                            " << std::endl;
    std::cout << "Git Hash: "<< REPO_VERSION << " (core) " << std::endl;
    if( !std::string(EXTERNAL_REPO_VERSION).empty() )
      std::cout << "          " << EXTERNAL_REPO_VERSION  << " (mod) "<< std::endl;
    if( !std::string(INTERNAL_REPO_VERSION).empty() )
      std::cout << "          " << INTERNAL_REPO_VERSION << " (internal) " << std::endl;
  }
}

void DisplaySplash(int rank)
{
  if( rank==0 )
  {
    std::cout << "                            " << std::endl;
    std::cout << "############################" << std::endl;
    std::cout << "                            " << std::endl;
    std::cout << "     GEOS Version 3.0.1     " << std::endl;
    std::cout << "                            " << std::endl;
    std::cout << "############################" << std::endl;
    std::cout << "                            " << std::endl;
  }
}

void DisplayVersionHistory(int rank)
{
  if( rank==0 )
  {

    // I am a bad man
#ifdef INCLUDE_VERSION_HISTORY
#include "VersionHistory.temp.h"
    std::cout << std::endl;
#endif

  }

}

std::string GetRepoVersionString(){
  std::string rv = std::string(REPO_VERSION) + " (core) ";
  if( !std::string(EXTERNAL_REPO_VERSION).empty() )
    rv += std::string(EXTERNAL_REPO_VERSION)  + " (external) ";
  if( !std::string(INTERNAL_REPO_VERSION).empty() )
    rv += std::string(INTERNAL_REPO_VERSION) + " (internal) ";

  return rv;
}
