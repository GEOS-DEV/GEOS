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
