/*
 * SidreWrapper.cpp
 *
 *  Created on: Jul 22, 2016
 *      Author: rrsettgast
 */
#include <string>
#if ATK_FOUND
#include "SidreWrapper.hpp"
#include <mpi.h>
#include "spio/IOManager.hpp"
#endif

namespace geosx
{
namespace dataRepository
{

SidreWrapper::SidreWrapper()
{}

SidreWrapper::~SidreWrapper()
{}

axom::sidre::DataStore& SidreWrapper::dataStore()
{
#if ATK_FOUND
  static axom::sidre::DataStore datastore;
  return datastore;
#endif
}

/* Write out a restart file. */
void SidreWrapper::writeTree(int num_files, const std::string & path, const std::string & protocol, MPI_Comm comm) 
{
#if ATK_FOUND
  axom::spio::IOManager ioManager(comm);
  ioManager.write(SidreWrapper::dataStore().getRoot(), num_files, path, protocol);
#endif
}

void SidreWrapper::reconstructTree(const std::string & root_path, const std::string & protocol, MPI_Comm comm) 
{
#if ATK_FOUND
  if (!SidreWrapper::dataStore().hasAttribute("__sizedFromParent__"))
  {
    SidreWrapper::dataStore().createAttributeScalar("__sizedFromParent__", -1);
  }

  axom::spio::IOManager ioManager(comm);
  ioManager.read(SidreWrapper::dataStore().getRoot(), root_path, protocol);
#endif
}

/* Load sidre external data. */
void SidreWrapper::loadExternalData(const std::string & root_path, MPI_Comm comm)
{
#if ATK_FOUND
  axom::spio::IOManager ioManager(comm);
  ioManager.loadExternalData(SidreWrapper::dataStore().getRoot(), root_path);
#endif
}

} /* end namespace dataRepository */
} /* end namespace geosx */

