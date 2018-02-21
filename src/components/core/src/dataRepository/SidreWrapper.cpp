/*
 * SidreWrapper.cpp
 *
 *  Created on: Jul 22, 2016
 *      Author: rrsettgast
 */
#include "SidreWrapper.hpp"
#include "dataRepository/Buffer.hpp"

#ifdef USE_ATK
#include "spio/IOManager.hpp"
#include "slic/slic.hpp"
#endif

#include <string>
#include <cstdio>
#include <mpi.h>


namespace geosx
{
namespace dataRepository
{

#ifdef USE_ATK
using namespace axom::sidre;
#endif

SidreWrapper::SidreWrapper()
{}

SidreWrapper::~SidreWrapper()
{}

#ifdef USE_ATK
DataStore& SidreWrapper::dataStore()
{
  static DataStore datastore;
  return datastore;
}
#endif


/* Write out a restart file. */
void SidreWrapper::writeTree(int num_files, const std::string & path, const std::string & protocol, MPI_Comm comm) 
{
#ifdef USE_ATK
  axom::spio::IOManager ioManager(comm);
  ioManager.write(SidreWrapper::dataStore().getRoot(), num_files, path, protocol);
#endif
}

void SidreWrapper::reconstructTree(const std::string & root_path, const std::string & protocol, MPI_Comm comm) 
{
#ifdef USE_ATK
  if (!SidreWrapper::dataStore().hasAttribute("__sizedFromParent__"))
  {
    SidreWrapper::dataStore().createAttributeScalar("__sizedFromParent__", -1);
  }
  
  axom::spio::IOManager ioManager(comm);
  ioManager.read(SidreWrapper::dataStore().getRoot(), root_path, protocol);
#endif
}


bool checkSidreTree(Group * group) 
{
#ifdef USE_ATK
  bool valid = true;
  for (int i = group->getFirstValidViewIndex(); i != InvalidIndex; i = group->getNextValidViewIndex(i)) 
  {
    View * view = group->getView(i);
    if (view->isExternal() && view->getVoidPtr() == nullptr) {
      SLIC_WARNING("Pointer should be valid: " << view->getPathName() << std::endl);
      valid = false;
    }
  }

  for (int i = group->getFirstValidGroupIndex(); i != InvalidIndex; i = group->getNextValidGroupIndex(i)) 
  {
    Group * childGroup = group->getGroup(i);
    valid &= checkSidreTree(childGroup);
  }

  return valid;

#else
  return false;
#endif
}

/* Load sidre external data. */
void SidreWrapper::loadExternalData(const std::string & root_path, MPI_Comm comm)
{
#ifdef USE_ATK
  if (!checkSidreTree(SidreWrapper::dataStore().getRoot())) {
    SLIC_ERROR("Tree not valid");
  }

  axom::spio::IOManager ioManager(comm);
  ioManager.loadExternalData(SidreWrapper::dataStore().getRoot(), root_path);
#endif
}

} /* end namespace dataRepository */
} /* end namespace geosx */

