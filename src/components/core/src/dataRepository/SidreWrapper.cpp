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
 * SidreWrapper.cpp
 *
 *  Created on: Jul 22, 2016
 *      Author: rrsettgast
 */
#include "SidreWrapper.hpp"
#include "dataRepository/Buffer.hpp"

#ifdef USE_ATK
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
  axom::sidre::IOManager ioManager(comm);
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
  
  axom::sidre::IOManager ioManager(comm);
  ioManager.read(SidreWrapper::dataStore().getRoot(), root_path, protocol);
#endif
}


/* Load sidre external data. */
void SidreWrapper::loadExternalData(const std::string & root_path, MPI_Comm comm)
{
#ifdef USE_ATK
  axom::sidre::IOManager ioManager(comm);
  ioManager.loadExternalData(SidreWrapper::dataStore().getRoot(), root_path);
#endif
}

} /* end namespace dataRepository */
} /* end namespace geosx */

