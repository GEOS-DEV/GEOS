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
 * SidreWrapper.hpp
 *
 *  Created on: Jul 22, 2016
 *      Author: rrsettgast
 */

#ifndef COMPONENTS_CORE_SRC_DATAREPOSITORY_SIDREWRAPPER_HPP_
#define COMPONENTS_CORE_SRC_DATAREPOSITORY_SIDREWRAPPER_HPP_

#include "common/GeosxConfig.hpp"
#include <string>

#ifdef USE_ATK
#include "sidre/DataStore.hpp"
#include "sidre/IOManager.hpp"
#endif
#include <mpi.h>

namespace geosx
{
namespace dataRepository
{

class SidreWrapper
{
public:
  SidreWrapper();
  ~SidreWrapper();
  
#ifdef USE_ATK
  static axom::sidre::DataStore& dataStore();
#endif

  static void writeTree(int num_files, const std::string & path, const std::string & protocol, MPI_Comm comm);

  static void reconstructTree(const std::string & root_path, const std::string & protocol, MPI_Comm comm);

  static void loadExternalData(const std::string & root_path, MPI_Comm comm);

private:

};

} /* namespace dataRepository */
} /* namespace geosx */

#endif /* COMPONENTS_CORE_SRC_DATAREPOSITORY_SIDREWRAPPER_HPP_ */
