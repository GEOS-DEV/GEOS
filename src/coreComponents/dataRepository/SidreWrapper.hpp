/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

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

#ifdef GEOSX_USE_ATK
#include "axom/sidre/core/sidre.hpp"
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

#ifdef GEOSX_USE_ATK
  static axom::sidre::DataStore & dataStore();
#endif

  static void writeTree( int num_files, const std::string & path, const std::string & protocol, MPI_Comm comm );

  static void reconstructTree( const std::string & root_path, const std::string & protocol, MPI_Comm comm );

  static void loadExternalData( const std::string & root_path, MPI_Comm comm );

private:

};

} /* namespace dataRepository */
} /* namespace geosx */

#endif /* COMPONENTS_CORE_SRC_DATAREPOSITORY_SIDREWRAPPER_HPP_ */
