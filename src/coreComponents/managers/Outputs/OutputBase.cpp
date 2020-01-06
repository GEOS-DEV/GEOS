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

/**
 * @file OutputBase.cpp
 */

#include "OutputBase.hpp"
#include "mpiCommunications/MpiWrapper.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

OutputBase::OutputBase( std::string const & name,
                        Group * const parent ):
  ExecutableGroup( name, parent),
  m_slaveDirectory(),
  m_parallelThreads(1)
{
  setInputFlags(InputFlags::OPTIONAL_NONUNIQUE);

  registerWrapper(viewKeysStruct::slaveDirectoryString, &m_slaveDirectory, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("slave directory path");

  registerWrapper(viewKeysStruct::parallelThreadsString, &m_parallelThreads, false )->
    setApplyDefaultValue(1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Number of plot files.");

}

OutputBase::~OutputBase()
{}

OutputBase::CatalogInterface::CatalogType& OutputBase::GetCatalog()
{
  static OutputBase::CatalogInterface::CatalogType catalog;
  return catalog;
}


void OutputBase::InitializePreSubGroups( Group * const GEOSX_UNUSED_ARG( group ) )
{
  // This command doesn't seem to work anymore
  // SetupDirectoryStructure();
}


void OutputBase::SetupDirectoryStructure()
{
  string slaveDirectory = m_slaveDirectory;

  int const rank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX);
  if (rank  == 0)
  {
    if (!slaveDirectory.empty())
    {
      string cmd = "mkdir -p " + slaveDirectory;
      int ret = std::system(cmd.c_str());
      GEOSX_ERROR_IF(ret != 0, "Command '" << cmd << "' exited with code " << std::to_string(ret));
    }
  }
}


} /* namespace geosx */
