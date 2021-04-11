/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file RestartOutput.cpp
 */

#include "RestartOutput.hpp"
#include "fileIO/silo/SiloFile.hpp"

namespace geosx
{

using namespace dataRepository;

RestartOutput::RestartOutput( string const & name,
                              Group * const parent ):
  OutputBase( name, parent )
{}

RestartOutput::~RestartOutput()
{}

bool RestartOutput::execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                             real64 const GEOSX_UNUSED_PARAM( dt ),
                             integer const cycleNumber,
                             integer const GEOSX_UNUSED_PARAM( eventCounter ),
                             real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                             DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  GEOSX_MARK_FUNCTION;

  Group & rootGroup = this->getGroupByPath( "/Problem" );

  // Ignoring the eventProgress indicator for now to be compliant with the integrated test repo
  // integer const eventProgressPercent = static_cast<integer const>(eventProgress * 100.0);
  char fileName[200] = {0};
  sprintf( fileName, "%s_%s_%09d", getFileNameRoot().c_str(), "restart", cycleNumber );

  rootGroup.prepareToWrite();
  writeTree( joinPath( OutputBase::getOutputDirectory(), fileName ), *(rootGroup.getConduitNode().parent()) );
  rootGroup.finishWriting();

  return false;
}


REGISTER_CATALOG_ENTRY( OutputBase, RestartOutput, string const &, Group * const )
} /* namespace geosx */
