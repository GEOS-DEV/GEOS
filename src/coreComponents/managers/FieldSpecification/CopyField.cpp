/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "CopyField.hpp"
#include "managers/DomainPartition.hpp"
#include "FieldSpecificationBase.hpp"

namespace geosx
{

using namespace dataRepository;

CopyField::CopyField( std::string const & name,
                      Group * const parent ):
  TaskBase( name, parent )
{
  registerWrapper( viewKeyStruct::fromString, &m_from )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the field to copy." );

  registerWrapper( viewKeyStruct::toString, &m_to )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the field that will be created and filled." );

  registerWrapper( viewKeyStruct::targetRegionsString, &m_targetRegions )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Regions on which the copy will be done." );
}

CopyField::~CopyField()
{}

void CopyField::Execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                         real64 const GEOSX_UNUSED_PARAM( dt ),
                         integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                         integer const GEOSX_UNUSED_PARAM( eventCounter ),
                         real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                         Group * GEOSX_UNUSED_PARAM( domain ) )
{
  GEOSX_ERROR("Not implemented");
}

REGISTER_CATALOG_ENTRY( TaskBase, CopyField, std::string const &, Group * const )

} /* namespace */
