/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

/*
 * NewComponent.cpp
 *
 *  Created on: Jun 8, 2016
 *      Author: settgast
 */

#include "NewComponent.hpp"

namespace geosx
{

NewComponent::NewComponent( std::string const & name,
                            ManagedGroup * const parent ):
    SolverBase(name,parent)
{

}

NewComponent::~NewComponent()
{}



void NewComponent::InitializePreSubGroups( ManagedGroup * const )
{}

void NewComponent::InitializePostSubGroups( ManagedGroup * const )
{}



void NewComponent::ReadXML_PostProcess( )
{}




real64 NewComponent::SolverStep( real64 const & /*time_n*/,
               real64 const & /*dt*/,
               integer const /*cycleNumber*/,
               DomainPartition * /*domain*/ )
{
  return 0;
}

REGISTER_CATALOG_ENTRY( SolverBase, NewComponent, std::string const &, dataRepository::ManagedGroup * const )

} /* namespace geosx */
