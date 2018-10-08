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
 * GEOSXMeshBody.hpp
 *
 *  Created on: Oct 10, 2018
 *      Author: Antoine Mazuyer
 */

#include "GEOSXMeshBody.hpp"

namespace geosx {
using namespace dataRepository;
GEOSXMeshBody::GEOSXMeshBody( const std::string& name,
                              ManagedGroup * const parent ) :
    MeshBody(name, parent)
{
}

REGISTER_CATALOG_ENTRY(MeshBody, GEOSXMeshBody, std::string const &, ManagedGroup * const)
} // namespace
