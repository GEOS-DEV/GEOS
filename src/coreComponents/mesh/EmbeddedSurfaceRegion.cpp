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

/**
 * @file EmbeddedSurfaceRegion.cpp
 */

#include "EmbeddedSurfaceRegion.hpp"

#include "EdgeManager.hpp"

namespace geosx
{
using namespace dataRepository;

EmbeddedSurfaceRegion::EmbeddedSurfaceRegion( string const & name, Group * const parent ):
  ElementRegionBase( name, parent )
{
  this->GetGroup(viewKeyStruct::elementSubRegions)->RegisterGroup<EmbeddedSurfaceSubRegion>("default");
}

EmbeddedSurfaceRegion::~EmbeddedSurfaceRegion()
{}

/*

localIndex EmbeddedSurfaceRegion::AddToFractureMesh( EdgeManager * const edgeManager,
                                                 FaceManager const * const faceManager,
                                                 array1d< array1d<localIndex> > const & originalFaceToEdges,
                                                 string const & subRegionName,
                                                 localIndex const faceIndices[2]  )
{
  localIndex rval = -1;
  // Matteo: has to be filled in
  return rval;
}
*/

REGISTER_CATALOG_ENTRY( ObjectManagerBase, EmbeddedSurfaceRegion, std::string const &, Group * const )

} /* namespace geosx */
