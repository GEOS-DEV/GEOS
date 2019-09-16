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

#include  "WellGeneratorBase.hpp"

#include "PerforationData.hpp"
#include "Perforation.hpp"

namespace geosx
{

using namespace dataRepository;

WellGeneratorBase::WellGeneratorBase( string const & name, Group * const parent ):
  MeshGeneratorBase( name, parent ),
  m_numElemsPerSegment(0),
  m_crossSectionArea(0),
  m_wellRegionName(""),
  m_wellControlsName(""),
  m_meshBodyName(""),
  m_numElems(0),
  m_numNodesPerElem(2),
  m_numNodes(0),
  m_numPerforations(0),
  m_nDims(3),
  m_polylineHeadNodeId(-1)
{
}

WellGeneratorBase::~WellGeneratorBase()
{
}

void WellGeneratorBase::PostProcessInput()
{
}

Group * WellGeneratorBase::CreateChild( string const & childKey, string const & childName )
{
  if ( childKey == keys::perforation )
  {
    ++m_numPerforations;

    // keep track of the perforations that have been added
    m_perforationList.push_back( childName );

    return RegisterGroup<Perforation>( childName );
  }
  else
  {
    GEOS_ERROR( "Unrecognized node: " << childKey );
  }
  return nullptr;
}

} // namespace
