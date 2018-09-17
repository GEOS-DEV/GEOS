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
 * @file Perforation.cpp
 *
 */

#include "Perforation.hpp"

namespace geosx
{


Perforation::Perforation(string const & name, dataRepository::ManagedGroup * const parent)
  : ManagedGroup(name, parent),
    m_location(),
    m_transmissibility()
{
  RegisterViewWrapper( viewKeys.location.Key(), &m_location, false );
  RegisterViewWrapper( viewKeys.transmissibility.Key(), &m_transmissibility, false );
  RegisterViewWrapper( viewKeys.segmentName.Key(), &m_segmentName, false );
}

Perforation::~Perforation()
{

}

void Perforation::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->AllocateChildNode( viewKeys.location.Key(),
                              viewKeys.location.Key(),
                              -1,
                              "R1Tensor",
                              "R1Tensor",
                              "Perforation physical coordinates",
                              "Perforation physical coordinates",
                              "REQUIRED",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.transmissibility.Key(),
                              viewKeys.transmissibility.Key(),
                              -1,
                              "real64",
                              "real64",
                              "Connection transmissibility",
                              "Connection transmissibility",
                              "-1",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.segmentName.Key(),
                              viewKeys.segmentName.Key(),
                              -1,
                              "string",
                              "string",
                              "Well segment name",
                              "Well segment name",
                              "0",
                              "",
                              0,
                              1,
                              0 );
}


} //namespace geosx