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

#include "dataRepository/InputFlags.hpp"

namespace geosx
{

using namespace dataRepository;

Perforation::Perforation(string const & name, ManagedGroup * const parent)
  : ManagedGroup(name, parent),
    m_location(),
    m_transmissibility()
{
  RegisterViewWrapper( viewKeysPerforation.location.Key(), &m_location, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Perforation physical coordinates");

  RegisterViewWrapper( viewKeysPerforation.transmissibility.Key(), &m_transmissibility, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Connection transmissibility");

  RegisterViewWrapper( viewKeysPerforation.segmentName.Key(), &m_segmentName, false )->
    setDefaultValue("0")->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Well segment name (can be omitted for single-segment wells");
}

Perforation::~Perforation()
{

}


} //namespace geosx
