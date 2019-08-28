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
    m_distanceFromHead(0),
    m_transmissibility(0)
{
  RegisterViewWrapper( viewKeyStruct::distanceFromHeadString, &m_distanceFromHead, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Perforation linear distance from well head");

  RegisterViewWrapper( viewKeyStruct::transmissibilityString, &m_transmissibility, false )->
    setApplyDefaultValue(-1.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Perforation transmissibility");

}

void Perforation::PostProcessInput()
{
  GEOS_ERROR_IF( m_distanceFromHead <= 0,
                 "Invalid distance well head to perforation " << getName() );

}

Perforation::~Perforation()
{
}


} //namespace geosx
