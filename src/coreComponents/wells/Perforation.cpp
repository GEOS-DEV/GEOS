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
    m_transmissibility(0),
    m_wellElementName("")
{
  RegisterViewWrapper( viewKeyStruct::locationString, &m_location, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Perforation physical coordinates");

  RegisterViewWrapper( viewKeyStruct::transmissibilityString, &m_transmissibility, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Perforation transmissibility");

  RegisterViewWrapper( viewKeyStruct::wellElementNameString, &m_wellElementName, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Segment name");
}

void Perforation::PostProcessInput()
{  
  GEOS_ERROR_IF( m_wellElementName.empty(), 
                 "Invalid segment name in segment " << getName() );
  
  GEOS_ERROR_IF( m_transmissibility <= 0, 
                 "Negative transmissibility value in segment " << getName() );
}


Perforation::~Perforation()
{

}


} //namespace geosx
