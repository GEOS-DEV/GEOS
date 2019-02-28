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
 * @file Connection.cpp
 *
 */

#include "Connection.hpp"

#include "dataRepository/InputFlags.hpp"

namespace geosx
{

using namespace dataRepository;

Connection::Connection(string const & name, ManagedGroup * const parent)
  : ManagedGroup(name, parent),
    m_nextWellElementName(""),
    m_prevWellElementName("")
{
  RegisterViewWrapper( viewKeyStruct::nextWellElementNameString, &m_nextWellElementName, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Next well element name");

  RegisterViewWrapper( viewKeyStruct::prevWellElementNameString, &m_prevWellElementName, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Previous well element name");

  RegisterViewWrapper( viewKeyStruct::areaString, &m_area, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Area of the cross section");

}

Connection::~Connection()
{

}


} //namespace geosx
