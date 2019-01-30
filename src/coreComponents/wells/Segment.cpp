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
 * @file Segment.cpp
 *
 */

#include "Segment.hpp"

#include "dataRepository/InputFlags.hpp"

namespace geosx
{

using namespace dataRepository;

Segment::Segment(string const & name, ManagedGroup * const parent)
  : ManagedGroup(name, parent)
{
  RegisterViewWrapper( viewKeysSegment.segmentName.Key(), &m_segmentName, false )->
    setDefaultValue("0")->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Well segment name (can be omitted for single-segment wells");
}

Segment::~Segment()
{

}


} //namespace geosx
