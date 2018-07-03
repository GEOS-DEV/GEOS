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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file ElementRegion.cpp
 * @author settgast1
 * @date Dec 2, 2010
 */

#include "ElementRegion.h"

ElementRegion::ElementRegion():
  ObjectManager(),
  m_numElements(ObjectManager::m_numObjects),
  xi(),
  m_dNdXi(),
  m_dNdX(),
  detJ()
{
  // TODO Auto-generated constructor stub

}

ElementRegion::~ElementRegion()
{
  // TODO Auto-generated destructor stub
}
