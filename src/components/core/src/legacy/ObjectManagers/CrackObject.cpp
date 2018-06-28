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

/*
 * CrackSurface.cpp
 *
 *  Created on: Nov 4, 2014
 *      Author: annavarapusr1, rrsettgast
 */

#include "CrackObject.h"



CrackObjectManager::CrackObjectManager():
  ObjectDataStructureBaseT( ObjectDataStructureBaseT::CrackSurfaceManager ),
  m_crackSurface1ToCutFaces(m_VariableOneToManyMaps["crackSurface1ToCutFaces"]),
  m_crackSurface2ToCutFaces(m_VariableOneToManyMaps["crackSurface2ToCutFaces"]),
  m_crackToElements(m_VariableOneToManyMaps["crackToElements"])
{
//  // For now use this definition to define a crack
//  this->AddKeylessDataField<R1Tensor>("crackPositionStart");
//  this->AddKeylessDataField<R1Tensor>("crackPositionEnd");
//  m_crackFaceManager.resize(1);

}

CrackObjectManager::~CrackObjectManager()
{
  // TODO Auto-generated destructor stub
}
