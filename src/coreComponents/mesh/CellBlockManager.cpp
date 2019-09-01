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

/*
 * ElementManagerT.cpp
 *
 *  Created on: Sep 14, 2010
 *      Author: settgast1
 */

#include "CellBlockManager.hpp"

#include "FaceManager.hpp"
//#include "legacy/IO/BinStream.h"
#include <map>
#include <vector>
//#include "legacy/Constitutive/Material/MaterialFactory.h"
//#include "legacy/ArrayT/ArrayT.h"

namespace geosx
{
using namespace dataRepository;

CellBlockManager::CellBlockManager(  string const & name, Group * const parent ):
  ObjectManagerBase(name,parent)
{
  this->RegisterGroup<Group>(keys::cellBlocks);
}

CellBlockManager::~CellBlockManager()
{
  // TODO Auto-generated destructor stub
}

void CellBlockManager::resize( integer_array const & numElements,
                               string_array const & regionNames,
                               string_array const & elementTypes )
{
  localIndex const numRegions = integer_conversion<localIndex>(regionNames.size());
//  ManagedGroup * elementRegions = this->GetGroup(keys::cellBlocks);
  for( localIndex reg=0 ; reg<numRegions ; ++reg )
  {
    CellBlock * elemRegion = this->GetRegion( regionNames[reg] );
    elemRegion->resize(numElements[reg]);
  }
}


//CellBlock & CellBlockManager::CreateRegion( string const & regionName,
//                                             string const & elementType,
//                                             integer const & numElements )
//{
////  ElementRegion * elemRegion = elementRegions.RegisterGroup( regionNames );
////  elemRegion->resize(numElements);
//}

Group * CellBlockManager::CreateChild( string const & childKey, string const & childName )
{
  return nullptr;
}



REGISTER_CATALOG_ENTRY( ObjectManagerBase, CellBlockManager, string const &, Group * const )
}
