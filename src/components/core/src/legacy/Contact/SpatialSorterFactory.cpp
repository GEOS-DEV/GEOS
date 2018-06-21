// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
/**
 * @file SpatialSorterFactory.cpp
 * @author Scott Johnson
 * @date Oct 6, 2013
 */

#include "SpatialSorterFactory.h"

namespace SpatialSorting
{
SpatialSorterCatalogueType& SpatialSorterFactory::GetSpatialSorterCatalogue()
{
  static SpatialSorterCatalogueType theCatalogue;
  return theCatalogue;
}

void SpatialSorterFactory::GetSpatialSorterNames(std::vector<std::string>& nameList)
{
  for (SpatialSorterCatalogueType::const_iterator it = GetSpatialSorterCatalogue().begin() ;
       it != GetSpatialSorterCatalogue().end() ; ++it)
  {
    nameList.push_back(it->first);
  }
  ;
}

SpatialSorterBase* SpatialSorterFactory::NewSpatialSorter(const std::string& spatialSorterName)
{

  SpatialSorterInitializer* spatialSorterInitializer = GetSpatialSorterCatalogue()[spatialSorterName];
  SpatialSorterBase *theNewSpatialSorter = NULL;

  if (!spatialSorterInitializer)
  {
    std::string msg = "ERROR: Could not create unrecognized spatialSorter ";
    msg = msg + spatialSorterName;
    throw GPException(msg);
  }
  else
  {
    theNewSpatialSorter = spatialSorterInitializer->InitializeSpatialSorter();
  }

  return theNewSpatialSorter;
}
}
