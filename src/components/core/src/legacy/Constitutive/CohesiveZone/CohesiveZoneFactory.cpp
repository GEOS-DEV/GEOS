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
 * @file CohesiveZoneFactory.cpp
 * @author Scott Johnson
 * @date Oct 6, 2013
 */

#include "CohesiveZoneFactory.h"

CohesiveZoneCatalogueType& CohesiveZoneFactory::GetCohesiveZoneCatalogue()
{
  static CohesiveZoneCatalogueType theCatalogue;
  return theCatalogue;
}

void CohesiveZoneFactory::GetCohesiveZoneNames(std::vector<std::string>& nameList)
{
  for (CohesiveZoneCatalogueType::const_iterator it = GetCohesiveZoneCatalogue().begin() ;
       it != GetCohesiveZoneCatalogue().end() ; ++it)
  {
    nameList.push_back(it->first);
  }
  ;
}

#if USECPP11==1
std::unique_ptr<CohesiveZoneBase>
#else
CohesiveZoneBase*
#endif
CohesiveZoneFactory::NewCohesiveZone(const std::string& cohesiveZoneName,
                                     TICPP::HierarchicalDataNode* hdn)
{

  CohesiveZoneInitializer* cohesiveZoneInitializer = GetCohesiveZoneCatalogue()[cohesiveZoneName];
  CohesiveZoneBase *theNewCohesiveZone = NULL;

  if (!cohesiveZoneInitializer)
  {
    std::string msg = "ERROR: Could not create unrecognized cohesive zone ";
    msg = msg + cohesiveZoneName;
    throw GPException(msg);
  }
  else
  {
    theNewCohesiveZone = cohesiveZoneInitializer->InitializeCohesiveZone(hdn);
  }

  return theNewCohesiveZone;
}
