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
 * @file FractunatorFactory.cpp
 * @author Scott Johnson
 * @date Oct 26, 2013
 */

#include "FractunatorFactory.h"

FractunatorCatalogueType& FractunatorFactory::GetFractunatorCatalogue()
{
  static FractunatorCatalogueType theCatalogue;
  return theCatalogue;
}

void FractunatorFactory::GetFractunatorNames(std::vector<std::string>& nameList)
{
  for (FractunatorCatalogueType::const_iterator it = GetFractunatorCatalogue().begin() ;
       it != GetFractunatorCatalogue().end() ; ++it)
  {
    nameList.push_back(it->first);
  }
  ;
}

FractunatorBase* FractunatorFactory::NewFractunator(const std::string& fractunatorName,
                                                    TICPP::HierarchicalDataNode* hdn)
{

  FractunatorInitializer* fractunatorInitializer = GetFractunatorCatalogue()[fractunatorName];
  FractunatorBase *theNewFractunator = NULL;

  if (!fractunatorInitializer)
  {
    std::string msg = "ERROR: Could not create unrecognized fractunator ";
    msg = msg + fractunatorName;
    throw GPException(msg);
  }
  else
  {
    theNewFractunator = fractunatorInitializer->InitializeFractunator(hdn);
  }

  return theNewFractunator;
}
