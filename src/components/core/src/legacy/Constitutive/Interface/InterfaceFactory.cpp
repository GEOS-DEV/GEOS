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
 * @file InterfaceFactory.cpp
 * @author Scott Johnson
 * @date Oct 6, 2013
 */

#include "InterfaceFactory.h"

InterfaceCatalogueType& InterfaceFactory::GetInterfaceCatalogue()
{
  static InterfaceCatalogueType theCatalogue;
  return theCatalogue;
}

void InterfaceFactory::GetInterfaceNames(std::vector<std::string>& nameList)
{
  for (InterfaceCatalogueType::const_iterator it = GetInterfaceCatalogue().begin() ;
       it != GetInterfaceCatalogue().end() ; ++it)
  {
    nameList.push_back(it->first);
  }
  ;
}

InterfaceBase* InterfaceFactory::NewInterface(const std::string& interfaceName,
                                              TICPP::HierarchicalDataNode* hdn)
{

  InterfaceInitializer* interfaceInitializer = GetInterfaceCatalogue()[interfaceName];
  InterfaceBase *theNewInterface = NULL;

  if (!interfaceInitializer)
  {
    std::string msg = "ERROR: Could not create unrecognized interface ";
    msg = msg + interfaceName;
    throw GPException(msg);
  }
  else
  {
    theNewInterface = interfaceInitializer->InitializeInterface(hdn);
  }

  return theNewInterface;
}
