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
 * @file MaterialFactory.cpp
 * @author Scott Johnson
 * @date Oct 6, 2013
 */

#include "MaterialFactory.h"

MaterialCatalogueType& MaterialFactory::GetMaterialCatalogue()
{
  static MaterialCatalogueType theCatalogue;
  return theCatalogue;
}

void MaterialFactory::GetMaterialNames(std::vector<std::string>& nameList)
{
  for (MaterialCatalogueType::const_iterator it = GetMaterialCatalogue().begin() ;
       it != GetMaterialCatalogue().end() ; ++it)
  {
    nameList.push_back(it->first);
  }
  ;
}

#if USECPP11==1
std::unique_ptr<MaterialBase>
#else
MaterialBase*
#endif
MaterialFactory::NewMaterial(const std::string& materialName,
                             TICPP::HierarchicalDataNode* hdn)
{
  unsigned pos = materialName.find("Material"); // position of "Material" in
                                                // materialName
  const std::string str = materialName.substr(0, pos); // get up to "Material"
                                                       // // SW- why not use
                                                       // full name?

  MaterialInitializer* materialInitializer = GetMaterialCatalogue()[str];
  MaterialBase *theNewMaterial = NULL;

  if (!materialInitializer)
  {
    std::string msg = "ERROR: Could not create unrecognized material ";
    msg = msg + str;
    throw GPException(msg);
  }
  else
  {
    theNewMaterial = materialInitializer->InitializeMaterial(hdn);
  }

  return theNewMaterial;


//#if USECPP11==1
//  return
// std::unique_ptr<MaterialBase>(materialInitializer->InitializeMaterial(hdn));
//#else
//  return materialInitializer->InitializeMaterial(hdn);
//#endif
}
