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
 * @file ParallelPlateFluidModelBase.cpp
 * @author walsh24
 * @date Feb 28, 2011
 */

#include "ParallelPlateFluidModelBase.h"


////////////////////////////////////////////////
// Parallel Plate Fluid Model Base Class

/* empty */


////////////////////////////////////////////////
// Parallel Plate Fluid Model Factory

ParallelPlateFluidModelCatalogueType& ParallelPlateFluidModelFactory::GetParallelPlateFluidModelCatalogue()
{
  static ParallelPlateFluidModelCatalogueType theCatalogue;
  return theCatalogue;
}

void ParallelPlateFluidModelFactory::GetParallelPlateFluidModelNames(std::vector<std::string>& nameList)
{
  for (ParallelPlateFluidModelCatalogueType::const_iterator it = GetParallelPlateFluidModelCatalogue().begin() ;
       it != GetParallelPlateFluidModelCatalogue().end() ; ++it)
  {
    nameList.push_back(it->first);
  }
  ;
}

ParallelPlateFluidModelBase* ParallelPlateFluidModelFactory::NewParallelPlateFluidModel(const std::string& parallelPlateFluidModelName,
                                                                                        TICPP::HierarchicalDataNode* const hdn)
{

  ParallelPlateFluidModelInitializer* parallelPlateFluidModelInitializer = GetParallelPlateFluidModelCatalogue()[parallelPlateFluidModelName];
  ParallelPlateFluidModelBase *theNewParallelPlateFluidModel = NULL;

  if (!parallelPlateFluidModelInitializer)
  {
    std::string msg = "ERROR: Could not create unrecognized parallelPlateFluidModel ";
    msg = msg + parallelPlateFluidModelName;
    throw GPException(msg);
  }
  else
  {
    theNewParallelPlateFluidModel = parallelPlateFluidModelInitializer->InitializeParallelPlateFluidModel(hdn);
  }

  return theNewParallelPlateFluidModel;
}
