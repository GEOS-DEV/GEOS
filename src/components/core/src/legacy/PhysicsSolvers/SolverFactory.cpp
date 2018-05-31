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
 * @file SolverFactory.cpp
 * @author walsh24
 * @date Feb 28, 2011
 */

#include "SolverFactory.h"

SolverCatalogueType& SolverFactory::GetSolverCatalogue()
{
  static SolverCatalogueType theCatalogue;
  return theCatalogue;
}

void SolverFactory::GetSolverNames(std::vector<std::string>& nameList)
{
  for (SolverCatalogueType::const_iterator it = GetSolverCatalogue().begin() ;
       it != GetSolverCatalogue().end() ; ++it)
  {
    nameList.push_back(it->first);
  }
  ;
}

SolverBase* SolverFactory::NewSolver(const std::string& solverName,
                                     TICPP::HierarchicalDataNode* const hdn,
                                     ProblemManagerT* const pm)
{

  SolverInitializer* solverInitializer = GetSolverCatalogue()[solverName];
  SolverBase *theNewSolver = NULL;

  if (!solverInitializer)
  {
    std::string msg = "ERROR: Could not create unrecognized solver ";
    msg = msg + solverName;
    throw GPException(msg);
  }
  else
  {
    theNewSolver = solverInitializer->InitializeSolver(hdn, pm);
  }

  return theNewSolver;
}
