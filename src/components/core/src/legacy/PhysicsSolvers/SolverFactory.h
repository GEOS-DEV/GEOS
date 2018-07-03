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
 * @file SolverFactory.h
 * @author walsh24
 * @date Feb 28, 2011
 */

#ifndef SOLVERFACTORY_H_
#define SOLVERFACTORY_H_

#include "SolverBase.h"
#include "Utilities/StringUtilities.h"

#include <map>
#include <string>
#include <vector>

//////////////////////////

// Solver Factory
//
// Consists of the following parts:
//   * The function to generate new solver pointers: "newSolver"
//   * A base class to derive the functions to generate solver pointers:
// "SolverInitializer"
//   * A String-to-Solver-Intializer map hidden behind the getSolverCatalogue
// function
//   * A template to create solver initializers: "SolverRegistrator"
//   * A compiler directive to simplify autoregistration: "REGISTER_SOLVER"
//
// Most solvers will only need to use one or two of the parts:
//   * To register a new solver in the factory: REGISTER_SOLVER( SolverClassName
// )
//   * To load a solver pointer from the factory:       SolverBase* aSolverPtr =
// SolverFactory::NewSolver(solverString, args );

/// Base class to generate new Solver pointers
class SolverInitializer
{
public:
  virtual SolverBase* InitializeSolver(TICPP::HierarchicalDataNode* const hdn, ProblemManagerT* const pm) = 0;

  virtual ~SolverInitializer()
  {}
};

typedef std::map<std::string, SolverInitializer*> SolverCatalogueType;

class SolverFactory
{
public:
  /// The Solver Factory.
  static SolverBase* NewSolver(const std::string& solverName,
                               TICPP::HierarchicalDataNode* const hdn,
                               ProblemManagerT* const pm);

  /// Interface to the Solver name -> Solver initializer map
  static SolverCatalogueType& GetSolverCatalogue();

  /// Return a list of supported solver names
  static void GetSolverNames(std::vector<std::string>& nameList);
};

/// Template for creating classes derived from SolverInitializer
template<class SolverType>
class SolverRegistrator : public SolverInitializer
{
public:
  SolverRegistrator(void)
  {
    const std::string solverName(SolverType::SolverName());
    SolverFactory::GetSolverCatalogue()[solverName] = this;
  }

  SolverBase* InitializeSolver(TICPP::HierarchicalDataNode* const hdn, ProblemManagerT* const pm)
  {
    if(!hdn)
      throw GPException("Need to specify a valid HierarchicalDataNode to InitializeSolver");

    const std::string name(hdn->GetAttributeString("name"));

    if(name.empty())
      throw GPException("Need to specify a valid name to InitializeSolver");

    SolverBase* ret = new SolverType(name, pm);
    ret->ReadXML(hdn);
    return ret;
  }
};

/// Compiler directive to simplify autoregistration
#define REGISTER_SOLVER( ClassName ) namespace { SolverRegistrator<ClassName> reg_ ## ClassName; }

#endif
