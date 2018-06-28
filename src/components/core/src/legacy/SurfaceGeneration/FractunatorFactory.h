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
 * @file FractunatorFactory.h
 * @author Scott Johnson
 * @date Oct 26, 2013
 */

#ifndef FRACTUNATORFACTORY_H_
#define FRACTUNATORFACTORY_H_

#include "FractunatorBase.h"
#include "Utilities/StringUtilities.h"

#include <map>
#include <string>
#include <vector>

//////////////////////////

// Fractunator Factory
//
// Consists of the following parts:
//   * The function to generate new fractunator pointers: "NewFractunator"
//   * A base class to derive the functions to generate fractunator pointers:
// "FractunatorInitializer"
//   * A String-to-Fractunator-Intializer map hidden behind the
// GetFractunatorCatalogue function
//   * A template to create fractunator initializers: "FractunatorRegistrator"
//   * A compiler directive to simplify autoregistration: "REGISTER_FRACTUNATOR"
//
// Most fractunators will only need to use one or two of the parts:
//   * To register a new fractunator in the factory: REGISTER_FRACTUNATOR(
// FractunatorClassName )
//   * To load a fractunator pointer from the factory:       FractunatorBase*
// aFractunatorPtr = FractunatorFactory::NewFractunator(fractunatorString, args
// );

/// Base class to generate new Fractunator pointers
class FractunatorInitializer
{
public:
  virtual FractunatorBase* InitializeFractunator(TICPP::HierarchicalDataNode* hdn) = 0;

  virtual ~FractunatorInitializer()
  {}
};

typedef std::map<std::string, FractunatorInitializer*> FractunatorCatalogueType;

class FractunatorFactory
{
public:
  /// The Fractunator Factory.
  static FractunatorBase* NewFractunator(const std::string& fractunatorName,
                                         TICPP::HierarchicalDataNode* hdn = 0);

  /// Fractunator to the Fractunator name -> Fractunator initializer map
  static FractunatorCatalogueType& GetFractunatorCatalogue();

  /// Return a list of supported fractunator names
  static void GetFractunatorNames(std::vector<std::string>& nameList);
};

/// Template for creating classes derived from FractunatorInitializer
template<class FractunatorType>
class FractunatorRegistrator : public FractunatorInitializer
{
public:
  FractunatorRegistrator(void)
  {
    FractunatorFactory::GetFractunatorCatalogue()[FractunatorType::FractunatorName()] = this;
  }

  FractunatorBase* InitializeFractunator(TICPP::HierarchicalDataNode* hdn)
  {
    FractunatorBase* tmp = new FractunatorType();
    if(hdn)
      tmp->ReadXML(*hdn);
    return tmp;
  }
};

/// Compiler directive to simplify autoregistration
#define REGISTER_FRACTUNATOR( ClassName ) namespace { FractunatorRegistrator<ClassName> reg_ ## ClassName; }

#endif
