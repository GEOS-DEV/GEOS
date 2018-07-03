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
 * @file SpatialSorterFactory.h
 * @author Scott Johnson
 * @date Oct 6, 2013
 */

#ifndef SPATIALSORTERFACTORY_H_
#define SPATIALSORTERFACTORY_H_

#include "SpatialSorterBase.h"
#include "Utilities/StringUtilities.h"

#include <map>
#include <string>
#include <vector>

//////////////////////////

// SpatialSorter Factory
//
// Consists of the following parts:
//   * The function to generate new spatialSorter pointers: "NewSpatialSorter"
//   * A base class to derive the functions to generate spatialSorter pointers:
// "SpatialSorterInitializer"
//   * A String-to-SpatialSorter-Intializer map hidden behind the
// GetSpatialSorterCatalogue function
//   * A template to create spatialSorter initializers:
// "SpatialSorterRegistrator"
//   * A compiler directive to simplify autoregistration:
// "REGISTER_SPATIALSORTER"
//
// Most spatialSorters will only need to use one or two of the parts:
//   * To register a new spatialSorter in the factory: REGISTER_SPATIALSORTER(
// SpatialSorterClassName )
//   * To load a spatialSorter pointer from the factory:
//       SpatialSorterBase* aSpatialSorterPtr =
// SpatialSorterFactory::NewSpatialSorter(spatialSorterString, args );

/// Base class to generate new SpatialSorter pointers

namespace SpatialSorting
{
class SpatialSorterInitializer
{
public:
  virtual SpatialSorterBase* InitializeSpatialSorter() = 0;

  virtual ~SpatialSorterInitializer()
  {}
};

typedef std::map<std::string, SpatialSorterInitializer*> SpatialSorterCatalogueType;

class SpatialSorterFactory
{
public:
  /// The SpatialSorter Factory.
  static SpatialSorterBase* NewSpatialSorter(const std::string& spatialSorterName);

  /// Interface to the SpatialSorter name -> SpatialSorter initializer map
  static SpatialSorterCatalogueType& GetSpatialSorterCatalogue();

  /// Return a list of supported spatialSorter names
  static void GetSpatialSorterNames(std::vector<std::string>& nameList);
};

/// Template for creating classes derived from SpatialSorterInitializer
template<class SpatialSorterType>
class SpatialSorterRegistrator : public SpatialSorterInitializer
{
public:
  SpatialSorterRegistrator(void)
  {
    SpatialSorterFactory::GetSpatialSorterCatalogue()[SpatialSorterType::SpatialSorterName()] = this;
  }

  SpatialSorterBase* InitializeSpatialSorter()
  {
    return new SpatialSorterType();
  }
};
}

/// Compiler directive to simplify autoregistration
#define REGISTER_SPATIALSORTER( ClassName ) namespace SpatialSorting { SpatialSorterRegistrator<ClassName> reg_ ## ClassName; }

#endif
