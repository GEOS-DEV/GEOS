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
 * @file CohesiveZoneFactory.h
 * @author Scott Johnson
 * @date Oct 6, 2013
 */

#ifndef COHESIVEZONEFACTORY_H_
#define COHESIVEZONEFACTORY_H_

#include "CohesiveZoneBase.h"
#include "Utilities/StringUtilities.h"

#include <map>
#include <string>
#include <vector>

//////////////////////////

// CohesiveZone Factory
//
// Consists of the following parts:
//   * The function to generate new cohesive zone pointers: "NewCohesiveZone"
//   * A base class to derive the functions to generate cohesive zone pointers:
// "CohesiveZoneInitializer"
//   * A String-to-CohesiveZone-Intializer map hidden behind the
// GetCohesiveZoneCatalogue function
//   * A template to create cohesive zone initializers:
// "CohesiveZoneRegistrator"
//   * A compiler directive to simplify autoregistration:
// "REGISTER_COHESIVEZONE"
//
// Most cohesive zones will only need to use one or two of the parts:
//   * To register a new cohesive zone in the factory: REGISTER_COHESIVEZONE(
// CohesiveZoneClassName )
//   * To load a cohesive zone pointer from the factory:       CohesiveZoneBase*
// aCohesiveZonePtr = CohesiveZoneFactory::NewCohesiveZone(cohesiveZoneString,
// args );

/// Base class to generate new CohesiveZone pointers
class CohesiveZoneInitializer
{
public:
  virtual
#if USECPP11==1
  std::unique_ptr<CohesiveZoneBase>
#else
  CohesiveZoneBase*
#endif
  InitializeCohesiveZone(TICPP::HierarchicalDataNode* hdn) = 0;

  virtual ~CohesiveZoneInitializer()
  {}
};

typedef std::map<std::string, CohesiveZoneInitializer*> CohesiveZoneCatalogueType;

class CohesiveZoneFactory
{
public:
  /// The CohesiveZone Factory.
  static
#if USECPP11==1
  std::unique_ptr<CohesiveZoneBase>
#else
  CohesiveZoneBase*
#endif
  NewCohesiveZone(const std::string& cohesiveZoneName,
                  TICPP::HierarchicalDataNode* hdn = 0);

  /// CohesiveZone to the CohesiveZone name -> CohesiveZone initializer map
  static CohesiveZoneCatalogueType& GetCohesiveZoneCatalogue();

  /// Return a list of supported cohesive zone names
  static void GetCohesiveZoneNames(std::vector<std::string>& nameList);
};

/// Template for creating classes derived from CohesiveZoneInitializer
template<class CohesiveZoneType>
class CohesiveZoneRegistrator : public CohesiveZoneInitializer
{
public:
  CohesiveZoneRegistrator(void)
  {
    CohesiveZoneFactory::GetCohesiveZoneCatalogue()[CohesiveZoneType::Name()] = this;
  }

#if USECPP11==1
  std::unique_ptr<CohesiveZoneBase>
#else
  CohesiveZoneBase*
#endif
  InitializeCohesiveZone(TICPP::HierarchicalDataNode* hdn)
  {

#if USECPP11==1
    std::unique_ptr<CohesiveZoneBase> tmp(new CohesiveZoneType());
#else
    CohesiveZoneBase* tmp = new CohesiveZoneType();
#endif

    if (hdn)
    {
      tmp->resize(0, 1);
      tmp->ReadXML(*hdn);
    }
    return tmp;
  }
};

/// Compiler directive to simplify autoregistration
#define REGISTER_COHESIVEZONE( ClassName ) namespace { CohesiveZoneRegistrator<ClassName> reg_ ## ClassName; }

#endif
