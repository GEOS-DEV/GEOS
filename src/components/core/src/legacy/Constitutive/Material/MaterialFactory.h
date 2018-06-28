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
 * @file MaterialFactory.h
 * @author Scott Johnson
 * @date Oct 6, 2013
 */

#ifndef MATERIALFACTORY_H_
#define MATERIALFACTORY_H_

#include "MaterialBase.h"
#include "legacy/Utilities/StringUtilities.h"

#include <map>
#include <string>
#include <vector>

//////////////////////////

// Material Factory
//
// Consists of the following parts:
//   * The function to generate new material pointers: "NewMaterial"
//   * A base class to derive the functions to generate material pointers:
// "MaterialInitializer"
//   * A String-to-Material-Intializer map hidden behind the
// GetMaterialCatalogue function
//   * A template to create material initializers: "MaterialRegistrator"
//   * A compiler directive to simplify autoregistration: "REGISTER_MATERIAL"
//
// Most materials will only need to use one or two of the parts:
//   * To register a new material in the factory: REGISTER_MATERIAL(
// MaterialClassName )
//   * To load a material pointer from the factory:       MaterialBase*
// aMaterialPtr = MaterialFactory::NewMaterial(materialString, args );

/// Base class to generate new Material pointers
class MaterialInitializer
{
public:
  virtual
#if USECPP11==1
  std::unique_ptr<MaterialBase>
#else
  MaterialBase*
#endif
  InitializeMaterial(TICPP::HierarchicalDataNode* hdn) = 0;

  virtual ~MaterialInitializer()
  {}
};

typedef std::map<std::string, MaterialInitializer*> MaterialCatalogueType;

class MaterialFactory
{
public:
  /// The Material Factory.
  static
#if USECPP11==1
  std::unique_ptr<MaterialBase>
#else
  MaterialBase*
#endif
  NewMaterial(const std::string& materialName,
              TICPP::HierarchicalDataNode* hdn = 0);

  /// Interface to the Material name -> Material initializer map
  static MaterialCatalogueType& GetMaterialCatalogue();

  /// Return a list of supported material names
  static void GetMaterialNames(std::vector<std::string>& nameList);
};

/// Template for creating classes derived from MaterialInitializer
/// But note that not all materials need to be registered in this way (eg.
// Geodyn material).
template<class MaterialType>
class MaterialRegistrator : public MaterialInitializer
{
public:
  MaterialRegistrator(void)
  {
    MaterialFactory::GetMaterialCatalogue()[MaterialType::Name()] = this;
  }

#if USECPP11==1
  std::unique_ptr<MaterialBase>
#else
  MaterialBase*
#endif
  InitializeMaterial(TICPP::HierarchicalDataNode* hdn)
  {

#if USECPP11==1
    std::unique_ptr<MaterialBase> tmp = std::unique_ptr<MaterialBase>(new MaterialType());
#else
    MaterialBase* tmp = new MaterialType();
#endif

    if(hdn)
    {
      tmp->resize(0,1);
      tmp->ReadXML(*hdn);
    }
    return tmp;
  }
};

/// Compiler directive to simplify autoregistration
#define REGISTER_MATERIAL( ClassName ) namespace { MaterialRegistrator<ClassName> reg_ ## ClassName; }

#endif
