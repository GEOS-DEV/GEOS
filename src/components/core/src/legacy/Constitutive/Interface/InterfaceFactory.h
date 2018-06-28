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
 * @file InterfaceFactory.h
 * @author Scott Johnson
 * @date Oct 6, 2013
 */

#ifndef INTERFACEFACTORY_H_
#define INTERFACEFACTORY_H_

#include "InterfaceBase.h"
#include "Utilities/StringUtilities.h"

#include <map>
#include <string>
#include <vector>

//////////////////////////

// Interface Factory
//
// Consists of the following parts:
//   * The function to generate new interface pointers: "NewInterface"
//   * A base class to derive the functions to generate interface pointers:
// "InterfaceInitializer"
//   * A String-to-Interface-Intializer map hidden behind the
// GetInterfaceCatalogue function
//   * A template to create interface initializers: "InterfaceRegistrator"
//   * A compiler directive to simplify autoregistration: "REGISTER_INTERFACE"
//
// Most interfaces will only need to use one or two of the parts:
//   * To register a new interface in the factory: REGISTER_INTERFACE(
// InterfaceClassName )
//   * To load a interface pointer from the factory:       InterfaceBase*
// aInterfacePtr = InterfaceFactory::NewInterface(interfaceString, args );

/// Base class to generate new Interface pointers
class InterfaceInitializer
{
public:
  virtual InterfaceBase* InitializeInterface(TICPP::HierarchicalDataNode* hdn) = 0;

  virtual ~InterfaceInitializer()
  {}
};

typedef std::map<std::string, InterfaceInitializer*> InterfaceCatalogueType;

class InterfaceFactory
{
public:
  /// The Interface Factory.
  static InterfaceBase* NewInterface(const std::string& interfaceName,
                                     TICPP::HierarchicalDataNode* hdn = 0);

  /// Interface to the Interface name -> Interface initializer map
  static InterfaceCatalogueType& GetInterfaceCatalogue();

  /// Return a list of supported interface names
  static void GetInterfaceNames(std::vector<std::string>& nameList);
};

/// Template for creating classes derived from InterfaceInitializer
template<class InterfaceType>
class InterfaceRegistrator : public InterfaceInitializer
{
public:
  InterfaceRegistrator(void)
  {
    InterfaceFactory::GetInterfaceCatalogue()[InterfaceType::Name()] = this;
  }

  InterfaceBase* InitializeInterface(TICPP::HierarchicalDataNode* hdn)
  {
    InterfaceBase* tmp = new InterfaceType();
    if(hdn)
    {
      tmp->resize(0,1);
      tmp->ReadXML(*hdn);
    }
    return tmp;
  }
};

/// Compiler directive to simplify autoregistration
#define REGISTER_INTERFACE( ClassName ) namespace { InterfaceRegistrator<ClassName> reg_ ## ClassName; }

#endif
