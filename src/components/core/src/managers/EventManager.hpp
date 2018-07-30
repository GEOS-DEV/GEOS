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
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/** A class for managing code events */

#ifndef SRC_COMPONENTS_CORE_SRC_EVENTMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_EVENTMANAGER_HPP_

#include "dataRepository/ManagedGroup.hpp"
#include "managers/Events/EventBase.hpp"


namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const Events("Events");
}
}

class EventManager : public dataRepository::ManagedGroup
{
public:
  EventManager( std::string const & name,
                ManagedGroup * const parent );

  virtual ~EventManager() override;

  virtual void FillDocumentationNode() override;

  virtual void CreateChild( string const & childKey, string const & childName ) override;

  void Run(dataRepository::ManagedGroup * domain);

  struct viewKeyStruct
  {
    dataRepository::ViewKey time = { "time" };
    dataRepository::ViewKey dt = { "dt" };
    dataRepository::ViewKey cycle = { "cycle" };
    dataRepository::ViewKey maxTime = { "maxTime" };
    dataRepository::ViewKey maxCycle = { "maxCycle" };
    dataRepository::ViewKey verbosity = { "verbosity" };
  } viewKeys;

  using CatalogInterface = cxx_utilities::CatalogInterface< EventBase, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_EVENTMANAGER_HPP_ */
