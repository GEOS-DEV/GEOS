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

/**
 * @class EventManager
 *
 * A class for managing code events.
 */
class EventManager : public dataRepository::ManagedGroup
{
public:
  /// Main constructor
  EventManager( std::string const & name,
                ManagedGroup * const parent );

  /// Destructor
  virtual ~EventManager() override;

  /// Documentation assignment
  virtual void FillDocumentationNode() override;

  /// XML post processing
  virtual void ReadXML_PostProcess() override;

  /// A method to add child events
  virtual void CreateChild( string const & childKey, string const & childName ) override;

  /**
   * The main execution loop for the code.  During each cycle, it will:
   *   - Calculate the event forecast (number of cycles until its expected execution)
   *   - Signal an event to prepare (forecast == 1)
   *   - Execute an event (forecast == 0)
   *   - Determine dt for the next cycle
   *   - Advance time, cycle, etc.
   */
  void Run(dataRepository::ManagedGroup * domain);

  struct viewKeyStruct
  {
    dataRepository::ViewKey time = { "time" };
    dataRepository::ViewKey dt = { "dt" };
    dataRepository::ViewKey cycle = { "cycle" };
    dataRepository::ViewKey maxTime = { "maxTime" };
    dataRepository::ViewKey maxCycle = { "maxCycle" };
    dataRepository::ViewKey verbosity = { "verbosity" };
    dataRepository::ViewKey currentSubEvent = { "currentSubEvent" };
    dataRepository::ViewKey currentMaxDt = { "currentMaxDt" };
  } viewKeys;

  /// Catalog interface
  using CatalogInterface = cxx_utilities::CatalogInterface< EventBase, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_EVENTMANAGER_HPP_ */
