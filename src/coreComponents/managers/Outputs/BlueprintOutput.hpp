/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

/**
 * @file BlueprintOutput.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_BLUEPRINTOUTPUT_HPP_
#define SRC_COMPONENTS_CORE_SRC_BLUEPRINTOUTPUT_HPP_

#include "OutputBase.hpp"
#include "fileIO/blueprint/Blueprint.hpp"

namespace geosx
{

/**
 * @class BlueprintOutput
 *
 * A class for creating blueprint-based outputs
 */
class BlueprintOutput : public OutputBase
{
public:
  /// Main constructor
  BlueprintOutput( std::string const & name,
              ManagedGroup * const parent );

  /// Destructor
  virtual ~BlueprintOutput() override;

  /// Catalog name interface
  static string CatalogName() { return "Blueprint"; }

  /// This method will be called by the event manager if triggered
  virtual void Execute( real64 const & time_n,
                        real64 const & dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const & eventProgress,
                        dataRepository::ManagedGroup * domain ) override;

  /// Write one final output as the code exits
  virtual void Cleanup( real64 const & time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const & eventProgress,
                        dataRepository::ManagedGroup * domain ) override
  {
    Execute(time_n, 0, cycleNumber, eventCounter, eventProgress, domain);
  }

  
  struct viewKeyStruct
  {
    dataRepository::ViewKey writeFEMFaces = { "writeFEMFaces" };
  } viewKeys;

};


} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_BLUEPRINTOUTPUT_HPP_ */
