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
 * @file RestartOutput.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_RESTARTOUTPUT_HPP_
#define SRC_COMPONENTS_CORE_SRC_RESTARTOUTPUT_HPP_

#include "OutputBase.hpp"


namespace geosx
{

/**
 * @class RestartOutput
 *
 * A class for creating restart-based outputs
 * Note: this class is redundant, and is expected to be retired
 */
class RestartOutput : public OutputBase
{
public:
  /// Main constructor
  RestartOutput( std::string const & name,
              ManagedGroup * const parent );

  /// Destructor
  virtual ~RestartOutput() override;

  /// Catalog name interface
  static string CatalogName() { return "Restart"; }

  /// Documentation assignment
  virtual void FillDocumentationNode() override;

  /// This method will be called by the event manager if triggered
  virtual void Execute( real64 const & time_n,
                        real64 const & dt,
                        const integer cycleNumber,
                        real64 const & eventPosition,
                        dataRepository::ManagedGroup * domain ) override;

  /// Write one final output as the code exits
  virtual void Cleanup( real64 const & time_n,
                        int const cycleNumber,
                        real64 const & eventPosition,
                        dataRepository::ManagedGroup * domain ) override { Execute(time_n, 0, cycleNumber, eventPosition, domain); }

  struct viewKeyStruct
  {
    dataRepository::ViewKey writeFEMFaces = { "writeFEMFaces" };
  } viewKeys;

};


} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_RESTARTOUTPUT_HPP_ */
