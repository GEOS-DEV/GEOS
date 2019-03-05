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
 * @file ChomboIO.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CHOMBOIO_HPP_
#define SRC_COMPONENTS_CORE_SRC_CHOMBOIO_HPP_

#include "OutputBase.hpp"
#include "fileIO/coupling/ChomboCoupler.hpp"

namespace geosx
{

/**
 * @class ChomboIO
 *
 * A class for coupling to CHOMBO
 */
class ChomboIO : public OutputBase
{
public:
  /// Main constructor
  ChomboIO(std::string const & name, ManagedGroup * const parent);

  /// Destructor
  virtual ~ChomboIO() final override;

  /// Catalog name interface
  static string CatalogName()
  { return "ChomboIO"; }

  /// This method will be called by the event manager if triggered
  virtual void Execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        dataRepository::ManagedGroup * const domain ) final override;

  /// Write one final output as the code exits
  virtual void Cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        dataRepository::ManagedGroup * const domain ) final override 
  { 
    m_waitForInput = 0;
    Execute(time_n, 0.0, cycleNumber, eventCounter, eventProgress, domain);
  }

  struct viewKeyStruct
  {
    static constexpr auto outputPathString = "outputPath";
    static constexpr auto beginCycleString = "beginCycle";
    static constexpr auto inputPathString = "inputPath";
    static constexpr auto waitForInputString = "waitForInput";
    static constexpr auto useChomboPressuresString = "useChomboPressures";

    dataRepository::ViewKey outputPath = { outputPathString };
    dataRepository::ViewKey beginCycle = { beginCycleString };
    dataRepository::ViewKey inputPath = { inputPathString };
    dataRepository::ViewKey waitForInput = { waitForInputString };
    dataRepository::ViewKey useChomboPressures = { useChomboPressuresString };
  } viewKeys;

private:
  ChomboCoupler* m_coupler;
  std::string m_outputPath;
  double m_beginCycle;
  std::string m_inputPath;
  integer m_waitForInput;
  integer m_useChomboPressures;
};


} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CHOMBOIO_HPP_ */
