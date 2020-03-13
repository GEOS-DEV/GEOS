/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ChomboIO.hpp
 */

#ifndef GEOSX_MANAGERS_OUTPUTS_CHOMBOIO_HPP_
#define GEOSX_MANAGERS_OUTPUTS_CHOMBOIO_HPP_

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
  ChomboIO( std::string const & name, Group * const parent );

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
                        dataRepository::Group * const domain ) final override;

  /// Write one final output as the code exits
  virtual void Cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        dataRepository::Group * const domain ) final override
  {
    m_waitForInput = 0;
    Execute( time_n, 0.0, cycleNumber, eventCounter, eventProgress, domain );
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
  ChomboCoupler * m_coupler;
  std::string m_outputPath;
  double m_beginCycle;
  std::string m_inputPath;
  integer m_waitForInput;
  integer m_useChomboPressures;
};


} /* namespace geosx */

#endif /* GEOSX_MANAGERS_OUTPUTS_CHOMBOIO_HPP_ */
