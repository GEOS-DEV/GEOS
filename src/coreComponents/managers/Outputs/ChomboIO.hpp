/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
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
class ChomboIO final : public OutputBase
{
public:
  /// @copydoc geosx::dataRepository::Group::Group( string const & name, Group * const parent )
  ChomboIO( string const & name, Group * const parent );

  /// Destructor
  virtual ~ChomboIO() override;

  /**
   * @brief Catalog name interface
   * @return This type's catalog name
   */
  static string catalogName()
  { return "ChomboIO"; }

  /**
   * @brief Writes out a Chombo plot file.
   * @copydetails EventBase::execute()
   */
  virtual void execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        dataRepository::Group * const domain ) override;

  /**
   * @brief Writes out a Chombo plot file at the end of the simulation.
   * @copydetails ExecutableGroup::cleanup()
   */
  virtual void cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        dataRepository::Group * const domain ) override
  {
    m_waitForInput = 0;
    execute( time_n, 0.0, cycleNumber, eventCounter, eventProgress, domain );
  }

  /// @cond DO_NOT_DOCUMENT
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
  /// @endcond

private:
  ChomboCoupler * m_coupler;
  string m_outputPath;
  double m_beginCycle;
  string m_inputPath;
  integer m_waitForInput;
  integer m_useChomboPressures;
};


} /* namespace geosx */

#endif /* GEOSX_MANAGERS_OUTPUTS_CHOMBOIO_HPP_ */
