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
  /**
   * @brief Main constructor
   * @param name The name of the object in the data repository
   * @param parent The parent of this object in the data repository
   **/
  ChomboIO( std::string const & name, Group * const parent );

  /// Destructor
  virtual ~ChomboIO() final override;

  /**
   * @brief Catalog name interface
   * @return This type's catalog name
   */
  static string CatalogName()
  { return "ChomboIO"; }

  /**
   * @brief Writes out a Chombo plot file.
   * @param time_n The current simulation time.
   * @param dt The current time step.
   * @param cycleNumber The current cycle.
   * @param eventCounter The event counter.
   * @param eventProgress The event progress.
   * @param domain The DomainPartition to write out up-casted to a Group.
   */
  virtual void Execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        dataRepository::Group * const domain ) final override;

  /**
   * @brief Writes out a Chombo plot file at the end of the simulation.
   * @param time_n The last simulation time.
   * @param cycleNumber The last cycle.
   * @param eventCounter The event counter.
   * @param eventProgress The event progress.
   * @param domain The DomainPartition to write out up-casted to a Group.
   */
  virtual void Cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        dataRepository::Group * const domain ) final override
  {
    m_waitForInput = 0;
    Execute( time_n, 0.0, cycleNumber, eventCounter, eventProgress, domain );
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
  std::string m_outputPath;
  double m_beginCycle;
  std::string m_inputPath;
  integer m_waitForInput;
  integer m_useChomboPressures;
};


} /* namespace geosx */

#endif /* GEOSX_MANAGERS_OUTPUTS_CHOMBOIO_HPP_ */
