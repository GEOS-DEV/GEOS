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
 * @file RestartOutput.hpp
 */

#ifndef GEOSX_MANAGERS_OUTPUTS_RESTARTOUTPUT_HPP_
#define GEOSX_MANAGERS_OUTPUTS_RESTARTOUTPUT_HPP_

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
  /**
   * @brief Main constructor.
   * @param name The name of the object in the data repository.
   * @param parent The parent of this object in the data repository.
   **/
  RestartOutput( std::string const & name,
                 Group * const parent );

  /// Destructor
  virtual ~RestartOutput() override;

  /**
   * @brief Catalog name interface
   * @return This type's catalog name
   */
  static string CatalogName() { return "Restart"; }

  /**
   * @brief Writes out a restart file.
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
                        dataRepository::Group * domain ) override;

  /**
   * @brief Write one final restart file as the code exits
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
                        dataRepository::Group * domain ) override
  {
    Execute( time_n, 0, cycleNumber, eventCounter, eventProgress, domain );
  }

  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    dataRepository::ViewKey writeFEMFaces = { "writeFEMFaces" };
  } viewKeys;
  /// @endcond
};


} /* namespace geosx */

#endif /* GEOSX_MANAGERS_OUTPUTS_RESTARTOUTPUT_HPP_ */
