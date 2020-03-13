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
  /// Main constructor
  RestartOutput( std::string const & name,
                 Group * const parent );

  /// Destructor
  virtual ~RestartOutput() override;

  /// Catalog name interface
  static string CatalogName() { return "Restart"; }

  /// This method will be called by the event manager if triggered
  virtual void Execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        dataRepository::Group * domain ) override;

  /// Write one final output as the code exits
  virtual void Cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        dataRepository::Group * domain ) override
  {
    Execute( time_n, 0, cycleNumber, eventCounter, eventProgress, domain );
  }

  struct viewKeyStruct
  {
    dataRepository::ViewKey writeFEMFaces = { "writeFEMFaces" };
  } viewKeys;

};


} /* namespace geosx */

#endif /* GEOSX_MANAGERS_OUTPUTS_RESTARTOUTPUT_HPP_ */
