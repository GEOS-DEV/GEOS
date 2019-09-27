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
 * @file BlueprintOutput.hpp
 */

#ifndef GEOSX_BLUEPRINTOUTPUT_HPP_
#define GEOSX_BLUEPRINTOUTPUT_HPP_

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
              Group * const parent );

  /// Destructor
  virtual ~BlueprintOutput() override;

  /// Catalog name interface
  static string CatalogName() { return "Blueprint"; }

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
    Execute(time_n, 0, cycleNumber, eventCounter, eventProgress, domain);
  }

  
  struct viewKeyStruct
  {
    dataRepository::ViewKey writeFEMFaces = { "writeFEMFaces" };
  } viewKeys;

};


} /* namespace geosx */

#endif /* GEOSX_BLUEPRINTOUTPUT_HPP_ */
