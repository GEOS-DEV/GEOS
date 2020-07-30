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
 * @file PythonOutput.hpp
 */

#ifndef GEOSX_MANAGERS_OUTPUTS_PYTHONOUTPUT_HPP_
#define GEOSX_MANAGERS_OUTPUTS_PYTHONOUTPUT_HPP_

#include "managers/Outputs/OutputBase.hpp"

namespace geosx
{

class PythonOutput : public OutputBase
{
public:

  PythonOutput( std::string const & name,
                Group * const parent ):
    OutputBase( name, parent )
  {}

  /**
   * @brief Destructor.
   */
  virtual ~PythonOutput() override
  {}

  /**
   * @brief Get the name used to register this object in an XML file.
   * @return The string "Python".
   */
  static string catalogName() { return "Python"; }

  /**
   * @brief Writes out a Blueprint plot file.
   * @copydetails EventBase::Execute()
   */
  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        dataRepository::Group * domain ) override
  {
    GEOSX_UNUSED_VAR( time_n );
    GEOSX_UNUSED_VAR( dt );
    GEOSX_UNUSED_VAR( cycleNumber );
    GEOSX_UNUSED_VAR( eventCounter );
    GEOSX_UNUSED_VAR( eventProgress );
    GEOSX_UNUSED_VAR( domain );
    return true;
  }
};


} // namespace geosx

#endif // GEOSX_MANAGERS_OUTPUTS_PYTHONOUTPUT_HPP_
