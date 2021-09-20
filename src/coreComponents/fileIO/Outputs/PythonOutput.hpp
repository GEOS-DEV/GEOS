/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PythonOutput.hpp
 */

#ifndef GEOSX_FILEIO_OUTPUTS_PYTHONOUTPUT_HPP_
#define GEOSX_FILEIO_OUTPUTS_PYTHONOUTPUT_HPP_

#include "fileIO/Outputs/OutputBase.hpp"

namespace geosx
{

/**
 * @class PythonOutput
 * @brief Performs no actual output but returns control to the Python interpreter.
 */
class PythonOutput : public OutputBase
{
public:

  /**
   * @brief Constructor.
   * @param name The name of this PythonOutput in the data repository.
   * @param parent The parent group.
   */
  PythonOutput( string const & name,
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
   * @brief Signals the EventManager to exit the loop early, and return to Python.
   * @note It is an error to use this event when not running through pygeosx.
   * @copydetails EventBase::execute()
   */
  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override
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

#endif // GEOSX_FILEIO_OUTPUTS_PYTHONOUTPUT_HPP_
