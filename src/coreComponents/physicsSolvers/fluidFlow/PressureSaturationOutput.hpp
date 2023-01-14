/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PressureSaturationOutput.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_PRESSURESATURATIONOUTPUT_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_PRESSURESATURATIONOUTPUT_HPP_

#include "events/tasks/TaskBase.hpp"

namespace geosx
{

class CompositionalMultiphaseBase;

/**
 * @class PressureSaturationOutput
 *
 */
class PressureSaturationOutput : public TaskBase
{
public:

  /**
   * @brief Constructor for the state reset class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  PressureSaturationOutput( const string & name,
                            Group * const parent );

  /// Destructor for the class
  ~PressureSaturationOutput() override;

  /// Accessor for the catalog name
  static string catalogName() { return "PressureSaturationOutput"; }

  /**
   * @defgroup Tasks Interface Functions
   *
   * This function implements the interface defined by the abstract TaskBase class
   */
  /**@{*/

  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /**@}*/

private:

  void postProcessInput() override;

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct
  {
    /// String for the flow solver name
    constexpr static char const * flowSolverNameString() { return "flowSolverName"; }
  };

  /// Name of the flow solver
  string m_flowSolverName;

  /// Pointer to the flow solver
  CompositionalMultiphaseBase * m_flowSolver;

};


} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_PRESSURESATURATIONOUTPUT_HPP_ */
