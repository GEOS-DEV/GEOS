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
 * @file SolidMechanicsPostEquilibrationStep.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSPOSTEQUILIBRATIONSTEP_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSPOSTEQUILIBRATIONSTEP_HPP_

#include "events/tasks/TaskBase.hpp"

namespace geosx
{

class SolidMechanicsLagrangianFEM;

/**
 * @class SolidMechanicsPostEquilibrationStep
 *
 * This class assumes that a equilibration step has been performed at the start of the simulation.
 * The role of the class is to implement a task that resets the total displacements to zero after the equilibration step.
 */
class SolidMechanicsPostEquilibrationStep : public TaskBase
{
public:

  /**
   * @brief Constructor for the post-equilibration step class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  SolidMechanicsPostEquilibrationStep( const string & name,
                                       Group * const parent );

  /// Destructor for the class
  ~SolidMechanicsPostEquilibrationStep() override;

  /// Accessor for the catalog name
  static string catalogName() { return "SolidMechanicsPostEquilibrationStep"; }

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

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct
  {
    /// String for the solid solver name
    constexpr static char const * solidSolverNameString() { return "solidSolverName"; }
  };

private:

  void postProcessInput() override;

  /// Name of the solid mechanics solver
  string m_solidSolverName;

  /// Pointer to the solid mechanics solver
  SolidMechanicsLagrangianFEM * m_solidSolver;

};


} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSPOSTEQUILIBRATIONSTEP_HPP_ */
