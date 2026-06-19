/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DLPhysicsSolverBase.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_DLPHYSICSSOLVERBASE_HPP_
#define GEOS_PHYSICSSOLVERS_DLPHYSICSSOLVERBASE_HPP_

#include "codingUtilities/traits.hpp"
#include "common/DataTypes.hpp"
#include "common/format/LogPart.hpp"
#include "dataRepository/ExecutableGroup.hpp"
#include "dataRepository/RestartFlags.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/utilities/LinearSolverResult.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "mesh/MeshBody.hpp"
#include "physicsSolvers/NonlinearSolverParameters.hpp"
#include "physicsSolvers/LinearSolverParameters.hpp"
#include "physicsSolvers/SolverStatistics.hpp"
#include "physicsSolvers/PhysicsSolverBase.hpp"
#include "DLSharedMemoryManager.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"

#include <limits>

namespace geos
{

class DomainPartition;

/**
 * @class DLPhysicsSolverBase
 * @brief Base class for all DL physics solvers
 *
 * This class provides the base interface for all DL physics solvers. It provides the basic
 * functionality for setting up and solving a system using DL, as well as the interface for
 * performing a timestep.
 */
class DLPhysicsSolverBase : public PhysicsSolverBase
{
public:
  /**
   * @brief Constructor for DLPhysicsSolverBase
   * @param name the name of this instantiation of DLPhysicsSolverBase
   * @param parent the parent group of this instantiation of DLPhysicsSolverBase
   */
  explicit DLPhysicsSolverBase( string const & name,
                                Group * const parent );

  /**
   * @brief Move constructor for DLPhysicsSolverBase
   */
  DLPhysicsSolverBase( DLPhysicsSolverBase && ) = default;

  /**
   * @brief Destructor for DLPhysicsSolverBase
   */
  virtual ~DLPhysicsSolverBase() override;

  /**
   * @brief Deleted constructor
   */
  DLPhysicsSolverBase() = delete;

  /**
   * @brief Deleted copy constructor
   */
  DLPhysicsSolverBase( DLPhysicsSolverBase const & ) = delete;

  /**
   * @brief Deleted copy assignment operator
   */
  DLPhysicsSolverBase & operator=( DLPhysicsSolverBase const & ) = delete;

  /**
   * @brief Deleted move assignment operator
   */
  DLPhysicsSolverBase & operator=( DLPhysicsSolverBase && ) = delete;

  /**
   * @brief Function for a nonlinear implicit integration step
   * @param time_n time at the beginning of the step
   * @param dt the perscribed timestep
   * @param cycleNumber the current cycle number
   * @param domain the domain object
   * @return return the timestep that was achieved during the step.
   *
   * This function implements a nonlinear newton method for implicit DL problems. It requires that the
   * other functions in the solver interface are implemented in the derived physics solver.
   */
  virtual real64 nonlinearImplicitStep( real64 const & time_n,
                                        real64 const & dt,
                                        integer const cycleNumber,
                                        DomainPartition & domain ) override;



  virtual void
  setupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               ParallelVector & rhs,
               ParallelVector & solution,
               bool const setSparsity = true ) override;

#if defined(GEOS_USE_PYGEOSX)
  PyTypeObject * getPythonType() const;
#endif

protected:
  /**
   * @brief Solve a nonlinear system using a DL approach
   * @param time_n the time at the beginning of the step
   * @param dt the desired timestep
   * @param cycleNumber the current cycle number
   * @param domain the domain partition
   * @return true if the nonlinear system was solved, false otherwise
   */
  virtual bool solveNonlinearSystemUsingDL( real64 const & time_n,
                                            real64 const & dt,
                                            integer const cycleNumber,
                                            DomainPartition & domain );

  virtual void initializePostInitialConditionsPostSubGroups() override;


  /**
   * @brief Populate m_dofXCoords/m_dofYCoords/m_dofZCoords for different dofs.
   *
   * @param dofManager  DofManager holding DOF layout and rank offset
   * @param domain      DomainPartition used to iterate mesh levels / element subregions
   * @param elemDofFieldKey  The field key under which element DOF numbers are stored (e.g. "singlePhaseVariables")
   */
  virtual void populateDofCoords( DofManager const & dofManager,
                                  DomainPartition & domain,
                                  string const & elemDofFieldKey );

  /**
   * @brief Populate m_strainTrace for different dofs.
   *
   * @param dofManager  DofManager holding DOF layout and rank offset
   * @param domain      DomainPartition used to iterate mesh levels / element subregions
   * @param elemDofFieldKey  The field key under which element DOF numbers are stored (e.g. "singlePhaseVariables")
   */
  virtual void populateDofStrainTrace( DofManager const & dofManager,
                                       DomainPartition & domain,
                                       string const & elemDofFieldKey );

  /**
   * @brief Populate m_prevSolution for different dofs.
   *
   * @param dofManager  DofManager holding DOF layout and rank offset
   * @param domain      DomainPartition used to iterate mesh levels / element subregions
   * @param elemDofFieldKey  The field key under which element DOF numbers are stored (e.g. "singlePhaseVariables")
   */
  virtual void populateDofPrevSolution( DofManager const & dofManager,
                                        DomainPartition & domain,
                                        string const & elemDofFieldKey );


  /**
   * @brief Initialize shared memories needed for DL simulations
   * @return void
   */
  virtual void initializeSharedMemories();

  /**
   * @brief Share the DL model inputs through shared memory
   * @return void
   */
  virtual void shareDLModelInputs( real64 const & time_n,
                                   real64 const & dt,
                                   integer const cycleNumber,
                                   DomainPartition & domain );

  /**
   * @brief Read the DL model outputs through shared memory
   * @return void
   */
  virtual void readDLModelOutputs( real64 const & time_n,
                                   real64 const & dt,
                                   integer const cycleNumber,
                                   DomainPartition & domain );

  DLSharedMemoryManager m_sharedMemoryManager;

  // Data vectors for DL solver
  ParallelVector m_dofXCoords;
  ParallelVector m_dofYCoords;
  ParallelVector m_dofZCoords;
  ParallelVector m_prevSolution;
  ParallelVector m_strainTrace;   //TODO: Consider Removing/Moving to a more suitable class. Like IFENNPoroMechanics


private:
};

} // namespace geos

#endif /* GEOS_PHYSICSSOLVERS_DLPHYSICSSOLVERBASE_HPP_ */
