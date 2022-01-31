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
 * @file DARTSSuperEngine.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_DARTSSUPERENGINE_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_DARTSSUPERENGINE_HPP_

#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "functions/MultivariableTableFunction.hpp"

namespace geosx
{

/**
 * @class DARTSSuperEngine
 *
 * A compositional multiphase thermal reactive solver based on DARTS super engine
 * Takes into account diffusion, capillarity, kinetic and equilibrium reactions, precipitation/dissolution
 * through OBL approach using only cell-centered variables
 * works with TPFA
 */
//START_SPHINX_INCLUDE_00
class DARTSSuperEngine : public FlowSolverBase
{
//END_SPHINX_INCLUDE_00
public:

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  DARTSSuperEngine( const string & name,
                    Group * const parent );

  /// deleted default constructor
  DARTSSuperEngine() = delete;

  /// deleted copy constructor
  DARTSSuperEngine( DARTSSuperEngine const & ) = delete;

  /// default move constructor
  DARTSSuperEngine( DARTSSuperEngine && ) = default;

  /// deleted assignment operator
  DARTSSuperEngine & operator=( DARTSSuperEngine const & ) = delete;

  /// deleted move operator
  DARTSSuperEngine & operator=( DARTSSuperEngine && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~DARTSSuperEngine() override = default;

//START_SPHINX_INCLUDE_01
  /**
   * @brief name of the solver in the object catalog
   * @return string that contains the catalog name to generate a new object through the object catalog.
   */
  static string catalogName() { return "DARTSSuperEngine"; }

  virtual void registerDataOnMesh( Group & meshBodies ) override;
//END_SPHINX_INCLUDE_01

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual real64
  solverStep( real64 const & time_n,
              real64 const & dt,
              integer const cycleNumber,
              DomainPartition & domain ) override;

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  virtual void
  assembleSystem( real64 const time_n,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  virtual void
  applyBoundaryConditions( real64 const time_n,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override;

  virtual void
  solveSystem( DofManager const & dofManager,
               ParallelMatrix & matrix,
               ParallelVector & rhs,
               ParallelVector & solution ) override;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual void
  implicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain ) override;

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual real64
  calculateResidualNorm( DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override;

  virtual real64
  scalingForSystemSolution( DomainPartition const & domain,
                            DofManager const & dofManager,
                            arrayView1d< real64 const > const & localSolution ) override;

  virtual bool
  checkSystemSolution( DomainPartition const & domain,
                       DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor ) override;

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override;

  virtual void updateState( DomainPartition & domain ) override final;


  /**
   * @brief Recompute operator values and derivatives from primary variables
   * @param dataGroup the group storing the required fields
   */
  void updateOBLOperators( ObjectManagerBase & dataGroup ) const;

  /**
   * @brief Get the number of fluid components (species)
   * @return the number of components
   */
  localIndex numFluidComponents() const { return m_numComponents; }

  /**
   * @brief Get the number of fluid phases
   * @return the number of phases
   */
  localIndex numFluidPhases() const { return m_numPhases; }

  /**
   * @brief assembles the accumulation and volume balance terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
   */
  void assembleAccumulationAndVolumeBalanceTerms( real64 const dt,
                                                  DomainPartition & domain,
                                                  DofManager const & dofManager,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs ) const;

  /**@}*/

  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {
    static constexpr char const * elemDofFieldString() { return "compositionalMolarVariables"; }

    // inputs

    static constexpr char const * numPhasesString() { return "numPhases"; }

    static constexpr char const * numComponentsString() { return "numComponents"; }

    static constexpr char const * enableEnergyBalanceString() { return "enableEnergyBalance"; }

    static constexpr char const * componentNamesString() { return "componentNames"; }

    static constexpr char const * phaseNamesString() { return "phaseNames"; }

    static constexpr char const * OBLOperatorsTableFileString() { return "OBLOperatorsTableFile"; }

    static constexpr char const * computeCFLNumbersString() { return "computeCFLNumbers"; }

    static constexpr char const * maxCompFracChangeString() { return "maxCompFractionChange"; }

    static constexpr char const * allowLocalCompFractionChoppingString() { return "allowLocalCompFractionChopping"; }
  };

  /**
   * @brief Backup current values of all constitutive fields that participate in the accumulation term
   * @param domain the domain containing the mesh and fields
   */
  void backupFields( MeshLevel & mesh ) const;

  /**
   * @brief Function to perform the Application of Dirichlet type BC's
   * @param time current time
   * @param dt time step
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param localMatrix local system matrix
   * @param localRhs local system right-hand side vector
   */
  void applyDirichletBC( real64 const time,
                         real64 const dt,
                         DofManager const & dofManager,
                         DomainPartition & domain,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs ) const;

  /**
   * @brief Apply source flux boundary conditions to the system
   * @param time current time
   * @param dt time step
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param localMatrix local system matrix
   * @param localRhs local system right-hand side vector
   */
  void applySourceFluxBC( real64 const time,
                          real64 const dt,
                          DofManager const & dofManager,
                          DomainPartition & domain,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs ) const;

  /**
   * @brief Sets all the negative component densities (if any) to zero.
   * @param domain the physical domain object
   */
  void chopNegativeDensities( DomainPartition & domain );

  virtual void initializePostInitialConditionsPreSubGroups() override;


  /**
   * @brief Compute the largest CFL number in the domain
   * @param dt the time step size
   * @param domain the domain containing the mesh and fields
   */
  void
  computeCFLNumbers( real64 const & dt, DomainPartition & domain );


  virtual void initializePreSubGroups() override;


protected:

  virtual void postProcessInput() override;

  /// the max number of fluid phases
  integer m_numPhases;

  /// the number of fluid components
  integer m_numComponents;

  /// list of component names
  array1d< string > m_componentNames;

  /// list of phase names names
  array1d< string > m_phaseNames;

  /// the number of OBL operators
  integer m_numOBLOperators;

  /// OBL operators table file (if OBL physics becomes consitutive, multiple regions will be supported )
  string m_OBLOperatorsTableFile;

  /// OBL operators table function tabulated vs all primary variables
  MultivariableTableFunction const * m_OBLOperatorsTable;

  /// flag indicating whether energy balance will be enabled or not
  integer m_enableEnergyBalance;

  /// flag indicating whether CFL numbers will be computed or not
  integer m_computeCFLNumbers;

  /// maximum (absolute) change in a component fraction between two Newton iterations
  real64 m_maxCompFracChange;

  /// minimum value of the scaling factor obtained by enforcing maxCompFracChange
  real64 m_minScalingFactor;

  /// flag indicating whether local (cell-wise) chopping of negative compositions is allowed
  integer m_allowCompFracChopping;
};


} // namespace geosx


#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_DARTSSUPERENGINE_HPP_
