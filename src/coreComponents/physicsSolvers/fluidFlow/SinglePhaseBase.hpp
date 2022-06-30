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
 * @file SinglePhaseBase.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASE_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASE_HPP_

#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/ThermalSinglePhaseBaseKernels.hpp"

namespace geosx
{

namespace constitutive
{

class ConstitutiveBase;

}

/**
 * @class SinglePhaseBase
 *
 * base class to perform a single phase finite volume solve.
 */
class SinglePhaseBase : public FlowSolverBase
{
public:
  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  SinglePhaseBase( const string & name,
                   Group * const parent );


  /// deleted default constructor
  SinglePhaseBase() = delete;

  /// deleted copy constructor
  SinglePhaseBase( SinglePhaseBase const & ) = delete;

  /// default move constructor
  SinglePhaseBase( SinglePhaseBase && ) = default;

  /// deleted assignment operator
  SinglePhaseBase & operator=( SinglePhaseBase const & ) = delete;

  /// deleted move operator
  SinglePhaseBase & operator=( SinglePhaseBase && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseBase() override = default;

  virtual void registerDataOnMesh( Group & meshBodies ) override;

  virtual void setupSystem( DomainPartition & domain,
                            DofManager & dofManager,
                            CRSMatrix< real64, globalIndex > & localMatrix,
                            ParallelVector & rhs,
                            ParallelVector & solution,
                            bool const setSparsity = true ) override;

  virtual real64
  solverStep( real64 const & time_n,
              real64 const & dt,
              integer const cycleNumber,
              DomainPartition & domain ) override;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  ///@{

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
  solveLinearSystem( DofManager const & dofManager,
                     ParallelMatrix & matrix,
                     ParallelVector & rhs,
                     ParallelVector & solution ) override;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual void
  implicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain ) override;

  ///@}

  /**
   * @brief assembles the accumulation terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
   */
  void assembleAccumulationTerms( DomainPartition & domain,
                                  DofManager const & dofManager,
                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                  arrayView1d< real64 > const & localRhs );

  /**
   * @brief assembles the flux terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
   */
  virtual void
  assembleFluxTerms( real64 const time_n,
                     real64 const dt,
                     DomainPartition const & domain,
                     DofManager const & dofManager,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs ) = 0;

  /**
   * @brief assembles the flux terms for all cells for the poroelastic case
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
   * @param jumpDofKey dofKey of the displacement jump
   */
  virtual void
  assemblePoroelasticFluxTerms( real64 const time_n,
                                real64 const dt,
                                DomainPartition const & domain,
                                DofManager const & dofManager,
                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                arrayView1d< real64 > const & localRhs,
                                string const & jumpDofKey ) = 0;

  /**
   * @brief assembles the flux terms for all cells for the hydrofracture case
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
   * @param dR_dAper
   */
  virtual void
  assembleHydrofracFluxTerms( real64 const time_n,
                              real64 const dt,
                              DomainPartition const & domain,
                              DofManager const & dofManager,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs,
                              CRSMatrixView< real64, localIndex const > const & dR_dAper ) = 0;

  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {
    static constexpr char const * elemDofFieldString() { return "primaryVariables"; }

    // inputs
    static constexpr char const * inputTemperatureString() { return "temperature"; }
    static constexpr char const * thermalConductivityNamesString() { return "thermalConductivityNames"; }
    static constexpr char const * solidInternalEnergyNamesString() { return "solidInternalEnergyNames"; }
  };

  /**
   * @brief Function to perform the Application of Dirichlet type BC's
   * @param time current time
   * @param dt time step
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param localMatrix local system matrix
   * @param localRhs local system right-hand side vector
   */
  void
  applyDirichletBC( real64 const time_n,
                    real64 const dt,
                    DomainPartition & domain,
                    DofManager const & dofManager,
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
  void
  applySourceFluxBC( real64 const time_n,
                     real64 const dt,
                     DomainPartition & domain,
                     DofManager const & dofManager,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs ) const;

  /**
   * @brief Apply aquifer boundary conditions to the system
   * @param time current time
   * @param dt time step
   * @param domain the domain
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix local system matrix
   * @param localRhs local system right-hand side vector
   */
  virtual void
  applyAquiferBC( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) const = 0;

  virtual void
  updateState ( DomainPartition & domain ) override final;

  /**
   * @brief Function to update all constitutive state and dependent variables
   * @param dataGroup group that contains the fields
   */
  void
  updateFluidState( ObjectManagerBase & subRegion ) const;


  /**
   * @brief Function to update all constitutive models
   * @param dataGroup group that contains the fields
   */
  virtual void
  updateFluidModel( ObjectManagerBase & dataGroup ) const;

  /**
   * @brief Update all relevant solid internal energy models using current values of temperature
   * @param dataGroup the group storing the required fields
   */
  void updateSolidInternalEnergyModel( ObjectManagerBase & dataGroup ) const;

  /**
   * @brief Function to update fluid mobility
   * @param dataGroup group that contains the fields
   */
  void
  updateMobility( ObjectManagerBase & dataGroup ) const;

  /**
   * @brief Setup stored views into domain data for the current step
   */

  virtual void initializePreSubGroups() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

  /**
   * @brief Compute the hydrostatic equilibrium using the compositions and temperature input tables
   */
  void computeHydrostaticEquilibrium();

protected:

  /**
   * @brief Checks constitutive models for consistency
   * @param[in] domain the domain partition
   */
  virtual void
  validateConstitutiveModels( DomainPartition & domain ) const;

  /**
   * @brief Initialize the aquifer boundary condition (gravity vector, water phase index)
   */
  void
  initializeAquiferBC() const;

  virtual void setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const override;


  /**
   * @brief Structure holding views into fluid properties used by the base solver.
   */
  struct FluidPropViews
  {
    arrayView2d< real64 const > const dens;             ///< density
    arrayView2d< real64 const > const dDens_dPres;      ///< derivative of density w.r.t. pressure
    arrayView2d< real64 const > const visc;             ///< viscosity
    arrayView2d< real64 const > const dVisc_dPres;      ///< derivative of viscosity w.r.t. pressure
    real64 const defaultDensity;                     ///< default density to use for new elements
    real64 const defaultViscosity;                    ///< default vi to use for new elements
  };

  /**
   * @brief Structure holding views into thermal fluid properties used by the base solver.
   */
  struct ThermalFluidPropViews
  {
    arrayView2d< real64 const > const dDens_dTemp;      ///< derivative of density w.r.t. temperature
    arrayView2d< real64 const > const dVisc_dTemp;      ///< derivative of viscosity w.r.t. temperature
  };

  /**
   * @brief Extract properties from a fluid.
   * @param fluid base reference to the fluid object
   * @return structure with property views
   *
   * This function enables derived solvers to substitute SingleFluidBase for a different,
   * unrelated fluid class, and customize property extraction. For example, it is used by
   * SinglePhaseProppantBase to allow using SlurryFluidBase, which does not inherit from
   * SingleFluidBase currently (but this design may need to be revised).
   */
  virtual FluidPropViews getFluidProperties( constitutive::ConstitutiveBase const & fluid ) const;

  virtual ThermalFluidPropViews getThermalFluidProperties( constitutive::ConstitutiveBase const & fluid ) const;

  /// the input temperature
  real64 m_inputTemperature;

private:
  virtual void setConstitutiveNames( ElementSubRegionBase & subRegion ) const override;

};

} /* namespace geosx */

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASE_HPP_
