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
 * @file SinglePhaseWell.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELL_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELL_HPP_

#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"

#include "constitutive/fluid/singlefluid/SingleFluidLayouts.hpp"

namespace geos
{

namespace dataRepository
{
class Group;
}

namespace constitutive
{
class SingleFluidBase;
}
class WellElementSubRegion;

/**
 * @class SinglePhaseWell
 *
 * A single-phase well solver
 */
class SinglePhaseWell : public WellControls
{
public:

  using DerivOffset = constitutive::singlefluid::DerivativeOffsetC< 1 >;
  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  SinglePhaseWell( const string & name,
                   Group * const parent );

  /// deleted default constructor
  SinglePhaseWell() = delete;

  /// deleted copy constructor
  SinglePhaseWell( SinglePhaseWell const & ) = delete;

  /// default move constructor
  SinglePhaseWell( SinglePhaseWell && ) = delete;

  /// deleted assignment operator
  SinglePhaseWell & operator=( SinglePhaseWell const & ) = delete;

  /// deleted move operator
  SinglePhaseWell & operator=( SinglePhaseWell && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseWell() override = default;

  void registerWellDataOnMesh( WellElementSubRegion & subRegion ) override;

  /**
   * @defgroup WellManager Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   * The "Well" versions apply to individual well subRegions, whereas the others apply to all wells
   */
  /**@{*/
  /**
   *   * @brief Initialize well for the beginning of a simulation or restart
   *   @param domain the domain
   *   @param mesh the mesh level
   *   @param subRegion the well subRegion
   *  @param time_n the current time
   */
  virtual void initializeWell( DomainPartition & domain, MeshLevel & mesh, WellElementSubRegion & subRegion, real64 const & time_n )override;

  virtual void initializeWellPostInitialConditionsPreSubGroups( WellElementSubRegion & subRegion )override;

  virtual bool isCompositional() const override { return false; }
  /**
   * @copydoc WellControls::assembleWellAccumulationTerms()
   */
  virtual void assembleWellAccumulationTerms( real64 const & time,
                                              real64 const & dt,
                                              WellElementSubRegion & subRegion,
                                              DofManager const & dofManager,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs ) override;
  /**
   * @copydoc WellControls::assembleWellConstraintTerms()
   */
  virtual void assembleWellPressureRelations( real64 const & time_n,
                                              real64 const & dt,
                                              WellElementSubRegion const & subRegion,
                                              DofManager const & dofManager,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs ) override;

  /**
   * @copydoc WellControls::assembleWellConstraintTerms()
   */
  virtual void assembleWellConstraintTerms( real64 const & time_n,
                                            real64 const & dt,
                                            WellElementSubRegion const & subRegion,
                                            DofManager const & dofManager,
                                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                            arrayView1d< real64 > const & localRhs ) override;

  /**
   * @copydoc WellControls::computeWellPerforationRates()
   */
  virtual void computeWellPerforationRates( real64 const & time_n,
                                            real64 const & GEOS_UNUSED_PARAM( dt ),
                                            ElementRegionManager & elemManager,
                                            WellElementSubRegion & subRegion ) override;

  /**
   * @copydoc WellControls::assembleFluxTerms()
   */
  virtual void assembleWellFluxTerms( real64 const & time,
                                      real64 const & dt,
                                      WellElementSubRegion const & subRegion,
                                      DofManager const & dofManager,
                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                      arrayView1d< real64 > const & localRhs ) override;
  /**@}*/

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   * The "Well" versions apply to individual well subRegions, whereas the others apply to all wells
   */
  /**@{*/

  virtual array1d< real64 >
  calculateLocalWellResidualNorm( real64 const & time_n,
                                  real64 const & dt,
                                  NonlinearSolverParameters const & nonlinearSolverParameters,
                                  WellElementSubRegion const & subRegion,
                                  DofManager const & dofManager,
                                  arrayView1d< real64 const > const & localRhs )override;


  virtual real64
  calculateWellResidualNorm( real64 const & time_n,
                             real64 const & dt,
                             NonlinearSolverParameters const & nonlinearSolverParameters,
                             WellElementSubRegion const & subRegion,
                             DofManager const & dofManager,
                             arrayView1d< real64 const > const & localRhs ) override;

  virtual real64 scalingForWellSystemSolution( WellElementSubRegion & subRegion,
                                               DofManager const & dofManager,
                                               arrayView1d< real64 const > const & localSolution ) override;
  /**
   * @copydoc WellControls::checkSystemSolution()
   */


  virtual bool
  checkWellSystemSolution( WellElementSubRegion & subRegion,
                           DofManager const & dofManager,
                           arrayView1d< real64 const > const & localSolution,
                           real64 const scalingFactor,
                           real64 & minPressure,
                           real64 & minDensity,
                           real64 & minTotalDensity,
                           ElementsReporterBuffer & negPressureIds,
                           ElementsReporterBuffer & negDensityIds,
                           ElementsReporterBuffer & negTotalDensityIds ) override;
  /**
   * @copydoc WellControls::applyWellSystemSolution()
   */
  virtual void
  applyWellSystemSolution( DofManager const & dofManager,
                           arrayView1d< real64 const > const & localSolution,
                           real64 const scalingFactor,
                           real64 const dt,
                           DomainPartition & domain,
                           MeshLevel & mesh,
                           WellElementSubRegion & subRegion ) override;

  virtual void applyWellBoundaryConditions ( real64 const GEOS_UNUSED_PARAM( time_n ),
                                             real64 const GEOS_UNUSED_PARAM( dt ),
                                             ElementRegionManager & GEOS_UNUSED_PARAM( elemManager ),
                                             WellElementSubRegion & GEOS_UNUSED_PARAM( subRegion ),
                                             DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                             arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ),
                                             CRSMatrixView< real64, globalIndex const > const & GEOS_UNUSED_PARAM( localMatrix ) )override {};

  virtual void resetStateToBeginningOfStep( ElementRegionManager const & elemManager, WellElementSubRegion & subRegion ) override;

  virtual void implicitStepSetup( real64 const & time_n,
                                  real64 const & GEOS_UNUSED_PARAM( dt ),
                                  ElementRegionManager & elemManager,
                                  WellElementSubRegion & subRegion )override;

  virtual void
  implicitStepComplete( real64 const & time,
                        real64 const & dt,
                        WellElementSubRegion const & subRegion ) override;

  virtual void printRates( real64 const & time_n,
                           real64 const & dt,
                           WellElementSubRegion const & subRegion ) override;

  virtual real64 updateSubRegionState( ElementRegionManager const & elemManager, WellElementSubRegion & subRegion ) override;

  /**@}*/

  virtual string wellElementDofName() const override { return viewKeyStruct::dofFieldString(); }

  virtual string resElementDofName() const override;

  virtual localIndex numFluidComponents() const override { return 0; }

  virtual localIndex numFluidPhases() const override { return 1; }

  /**
   * @brief Recompute the volumetric rate that are used in the well constraints
   * @param elemManager the well region manager
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  virtual void calculateReferenceElementRates( WellElementSubRegion & subRegion );

  /**
   * @brief Recompute the BHP pressure that is used in the well constraints
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  virtual void updateBHPForConstraint( WellElementSubRegion & subRegion );

  /**
   * @brief Update fluid constitutive model state
   * @param subRegion the well subRegion containing the well elements and their associated fields
   */
  virtual void updateFluidModel( WellElementSubRegion & subRegion ) const;
  /**
   * @brief Update separator model state
   * @param elemManager the element region manager
   * @param subRegion the well subRegion containing the separator
   */
  void updateSeparator( ElementRegionManager const & elemManager, WellElementSubRegion & subRegion );

  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive
   * models)

   * @param
   * @param subRegion the well subRegion containing the well elements and their associated
   */
  virtual real64 updateWellState( ElementRegionManager const & elemManager, WellElementSubRegion & subRegion ) override;


  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive
   * models)
   * @param elemManager the element region manager
   * @param subRegion the well subRegion containing the well elements and their associated fields
   */
  /*
   * @brief apply a special treatment to the wells that are shut

   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  void shutDownWell( WellElementSubRegion & subRegion,
                     DofManager const & dofManager,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs );
  struct viewKeyStruct : WellControls::viewKeyStruct
  {
    static constexpr char const * dofFieldString() { return "wellVars"; }


  };

protected:

  virtual void initializePostInitialConditionsPreSubGroups() override;

  void saveState( WellElementSubRegion & subRegion );

  virtual void postRestartInitialization( )override;
  /// flag if negative pressure is allowed
  integer m_allowNegativePressure;

private:

  virtual void setConstitutiveNames( ElementSubRegionBase & subRegion ) const override;


  /**
   * @brief Make sure that the well constraints are compatible
   * @param time_n the time at the beginning of the time step
   * @param dt the time step dt
   * @param subRegion the well subRegion
   */
  virtual void validateWellConstraints( real64 const & time_n,
                                        real64 const & dt,
                                        WellElementSubRegion const & subRegion
                                        ) override;



  /**
   * @brief Create well separator
   */
  virtual void createSeparator( WellElementSubRegion & subRegion ) override;

};

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELL_HPP_
