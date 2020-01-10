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
 * @file SinglePhaseFlow.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FINITEVOLUME_SINGLEPHASEFLOW_HPP_
#define GEOSX_PHYSICSSOLVERS_FINITEVOLUME_SINGLEPHASEFLOW_HPP_

#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "constitutive/fluid/CompressibleSinglePhaseFluid.hpp"
#include "constitutive/solid/solidSelector.hpp"

namespace geosx
{

namespace dataRepository
{
class Group;
}
class FieldSpecificationBase;

class FiniteElementBase;

class DomainPartition;

/**
 * @class SinglePhaseFlow
 *
 * class to perform a single phase finite volume solve.
 */
class SinglePhaseFlow : public FlowSolverBase
{
public:
  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  SinglePhaseFlow( const std::string & name,
                   Group * const parent );


  /// deleted default constructor
  SinglePhaseFlow() = delete;

  /// deleted copy constructor
  SinglePhaseFlow( SinglePhaseFlow const & ) = delete;

  /// default move constructor
  SinglePhaseFlow( SinglePhaseFlow && ) = default;

  /// deleted assignment operator
  SinglePhaseFlow & operator=( SinglePhaseFlow const & ) = delete;

  /// deleted move operator
  SinglePhaseFlow & operator=( SinglePhaseFlow && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseFlow() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName()
  { return "SinglePhaseFlow"; }

  virtual void RegisterDataOnMesh( Group * const MeshBodies ) override;

  virtual void SetInitialTimeStep( Group * const domain ) override;

  virtual real64
  SolverStep( real64 const & time_n,
              real64 const & dt,
              integer const cycleNumber,
              DomainPartition * domain ) override;

  virtual real64
  ExplicitStep( real64 const & time_n,
                real64 const & dt,
                integer const cycleNumber,
                DomainPartition * domain ) override;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void
  ExplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition * const domain) override;

  virtual void
  ImplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition * const domain,
                     DofManager & dofManager,
                     ParallelMatrix & matrix,
                     ParallelVector & rhs,
                     ParallelVector & solution ) override;

  virtual void
  SetupDofs( DomainPartition const * const domain,
             DofManager & dofManager ) const override;

  virtual void
  AssembleSystem( real64 const time_n,
                  real64 const dt,
                  DomainPartition * const domain,
                  DofManager const & dofManager,
                  ParallelMatrix & matrix,
                  ParallelVector & rhs ) override;

  virtual void
  ApplyBoundaryConditions( real64 const time_n,
                           real64 const dt,
                           DomainPartition * const domain,
                           DofManager const & dofManager,
                           ParallelMatrix & matrix,
                           ParallelVector & rhs ) override;

  virtual real64
  CalculateResidualNorm( DomainPartition const * const domain,
                         DofManager const & dofManager,
                         ParallelVector const & rhs ) override;

  virtual void
  SolveSystem( DofManager const & dofManager,
               ParallelMatrix & matrix,
               ParallelVector & rhs,
               ParallelVector & solution ) override;

  virtual void
  ApplySystemSolution( DofManager const & dofManager,
                       ParallelVector const & solution,
                       real64 const scalingFactor,
                       DomainPartition * const domain ) override;

  virtual void
  ResetStateToBeginningOfStep( DomainPartition * const domain ) override;

  virtual void
  ImplicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition * const domain ) override;

  template< bool ISPORO >
  void AccumulationLaunch( localIndex const er,
                           localIndex const esr,
                           CellElementSubRegion const * const subRegion,
                           DofManager const * const dofManager,
                           ParallelMatrix * const matrix,
                           ParallelVector * const rhs );

  template< bool ISPORO >
  void AccumulationLaunch( localIndex const er,
                           localIndex const esr,
                           FaceElementSubRegion const * const subRegion,
                           DofManager const * const dofManager,
                           ParallelMatrix * const matrix,
                           ParallelVector * const rhs );


  /**
   * @brief assembles the accumulation terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  template< bool ISPORO >
  void AssembleAccumulationTerms( DomainPartition const * const domain,
                                  DofManager const * const dofManager,
                                  ParallelMatrix * const matrix,
                                  ParallelVector * const rhs );

  /**
   * @brief assembles the flux terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  void AssembleFluxTerms( real64 const time_n,
                          real64 const dt,
                          DomainPartition const * const domain,
                          DofManager const * const dofManager,
                          ParallelMatrix * const matrix,
                          ParallelVector * const rhs );

  /**
   * @brief assembles the flux terms for all cells in the explicit solver
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   */
  void AssembleFluxTermsExplicit( real64 const time_n,
                                  real64 const dt,
                                  DomainPartition * domain,
                                  DofManager const * const dofManager);

  /**
   * @brief assemble the flux terms for all cells and then apply flux BC in the explicit solver
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   */
  void CalculateAndApplyMassFlux( real64 const time_n,
                                  real64 const dt,
								  DomainPartition * const domain ) override;

  /**@}*/

  /**
   * @brief Function to update all constitutive state and dependent variables
   * @param dataGroup group that contains the fields
   */
  void UpdateState( Group * dataGroup ) const;

  /**
   * @brief Function to update all constitutive models
   * @param dataGroup group that contains the fields
   */
  void UpdateFluidModel( Group * const dataGroup ) const;

  /**
   * @brief Function to update all constitutive state and dependent variables in the explicit solver
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   */
  void UpdateEOS( real64 const time_n,
				  real64 const dt,
                  DomainPartition * const domain ) override;

  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {
    static constexpr auto facePressureString = "facePressure";

    // intermediate fields
    static constexpr auto mobilityString = "mobility";
    static constexpr auto dMobility_dPressureString = "dMobility_dPressure";

    static constexpr auto massString = "mass";

    // face fields
    static constexpr auto faceDensityString = "faceDensity";
    static constexpr auto faceViscosityString = "faceViscosity";
    static constexpr auto faceMobilityString = "faceMobility";

  } viewKeysSinglePhaseFlow;

  viewKeyStruct & viewKeys()
  { return viewKeysSinglePhaseFlow; }

  viewKeyStruct const & viewKeys() const
  { return viewKeysSinglePhaseFlow; }

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {
  } groupKeysSinglePhaseFlow;

  groupKeyStruct & groupKeys()
  { return groupKeysSinglePhaseFlow; }

  groupKeyStruct const & groupKeys() const
  { return groupKeysSinglePhaseFlow; }

  /**
   * @brief Setup stored views into domain data for the current step
   */
  void ResetViews( DomainPartition * const domain ) override;

protected:

  virtual void PostProcessInput() override;

  virtual void InitializePostInitialConditions_PreSubGroups( dataRepository::Group * const rootGroup ) override;

private:


  /**
   * @brief Function to perform the application of Dirichlet BCs on faces
   * @param time_n current time
   * @param dt time step
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  void ApplyFaceDirichletBC_implicit( real64 const time_n,
                                      real64 const dt,
                                      DofManager const * const dofManager,
                                      DomainPartition * const domain,
                                      ParallelMatrix * const matrix,
                                      ParallelVector * const rhs );

  /**
   * @brief Function to update all constitutive models
   * @param dataGroup group that contains the fields
   */
  void UpdateSolidModel( Group * const dataGroup ) const;

  /**
   * @brief Function to update fluid mobility
   * @param dataGroup group that contains the fields
   */
  void UpdateMobility( Group * const dataGroup ) const;

  /// views into primary variable fields

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_pressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_deltaPressure;

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_deltaVolume;

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_mass;

  /// views into intermediate fields

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_porosity;

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_mobility;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_dMobility_dPres;

  /// views into backup fields

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_porosityOld;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_densityOld;

  /// views into material fields

  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_pvMult;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_dPvMult_dPres;

  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_density;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_dDens_dPres;

  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_viscosity;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_dVisc_dPres;

  /// coupling with mechanics

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_totalMeanStressOld;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_totalMeanStress;

  ElementRegionManager::MaterialViewAccessor<arrayView1d<real64>> m_bulkModulus;
  ElementRegionManager::MaterialViewAccessor<real64> m_biotCoefficient;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_totalCompressibility;
};


} /* namespace geosx */

#endif //GEOSX_PHYSICSSOLVERS_FINITEVOLUME_SINGLEPHASEFLOW_HPP_
