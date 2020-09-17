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
 * @file SinglePhaseBase.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASE_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASE_HPP_

#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseKernels.hpp"

namespace geosx
{


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
  SinglePhaseBase( const std::string & name,
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

  virtual void RegisterDataOnMesh( Group * const MeshBodies ) override;

  virtual void SetupSystem( DomainPartition * const domain,
                            DofManager & dofManager,
                            ParallelMatrix & matrix,
                            ParallelVector & rhs,
                            ParallelVector & solution ) override;

  virtual real64
  SolverStep( real64 const & time_n,
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
  ImplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition * const domain,
                     DofManager & dofManager,
                     ParallelMatrix & matrix,
                     ParallelVector & rhs,
                     ParallelVector & solution ) override;

  virtual void
  AssembleSystem( real64 const time_n,
                  real64 const dt,
                  DomainPartition * const domain,
                  DofManager const & dofManager,
                  ParallelMatrix & matrix,
                  ParallelVector & rhs ) override;


  virtual void
  SolveSystem( DofManager const & dofManager,
               ParallelMatrix & matrix,
               ParallelVector & rhs,
               ParallelVector & solution ) override;

  virtual void
  ResetStateToBeginningOfStep( DomainPartition * const domain ) override;

  virtual void
  ImplicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition * const domain ) override;

  template< bool ISPORO >
  void AccumulationLaunch( localIndex const er,
                           localIndex const esr,
                           CellElementSubRegion const & subRegion,
                           DofManager const & dofManager,
                           ParallelMatrix * const matrix,
                           ParallelVector * const rhs,
			   DomainPartition const * const domain);

  template< bool ISPORO >
  void AccumulationLaunch( localIndex const er,
                           localIndex const esr,
                           FaceElementSubRegion const & subRegion,
                           DofManager const & dofManager,
                           ParallelMatrix * const matrix,
                           ParallelVector * const rhs,
			   DomainPartition const * const domain);


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
  virtual void
  AssembleFluxTerms( real64 const time_n,
                     real64 const dt,
                     DomainPartition const * const domain,
                     DofManager const * const dofManager,
                     ParallelMatrix * const matrix,
                     ParallelVector * const rhs ) = 0;



  /**@}*/

  /**
   * @brief Function to update all constitutive state and dependent variables
   * @param dataGroup group that contains the fields
   */
  virtual void UpdateState( Group & dataGroup, localIndex const targetIndex ) const;

  /**
   * @brief Function to update all constitutive models
   * @param dataGroup group that contains the fields
   */
  virtual void UpdateFluidModel( Group & dataGroup, localIndex const targetIndex ) const;

  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {
    static constexpr auto boundaryFacePressureString = "boundaryFacePressure";

    // intermediate fields
    static constexpr auto mobilityString = "mobility";
    static constexpr auto dMobility_dPressureString = "dMobility_dPressure";
    static constexpr auto porosityString = "porosity";

    // face fields
    static constexpr auto boundaryFaceDensityString = "boundaryFaceDensity";
    static constexpr auto boundaryFaceViscosityString = "boundaryFaceViscosity";
    static constexpr auto boundaryFaceMobilityString = "boundaryFaceMobility";

    //backup fields
    static constexpr auto porosityOldString = "porosityOld";
    static constexpr auto densityOldString = "densityOld";
    static constexpr auto transTMultString = "transTMult";
    static constexpr auto poroMultString = "poroMult";

  } viewKeysSinglePhaseBase;

  viewKeyStruct & viewKeys()
  { return viewKeysSinglePhaseBase; }

  viewKeyStruct const & viewKeys() const
  { return viewKeysSinglePhaseBase; }

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {} groupKeysSinglePhaseBase;

  groupKeyStruct & groupKeys()
  { return groupKeysSinglePhaseBase; }

  groupKeyStruct const & groupKeys() const
  { return groupKeysSinglePhaseBase; }

  /**
   * @brief Setup stored views into domain data for the current step
   */
  virtual void ResetViews( DomainPartition * const domain ) override;

protected:

  virtual void ValidateFluidModels( DomainPartition const & domain ) const;

  virtual void InitializePreSubGroups( Group * const rootGroup ) override;

  virtual void InitializePostInitialConditions_PreSubGroups( dataRepository::Group * const rootGroup ) override;

  /**
   * @brief Function to update all constitutive models
   * @param dataGroup group that contains the fields
   */
  void UpdateSolidModel( Group & dataGroup, localIndex const targetIndex ) const;

  /**
   * @brief Function to update fluid mobility
   * @param dataGroup group that contains the fields
   */
  template< class FLUIDBASE >
  void UpdateMobility( Group & dataGroup, localIndex const targetIndex ) const;


  /// views into primary variable fields

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_pressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_deltaPressure;

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_deltaVolume;

  /// views into intermediate fields

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_porosity;

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_mobility;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_dMobility_dPres;

  /// views into backup fields

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_porosityOld;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_densityOld;

  /// views into material fields

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_pvMult;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_dPvMult_dPres;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_density;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_dDens_dPres;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_viscosity;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_dVisc_dPres;

  /// coupling with mechanics

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_totalMeanStressOld;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_totalMeanStress;

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_bulkModulus;
  ElementRegionManager::ElementViewAccessor< real64 > m_biotCoefficient;

  // coupling with proppant transport

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_poroMultiplier;

  ElementRegionManager::ElementViewAccessor< arrayView1d< R1Tensor > > m_transTMultiplier;

private:

  virtual void ResetViewsPrivate( ElementRegionManager * const elemManager );

};



template< class FLUIDBASE >
void SinglePhaseBase::UpdateMobility( Group & dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  // output

  arrayView1d< real64 > const & mob =
    dataGroup.getReference< array1d< real64 > >( viewKeyStruct::mobilityString );

  arrayView1d< real64 > const & dMob_dPres =
    dataGroup.getReference< array1d< real64 > >( viewKeyStruct::dMobility_dPressureString );

  // input

  FLUIDBASE & fluid = GetConstitutiveModel< FLUIDBASE >( dataGroup, m_fluidModelNames[targetIndex] );

  arrayView2d< real64 const > const & dens =
    fluid.template getReference< array2d< real64 > >( FLUIDBASE::viewKeyStruct::densityString );

  arrayView2d< real64 const > const & dDens_dPres =
    fluid.template getReference< array2d< real64 > >( FLUIDBASE::viewKeyStruct::dDens_dPresString );

  arrayView2d< real64 const > const & visc =
    fluid.template getReference< array2d< real64 > >( FLUIDBASE::viewKeyStruct::viscosityString );

  arrayView2d< real64 const > const & dVisc_dPres =
    fluid.template getReference< array2d< real64 > >( FLUIDBASE::viewKeyStruct::dVisc_dPresString );

  SinglePhaseBaseKernels::MobilityKernel::Launch( dataGroup.size(),
                                                  dens,
                                                  dDens_dPres,
                                                  visc,
                                                  dVisc_dPres,
                                                  mob,
                                                  dMob_dPres );

}

} /* namespace geosx */

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASE_HPP_
