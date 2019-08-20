/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file SinglePhaseFlow.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_SINGLEPHASEFLOW_HPP_
#define SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_SINGLEPHASEFLOW_HPP_

#include "physicsSolvers/FiniteVolume/FlowSolverBase.hpp"
#include "constitutive/Fluid/SingleFluidBase.hpp"

namespace geosx
{

namespace dataRepository
{
class ManagedGroup;
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
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  SinglePhaseFlow( const std::string & name,
                   ManagedGroup * const parent );


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

  virtual void RegisterDataOnMesh( ManagedGroup * const MeshBodies ) override;

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
                           ParallelMatrix * const matrix,
                           ParallelVector * const rhs );

  template< bool ISPORO >
  void AccumulationLaunch( localIndex const er,
                           localIndex const esr,
                           FaceElementSubRegion const * const subRegion,
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



  /**@}*/


  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {
    // dof numbering
    static constexpr auto blockLocalDofNumberString = "blockLocalDofNumber_SinglePhaseFlow" ;

    // primary solution field
    static constexpr auto pressureString = "pressure";
    static constexpr auto deltaPressureString = "deltaPressure";
    static constexpr auto facePressureString = "facePressure";

    static constexpr auto deltaVolumeString = "deltaVolume";

    // intermediate fields
    static constexpr auto mobilityString = "mobility";
    static constexpr auto dMobility_dPressureString = "dMobility_dPressure";
    static constexpr auto porosityString = "porosity";

    // face fields
    static constexpr auto faceDensityString = "faceDensity";
    static constexpr auto faceViscosityString = "faceViscosity";
    static constexpr auto faceMobilityString = "faceMobility";

    //backup fields
    static constexpr auto porosityOldString = "porosityOld";
    static constexpr auto densityOldString = "densityOld";

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

private:

  void SetupSystem( DomainPartition * const domain,
                    DofManager & dofManager,
                    ParallelMatrix & matrix,
                    ParallelVector & rhs,
                    ParallelVector & solution );

  /**
   * @brief set the sparsity pattern for the linear system
   * @param domain the domain partition
   * @param sparsity the sparsity pattern matrix
   */
  void SetSparsityPattern( DomainPartition const * const domain,
                           ParallelMatrix * const matrix ) const;

  /**
   * @brief sets the dof indices for this solver
   * @param meshLevel the mesh object (single level only)
   * @param numLocalRows the number of local rows on this partition
   * @param numGlobalRows the number of global rows in the problem
   * @param offset the DOF offset for this solver in the case of a non-block system
   *
   * This function sets the number of global rows, and sets the dof numbers for
   * this solver. dof numbers are referred to trilinosIndices currently.
   */
  void SetNumRowsAndTrilinosIndices( MeshLevel * const meshLevel,
                                     localIndex & numLocalRows,
                                     globalIndex & numGlobalRows,
                                     localIndex offset ) const;

protected:
  virtual void InitializePostInitialConditions_PreSubGroups( dataRepository::ManagedGroup * const rootGroup ) override;

private:

  /**
   * @brief Setup stored views into domain data for the current step
   */
  void ResetViews( DomainPartition * const domain ) override;

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
  void UpdateFluidModel( ManagedGroup * const dataGroup ) const;

  /**
   * @brief Function to update all constitutive models
   * @param dataGroup group that contains the fields
   */
  void UpdateSolidModel( ManagedGroup * const dataGroup ) const;

  /**
   * @brief Function to update fluid mobility
   * @param dataGroup group that contains the fields
   */
  void UpdateMobility( ManagedGroup * const dataGroup ) const;

  /**
   * @brief Function to update all constitutive state and dependent variables
   * @param dataGroup group that contains the fields
   */
  void UpdateState( ManagedGroup * dataGroup ) const;


  /// views into primary variable fields

  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> m_dofNumber; // TODO will move to DofManager

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_pressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_deltaPressure;

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_deltaVolume;

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
};


} /* namespace geosx */

#endif //SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_SINGLEPHASEFLOW_HPP_
