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
 * @file SolidMechanicsMixedVEM.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSMIXEDVEM_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSMIXEDVEM_HPP_

#include "physicsSolvers/PhysicsSolverBase.hpp"
#include "mixedVEM/MixedVEMTypes.hpp"

namespace geos
{

/**
 * @class SolidMechanicsMixedVEM
 *
 * Mixed VEM of Dassi, Lovadina, Visinoni (2020) for linear elasticity, unknowns in Sigma_h x U_h.
 * hybridization = 0 solves (19); hybridization = 1 solves H lambda = h and recovers (sigma_E, u_E).
 * A prescribed traction is essential, a prescribed displacement natural (Remark 1).
 */
class SolidMechanicsMixedVEM : public PhysicsSolverBase
{
public:

  /**
   * @brief Constructor.
   * @param name the name of this solver
   * @param parent the parent group of this solver
   */
  SolidMechanicsMixedVEM( string const & name, dataRepository::Group * const parent );

  /// Deleted default constructor.
  SolidMechanicsMixedVEM() = delete;

  /// Deleted copy constructor.
  SolidMechanicsMixedVEM( SolidMechanicsMixedVEM const & ) = delete;

  /// Default move constructor.
  SolidMechanicsMixedVEM( SolidMechanicsMixedVEM && ) = default;

  /// Deleted copy assignment operator.
  SolidMechanicsMixedVEM & operator=( SolidMechanicsMixedVEM const & ) = delete;

  /// Deleted move assignment operator.
  SolidMechanicsMixedVEM & operator=( SolidMechanicsMixedVEM && ) = delete;

  virtual ~SolidMechanicsMixedVEM() override = default;

  /**
   * @brief Static factory catalog name.
   * @return the catalog name
   */
  static string catalogName() { return "SolidMechanicsMixedVEM"; }

  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  virtual string getCatalogName() const override { return catalogName(); }

  virtual void registerDataOnMesh( dataRepository::Group & meshBodies ) override final;

  virtual void setupDofs( DomainPartition const & domain,
                          DofManager & dofManager ) const override final;

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override final;

  virtual void
  assembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override final;

  virtual void
  applyBoundaryConditions( real64 const time,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override final;

  virtual real64
  calculateResidualNorm( real64 const & time,
                         real64 const & dt,
                         DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override final;

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       real64 const dt,
                       DomainPartition & domain ) override final;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override final;

  virtual void
  updateState( DomainPartition & domain ) override final;

  virtual void
  implicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain ) override final;

  virtual std::unique_ptr< PreconditionerBase< LAInterface > >
  createPreconditioner( DomainPartition & domain ) const override final;

  virtual real64
  solverStep( real64 const & time_n,
              real64 const & dt,
              integer const cycleNumber,
              DomainPartition & domain ) override final;

  /**
   * @brief View keys.
   */
  struct viewKeyStruct : public PhysicsSolverBase::viewKeyStruct
  {
    /// @return The key for the solid material names
    static constexpr char const * solidMaterialNamesString() { return "solidMaterialNames"; }
  };

  /**
   * @brief @return whether the hybridized form is in use.
   */
  bool useHybridization() const { return m_useHybridization; }

protected:

  virtual void postInputInitialization() override final;

  virtual void initializePreSubGroups() override final;

  virtual void setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const override final;

private:

  /**
   * @brief Mark the boundary role of every face and evaluate the prescribed data.
   * @param time the current time
   * @param mesh the mesh level
   */
  void classifyFaces( real64 const time, MeshLevel & mesh ) const;

  /**
   * @brief Recover the element unknowns if needed and fill the writable cell fields.
   * @param domain the domain partition
   */
  void computeCellFields( DomainPartition & domain ) const;

  /**
   * @brief Build the near null space of the interface operator.
   * @param domain the domain partition
   *
   * The multiplier traces of the six global rigid body motions span ker H.
   */
  void computeNearNullSpace( DomainPartition & domain ) const;

  /// hybridization flag, read from the MixedVEMDiscretization named by m_discretizationName
  bool m_useHybridization;

  /// length h of s_E, equation (15)
  mixedVEM::StabilizationLength m_stabilizationLength;

  /// multiplier traces of the six global rigid body motions, the near null space of H
  mutable array1d< ParallelVector > m_nearNullSpace;

};

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSMIXEDVEM_HPP_
