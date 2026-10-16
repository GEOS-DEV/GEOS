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
 * @file SolidMechanicsLagrangeContact.hpp
 *
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSLAGRANGECONTACT_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSLAGRANGECONTACT_HPP_

#include "physicsSolvers/solidMechanics/contact/ContactSolverBase.hpp"

namespace geos
{

class NumericalMethodsManager;

class SolidMechanicsLagrangeContact : public ContactSolverBase
{
public:
  SolidMechanicsLagrangeContact( const string & name,
                                 Group * const parent );

  ~SolidMechanicsLagrangeContact() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName()
  {
    return "SolidMechanicsLagrangeContact";
  }
  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  virtual void initializePreSubGroups() override;

  virtual void registerDataOnMesh( Group & MeshBodies ) override final;

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void setSparsityPattern( DomainPartition & domain,
                                   DofManager & dofManager,
                                   CRSMatrix< real64, globalIndex > & localMatrix,
                                   SparsityPattern< globalIndex > & pattern ) override;

  virtual std::unique_ptr< PreconditionerBase< LAInterface > >
  createPreconditioner( DomainPartition & domain ) const override;

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override final;

  virtual void
  implicitStepComplete( real64 const & time_n,
                        real64 const & dt,
                        DomainPartition & domain ) override final;

  virtual void
  assembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  virtual real64
  calculateResidualNorm( real64 const & time,
                         real64 const & dt,
                         DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override;

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       real64 const dt,
                       DomainPartition & domain ) override;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  void updateState( DomainPartition & domain ) override final;

  void assembleContact( real64 const time,
                        real64 const dt,
                        DomainPartition & domain,
                        DofManager const & dofManager,
                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                        arrayView1d< real64 > const & localRhs );

  void assembleForceResidualDerivativeWrtTraction( MeshLevel const & mesh,
                                                   string_array const & regionNames,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs );

  void assembleTractionResidualDerivativeWrtDisplacementAndTraction( MeshLevel const & mesh,
                                                                     string_array const & regionNames,
                                                                     DofManager const & dofManager,
                                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                     arrayView1d< real64 > const & localRhs );

  void assembleForceResidualPressureContribution( MeshLevel const & mesh,
                                                  string_array const & regionNames,
                                                  DofManager const & dofManager,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs );

  void assembleStabilization( MeshLevel const & mesh,
                              NumericalMethodsManager const & numericalMethodManager,
                              DofManager const & dofManager,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs );

  bool resetConfigurationToDefault( DomainPartition & domain ) const override final;

  bool updateConfiguration( DomainPartition & domain,
                            integer const configurationLoopIter ) override final;

  bool isFractureAllInStickCondition( DomainPartition const & domain ) const;

  void computeRotationMatrices( DomainPartition & domain ) const;

  void computeTolerances( DomainPartition & domain ) const;

  real64 const machinePrecision = std::numeric_limits< real64 >::epsilon();

  string getStabilizationName() const { return m_stabilizationName; }
  bool hasStabilization() const { return true;}

protected:

  real64 calculateContactResidualNorm( DomainPartition const & domain,
                                       DofManager const & dofManager,
                                       arrayView1d< real64 const > const & localRhs );

  virtual void postInputInitialization() override final;

  void setMGRStrategy();

private:
  string m_stabilizationName;

  real64 const m_slidingCheckTolerance = 0.05;

  real64 m_stabilizationScalingCoefficient = 1.0;

  void computeFaceDisplacementJump( DomainPartition & domain );

  struct viewKeyStruct : ContactSolverBase::viewKeyStruct
  {
    constexpr static char const * stabilizationNameString() { return "stabilizationName"; }

    constexpr static char const * rotationMatrixString() { return "rotationMatrix"; }

    constexpr static char const * normalDisplacementToleranceString() { return "normalDisplacementTolerance"; }

    constexpr static char const * normalTractionToleranceString() { return "normalTractionTolerance"; }

    constexpr static char const * slidingToleranceString() { return "slidingTolerance"; }

    constexpr static char const * transMultiplierString() { return "penaltyStiffnessTransMultiplier"; }

    constexpr static char const * stabilizationScalingCoefficientString() { return "stabilizationScalingCoefficient"; }

    constexpr static char const * useLocalYieldAccelerationString() { return "useLocalYieldAcceleration"; }

    constexpr static char const * localYieldAccelerationBufferString() { return "localYieldAccelerationBuffer"; }
  };

  /// Member variables and functions needed for yield acceleration. Naming convention follows ( Jiang & Tchelepi, 2019 )

  integer m_useLocalYieldAcceleration; // flag for applying acceleration to yield
  real64 m_localYieldAccelerationBuffer; // buffer to control the acceleration

  void tryLocalYieldAcceleration( FaceElementSubRegion & subRegion,
                                  integer configurationLoopIter,
                                  localIndex kfe,
                                  real64 currentTau_unscaled,
                                  real64 limitTau,
                                  real64 currentTau,
                                  integer & fractureState );

};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSLAGRANGECONTACT_HPP_ */
