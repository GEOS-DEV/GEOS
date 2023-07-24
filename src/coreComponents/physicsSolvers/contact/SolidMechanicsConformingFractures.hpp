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
 * @file SolidMechanicsConformingFractures.hpp
 *
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSCONFORMINGFRACTURES_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSCONFORMINGFRACTURES_HPP_

#include "physicsSolvers/contact/ContactSolverBase.hpp"

namespace geos
{

class SolidMechanicsLagrangianFEM;

class SolidMechanicsConformingFractures : public ContactSolverBase
{
public:

  SolidMechanicsConformingFractures( const string & name,
                                     Group * const parent );

  ~SolidMechanicsConformingFractures() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName()
  {
    return "SolidMechanicsConformingFractures";
  }

  /// String used to form the solverName used to register single-physics solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "SolidMechanicsConformingFractures"; }

  virtual void initializePreSubGroups() override;

  virtual void registerDataOnMesh( Group & MeshBodies ) override final;

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  setupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               ParallelVector & rhs,
               ParallelVector & solution,
               bool const setSparsity = true ) override;

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
                       DomainPartition & domain ) override;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual real64
  setNextDt( real64 const & currentDt,
             DomainPartition & domain ) override;

  void updateState( DomainPartition & domain ) override final;

  
  void assembleForceResidualDerivativeWrtTraction( MeshLevel const & mesh,
                                                   arrayView1d< string const > const & regionNames,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs );

  void assembleTractionResidualDerivativeWrtDisplacementAndTraction( MeshLevel const & mesh,
                                                       arrayView1d< string const > const & regionNames,
                                                       DofManager const & dofManager,
                                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                       arrayView1d< real64 > const & localRhs );
                                                       
  string const & getContactRelationName() const { return m_contactRelationName; }

  bool resetConfigurationToDefault( DomainPartition & domain ) const override final;
  
  bool updateConfiguration( DomainPartition & domain ) override final;

  enum ContactEnforcementMethod : integer
  {
    Penalty = 0,                 ///< penalty-based contact enforcement  
    NodalLagrangeMultiplier = 1, ///< node-based Lagrange multipliers (node-to-node contact)
    FaceLagrangeMultiplier = 2   ///< face-based Lagrange multipliers (currently implemented in LagrangianContactSolver)
  };

protected:
  virtual void postProcessInput() override final;

  virtual void assembleNodalLagrangeMultiplierContact(DomainPartition & domain,
                                                      DofManager const & dofManager,
                                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                      arrayView1d< real64 > const & localRhs);


private:

  localIndex m_contactRelationFullIndex;

  real64 const m_slidingCheckTolerance = 0.05;

  real64 m_initialResidual[3] = {0.0, 0.0, 0.0};

  void computeNodalDisplacementJump( DomainPartition & domain ) const;
  
  void computeRotationMatrices( DomainPartition & domain ) const;

  void computeFaceNodalArea( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                             ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                             localIndex const kf0,
                             array1d< real64 > & nodalArea ) const;

  void computeNodalTractionDOFs( DomainPartition const & domain, DofManager & dofManager ) const;

  virtual void setConstitutiveNames( ElementSubRegionBase & subRegion ) const override;

  struct viewKeyStruct : ContactSolverBase::viewKeyStruct
  {
    constexpr static char const * rotationMatrixString() { return "rotationMatrix"; }
    constexpr static char const * contactEnforcementMethodString() { return "contactEnforcementMethod"; }

    constexpr static char const * normalDisplacementToleranceString() { return "normalDisplacementTolerance"; }
    constexpr static char const * normalTractionToleranceString() { return "normalTractionTolerance"; }
    constexpr static char const * slidingToleranceString() { return "slidingTolerance"; }
  };

  ContactEnforcementMethod m_contactEnforcementMethod;

};

ENUM_STRINGS( SolidMechanicsConformingFractures::ContactEnforcementMethod, "Penalty", "NodalLagrangeMultiplier", "FaceLagrangeMultiplier" );


} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSCONFORMINGFRACTURES_HPP_ */