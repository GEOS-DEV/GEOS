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

/*
 *  @file ContactSolverBase.hpp
 *
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_CONTACTSOLVERBASE_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_CONTACTSOLVERBASE_HPP_

#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/contact/ContactFields.hpp"

namespace geos
{

class ContactSolverBase : public SolidMechanicsLagrangianFEM
{
public:
  ContactSolverBase( const string & name,
                     Group * const parent );

  ~ContactSolverBase() override = default;

  virtual void registerDataOnMesh( dataRepository::Group & meshBodies ) override;

  virtual real64
  explicitStep( real64 const & time_n,
                real64 const & dt,
                integer const cycleNumber,
                DomainPartition & domain ) override final;

  virtual void
  applyBoundaryConditions( real64 const time,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override;

  string const & getContactRelationName() const { return m_contactRelationName; }

  string const & getFractureRegionName() const { return m_fractureRegionName; }

  void outputConfigurationStatistics( DomainPartition const & domain ) const override final;

//  SolidMechanicsLagrangianFEM * getSolidSolver() { return m_solidSolver; }

  void setSolidSolverDofFlags( bool const flag ) { m_setupSolidSolverDofs = flag; }

  void synchronizeFractureState( DomainPartition & domain ) const;

protected:

  virtual void postProcessInput() override;

  void computeFractureStateStatistics( MeshLevel const & mesh,
                                       globalIndex & numStick,
                                       globalIndex & numSlip,
                                       globalIndex & numOpen ) const;

  GEOS_HOST_DEVICE
  inline
  static bool compareFractureStates( integer const state0,
                                     integer const state1 )
  {
    return state0 == state1
           || ( state0 == fields::contact::FractureState::NewSlip && state1 == fields::contact::FractureState::Slip )
           || ( state0 == fields::contact::FractureState::Slip && state1 == fields::contact::FractureState::NewSlip );
  }

//  /// Solid mechanics solver name
//  string m_solidSolverName;

  /// fracture region name
  string m_fractureRegionName;

//  /// pointer to the solid mechanics solver
//  SolidMechanicsLagrangianFEM * m_solidSolver;

  /// contact relation name string
  string m_contactRelationName;

  ///
  bool m_setupSolidSolverDofs;

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
//    constexpr static char const * solidSolverNameString() { return "solidSolverName"; }

    constexpr static char const * contactRelationNameString() { return "contactRelationName"; }

    constexpr static char const * fractureRegionNameString() { return "fractureRegionName"; }

    constexpr static char const * fractureStateString() { return "fractureState"; }

    constexpr static char const * oldFractureStateString() { return "oldFractureState"; }

    constexpr static char const * initialFractureStateString() { return "initialFractureState"; }
  };
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_CONTACT_CONTACTSOLVERBASE_HPP_ */
