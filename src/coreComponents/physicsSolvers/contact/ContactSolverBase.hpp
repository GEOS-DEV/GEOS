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

#ifndef GEOSX_PHYSICSSOLVERS_CONTACT_CONTACTSOLVERBASE_HPP_
#define GEOSX_PHYSICSSOLVERS_CONTACT_CONTACTSOLVERBASE_HPP_

#include "physicsSolvers/SolverBase.hpp"

namespace geosx
{
class SolidMechanicsLagrangianFEM;

class ContactSolverBase : public SolverBase
{
public:
  ContactSolverBase( const string & name,
                     Group * const parent );

  ~ContactSolverBase() override = default;

  virtual void registerDataOnMesh( dataRepository::Group & meshBodies ) override;

  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             int const cycleNumber,
                             DomainPartition & domain ) override;

  virtual real64
  explicitStep( real64 const & time_n,
                real64 const & dt,
                integer const cycleNumber,
                DomainPartition & domain ) override final;

  // virtual bool updateConfiguration( DomainPartition & domain ) override;

  virtual void
  applyBoundaryConditions( real64 const time,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override;

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static char const * solidSolverNameString() { return "solidSolverName"; }

    constexpr static char const * contactRelationNameString() { return "contactRelationName"; }

    constexpr static char const * fractureRegionNameString() { return "fractureRegionName"; }

    constexpr static char const * fractureStateString() { return "fractureState"; }

    constexpr static char const * oldFractureStateString() { return "oldFractureState"; }
  };

  string const & getContactRelationName() const { return m_contactRelationName; }

  string const & getFractureRegionName() const { return m_fractureRegionName; }

  void outputConfigurationStatistics( DomainPartition const & domain ) const override final;

  /**
   * @struct FractureState
   *
   * A struct for the fracture states
   */
  struct FractureState
  {
    enum State : integer
    {
      Stick = 0, ///< element is closed: no jump across the discontinuity.
      Slip = 1, ///< element is sliding: no normal jump across the discontinuity, but sliding is allowed.
      NewSlip = 2, ///< element just starts sliding: no normal jump across the discontinuity, but sliding is allowed.
      Open = 3 ///< element is open: no constraints are imposed.
    };
  };

protected:

  virtual void postProcessInput() override;

  void computeFractureStateStatistics( MeshLevel const & mesh,
                                       globalIndex & numStick,
                                       globalIndex & numSlip,
                                       globalIndex & numOpen ) const;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static bool compareFractureStates( integer const state0,
                                     integer const state1 )
  {
    return state0 == state1
           || ( state0 == FractureState::NewSlip && state1 == FractureState::Slip )
           || ( state0 == FractureState::Slip && state1 == FractureState::NewSlip );
  }

  void initializeFractureState( SurfaceElementSubRegion & subRegion,
                                string const & fieldName ) const;

  void synchronizeFractureState( DomainPartition & domain ) const;

  /// Solid mechanics solver name
  string m_solidSolverName;

  /// fracture region name
  string m_fractureRegionName;

  /// pointer to the solid mechanics solver
  SolidMechanicsLagrangianFEM * m_solidSolver;

  /// contact relation name string
  string m_contactRelationName;
};

ENUM_STRINGS( ContactSolverBase::FractureState::State, "stick", "slip", "new_slip", "open" );


} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_CONTACT_CONTACTSOLVERBASE_HPP_ */
