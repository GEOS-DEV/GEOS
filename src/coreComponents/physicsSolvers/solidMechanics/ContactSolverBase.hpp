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
 * ContactSolverBase.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_CONTACTSOLVERBASE_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_CONTACTSOLVERBASE_HPP_

#include "physicsSolvers/SolverBase.hpp"

namespace geosx
{
using namespace constitutive;

class SolidMechanicsLagrangianFEM;

class ContactSolverBase : public SolverBase
{
public:
  ContactSolverBase( const string & name,
                     Group * const parent );

  ~ContactSolverBase() override;

  virtual void registerDataOnMesh( dataRepository::Group & meshBodies ) override final;

  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             int const cycleNumber,
                             DomainPartition & domain ) override;

  virtual bool updateConfiguration( DomainPartition & domain ) override;

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static char const * solidSolverNameString() { return "solidSolverName"; }

    constexpr static char const * contactRelationNameString() { return "contactRelationName"; }

    constexpr static char const * dispJumpString() { return "displacementJump"; }

    constexpr static char const * deltaDispJumpString() { return "deltaDisplacementJump"; }

    constexpr static char const * oldDispJumpString()  { return "oldDisplacementJump"; }

    constexpr static char const * fractureRegionNameString() { return "fractureRegionName"; }

    constexpr static char const * tractionString() { return "traction"; }

    constexpr static char const * fractureStateString() { return "fractureState"; }
  
    constexpr static char const * oldFractureStateString() { return "oldFractureState"; }
  };

  string const & getContactRelationName() const { return m_contactRelationName; }

  string const & getFractureRegionName() const { return m_fractureRegionName; }

protected:

  virtual void initializePostInitialConditionsPreSubGroups() override;

  virtual void postProcessInput() override;

  /**
   * @struct FractureState
   *
   * A struct for the fracture states
   */
  struct FractureState
  {
    static constexpr integer STICK = 0;    ///< element is closed: no jump across the discontinuity.
    static constexpr integer SLIP = 1;     ///< element is sliding: no normal jump across the discontinuity, but sliding is allowed.
    static constexpr integer NEW_SLIP = 2; ///< element just starts sliding: no normal jump across the discontinuity, but sliding is
                                           ///< allowed.
    static constexpr integer OPEN = 3;     ///< element is open: no constraints are imposed.
  };

  string fractureStateToString( integer const & state ) const
  {
    string stringState;
    switch( state )
    {
      case FractureState::STICK:
      {
        stringState = "stick";
        break;
      }
      case FractureState::SLIP:
      {
        stringState = "slip";
        break;
      }
      case FractureState::NEW_SLIP:
      {
        stringState = "new_slip";
        break;
      }
      case FractureState::OPEN:
      {
        stringState = "open";
        break;
      }
    }
    return stringState;
  }

  /// Solid mechanics solver name
  string m_solidSolverName;

  /// fracture region name
  string m_fractureRegionName;

  /// pointer to the solid mechanics solver
  SolidMechanicsLagrangianFEM * m_solidSolver;

  /// contact relation name string
  string m_contactRelationName;
};


} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_CONTACTSOLVERBASE_HPP_ */
