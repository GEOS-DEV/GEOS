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

/*
 * @file WellPhaseRateConstraints.hpp
 */


#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLPHASERATECONSTRAINTS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLPHASERATECONSTRAINTS_HPP

#include "common/format/EnumStrings.hpp"
#include "dataRepository/Group.hpp"
#include "functions/TableFunction.hpp"
#include "WellConstraintsBase.hpp"

namespace geos
{
namespace dataRepository
{
namespace keys
{
static constexpr auto phaseProductionConstraint = "PhaseProductionConstraint";
static constexpr auto phaseInjectionConstraint  =  "PhaseInjectionConstraint";
}
}

namespace constraintViewStruct
{

struct phaseConstraintKey
{
  /// String key for the well target phase rate
  static constexpr char const * phaseRateString() { return "phaseRate"; }
  /// String key for the well target phase name
  static constexpr char const * phaseNameString() { return "phaseName"; }
};

}
template< typename T >
localIndex getPhaseIndex( T const & fluidModel, std::string const & inputPhase )
{
  localIndex phaseIndex=-1;
  // Find target phase index for phase rate constraint
  for( integer ip = 0; ip < fluidModel.numFluidPhases(); ++ip )
  {
    if( fluidModel.phaseNames()[ip] == inputPhase )
    {
      phaseIndex = ip;
    }
  }
  return phaseIndex;
}

/**
 * @class PhaseProductionConstraint
 * @brief This class describes a phase rate constraint used to control a production well.
 */

class PhaseProductionConstraint : public WellConstraintBase
{
public:


  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for WellControls Objects.
   * @param[in] name the name of this instantiation of WellControls in the repository
   * @param[in] parent the parent group of this instantiation of WellControls
   */
  explicit PhaseProductionConstraint( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~PhaseProductionConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  PhaseProductionConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  PhaseProductionConstraint( PhaseProductionConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  PhaseProductionConstraint( PhaseProductionConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  PhaseProductionConstraint & operator=( PhaseProductionConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  PhaseProductionConstraint & operator=( PhaseProductionConstraint && ) = delete;

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  // Temp interface - tjb
  virtual WellControls::Control getControl() const override { return WellControls::Control::PHASEVOLRATE; };
  /**
   * @brief Get name of constraint
   * @return constraint key
   */
  virtual std::string getConstraintKey( ) const override { return dataRepository::keys::phaseProductionConstraint; };


  virtual bool checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const override;
  ///@}

  // Phase constraint defintion keys
  constraintViewStruct::phaseConstraintKey viewKeysPhaseConstraint;

  // Surface condition definition keyes
  constraintViewStruct::surfaceConditionsKey viewKeysSurfaceCondtions;

  /**
   * @brief Validate phase type is consistent with fluidmodel
   */
  template< typename T > void validatePhaseType( T const & fluidModel );
protected:

  virtual void postInputInitialization() override;

private:
  /// Name of the targeted phase
  string m_phaseName;

  /// Index of the target phase, used to impose the phase rate constraint
  localIndex m_phaseIndex;

  /// Flag to decide whether rates are controlled at rates or surface conditions
  integer m_useSurfaceConditions;

  /// Surface pressure
  real64 m_surfacePres;

  /// Surface temperature
  real64 m_surfaceTemp;
};

template< typename T >
void PhaseProductionConstraint::validatePhaseType( T const & fluidModel )
{
  // Find target phase index for phase rate constraint
  m_phaseIndex = getPhaseIndex( fluidModel, getReference< string >( constraintViewStruct::phaseConstraintKey::phaseNameString()));

  GEOS_THROW_IF( m_phaseIndex == -1,
                 "PhaseProductionConstraint " << getReference< string >( constraintViewStruct::phaseConstraintKey::phaseNameString())   <<
                 ": Invalid phase type for simulation fluid model",
                 InputError );
}

/**
 * @class PhaseInjectionConstraint
 * @brief This class describes a phase rate constraint used to control a injection well.
 */

class PhaseInjectionConstraint : public WellConstraintBase
{
public:


  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for WellControls Objects.
   * @param[in] name the name of this instantiation of WellControls in the repository
   * @param[in] parent the parent group of this instantiation of WellControls
   */
  explicit PhaseInjectionConstraint( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~PhaseInjectionConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  PhaseInjectionConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  PhaseInjectionConstraint( PhaseInjectionConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  PhaseInjectionConstraint( PhaseInjectionConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  PhaseInjectionConstraint & operator=( PhaseInjectionConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  PhaseInjectionConstraint & operator=( PhaseInjectionConstraint && ) = delete;

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  // Temp interface - tjb
  virtual WellControls::Control getControl() const override { return WellControls::Control::PHASEVOLRATE; };

  /**
   * @brief Get name of constraint
   * @return constraint key
   */
  virtual std::string getConstraintKey( ) const override { return dataRepository::keys::phaseInjectionConstraint; };
  ///@}

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */

  // Phase constraint defintion keys
  constraintViewStruct::phaseConstraintKey viewKeysPhaseConstraint;

  // Surface condition definition keys
  constraintViewStruct::surfaceConditionsKey viewKeysSurfaceCondtions;

  // Injection stream definition keys
  constraintViewStruct::injectionStreamKey viewKeysInjectionStream;

  /**
   * @brief Validate phase type is consistent with fluidmodel
   */
  template< typename T > void validatePhaseType( T const & fluidModel );
  ///@}

  virtual bool checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const override;
protected:

  virtual void postInputInitialization() override;

private:
  /// Name of the targeted phase
  string m_phaseName;

  /// Index of the target phase, used to impose the phase rate constraint
  localIndex m_phaseIndex;

  /// Vector with global component fractions at the injector
  array1d< real64 > m_injectionStream;

  /// Temperature at the injector
  real64 m_injectionTemperature;

  /// Flag to decide whether rates are controlled at rates or surface conditions
  integer m_useSurfaceConditions;

  /// Surface pressure
  real64 m_surfacePres;

  /// Surface temperature
  real64 m_surfaceTemp;


};


template< typename T >
void PhaseInjectionConstraint::validatePhaseType( T const & fluidModel )
{
  // Find target phase index for phase rate constraint
  m_phaseIndex = getPhaseIndex( fluidModel, getReference< string >( constraintViewStruct::phaseConstraintKey::phaseNameString()));

  GEOS_THROW_IF( m_phaseIndex == -1,
                 "PhaseInjectionConstraint " << getReference< string >( constraintViewStruct::phaseConstraintKey::phaseNameString())   <<
                 ": Invalid phase type for simulation fluid model",
                 InputError );
}

} //namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONSTRAINT_HPP
