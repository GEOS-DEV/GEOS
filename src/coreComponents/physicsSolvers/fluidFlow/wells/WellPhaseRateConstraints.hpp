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
#include "WellConstants.hpp"

namespace geos
{


template< typename T >
localIndex getPhaseIndexFromFluidModel( T const & fluidModel, std::string const & inputPhase )
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
};

/**
 * @class PhaseConstraint
 * @brief This class describes a phase rate constraint used to control a production well.
 */

class PhaseConstraint : public WellConstraintBase
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
  explicit PhaseConstraint( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~PhaseConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  PhaseConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  PhaseConstraint( PhaseConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  PhaseConstraint( PhaseConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  PhaseConstraint & operator=( PhaseConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  PhaseConstraint & operator=( PhaseConstraint && ) = delete;

  ///@}



  /**
   * @name Getters / Setters
   */
  ///@{

  // Temp interface - tjb
  virtual ConstraintTypeId getControl() const override { return ConstraintTypeId::PHASEVOLRATE; };

  /**
   * @brief Get the target phase name
   * @return the target phase name
   */
  const string & getPhaseName() const { return m_phaseName; }

  /**
   * @brief Get the target phase index
   * @return the target phase index
   */
  const localIndex & getPhaseIndex() const { return m_phaseIndex; }

  ///@}

  struct viewKeyStruct
  {
    /// String key for the well target phase rate
    static constexpr char const * phaseRateString() { return "phaseRate"; }
    /// String key for the well target phase name
    static constexpr char const * phaseNameString() { return "phaseName"; }
  };

protected:

  virtual void postInputInitialization() override;

protected:
  /// Name of the targeted phase
  string m_phaseName;

  /// Index of the target phase, used to impose the phase rate constraint
  localIndex m_phaseIndex;

};


/**
 * @class PhaseProductionConstraint
 * @brief This class describes a phase rate constraint used to control a production well.
 */

class PhaseProductionConstraint : public PhaseConstraint
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

  /**
   * @brief Get name of constraint
   * @return constraint key
   */
  virtual std::string getConstraintKey( ) const override { return "PhaseProductionConstraint"; };


  virtual bool checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const override;
  ///@}

  /**
   * @brief Validate phase type is consistent with fluidmodel
   */
  template< typename T > void validatePhaseType( T const & fluidModel );
protected:

  virtual void postInputInitialization() override;

};

template< typename T >
void PhaseProductionConstraint::validatePhaseType( T const & fluidModel )
{
  // Find target phase index for phase rate constraint
  m_phaseIndex = getPhaseIndexFromFluidModel( fluidModel, getReference< string >( PhaseConstraint::viewKeyStruct::phaseNameString()));

  GEOS_THROW_IF( m_phaseIndex == -1,
                 "PhaseProductionConstraint " << getReference< string >( PhaseConstraint::viewKeyStruct::phaseNameString())   <<
                 ": Invalid phase type for simulation fluid model",
                 InputError );
}

/**
 * @class PhaseInjectionConstraint
 * @brief This class describes a phase rate constraint used to control a injection well.
 */

class PhaseInjectionConstraint : public PhaseConstraint
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

  /**
   * @brief Get name of constraint
   * @return constraint key
   */
  virtual std::string getConstraintKey( ) const override { return "PhaseInjectionConstraint"; };
  ///@}


  /**
   * @brief Validate phase type is consistent with fluidmodel
   */
  template< typename T > void validatePhaseType( T const & fluidModel );
  ///@}

  virtual bool checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const override;
protected:

  virtual void postInputInitialization() override;

private:

  /// Vector with global component fractions at the injector
  array1d< real64 > m_injectionStream;

  /// Temperature at the injector
  real64 m_injectionTemperature;

};


template< typename T >
void PhaseInjectionConstraint::validatePhaseType( T const & fluidModel )
{
  // Find target phase index for phase rate constraint
  m_phaseIndex = getPhaseIndexFromFluidModel( fluidModel, getReference< string >( PhaseConstraint::viewKeyStruct::phaseNameString()));

  GEOS_THROW_IF( m_phaseIndex == -1,
                 "PhaseInjectionConstraint " << getReference< string >( PhaseConstraint::viewKeyStruct::phaseNameString())   <<
                 ": Invalid phase type for simulation fluid model",
                 InputError );
}

/**
 * @class PhaseInjectionConstraint
 * @brief This class describes a phase rate constraint used to control a injection well.
 */
template< typename WellConstraintType >
class PhaseConstraint1 : public WellConstraintType
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
  explicit PhaseConstraint1( string const & name, dataRepository::Group * const parent );

  /**
   * @brief Default destructor.
   */
  ~PhaseConstraint1() override;

  /**
   * @brief Deleted default constructor.
   */
  PhaseConstraint1() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  PhaseConstraint1( PhaseConstraint1 const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  PhaseConstraint1( PhaseConstraint1 && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  PhaseConstraint1 & operator=( PhaseConstraint1 const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  PhaseConstraint1 & operator=( PhaseConstraint1 && ) = delete;

  ///@}

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new MultiphasePoromechanics object through the object catalog.
   */
  static string catalogName()
  {
    if constexpr ( std::is_same_v< WellConstraintType, InjectionConstraint > )    // special case
    {
      return "PhaseInjectionConstraint";
    }
    else   // default
    {
      return "PhaseProductionConstraint";
    }
  }

  /**
   * @name Getters / Setters
   */
  ///@{

  // Temp interface - tjb
  virtual ConstraintTypeId getControl() const override { return ConstraintTypeId::PHASEVOLRATE; };

  /**
   * @brief Get the target phase name
   * @return the target phase name
   */
  const string & getPhaseName() const { return m_phaseName; }

  /**
   * @brief Get the target phase index
   * @return the target phase index
   */
  const localIndex & getPhaseIndex() const { return m_phaseIndex; }

  ///@}

  struct viewKeyStruct
  {
    /// String key for the well target phase rate
    static constexpr char const * phaseRateString() { return "phaseRate"; }
    /// String key for the well target phase name
    static constexpr char const * phaseNameString() { return "phaseName"; }
  };

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

};

template< typename WellConstraintType >
template< typename T >
void PhaseConstraint1< WellConstraintType >::validatePhaseType( T const & fluidModel )
{
  // Find target phase index for phase rate constraint
  m_phaseIndex = getPhaseIndexFromFluidModel( fluidModel, this->template getReference< string >( viewKeyStruct::phaseNameString()));

  GEOS_THROW_IF( m_phaseIndex == -1,
                 "PhaseConstraint " << this->template getReference< string >( viewKeyStruct::phaseNameString())   <<
                 ": Invalid phase type for simulation fluid model",
                 InputError );
}

} //namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONSTRAINT_HPP
