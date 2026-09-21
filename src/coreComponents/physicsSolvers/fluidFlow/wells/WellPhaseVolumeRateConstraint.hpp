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
 * @file WellPhaseVolumeRateConstraint.hpp
 */


#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLPHASEVOLUMERATECONSTRAINT_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLPHASEVOLUMERATECONSTRAINT_HPP

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
}

/**
 * @class PhaseVolumeRateConstraint
 * @brief This class describes a phase rate constraint used to control a well of WellConstraintType type (Injection or Production).
 */

class PhaseVolumeRateConstraint : public WellConstraintBase
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
  explicit PhaseVolumeRateConstraint( string const & name, dataRepository::Group * const parent );

  /**
   * @brief Default destructor.
   */
  ~PhaseVolumeRateConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  PhaseVolumeRateConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  PhaseVolumeRateConstraint( PhaseVolumeRateConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  PhaseVolumeRateConstraint( PhaseVolumeRateConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  PhaseVolumeRateConstraint & operator=( PhaseVolumeRateConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  PhaseVolumeRateConstraint & operator=( PhaseVolumeRateConstraint && ) = delete;

  ///@}

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new Constraint object through the object catalog.
   */
  static string catalogName()
  {
    return "PhaseVolumeRateConstraint";
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
   * @brief Set the target phase name
   * @param[in] phaseName the target phase name
   */
  void setPhaseName( const string & phaseName ) { m_phaseName = phaseName; }

  /**
   * @brief Get the target phase index
   * @return the target phase index
   */
  const localIndex & getPhaseIndex() const { return m_phaseIndex; }

  /**
   * @brief Set the target phase index
   * @param[in] phaseIndex the target phase index
   */
  void setPhaseIndex( const localIndex & phaseIndex ) { m_phaseIndex = phaseIndex; }
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

template< typename T >
void PhaseVolumeRateConstraint::validatePhaseType( T const & fluidModel )
{
  // Find target phase index for phase rate constraint
  m_phaseIndex = getPhaseIndexFromFluidModel( fluidModel, this->template getReference< string >( viewKeyStruct::phaseNameString()));

  GEOS_THROW_IF( m_phaseIndex == -1,
                 "PhaseVolumeRateConstraint " << this->template getReference< string >( viewKeyStruct::phaseNameString())   <<
                 ": Invalid phase type for simulation fluid model",
                 InputError );
}

} //namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLPHASEVOLUMERATECONSTRAINT_HPP
