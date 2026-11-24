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
 * @file WellLiquidRateConstraint.hpp
 */


#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLLIQUIDRATECONSTRAINT_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLLIQUIDRATECONSTRAINT_HPP

#include "common/format/EnumStrings.hpp"

#include "functions/TableFunction.hpp"
#include "WellConstraintsBase.hpp"

namespace geos
{


/**
 * @class LiquidRateConstraint
 * @brief This class describes a Liquid rate constraint used to control of type WellConstraintType
 */


class LiquidRateConstraint : public WellConstraintBase
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
  explicit LiquidRateConstraint( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~LiquidRateConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  LiquidRateConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  LiquidRateConstraint( LiquidRateConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  LiquidRateConstraint( LiquidRateConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  LiquidRateConstraint & operator=( LiquidRateConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  LiquidRateConstraint & operator=( LiquidRateConstraint && ) = delete;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new Constraint object through the object catalog.
   */
  static string catalogName()
  {
    return "LiquidRateConstraint";
  }
  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{
  /**
   * @brief Get the target phase name
   * @return the target phase name
   */
  const string_array & getPhaseNames() const { return m_phaseNames; }

  /**
   * @brief Set phases associated with liquid constraint
   * @param array of phase names
   */
  void setPhaseNames( const string_array & phaseNames )  { m_phaseNames=phaseNames; }

  /**
   * @brief Get the phase indices
   * @return array of phase indices
   */
  const array1d< integer > & getPhaseIndices() const { return m_phaseIndices; }

  ///@}
  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
    /// String key for the liquid rate
    static constexpr char const * liquidRateString() { return "liquidRate"; }
    /// String key for the phases names
    static constexpr char const * phaseNamesString() { return "phaseNames"; }
  };

  // Temp interface - tjb
  virtual ConstraintTypeId getControl() const override { return ConstraintTypeId::LIQUIDRATE; };

  virtual bool checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const override;

  /**
   * @brief Validate Liquid type is consistent with fluidmodel
   */
  template< typename T >
  void  validateLiquidType( T const & fluidModel );

protected:

  virtual void postInputInitialization() override;

protected:

  /// Name of the targeted phase
  string_array m_phaseNames;
  ///Indices of the phases defining the fluid
  array1d< integer > m_phaseIndices;

};

template< typename T >
void LiquidRateConstraint::validateLiquidType( T const & fluidModel )
{
  m_phaseIndices.resize( m_phaseNames.size());
  for( size_t ip =0; ip<m_phaseNames.size(); ip++ )
  {
    integer phaseIndex = fluidModel.getPhaseIndex( m_phaseNames[ip] );
    GEOS_THROW_IF( phaseIndex == -1,
                   "LiquidProductionConstraint " <<  viewKeyStruct::liquidRateString()    <<
                   ": Invalid Liquid type for simulation fluid model " << m_phaseNames[ip],
                   InputError );
    m_phaseIndices[ip]=phaseIndex;
  }
}

} //namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLLiquidRateConstraint_HPP
