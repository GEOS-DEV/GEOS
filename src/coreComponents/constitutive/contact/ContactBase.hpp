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
 * @file ContactBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_CONTACT_CONTACTBASE_HPP_
#define GEOSX_CONSTITUTIVE_CONTACT_CONTACTBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"

namespace geosx
{

namespace constitutive
{

/**
 * @class ContactBaseUpdates
 *
 * This class is used for in-kernel contact relation updates
 */
class ContactBaseUpdates
{
public:

  /// Default copy constructor
  ContactBaseUpdates( ContactBaseUpdates const & ) = default;

  /// Default move constructor
  ContactBaseUpdates( ContactBaseUpdates && ) = default;

  /// Deleted default constructor
  ContactBaseUpdates() = default;

  /// Deleted copy assignment operator
  ContactBaseUpdates & operator=( ContactBaseUpdates const & ) = delete;

  /// Deleted move assignment operator
  ContactBaseUpdates & operator=( ContactBaseUpdates && ) =  delete;

  /**
   * @brief Evaluate the effective aperture, and its derivative wrt aperture
   * @param[in] aperture the model aperture/gap
   * @param[out] dEffectiveAperture_dAperture the derivative of the effective aperture wrt aperture
   * @return An effective physical aperture that is always > 0
   */
  GEOSX_HOST_DEVICE
  inline
  virtual real64 computeEffectiveAperture( real64 const aperture,
                                           real64 & dEffectiveAperture_dAperture ) const
  { GEOSX_UNUSED_VAR( aperture, dEffectiveAperture_dAperture ); return 0.0; }


  /**
   * @brief Evaluate the traction vector and its derivatives wrt to pressure and jump
   * @param[in] dispJump the displacement jump
   * @param[in] tractionVector the traction vector
   * @param[out] dTractionVector_dJump the derivative of the traction vector wrt displacement jump
   */
  GEOSX_HOST_DEVICE
  inline
  virtual void computeTraction( arraySlice1d< real64 const > const & dispJump,
                                arraySlice1d< real64 > const & tractionVector,
                                arraySlice2d< real64 > const & dTractionVector_dJump ) const
  { GEOSX_UNUSED_VAR( dispJump, tractionVector, dTractionVector_dJump ); }

  /**
   * @brief Update the traction with the pressure term
   * @param[in] pressure the pressure term
   * @param[in] isOpen a flag specifying whether the fracture is open or closed
   * @param[inout] traction the current tractionVector
   * @param[out] dTraction_dPressure the derivative of the fist component of traction wrt pressure
   * @return the updated traction
   */
  GEOSX_HOST_DEVICE
  inline
  virtual void addPressureToTraction( real64 const & pressure,
                                      bool const isOpen,
                                      arraySlice1d< real64 >const & tractionVector,
                                      real64 & dTraction_dPressure ) const
  { GEOSX_UNUSED_VAR( pressure, isOpen, tractionVector, dTraction_dPressure ); }

  /**
   * @brief Evaluate the limit tangential traction norm and return the derivative wrt normal traction
   * @param[in] normalTraction the normal traction
   * @param[out] dLimitTangentialTractionNorm_dTraction the derivative of the limit tangential traction norm wrt normal traction
   * @return the limit tangential traction norm
   */
  GEOSX_HOST_DEVICE
  inline
  virtual real64 computeLimitTangentialTractionNorm( real64 const & normalTraction,
                                                     real64 & dLimitTangentialTractionNorm_dTraction ) const
  { GEOSX_UNUSED_VAR( normalTraction, dLimitTangentialTractionNorm_dTraction ); return 0.0; }

};


/**
 * @class ContactBase
 *
 * This class serves as the interface for implementing contact enforcement constitutive relations.
 * This does not include the actual enforcement algorithm, but only the constitutive relations that
 * govern the behavior of the contact. So things like penalty, or friction, or kinematic constraint.
 */
class ContactBase : public ConstitutiveBase
{
public:

  /**
   * @brief The standard data repository constructor
   * @param name The name of the relation in the data repository
   * @param parent The name of the parent Group that holds this relation object.
   */
  ContactBase( string const & name,
               Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~ContactBase() override;

  /**
   * @brief accessor for penalty stiffness
   * @return the stiffness
   */
  real64 stiffness() const { return m_penaltyStiffness; }

  /**
   * @struct Structure to hold scoped key names
   */
  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    /// string/key for penalty stiffness
    static constexpr char const * penaltyStiffnessString() { return "penaltyStiffness"; }
  };

protected:

  /// The value of penalty to penetration
  real64 m_penaltyStiffness;

};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_CONTACT_CONTACTBASE_HPP_ */
