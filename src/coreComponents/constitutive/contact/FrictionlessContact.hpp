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
 * @file FrictionlessContact.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_CONTACT_FRICTIONLESSCONTACT_HPP_
#define GEOSX_CONSTITUTIVE_CONTACT_FRICTIONLESSCONTACT_HPP_

#include "constitutive/contact/ContactBase.hpp"

namespace geosx
{

namespace constitutive
{

/**
 * @class FrictionlessContactUpdates
 *
 * This class is used for in-kernel contact relation updates
 */
class FrictionlessContactUpdates : public ContactBaseUpdates
{
public:

  FrictionlessContactUpdates( real64 const & penaltyStiffness,
                              real64 const & shearStiffness,
                              TableFunction const & apertureTable )
    : ContactBaseUpdates( penaltyStiffness, shearStiffness, apertureTable )
  {}

  /// Default copy constructor
  FrictionlessContactUpdates( FrictionlessContactUpdates const & ) = default;

  /// Default move constructor
  FrictionlessContactUpdates( FrictionlessContactUpdates && ) = default;

  /// Deleted default constructor
  FrictionlessContactUpdates() = default;

  /// Deleted copy assignment operator
  FrictionlessContactUpdates & operator=( FrictionlessContactUpdates const & ) = delete;

  /// Deleted move assignment operator
  FrictionlessContactUpdates & operator=( FrictionlessContactUpdates && ) =  delete;

  /**
   * @brief Evaluate the traction vector and its derivatives wrt to pressure and jump
   * @param[in] dispJump the displacement jump
   * @param[in] tractionVector the traction vector
   * @param[out] dTractionVector_dJump the derivative of the traction vector wrt displacement jump
   */
  GEOSX_HOST_DEVICE
  inline
  virtual void computeTraction( localIndex const k,
                                arraySlice1d< real64 const > const & oldDispJump,
                                arraySlice1d< real64 const > const & dispJump,
                                arraySlice1d< real64 > const & tractionVector,
                                arraySlice2d< real64 > const & dTractionVector_dJump ) const override final;

  /**
   * @brief Evaluate the limit tangential traction norm and return the derivative wrt normal traction
   * @param[in] normalTraction the normal traction
   * @param[out] dLimitTangentialTractionNorm_dTraction the derivative of the limit tangential traction norm wrt normal traction
   * @return the limit tangential traction norm
   */
  GEOSX_HOST_DEVICE
  inline
  virtual real64 computeLimitTangentialTractionNorm( real64 const & normalTraction,
                                                     real64 & dLimitTangentialTractionNorm_dTraction ) const override final
  { GEOSX_UNUSED_VAR( normalTraction, dLimitTangentialTractionNorm_dTraction ); return 0.0; }

private:
};


/**
 * @class FrictionlessContact
 *
 * This class serves as the interface for implementing contact enforcement constitutive relations.
 * This does not include the actual enforcement algorithm, but only the constitutive relations that
 * govern the behavior of the contact. So things like penalty, or friction, or kinematic constraint.
 */
class FrictionlessContact : public ContactBase
{
public:

  /**
   * @brief The standard data repository constructor
   * @param name The name of the relation in the data repository
   * @param parent The name of the parent Group that holds this relation object.
   */
  FrictionlessContact( string const & name,
                       Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~FrictionlessContact() override;

  /**
   * @return A string that is used to register/lookup this class in the registry
   */
  static string catalogName() { return "FrictionlessContact"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = FrictionlessContactUpdates;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  /**
   * @struct Structure to hold scoped key names
   */
  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {};

protected:

};

GEOSX_HOST_DEVICE
void FrictionlessContactUpdates::computeTraction( localIndex const k,
                                                  arraySlice1d< real64 const > const & oldDispJump,
                                                  arraySlice1d< real64 const > const & dispJump,
                                                  arraySlice1d< real64 > const & tractionVector,
                                                  arraySlice2d< real64 > const & dTractionVector_dJump ) const
{
  GEOSX_UNUSED_VAR( k, oldDispJump );

  tractionVector[0] = dispJump[0] >= 0 ? 0.0 : m_penaltyStiffness * dispJump[0];
  tractionVector[1] = 0.0;
  tractionVector[2] = 0.0;

  LvArray::forValuesInSlice( dTractionVector_dJump, []( real64 & val ){ val = 0.0; } );
  dTractionVector_dJump( 0, 0 ) = dispJump[0] >=0 ? 0.0 : m_penaltyStiffness;
}


} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_CONTACT_FRICTIONLESSCONTACT_HPP_ */
