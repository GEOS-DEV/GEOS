/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 *  @file RateAndStateFriction.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CONTACT_RATEANDSTATEFRICTION_HPP_
#define GEOS_CONSTITUTIVE_CONTACT_RATEANDSTATEFRICTION_HPP_

#include "FrictionBase.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class RateAndStateFrictionUpdates
 *
 * This class is used for in-kernel contact relation updates
 */
class RateAndStateFrictionUpdates : public FrictionBaseUpdates
{
public:
  RateAndStateFrictionUpdates( real64 const displacementJumpThreshold,
                               arrayView1d< real64 > const frictionCoefficient )
    : FrictionBaseUpdates( displacementJumpThreshold ),
    m_frictionCoefficient( frictionCoefficient )
  {}

  /// Default copy constructor
  RateAndStateFrictionUpdates( RateAndStateFrictionUpdates const & ) = default;

  /// Default move constructor
  RateAndStateFrictionUpdates( RateAndStateFrictionUpdates && ) = default;

  /// Deleted default constructor
  RateAndStateFrictionUpdates() = delete;

  /// Deleted copy assignment operator
  RateAndStateFrictionUpdates & operator=( RateAndStateFrictionUpdates const & ) = delete;

  /// Deleted move assignment operator
  RateAndStateFrictionUpdates & operator=( RateAndStateFrictionUpdates && ) =  delete;

  GEOS_HOST_DEVICE
  inline
  virtual void updateFractureState( localIndex const k,
                                    arraySlice1d< real64 const > const & dispJump,
                                    arraySlice1d< real64 const > const & tractionVector,
                                    integer & fractureState ) const override final;

private:
  /// The friction coefficient 
  arrayView1d< real64 > m_frictionCoefficient;

  /// Rate and State coefficient a  
  arrayView1d< real64 const > m_a;
  
  /// Rate and State coefficient b 
  arrayView1d< real64 const > m_b;
  
  /// Rate and State reference velocity
  arrayView1d< real64 const > m_V0;
};


/**
 * @class RateAndStateFriction
 *
 * Class to provide a RateAndStateFriction friction model.
 */
class RateAndStateFriction : public FrictionBase
{
public:

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  RateAndStateFriction( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~RateAndStateFriction() override;

  static string catalogName() { return "RateAndStateFriction"; }

  virtual string getCatalogName() const override { return catalogName(); }

  ///@}

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override final;

  /**
   * @brief Const accessor for friction angle
   * @return A const reference to arrayView1d<real64 const> containing the
   *         friction coefficient (at every element).
   */
  real64 const & frictionCoefficient() const { return m_frictionCoefficient; }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = RateAndStateFrictionUpdates;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

private:

  virtual void postInputInitialization() override;

  /// The friction coefficient for each upper level dimension (i.e. cell) of *this
  array1d< real64 > m_frictionCoefficient;

  /// Rate and State coefficient a  
  array1d< real64 > m_a;
  
  /// Rate and State coefficient b 
  array1d< real64 > m_b;
  
  /// Rate and State reference velocity
  array1d< real64 > m_V0;

/**
 * @struct Set of "char const *" and keys for data specified in this class.
 */
  struct viewKeyStruct : public FrictionBase::viewKeyStruct
  {
    /// string/key for friction coefficient
    static constexpr char const * frictionCoefficientString() { return "frictionCoefficient"; }
  };

};

GEOS_HOST_DEVICE
inline void RateAndStateFrictionUpdates::updateFractureState( localIndex const k,
                                                              arraySlice1d< real64 const > const & dispJump,
                                                              arraySlice1d< real64 const > const & tractionVector,
                                                              integer & fractureState ) const
{
  using namespace fields::contact;

  if( dispJump[0] >  -m_displacementJumpThreshold )
  {
    fractureState = FractureState::Open;
  } else
  {
    fractureState = FractureState::Slip;
  }
}

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_CONTACT_RATEANDSTATEFRICTION_HPP_ */
