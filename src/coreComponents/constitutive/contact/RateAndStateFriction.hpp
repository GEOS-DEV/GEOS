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

  
class KernelWrapper : public FrictionBaseUpdates
{
public:
  KernelWrapper( real64 const displacementJumpThreshold,
                 arrayView1d< real64 > frictionCoefficient, 
                 arrayView1d< real64 const > a, 
                 arrayView1d< real64 const > b,
                 arrayView1d< real64 const > Dc,
                 arrayView1d< real64 const > V0,
                 arrayView1d< real64 const > mu0 )
    : FrictionBaseUpdates( displacementJumpThreshold ),
    m_frictionCoefficient( frictionCoefficient ),
    m_a( a ),
    m_b( b ),
    m_Dc( Dc ),
    m_V0( V0 ),
    m_mu0( mu0 )
  {}

  /// Default copy constructor
  KernelWrapper( KernelWrapper const & ) = default;

  /// Default move constructor
  KernelWrapper( KernelWrapper && ) = default;

  /// Deleted default constructor
  KernelWrapper() = delete;

  /// Deleted copy assignment operator
  KernelWrapper & operator=( KernelWrapper const & ) = delete;

  /// Deleted move assignment operator
  KernelWrapper & operator=( KernelWrapper && ) =  delete;

  GEOS_HOST_DEVICE
  inline
  virtual void updateFractureState( localIndex const k,
                                    arraySlice1d< real64 const > const & dispJump,
                                    arraySlice1d< real64 const > const & tractionVector,
                                    integer & fractureState ) const override final;

GEOS_HOST_DEVICE
inline std::tuple< real64, real64, real64 > 
computeShearTraction( real64 const normalTraction, 
                      real64 const slipRate, 
                      real64 const stateVariable ) const;

GEOS_HOST_DEVICE
inline real64 frictionCoefficient( localIndex const k, 
                                   real64 const slipRate, 
                                   real64 const stateVariable ) const;

GEOS_HOST_DEVICE
inline real64 dfrictionCoefficient_dSlipRate( localIndex const k, 
                                              real64 const slipRate, 
                                              real64 const stateVariable ) const;

GEOS_HOST_DEVICE
inline real64 dfrictionCoefficient_dStateVariable( localIndex const k, 
                                                   real64 const slipRate, 
                                                   real64 const stateVariable) const;                                                                                                              

GEOS_HOST_DEVICE
inline real64 dStateVariabledT( localIndex const k, 
                                real64 const slipRate, 
                                real64 const stateVariable) const;

GEOS_HOST_DEVICE
inline real64 dStateVariabledT_dStateVariable( localIndex const k, 
                                               real64 const slipRate, 
                                               real64 const stateVariable) const;

GEOS_HOST_DEVICE
inline real64 dStateVariabledT_dSlipRate( localIndex const k, 
                                          real64 const slipRate, 
                                          real64 const stateVariable) const;                        
private:
  /// The friction coefficient 
  arrayView1d< real64 > m_frictionCoefficient;

  /// Rate and State coefficient a  
  arrayView1d< real64 const > m_a;
  
  /// Rate and State coefficient b 
  arrayView1d< real64 const > m_b;

   /// Rate and State characteristic length
  arrayView1d< real64 const > m_Dc;
  
  /// Rate and State reference velocity
  arrayView1d< real64 const > m_V0;

  /// Rate and State reference friction coefficient
  arrayView1d< real64 const > m_mu0;
};


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

  /// Rate and State characteristic length
  array1d< real64 > m_Dc;
  
  /// Rate and State reference velocity
  array1d< real64 > m_V0;

  /// Rate and State reference friction coefficient
  array1d< real64 > m_mu0;

/**
 * @struct Set of "char const *" and keys for data specified in this class.
 */
  struct viewKeyStruct : public FrictionBase::viewKeyStruct
  {
    /// string/key for friction coefficient
    static constexpr char const * frictionCoefficientString() { return "frictionCoefficient"; }
    /// string/key for friction coefficient
    static constexpr char const * aCoefficientString() { return "a"; }
    /// string/key for friction coefficient
    static constexpr char const * bCoefficientString() { return "b"; }
    /// string/key for friction coefficient
    static constexpr char const * DcCoefficientString() { return "Dc"; }
    /// string/key for friction coefficient
    static constexpr char const * referenceVelocityString() { return "referenceVelocity"; }
    /// string/key for friction coefficient
    static constexpr char const * referenceFrictionCoefficientString() { return "referenceFrictionCoefficient"; }
  };

};

GEOS_HOST_DEVICE
inline void RateAndStateFriction::KernelWrapper::updateFractureState( localIndex const k,
                                                                      arraySlice1d< real64 const > const & dispJump,
                                                                      arraySlice1d< real64 const > const & tractionVector,
                                                                      integer & fractureState ) const
{

  GEOS_UNUSED_VAR(tractionVector, k);
  using namespace fields::contact;

  if( dispJump[0] >  -m_displacementJumpThreshold )
  {
    fractureState = FractureState::Open;
  } else
  {
    fractureState = FractureState::Slip;
  }
}

GEOS_HOST_DEVICE
inline std::tuple< real64, real64, real64 > 
RateAndStateFriction::KernelWrapper::computeShearTraction( real64 const normalTraction, 
                                                           real64 const slipRate, 
                                                           real64 const stateVariable ) const 
{
  GEOS_UNUSED_VAR(normalTraction, slipRate, stateVariable);
  real64 shearTraction = 0.0;
  real64 dTauFriction[2] = {0.0, 0.0}; 

  return std::make_tuple( shearTraction, dTauFriction[0], dTauFriction[1] );
}

GEOS_HOST_DEVICE
inline real64 RateAndStateFriction::KernelWrapper::frictionCoefficient( localIndex const k, 
                                                                        real64 const slipRate, 
                                                                        real64 const stateVariable ) const
{
  real64 const arg = ( slipRate / (2 * m_V0[k]) ) * LvArray::math::exp( stateVariable / m_a[k] );
  real64 const frictionCoefficient = m_a[k] * LvArray::math::asin( arg ); //TODO: change!! asinh
  
  m_frictionCoefficient[k] = frictionCoefficient;

  return frictionCoefficient;
}         

GEOS_HOST_DEVICE
inline real64 RateAndStateFriction::KernelWrapper::dfrictionCoefficient_dSlipRate( localIndex const k, 
                                                                                   real64 const slipRate, 
                                                                                   real64 const stateVariable ) const 
{

  real64 const inner_expression = ( slipRate / (2 * m_V0[k]) ) * LvArray::math::exp( stateVariable / m_a[k] ); 

  return (m_a[k] / LvArray::math::sqrt(1 + inner_expression*inner_expression )) * (1 / (2 * m_V0[k])) * LvArray::math::exp( stateVariable / m_a[k] );
}
    
GEOS_HOST_DEVICE
inline real64 RateAndStateFriction::KernelWrapper::dfrictionCoefficient_dStateVariable( localIndex const k, 
                                                   real64 const slipRate, 
                                                   real64 const stateVariable) const
{
        
  real64 const arg = (slipRate / (2 * m_V0[k])) * LvArray::math::exp(stateVariable / m_a[k]);      
        
  return ( slipRate / (2 * m_V0[k])) * LvArray::math::exp(stateVariable / m_a[k]) / LvArray::math::sqrt(arg * arg + 1);
}

GEOS_HOST_DEVICE
inline real64 RateAndStateFriction::KernelWrapper::dStateVariabledT( localIndex const k, 
                                                                     real64 const slipRate, 
                                                                     real64 const stateVariable) const
{
  real64 const mu = frictionCoefficient(k, slipRate, stateVariable);
        
  return -slipRate / m_Dc[k] * (mu - m_mu0[k] + (m_b[k] - m_a[k]) * LvArray::math::log(slipRate / m_V0[k])) ;
}

GEOS_HOST_DEVICE
inline real64 RateAndStateFriction::KernelWrapper::dStateVariabledT_dStateVariable( localIndex const k,
                                                                                    real64 const slipRate, 
                                                                                    real64 const stateVariable) const
{
  return -slipRate / m_Dc[k] * dfrictionCoefficient_dStateVariable( k, slipRate, stateVariable );
}

GEOS_HOST_DEVICE
inline real64 RateAndStateFriction::KernelWrapper::dStateVariabledT_dSlipRate( localIndex const k,
                                                                               real64 const slipRate, 
                                                                               real64 const stateVariable ) const
{        
  real64 const part1 = - 1.0 / m_Dc[k] * (frictionCoefficient(k, slipRate, stateVariable) - m_mu0[k] + (m_b[k] - m_a[k]) * LvArray::math::log(slipRate / m_V0[k]));
        
  real64 const part2 = - slipRate / m_Dc[k] * (dfrictionCoefficient_dSlipRate(k, slipRate, stateVariable) + (m_b[k] - m_a[k]) / slipRate);
        
  return part1 + part2;
}

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_CONTACT_RATEANDSTATEFRICTION_HPP_ */
