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
 * @file Damage.hpp
 * @brief This class overrides the linear elastic constitutive updates to account for a damage field
 *
 * In a phase-field for fracture model, the damage variable affects the Elasticity equation
 * with the degradation of the stresses. Instead of sigma = C : epsilon, we have sigma = g(d)*C:epsilon,
 * where g(d) is the degradation function. This degradation function can either be a quadratic one
 * (set LORENTZ 0) or a quasi-quadratic one (set LORENTZ 1). In general, the quadratic one will give you
 * brittle fracture behaviour. The quasi-quadratic one, combined with linear dissipation, will give you
 * cohesive fracture behaviour, with a user-defined critical stress. If you use quadratic dissipation in
 * your damage solver, set QUADRATIC_DISSIPATION to 1.
 *
 * References:
 *
 * Miehe, Christian; Hofacker, Martina; Welschinger, Fabian. A phase field model for rate-independent crack
 * propagation: Robust algorithmic implementation based on operator splits.
 * Computer Methods in Applied Mechianics and Engineering, v. 199, n. 45-48, p. 2765-2778, 2010
 *
 * Borden, Micheal J., et al. A phase-field description of dynamic brittle fracture.
 * Computer Methods in Applied Mechanics and Engineering, v. 217, p. 77-95, 2012
 *
 * Bourdin, Blaise; Francfort, Gille A.; Marigo, Jean-Jacques. The variational approach to fracture.
 * Journal of Elasticity, v. 91, n. 1-3, p. 5-148, 2008.
 *
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_DAMAGE_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_DAMAGE_HPP_

#include "constitutive/solid/SolidBase.hpp"

namespace geosx
{
namespace constitutive
{

// DAMAGE MODEL UPDATES
//
// NOTE: This model uses the m_newStress array to represent the stress in an
//       elastic, "undamaged" configuration.  We then scale the results
//       by the damage factor whenever the true stress is requested through an update
//       function.  The developer should be very cautious if accessing the stress
//       directly through an arrayView, as it does not represent the true stress.
//
// NOTE: This model is designed to work with phase field implementation where m_damage
//       is updated externally to the material model routine.  A modified implementation
//       would be required to use this directly within a SolidMechanics-only solver, to
//       internally update the damage variable and consistently linearize the system.

/**
 * @class DamageUpdates
 *
 * Class to provide material updates for damaged solids that may be
 * called from a kernel function. The intact solid behavior is that of UDPATE_BASE
 * and this class adds the functions needed to account for damage.
 *
 */
template< typename UPDATE_BASE >
class DamageUpdates : public UPDATE_BASE
{
public:
  template< typename ... PARAMS >
  DamageUpdates( arrayView2d< real64 > const & inputDamage,
                 arrayView2d< real64 > const & inputStrainEnergyDensity,
                 real64 const & inputLengthScale,
                 real64 const & inputCriticalFractureEnergy,
                 real64 const & inputcriticalStrainEnergy,
                 PARAMS && ... baseParams ):
    UPDATE_BASE( std::forward< PARAMS >( baseParams )... ),
    m_damage( inputDamage ),
    m_strainEnergyDensity( inputStrainEnergyDensity ),
    m_lengthScale( inputLengthScale ),
    m_criticalFractureEnergy( inputCriticalFractureEnergy ),
    m_criticalStrainEnergy( inputcriticalStrainEnergy )
  {}

  using DiscretizationOps = typename UPDATE_BASE::DiscretizationOps;

  using UPDATE_BASE::smallStrainNoStateUpdate;
  using UPDATE_BASE::smallStrainUpdate;
  using UPDATE_BASE::smallStrainNoStateUpdate_StressOnly;
  using UPDATE_BASE::smallStrainUpdate_StressOnly;
  using UPDATE_BASE::hypoUpdate;
  using UPDATE_BASE::hypoUpdate_StressOnly;
  using UPDATE_BASE::hyperUpdate;
  using UPDATE_BASE::hyperUpdate_StressOnly;
  using UPDATE_BASE::saveConvergedState;

  using UPDATE_BASE::m_disableInelasticity;


  /**
   * @brief quadratic degradation function, used to reduce the material stiffness when there is damage
   *
   * @param k the element number
   * @param q the quadrature point
   * @return the value of the degradation function at element k, quadrature point q
   */
  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  virtual real64 getDegradationValue( localIndex const k,
                                      localIndex const q ) const
  {
    return (1 - m_damage( k, q ))*(1 - m_damage( k, q ));
  }


  /**
   * @brief derivative of the degratation function
   *
   * @param d the damage value
   * @return the derivative of the degradation function at d
   */
  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  virtual real64 getDegradationDerivative( real64 const d ) const
  {
    return -2*(1 - d);
  }

  /**
   * @brief the second derivative of the degradation function
   *
   * @param d the damage value
   * @return the second derivative of the degradation at d
   */
  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  virtual real64 getDegradationSecondDerivative( real64 const d ) const
  {
    GEOSX_UNUSED_VAR( d );
    return 2.0;
  }

  ///modified smallStrainUpdate to account for the presence of damage
  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( & stress )[6],
                                  DiscretizationOps & stiffness ) const override
  {
    UPDATE_BASE::smallStrainUpdate( k, q, strainIncrement, stress, stiffness );

    if( m_disableInelasticity )
    {
      return;
    }

    real64 factor = getDegradationValue( k, q );
    LvArray::tensorOps::scale< 6 >( stress, factor );
    stiffness.scaleParams( factor );
  }


  /**
   * @defgroup Getter Functions
   *
   * These are getters to access members of the Damage class
   */
  /**@{*/

  // TODO: The code below assumes the strain energy density will never be
  //       evaluated in a non-converged / garbage configuration.
  // the history approach (passing the max of sed over time) is being used
  // to enforce damage irreversibilty

  GEOSX_HOST_DEVICE
  virtual real64 getStrainEnergyDensity( localIndex const k,
                                         localIndex const q ) const override
  {
    real64 const sed = SolidBaseUpdates::getStrainEnergyDensity( k, q );

    if( sed > m_strainEnergyDensity( k, q ) )
    {
      m_strainEnergyDensity( k, q ) = sed;
    }

    return m_strainEnergyDensity( k, q );
  }

  GEOSX_HOST_DEVICE
  real64 getRegularizationLength() const
  {
    return m_lengthScale;
  }

  GEOSX_HOST_DEVICE
  real64 getCriticalFractureEnergy() const
  {
    return m_criticalFractureEnergy;
  }

  GEOSX_HOST_DEVICE
  virtual real64 getEnergyThreshold() const
  {
    #if LORENTZ
    return m_criticalStrainEnergy;
    #else
    return 3*m_criticalFractureEnergy/(16 * m_lengthScale);
    #endif
  }

  /**@}*/
  ///an array with the values of damage at quadrature points
  arrayView2d< real64 > const m_damage;
  ///an array with the values of the strain energy density at quadrature points
  arrayView2d< real64 > const m_strainEnergyDensity;
  ///the phase-field regularization length
  real64 const m_lengthScale;
  ///the critical energy release rate (Gc)
  real64 const m_criticalFractureEnergy;
  ///the critical strain energy denstiy (Psi_c) - threshold for damage initiation
  real64 const m_criticalStrainEnergy;
};



class DamageBase : public SolidBase
{};

/**
 * @class Damage
 *
 * Template class to account for damage effects in a solid. The damage formulation is
 * that of the phase-field approach to fracture (which can also be interpreted as a
 * gradient damage model)
 *
 */
template< typename BASE >
class Damage : public BASE
{
public:

  /// @typedef Alias for LinearElasticIsotropicUpdates
  using KernelWrapper = DamageUpdates< typename BASE::KernelWrapper >;

  Damage( string const & name, dataRepository::Group * const parent );
  virtual ~Damage() override = default;

  static string catalogName() { return string( "Damage" ) + BASE::m_catalogNameString; }
  virtual string getCatalogName() const override { return catalogName(); }

  virtual void postProcessInput() override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;


  KernelWrapper createKernelUpdates() const
  {
    return BASE::template createDerivedKernelUpdates< KernelWrapper >( m_damage.toView(),
                                                                       m_strainEnergyDensity.toView(),
                                                                       m_lengthScale,
                                                                       m_criticalFractureEnergy,
                                                                       m_criticalStrainEnergy );
  }

  struct viewKeyStruct : public BASE::viewKeyStruct
  {
    static constexpr char const * damageString() { return "damage"; }
    static constexpr char const * strainEnergyDensityString() { return "strainEnergyDensity"; }
    /// string/key for regularization length
    static constexpr char const * lengthScaleString() { return "lengthScale"; }
    /// string/key for Gc
    static constexpr char const * criticalFractureEnergyString() { return "criticalFractureEnergy"; }
    /// string/key for sigma_c
    static constexpr char const * criticalStrainEnergyString() { return "criticalStrainEnergy"; }
  };


protected:
  ///an array with the values of damage at quadrature points
  array2d< real64 > m_damage;
  ///an array with the values of the strain energy density at quadrature points
  array2d< real64 > m_strainEnergyDensity;
  ///the phase-field regularization length
  real64 m_lengthScale;
  ///the critical energy release rate (Gc)
  real64 m_criticalFractureEnergy;
  ///the critical strain energy denstiy (Psi_c) - threshold for damage initiation
  real64 m_criticalStrainEnergy;
};

}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_DAMAGE_HPP_ */
