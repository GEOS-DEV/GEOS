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
 * @file DamageSpectral.hpp
 * @brief Overrides the SSLE constitutive updates to account for a damage varible and spectral split.
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_DAMAGESPECTRAL_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_DAMAGESPECTRAL_HPP_
#include "Damage.hpp"
#include "DamageSpectralUtilities.hpp"
#include "PropertyConversions.hpp"
#include "SolidBase.hpp"
#include "SolidModelDiscretizationOpsFullyAnisotroipic.hpp"

#define QUADRATIC_DISSIPATION 0

using namespace LvArray;

namespace geosx
{
namespace constitutive
{

/**
 * @class DamageSpectralUpdates
 *
 * @tparam UPDATE_BASE the underlying intact solid model
 *
 * this class implements the material updates for the case of a damage response
 * that is assymetric in tension and compression. The spectral decomposition of
 * the strain tensor is used to effect the assymetry.
 *
 * In this implementation, the fracture behavior of the material is assumed to be
 * cohesive, and the degradation function proposed by Lorentz et al. is used.
 *
 * References: (for the cohesive model)
 *
 * Lorentz, E., Cuvilliez, S., Kazymyrenko, K., 2011. Convergence of a gradient damage model toward a cohesive zone
 * model. Comptes Rendus Mcanique 339 (1), 20 – 26
 *
 * Lorentz, E., Cuvilliez, S., Kazymyrenko, K., 2012. Modelling large crack propagation: from gradient damage to
 * cohesive zone models. International Journal of Fracture 178 (1), 85–95
 *
 * Geelen, R., Liu, Y., Hu, T., Tupek, M., Dolbow, J., 2019. A phase-field formulation for dynamic cohesive fracture.
 * Computer Methods in Applied Mechanics and Engineering 348, 680-711
 *
 */
template< typename UPDATE_BASE >
class DamageSpectralUpdates : public DamageUpdates< UPDATE_BASE >
{
public:
  template< typename ... PARAMS >
  DamageSpectralUpdates( arrayView2d< real64 > const & inputDamage,
                         arrayView2d< real64 > const & inputStrainEnergyDensity,
                         real64 const & inputLengthScale,
                         real64 const & inputCriticalFractureEnergy,
                         real64 const & inputcriticalStrainEnergy,
                         PARAMS && ... baseParams ):
    DamageUpdates< UPDATE_BASE >( inputDamage, inputStrainEnergyDensity, inputLengthScale,
                                  inputCriticalFractureEnergy, inputcriticalStrainEnergy,
                                  std::forward< PARAMS >( baseParams )... )
  {}

  using DiscretizationOps = SolidModelDiscretizationOpsFullyAnisotroipic; // could maybe optimize, but general for now

  using DamageUpdates< UPDATE_BASE >::smallStrainUpdate;
  using DamageUpdates< UPDATE_BASE >::saveConvergedState;

  using DamageUpdates< UPDATE_BASE >::getDegradationValue;
  using DamageUpdates< UPDATE_BASE >::getDegradationDerivative;
  using DamageUpdates< UPDATE_BASE >::getDegradationSecondDerivative;
  using DamageUpdates< UPDATE_BASE >::getEnergyThreshold;

  using DamageUpdates< UPDATE_BASE >::m_strainEnergyDensity;
  using DamageUpdates< UPDATE_BASE >::m_criticalStrainEnergy;
  using DamageUpdates< UPDATE_BASE >::m_criticalFractureEnergy;
  using DamageUpdates< UPDATE_BASE >::m_lengthScale;
  using DamageUpdates< UPDATE_BASE >::m_damage;
  using DamageUpdates< UPDATE_BASE >::m_disableInelasticity;

  using UPDATE_BASE::m_bulkModulus;  // TODO: model below strongly assumes iso elasticity, templating not so useful
  using UPDATE_BASE::m_shearModulus;


  ///degradation function of the quasi-quadratic type, used to model a cohesive fracture behavior
  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  virtual real64 getDegradationValue( localIndex const k,
                                      localIndex const q ) const override
  {
    #if QUADRATIC_DISSIPATION
    real64 m = m_criticalFractureEnergy/(2*m_lengthScale*m_criticalStrainEnergy);
    #else
    real64 m = 3*m_criticalFractureEnergy/(8*m_lengthScale*m_criticalStrainEnergy);
    #endif
    real64 p = 1;
    return pow( 1 - m_damage( k, q ), 2 ) /( pow( 1 - m_damage( k, q ), 2 ) + m * m_damage( k, q ) * (1 + p*m_damage( k, q )) );
  }

  ///first derivative of the quasi-quadratic degradation function
  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  virtual real64 getDegradationDerivative( real64 const d ) const override
  {
    #if QUADRATIC_DISSIPATION
    real64 m = m_criticalFractureEnergy/(2*m_lengthScale*m_criticalStrainEnergy);
    #else
    real64 m = 3*m_criticalFractureEnergy/(8*m_lengthScale*m_criticalStrainEnergy);
    #endif
    real64 p = 1;
    return -m*(1 - d)*(1 + (2*p + 1)*d) / pow( pow( 1-d, 2 ) + m*d*(1+p*d), 2 );
  }

  ///second derivative of the quasi-quadratic degradation function
  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  virtual real64 getDegradationSecondDerivative( real64 const d ) const override
  {
    #if QUADRATIC_DISSIPATION
    real64 m = m_criticalFractureEnergy/(2*m_lengthScale*m_criticalStrainEnergy);
    #else
    real64 m = 3*m_criticalFractureEnergy/(8*m_lengthScale*m_criticalStrainEnergy);
    #endif
    real64 p = 1;
    return -2*m*( pow( d, 3 )*(2*m*p*p + m*p + 2*p + 1) + pow( d, 2 )*(-3*m*p*p -3*p) + d*(-3*m*p - 3) + (-m+p+2) )/pow( pow( 1-d, 2 ) + m*d*(1+p*d), 3 );
  }

  ///this function implements the modified stress update to account for the spectral split - it also computes the stiffness and active part
  /// of the strain energy
  /**
   * @brief performs the update of stresses and stiffness, using damage with spectral decomposition
   *
   * This function performs the updates of the mechanical state of a damaged solid, assuming the spectral decomposition
   * for assymetric damage behavior in tension and compression. From the strain increment computed in the solver, it
   * updates the stresses and stiffness of the material at each quadrature point.
   *
   * @note The active part of the strain energy density, stored in m_strainEnergyDensity is also updated in this step.
   * This active part depends on the strain decomposition, since only the positive part of strain contributes to damage
   * evolution.
   *
   * @param[in] k The element index.
   * @param[in] q The quadrature point index.
   * @param[in] strainIncrement Strain increment in Voight notation (linearized strain)
   * @param[out] stress New stress value (Cauchy stress)
   * @param[out] stiffness New stiffness value
   */
  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( & stress )[6],
                                  real64 ( & stiffness )[6][6] ) const override final
  {
    // perform elastic update for "undamaged" stress

    UPDATE_BASE::smallStrainUpdate( k, q, strainIncrement, stress, stiffness );  // elastic trial update

    if( m_disableInelasticity )
    {
      return;
    }

    // get undamaged elastic strain

    real64 strain[6];
    UPDATE_BASE::getElasticStrain( k, q, strain );

    strain[3] = strain[3]/2; // eigen-decomposition below does not use engineering strains
    strain[4] = strain[4]/2;
    strain[5] = strain[5]/2;

    real64 traceOfStrain = strain[0] + strain[1] + strain[2];

    real64 mu = m_shearModulus[k];
    real64 lambda = conversions::bulkModAndShearMod::toFirstLame( m_bulkModulus[k], mu );
    real64 damageFactor = getDegradationValue( k, q );

    // get eigenvalues and eigenvectors

    real64 eigenValues[3] = {};
    real64 eigenVectors[3][3] = {};
    LvArray::tensorOps::symEigenvectors< 3 >( eigenValues, eigenVectors, strain );

    // tranpose eigenVectors matrix

    real64 temp[3][3] = {};
    LvArray::tensorOps::transpose< 3, 3 >( temp, eigenVectors );
    LvArray::tensorOps::copy< 3, 3 >( eigenVectors, temp );

    // get trace+ and trace-

    real64 tracePlus = fmax( traceOfStrain, 0.0 );
    real64 traceMinus = fmin( traceOfStrain, 0.0 );

    // build symmetric matrices of positive and negative eigenvalues

    real64 eigenPlus[6] = {};
    real64 eigenMinus[6] = {};
    real64 Itensor[6] = {};

    for( int i = 0; i < 3; i++ )
    {
      Itensor[i] = 1;
      eigenPlus[i] = fmax( eigenValues[i], 0.0 );
      eigenMinus[i] = fmin( eigenValues[i], 0.0 );
    }

    real64 positivePartOfStrain[6] = {};
    real64 negativePartOfStrain[6] = {};
    LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( positivePartOfStrain, eigenVectors, eigenPlus );
    LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( negativePartOfStrain, eigenVectors, eigenMinus );

    // stress

    real64 positiveStress[6] = {};
    real64 negativeStress[6] = {};
    LvArray::tensorOps::scaledCopy< 6 >( positiveStress, Itensor, lambda*tracePlus );
    LvArray::tensorOps::scaledCopy< 6 >( negativeStress, Itensor, lambda*traceMinus );

    LvArray::tensorOps::scaledAdd< 6 >( positiveStress, positivePartOfStrain, 2*mu );
    LvArray::tensorOps::scaledAdd< 6 >( negativeStress, negativePartOfStrain, 2*mu );

    LvArray::tensorOps::copy< 6 >( stress, negativeStress );
    LvArray::tensorOps::scaledAdd< 6 >( stress, positiveStress, damageFactor );

    // stiffness

    real64 IxITensor[6][6] = {};
    for( int i=0; i < 3; i++ )
    {
      for( int j=0; j < 3; j++ )
      {
        IxITensor[i][j] = 1.0;
      }
    }

    real64 cPositive[6][6] = {};
    real64 positiveProjector[6][6] = {};
    real64 negativeProjector[6][6] = {};

    positiveProjectorTensor( eigenValues, eigenVectors, positiveProjector );
    negativeProjectorTensor( eigenValues, eigenVectors, negativeProjector );

    LvArray::tensorOps::scaledCopy< 6, 6 >( cPositive, IxITensor, lambda*heaviside( traceOfStrain ));
    LvArray::tensorOps::scaledCopy< 6, 6 >( stiffness, IxITensor, lambda*heaviside( -traceOfStrain ));

    LvArray::tensorOps::scale< 6, 6 >( positiveProjector, 2*mu );
    LvArray::tensorOps::scale< 6, 6 >( negativeProjector, 2*mu );

    LvArray::tensorOps::add< 6, 6 >( cPositive, positiveProjector );
    LvArray::tensorOps::add< 6, 6 >( stiffness, negativeProjector );

    LvArray::tensorOps::scale< 6, 6 >( cPositive, damageFactor );
    LvArray::tensorOps::add< 6, 6 >( stiffness, cPositive );

    // compute strain energy density

    real64 const sed = 0.5 * lambda * tracePlus * tracePlus + mu * doubleContraction( positivePartOfStrain, positivePartOfStrain );

    if( sed > m_strainEnergyDensity( k, q ) )
    {
      m_strainEnergyDensity( k, q ) = sed;
    }
  }

  ///alternatative form of the strain update if stiffness is passed as a DiscretizationOps
  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( & stress )[6],
                                  DiscretizationOps & stiffness ) const final
  {
    smallStrainUpdate( k, q, strainIncrement, stress, stiffness.m_c );
  }

  ///accessor for strain energy density - in this case, it only cares about the positive part
  GEOSX_HOST_DEVICE
  virtual real64 getStrainEnergyDensity( localIndex const k,
                                         localIndex const q ) const override final
  {
    return m_strainEnergyDensity( k, q );
  }

  /**
  * @brief Accessor for the energy threshold, which, in the cohesive formulation, 
  * is independent of the regularization length.
  * 
  * @return the critical energy threshold value.
  */
  GEOSX_HOST_DEVICE
  virtual real64 getEnergyThreshold() const override final
  {
    return m_criticalStrainEnergy;
  }

};

/**
 * @class DamageSpectral
 *
 * This class implements the changes to the update functions that are needed
 * to account for an assymetric damage response in tension and compression.
 * In this case, the split between tension and compression is effected through the
 * use of a spectral decomposition of the strain tensor. Only the positive part
 * of the strain tensor in this decomposition will be degraded. Also, only this
 * part will contribute to the active part of the strain energy that drives damage
 * evolution.
 *
 * Reference:
 *
 * Miehe, Christian; Hofacker, Martina; Welschinger, Fabian. A phase field model for rate-independent crack
 * propagation: Robust algorithmic implementation based on operator splits.
 * Computer Methods in Applied Mechianics and Engineering, v. 199, n. 45-48, p. 2765-2778, 2010
 *
 */
template< typename BASE >
class DamageSpectral : public Damage< BASE >
{
public:

  using KernelWrapper = DamageSpectralUpdates< typename BASE::KernelWrapper >;

  using Damage< BASE >::m_damage;
  using Damage< BASE >::m_strainEnergyDensity;
  using Damage< BASE >::m_criticalFractureEnergy;
  using Damage< BASE >::m_lengthScale;
  using Damage< BASE >::m_criticalStrainEnergy;

  DamageSpectral( string const & name, dataRepository::Group * const parent );
  virtual ~DamageSpectral() override;


  static string catalogName() { return string( "DamageSpectral" ) + BASE::m_catalogNameString; }
  virtual string getCatalogName() const override { return catalogName(); }


  KernelWrapper createKernelUpdates() const
  {
    return BASE::template createDerivedKernelUpdates< KernelWrapper >( m_damage.toView(),
                                                                       m_strainEnergyDensity.toView(),
                                                                       m_lengthScale,
                                                                       m_criticalFractureEnergy,
                                                                       m_criticalStrainEnergy );
  }

};


}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_DAMAGESPECTRAL_HPP_ */
