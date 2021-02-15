/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
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
#include "constitutive/solid/SolidBase.hpp"
#include "constitutive/solid/Damage.hpp"
#include "constitutive/solid/AuxiliaryFunctionsSpectral.hpp"
#include <cstdio>
#define QUADRATIC_DISSIPATION_SPECTRAL 0
#define LORENTZ_SPECTRAL 1

using namespace LvArray;

namespace geosx
{
namespace constitutive
{

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


  using DamageUpdates< UPDATE_BASE >::getStiffness;
  using DamageUpdates< UPDATE_BASE >::smallStrainNoState;
  using DamageUpdates< UPDATE_BASE >::smallStrain;
  using DamageUpdates< UPDATE_BASE >::hypoElastic;
  using DamageUpdates< UPDATE_BASE >::hyperElastic;
  using DamageUpdates< UPDATE_BASE >::m_damage;
  using DamageUpdates< UPDATE_BASE >::m_strainEnergyDensity;
  using DamageUpdates< UPDATE_BASE >::m_criticalStrainEnergy;
  using DamageUpdates< UPDATE_BASE >::m_criticalFractureEnergy;
  using DamageUpdates< UPDATE_BASE >::m_lengthScale;

  #if LORENTZ_SPECTRAL

  //Lorentz type Degradation Function

  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  virtual real64 getDegradationValue( localIndex const k,
                                      localIndex const q ) const override
  {
    //std::cout<<"Lorentz degradation"<<std::endl;
    #if QUADRATIC_DISSIPATION_SPECTRAL
    real64 m = m_criticalFractureEnergy/(2*m_lengthScale*m_criticalStrainEnergy);
    #else
    real64 m = 3*m_criticalFractureEnergy/(8*m_lengthScale*m_criticalStrainEnergy);
    #endif
    real64 p = 1;
    return pow( 1 - m_damage( k, q ), 2 ) /( pow( 1 - m_damage( k, q ), 2 ) + m * m_damage( k, q ) * (1 + p*m_damage( k, q )) );
  }

  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  virtual real64 getDegradationDerivative( real64 const d ) const override
  {
    //std::cout<<"Lorentz derivative"<<std::endl;
    #if QUADRATIC_DISSIPATION_SPECTRAL
    real64 m = m_criticalFractureEnergy/(2*m_lengthScale*m_criticalStrainEnergy);
    #else
    real64 m = 3*m_criticalFractureEnergy/(8*m_lengthScale*m_criticalStrainEnergy);
    #endif
    real64 p = 1;
    return -m*(1 - d)*(1 + (2*p + 1)*d) / pow( pow( 1-d, 2 ) + m*d*(1+p*d), 2 );
  }

  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  virtual real64 getDegradationSecondDerivative( real64 const d ) const override
  {
    //std::cout<<"Lorentz 2nd derivative"<<std::endl;
    #if QUADRATIC_DISSIPATION_SPECTRAL
    real64 m = m_criticalFractureEnergy/(2*m_lengthScale*m_criticalStrainEnergy);
    #else
    real64 m = 3*m_criticalFractureEnergy/(8*m_lengthScale*m_criticalStrainEnergy);
    #endif
    real64 p = 1;
    return -2*m*( pow( d, 3 )*(2*m*p*p + m*p + 2*p + 1) + pow( d, 2 )*(-3*m*p*p -3*p) + d*(-3*m*p - 3) + (-m+p+2) )/pow( pow( 1-d, 2 ) + m*d*(1+p*d), 3 );
  }

  #else
  //Quadratic Degradation

  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  real64 GetDegradationValue( localIndex const k,
                              localIndex const q ) const override
  {
    return (1 - m_damage( k, q ))*(1 - m_damage( k, q ));
  }

  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  real64 GetDegradationDerivative( real64 const d ) const override
  {
    return -2*(1 - d);
  }

  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  real64 GetDegradationSecondDerivative( real64 const d ) const override
  {
    GEOSX_UNUSED_VAR( d );
    return 2.0;
  }
  #endif

  //Modified getStiffness function to account for Spectral Decomposition of Stresses.
  GEOSX_HOST_DEVICE inline
  virtual void getStiffness( localIndex const k,
                             localIndex const q,
                             real64 (& c)[6][6] ) const override final
  {

    //Spectral Split
    UPDATE_BASE::getStiffness( k, q, c );
    real64 const damageFactor = getDegradationValue( k, q );
    real64 const K = UPDATE_BASE::getBulkModulus( k );
    real64 const mu = UPDATE_BASE::getShearModulus( k );
    real64 const lambda = K - 2*mu/3;
    //get strain tensor in voigt form
    real64 strain[6] = {};
    recoverStrainFromStress( this->m_stress[k][q].toSliceConst(), strain, K, mu );
    real64 traceOfStrain = strain[0] + strain[1] + strain[2];
    //get eigenvalues and eigenvectors
    real64 eigenValues[3] = {};
    real64 eigenVectors[3][3] = {};
    LvArray::tensorOps::symEigenvectors< 3 >( eigenValues, eigenVectors, strain );
    //tranpose eigenVectors matrix
    real64 temp[3][3] = {};
    LvArray::tensorOps::transpose< 3, 3 >( temp, eigenVectors );
    LvArray::tensorOps::copy< 3, 3 >( eigenVectors, temp );
    //construct 4th order IxI tensor
    real64 IxITensor[6][6] = {};
    for( int i=0; i < 3; i++ )
    {
      for( int j=0; j < 3; j++ )
      {
        IxITensor[i][j] = 1.0;
      }
    }

    //construct positive part
    real64 cPositive[6][6] = {};
    real64 positiveProjector[6][6] = {};
    PositiveProjectorTensor( eigenValues, eigenVectors, positiveProjector );
    LvArray::tensorOps::scaledCopy< 6, 6 >( cPositive, IxITensor, lambda*heaviside( traceOfStrain ));
    LvArray::tensorOps::scale< 6, 6 >( positiveProjector, 2*mu );
    LvArray::tensorOps::add< 6, 6 >( cPositive, positiveProjector );

    //construct negative part
    real64 negativeProjector[6][6] = {};
    NegativeProjectorTensor( eigenValues, eigenVectors, negativeProjector );
    LvArray::tensorOps::scaledCopy< 6, 6 >( c, IxITensor, lambda*heaviside( -traceOfStrain ));
    LvArray::tensorOps::scale< 6, 6 >( negativeProjector, 2*mu );
    LvArray::tensorOps::add< 6, 6 >( c, negativeProjector );
    //finish up
    LvArray::tensorOps::scale< 6, 6 >( cPositive, damageFactor );
    LvArray::tensorOps::add< 6, 6 >( c, cPositive );
  }

  //With Spectral Decomposition, only the Positive Part of the SED drives damage growth
  GEOSX_HOST_DEVICE
  virtual real64 calculateActiveStrainEnergyDensity( localIndex const k,
                                                     localIndex const q ) const override final
  {
    real64 const K = UPDATE_BASE::getBulkModulus( k );
    real64 const mu = UPDATE_BASE::getShearModulus( k );
    real64 const lambda = K - 2*mu/3;
    //get strain tensor in voigt form
    real64 strain[6] = {};
    recoverStrainFromStress( this->m_stress[k][q].toSliceConst(), strain, K, mu );
    real64 traceOfStrain = strain[0] + strain[1] + strain[2];
    //get eigenvalues and eigenvectors
    real64 eigenValues[3] = {0};
    real64 eigenVectors[3][3] = {};
    LvArray::tensorOps::symEigenvectors< 3 >( eigenValues, eigenVectors, strain );
    //transpose eigenValues matrix
    real64 temp[3][3] = {};
    LvArray::tensorOps::transpose< 3, 3 >( temp, eigenVectors );
    LvArray::tensorOps::copy< 3, 3 >( eigenVectors, temp );
    real64 tracePlus = fmax( traceOfStrain, 0.0 );
    //build symmetric matrices of positive and negative eigenvalues
    real64 eigenPlus[6] = {0};
    for( int i = 0; i < 3; i++ )
    {
      eigenPlus[i] = fmax( eigenValues[i], 0.0 );
    }
    real64 positivePartOfStrain[6] = {0};
    LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( positivePartOfStrain, eigenVectors, eigenPlus );

    real64 const sed = 0.5 * lambda * tracePlus * tracePlus + mu * doubleContraction( positivePartOfStrain, positivePartOfStrain );

    //enforce irreversibility using history field for the strain energy density
    if( sed > m_strainEnergyDensity( k, q ) )
    {
      m_strainEnergyDensity( k, q ) = sed;
    }
    return m_strainEnergyDensity( k, q );
  }

  //Modified getStress
  GEOSX_HOST_DEVICE
  virtual void getStress( localIndex const k,
                          localIndex const q,
                          real64 (& stress)[6] ) const override
  {

    //Spectral split
    real64 const damageFactor = getDegradationValue( k, q );
    real64 const K = UPDATE_BASE::getBulkModulus( k );
    real64 const mu = UPDATE_BASE::getShearModulus( k );
    real64 const lambda = K - 2*mu/3;
    //get strain tensor in voigt form
    real64 strain[6] = {};
    recoverStrainFromStress( this->m_stress[k][q].toSliceConst(), strain, K, mu );
    real64 traceOfStrain = strain[0] + strain[1] + strain[2];
    //get eigenvalues and eigenvectors
    real64 eigenValues[3] = {};
    real64 eigenVectors[3][3] = {};
    LvArray::tensorOps::symEigenvectors< 3 >( eigenValues, eigenVectors, strain );
    //transpose the eigenValues matrix
    real64 temp[3][3] = {};
    LvArray::tensorOps::transpose< 3, 3 >( temp, eigenVectors );
    LvArray::tensorOps::copy< 3, 3 >( eigenVectors, temp );
    //get trace+ and trace-
    real64 tracePlus = fmax( traceOfStrain, 0.0 );
    real64 traceMinus = fmin( traceOfStrain, 0.0 );
    //build symmetric matrices of positive and negative eigenvalues
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
    real64 positiveStress[6] = {};
    real64 negativeStress[6] = {};
    LvArray::tensorOps::scaledCopy< 6 >( positiveStress, Itensor, lambda*tracePlus );
    LvArray::tensorOps::scaledCopy< 6 >( negativeStress, Itensor, lambda*traceMinus );
    LvArray::tensorOps::scaledAdd< 6 >( positiveStress, positivePartOfStrain, 2*mu );
    LvArray::tensorOps::scaledAdd< 6 >( negativeStress, negativePartOfStrain, 2*mu );
    LvArray::tensorOps::copy< 6 >( stress, negativeStress );
    LvArray::tensorOps::scaledAdd< 6 >( stress, positiveStress, damageFactor );
  }

  //if we use the Linear Local Dissipation, we need an energy threshold
  GEOSX_HOST_DEVICE
  virtual real64 getEnergyThreshold() const override
  {
    #if LORENTZ_SPECTRAL
    return m_criticalStrainEnergy;
    #else
    return 3*m_criticalFractureEnergy/(16 * m_lengthScale);
    #endif
  }

};

template< typename BASE >
class DamageSpectral : public Damage< BASE >
{
public:

  /// @typedef Alias for LinearElasticIsotropicUpdates
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


  KernelWrapper createKernelUpdates()
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
