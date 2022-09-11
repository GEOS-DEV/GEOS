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
 * @file LohrenzBrayClarkViscosity.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_LOHRENZBRAYCLARKVISCOSITY_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_LOHRENZBRAYCLARKVISCOSITY_HPP_

// Compile-time option to use Wilke mixing rule instead.
#define USE_WILKE_MIXING true

#include "PVTCompositionalFunctionBase.hpp"

namespace geosx
{

namespace constitutive
{

namespace PVTProps
{

class LohrenzBrayClarkViscosityUpdate final : public PVTCompositionalFunctionBaseUpdate
{
public:

  LohrenzBrayClarkViscosityUpdate( arrayView1d< real64 const > const & componentMolarWeight,
                                   arrayView1d< real64 const > const & componentCriticalPressure,
                                   arrayView1d< real64 const > const & componentCriticalTemperature,
                                   arrayView1d< real64 const > const & componentCriticalVolume,
                                   arrayView1d< real64 const > const & componentAcentricFactor,
                                   arrayView1d< real64 const > const & componentVolumeShift,
                                   arrayView2d< real64 const > const & componentBinaryCoeff )
    : PVTCompositionalFunctionBaseUpdate( componentMolarWeight,
                                          componentCriticalPressure,
                                          componentCriticalTemperature,
                                          componentCriticalVolume,
                                          componentAcentricFactor,
                                          componentVolumeShift,
                                          componentBinaryCoeff )
  {}


  /// Compute phase viscosities (no derivatives)

  //template< int USD1, int USD2 >
  GEOSX_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                PhaseProp::SliceType const & phaseDensity,
                PhaseComp::SliceType const & phaseComposition,
                PhaseProp::SliceType const & phaseViscosity ) const
  {
    GEOSX_UNUSED_VAR( pressure ); // no-direct pressure dependence (instead through density)

    integer const numPhases = phaseComposition.value.size(0); 
    integer const numComponents = phaseComposition.value.size(1);

    // Estimate pure component properties at dilute-gas conditions (pressure near atmospheric) using
    // Stiel and Thodos [1961] correlation: https://doi.org/10.1002/aic.690070416
    // Dilute viscosity is solely temperature dependent
    // Units are converted so componentViscosity is in centipoise to match original reference

    array1d< real64 > componentDiluteViscosity( numComponents ); 
    array1d< real64 > dComponentDiluteViscosity_dTemperature( numComponents ); 

    computeComponentDiluteViscosity_StielThodos( numComponents,
                                                 temperature,
                                                 componentDiluteViscosity,
                                                 dComponentDiluteViscosity_dTemperature );

    // Estimate phase viscosity (in cp) at dilute gas conditions using either 
    // the Herning and Zipperer [1936] or Wilke [1950] mixture rule.  
    // Store in final phaseViscosity array.
    // Wilke is more expensive.  The classic LBC model uses Herning-Zipperer.

    #if USE_WILKE_MIXING
      computePhaseDiluteViscosity_Wilke( numComponents,
                                         numPhases,
                                         componentDiluteViscosity,
                                         dComponentDiluteViscosity_dTemperature,
                                         phaseComposition,
                                         phaseViscosity );

// FINITE DIFFERENCE TESTING (TO REMOVE)
/* 
                                         phaseComposition.value[1][0] += 0.0001;
                                         real64 tmp = phaseViscosity.value[1];

      computePhaseDiluteViscosity_Wilke( numComponents,
                                         numPhases,
                                         componentDiluteViscosity,
                                         dComponentDiluteViscosity_dTemperature,
                                         phaseComposition,
                                         phaseViscosity );
                                         
                                         phaseComposition.value[1][0] -= 0.0001;

    printf("LBC (%.3e) %.3e %.3e\n",phaseViscosity.value[1], phaseViscosity.derivs[1][2],(phaseViscosity.value[1]-tmp)/0.0001);
*/

    #else
      computePhaseDiluteViscosity_HerningZipperer( numComponents,
                                                   numPhases,
                                                   componentDiluteViscosity,
                                                   dComponentDiluteViscosity_dTemperature,
                                                   phaseComposition,
                                                   phaseViscosity );
    #endif

    // Estimate phase viscosity at given (P,T) conditions using LBC [1964] correlation.
    // This is an additional term added to the dilute gas estimate above.

    computePhaseViscosity_LohrenzBrayClark( numComponents,
                                            numPhases,
                                            phaseDensity,
                                            phaseComposition,
                                            phaseViscosity );


    // scale centipoise to pascal.seconds

    for( integer p=0; p<numPhases; ++p )
    {
      phaseViscosity.value[p] *= 1e-3; 
      for( integer k=0; k<numComponents; ++k)
      {
        phaseViscosity.derivs[p][k] *= 1e-3;
      }
    }
  }


  /// Compute "1/chi" parameter (inverse of the viscosity-reducing parameter) from [ST 1961, LBC 1964].
  /// Using units of (K, atm, amu).
  /// This version returns value and derivatives.
  GEOSX_HOST_DEVICE
  void inverseChiParameter( real64 const criticalTemperature, 
                            real64 const criticalPressure, 
                            real64 const molarWeight,
                            real64 & value,
                            real64 & derivP,
                            real64 & derivT,
                            real64 & derivM) const
  {
    real64 T  = pow( criticalTemperature, 1.0/6.0 );
    real64 dT = (1.0/6.0) * pow( criticalTemperature, -5.0/6.0 );

    real64 M  = sqrt( 1000*molarWeight ); // note: kg/mol to atomic mass units
    real64 dM = 5 * sqrt(10) / sqrt( molarWeight );

    real64 P  = pow( criticalPressure / 101325., 2.0/3.0 ); // note: pascal to atm conversion
    real64 dP = pow( 101325, -2.0/3.0 ) * pow( criticalPressure, -1.0/3.0 ) * 2.0/3.0;

    value  = M*P/T;
    derivP = M*dP/T;
    derivT = -M*P*dT/(T*T);
    derivM = dM*P/T;
  }


  /// Compute "1/chi" parameter (inverse of the viscosity-reducing parameter) from [ST 1961, LBC 1964].
  /// This version returns only the value, without derivatives.
  GEOSX_HOST_DEVICE
  real64 inverseChiParameter( real64 const criticalTemperature, 
                              real64 const criticalPressure, 
                              real64 const molarWeight) const
  {
    real64 value, discard;
    inverseChiParameter( criticalTemperature, criticalPressure, molarWeight, value, discard, discard, discard );
    return value;
  }


  /// Estimate pure component properties at dilute-gas conditions (pressure near atmospheric) using
  /// Stiel and Thodos [1961] correlation: https://doi.org/10.1002/aic.690070416.
  /// Dilute viscosity is solely temperature dependent.
  /// Units are converted so componentViscosity is in centipoise to match original reference.

  GEOSX_HOST_DEVICE
  void computeComponentDiluteViscosity_StielThodos( integer const numComponents,
                                                    real64 const temperature,
                                                    array1d<real64> & componentDiluteViscosity,
                                                    array1d<real64> & dComponentDiluteViscosity_dTemperature ) const
  {
    for(integer c=0; c<numComponents; ++c)
    {
      real64 reducedTemperature = temperature / m_componentCriticalTemperature[c];
      real64 inverseComponentChi = inverseChiParameter( m_componentCriticalTemperature[c], m_componentCriticalPressure[c], m_componentMolarWeight[c] );

      if( m_componentMolarWeight[c] < 2.1e-3) // hydrogen correlation, Stiel & Thodos, 1961, Eq. 12
      {
        componentDiluteViscosity[c] = 90.71e-5 * pow( 0.1375*temperature-1.67, 0.625 );
        dComponentDiluteViscosity_dTemperature[c] = 90.71e-5 * 0.625 * 0.1375 * pow( 0.1375*temperature-1.67, -0.375); 
      }
      else if( reducedTemperature <= 1.5 ) // nonpolar gas correlation at low temp, Eq. 9
      {
        componentDiluteViscosity[c] = 34e-5 * pow( reducedTemperature, 0.94 ) * inverseComponentChi;
        dComponentDiluteViscosity_dTemperature[c] = 34e-5 * 0.94 * pow( reducedTemperature, -0.06) * inverseComponentChi / m_componentCriticalTemperature[c];
      }
      else // nonpolar gas correlation at high temp, Eq. 10
      {
        componentDiluteViscosity[c] = 17.78e-5 * pow( 4.58*reducedTemperature-1.67, 0.625 ) * inverseComponentChi;
        dComponentDiluteViscosity_dTemperature[c] = 17.78e-5 * 4.58 * 0.625 * pow( 4.58*reducedTemperature-1.67, -0.375) * inverseComponentChi / m_componentCriticalTemperature[c];
      }
    }
  }


  /// Estimate phase viscosity (in cp) at dilute gas conditions using the 
  /// Herning and Zipperer [1936] mixture rule.  Store in phaseViscosity array.

  GEOSX_HOST_DEVICE
  void computePhaseDiluteViscosity_HerningZipperer( integer const numComponents,
                                                    integer const numPhases,
                                                    array1d< real64 > const & componentDiluteViscosity,
                                                    array1d< real64 > const & dComponentDiluteViscosity_dTemperature,
                                                    PhaseComp::SliceType const & phaseComposition,
                                                    PhaseProp::SliceType const & phaseViscosity ) const
  {
    for( integer p=0; p<numPhases; ++p )
    {
      real64 A = 0, dA_dT = 0;
      real64 B = 0; 

      for( integer c=0; c<numComponents; ++c)
      {
        A += phaseComposition.value[p][c] * sqrt(m_componentMolarWeight[c]) * componentDiluteViscosity[c];
        B += phaseComposition.value[p][c] * sqrt(m_componentMolarWeight[c]);
        dA_dT += phaseComposition.value[p][c] * sqrt(m_componentMolarWeight[c]) * dComponentDiluteViscosity_dTemperature[c];
      }

      phaseViscosity.value[p] = A/B;
      phaseViscosity.derivs[p][Deriv::dP] = 0;
      phaseViscosity.derivs[p][Deriv::dT] = dA_dT/B;

      for( integer c=0; c<numComponents; ++c)
      {
        phaseViscosity.derivs[p][Deriv::dC+c] = sqrt(m_componentMolarWeight[c]) * componentDiluteViscosity[c] / B - A * sqrt(m_componentMolarWeight[c]) / (B * B);
      }
    }
  }


  /// Estimate phase viscosity (in cp) at dilute gas conditions using the 
  /// Wilke [1950] mixture rule.  Store in phaseViscosity array.

  GEOSX_HOST_DEVICE
  void computePhaseDiluteViscosity_Wilke( integer const numComponents,
                                          integer const numPhases,
                                          array1d< real64 > const & componentDiluteViscosity,
                                          array1d< real64 > const & dComponentDiluteViscosity_dTemperature,
                                          PhaseComp::SliceType const & phaseComposition,
                                          PhaseProp::SliceType const & phaseViscosity ) const
  {
    // compute the "phi" interaction matrix (and its temperature derivatives)

    array2d< real64> phi( numComponents, numComponents );
    array2d< real64> dPhi_dT( numComponents, numComponents );

    phi.zero();
    dPhi_dT.zero();

    for( integer c=0; c<numComponents; ++c )
    {
      for( integer d=0; d<numComponents; ++d )
      {
        real64 weightRatio = m_componentMolarWeight[c] / m_componentMolarWeight[d];
        real64 viscosityRatio = componentDiluteViscosity[c] / componentDiluteViscosity[d];

        real64 dViscosityRatio_dT = dComponentDiluteViscosity_dTemperature[c] / componentDiluteViscosity[d] 
                                    - componentDiluteViscosity[c] * pow(componentDiluteViscosity[d],-2.0) * dComponentDiluteViscosity_dTemperature[d]; 

        real64 A = 1+sqrt(viscosityRatio)*pow(weightRatio,-0.25); 
        real64 dA_dT = 0.5*pow(weightRatio,-0.25)/sqrt(viscosityRatio)*dViscosityRatio_dT;

        real64 B = pow(A,2.0);
        real64 C = sqrt(8)*sqrt(1+weightRatio);

        phi[c][d] = B/C;
        dPhi_dT[c][d] = 2*A*dA_dT/C;
      } 
    }

    // compute phase viscosity via Wilke mixing rule

    for( integer p=0; p<numPhases; ++p )
    {
      phaseViscosity.value[p] = 0;
      phaseViscosity.derivs[p][Deriv::dP] = 0; // no pressure dependence
      phaseViscosity.derivs[p][Deriv::dT] = 0;

      for( integer c=0; c<numComponents; ++c )
      {
        phaseViscosity.derivs[p][Deriv::dC+c] = 0;
      } 

      for( integer c=0; c<numComponents; ++c )
      {
        real64 A = 0;
        real64 dA_dT = 0;

        for( integer d=0; d<numComponents; ++d )
        { 
          A += phi[c][d]*phaseComposition.value[p][d];
          dA_dT += dPhi_dT[c][d]*phaseComposition.value[p][d];
        } 

        phaseViscosity.value[p] += phaseComposition.value[p][c] * componentDiluteViscosity[c] / A;
        phaseViscosity.derivs[p][Deriv::dT] += phaseComposition.value[p][c] * dComponentDiluteViscosity_dTemperature[c] / A
                                               - phaseComposition.value[p][c] * componentDiluteViscosity[c] / (A*A) * dA_dT;

        phaseViscosity.derivs[p][Deriv::dC+c] += componentDiluteViscosity[c] / A; 
        for( integer d=0; d<numComponents; ++d )
        {  
          phaseViscosity.derivs[p][Deriv::dC+d] -= phaseComposition.value[p][c] * componentDiluteViscosity[c] / (A*A) * phi[c][d];
        } 
      }
    }
  }


  /// Compute additional phase viscosity contribution at given (P,T) conditions using the LBC [1964] correlation, 
  /// and add to phaseViscosity array.

  GEOSX_HOST_DEVICE
  void computePhaseViscosity_LohrenzBrayClark( integer const numComponents,
                                               integer const numPhases,
                                               PhaseProp::SliceType const & phaseDensity,
                                               PhaseComp::SliceType const & phaseComposition,
                                               PhaseProp::SliceType const & phaseViscosity ) const
  {
    for( integer p=0; p<numPhases; ++p )
    {
      // compute phase pseudo properties via Kay's mixing rule

      real64 phaseCriticalTemperature = 0;
      real64 phaseCriticalPressure = 0;
      real64 phaseCriticalVolume = 0;
      real64 phaseMolarWeight = 0;

      for( integer c=0; c<numComponents; ++c)
      {
        phaseCriticalTemperature += phaseComposition.value[p][c] * m_componentCriticalTemperature[c];
        phaseCriticalPressure    += phaseComposition.value[p][c] * m_componentCriticalPressure[c];
        phaseCriticalVolume      += phaseComposition.value[p][c] * m_componentCriticalVolume[c];
        phaseMolarWeight         += phaseComposition.value[p][c] * m_componentMolarWeight[c];
      }

      // compute LBC polynomial

      real64 reducedDensity = phaseDensity.value[p] * phaseCriticalVolume / phaseMolarWeight;
      real64 inversePhaseChi, dInversePhaseChi_dPc, dInversePhaseChi_dTc, dInversePhaseChi_dMw;

      inverseChiParameter( phaseCriticalTemperature, 
                           phaseCriticalPressure, 
                           phaseMolarWeight,
                           inversePhaseChi, 
                           dInversePhaseChi_dPc, 
                           dInversePhaseChi_dTc, 
                           dInversePhaseChi_dMw );

      real64 polynomialOne = 0.1023 
                             + 0.023364*reducedDensity 
                             + 0.058533*pow(reducedDensity,2) 
                             - 0.040758*pow(reducedDensity,3) 
                             + 0.0093324*pow(reducedDensity,4);

      real64 polynomialTwo = pow(polynomialOne,4)-1e-4;

      // add polynomial contribution to dilute term to get final phase viscosity

      phaseViscosity.value[p] += polynomialTwo * inversePhaseChi;

      // get derivatives, noting phaseDensity is a function of (pressure, temperature, composition)
      // and inversePhaseChi is a function of composition.
      // these derivatives are *added* to the ones already present from the dilute terms above.

      real64 dPolynomialOne_dReducedDensity = 0.023364 + 0.117066*reducedDensity - 0.122274*pow(reducedDensity,2) + 0.0373296*pow(reducedDensity,3);
      real64 dPolynomialTwo_dReducedDensity = 4 * pow(polynomialOne,3) * dPolynomialOne_dReducedDensity;
    
      real64 dViscosity_dDensity       = dPolynomialTwo_dReducedDensity * inversePhaseChi * phaseCriticalVolume / phaseMolarWeight;
      real64 dViscosity_dCriticalRatio = dPolynomialTwo_dReducedDensity * inversePhaseChi * phaseDensity.value[p];
 
      real64 dViscosity_dPc = polynomialTwo * dInversePhaseChi_dPc;
      real64 dViscosity_dTc = polynomialTwo * dInversePhaseChi_dTc;
      real64 dViscosity_dMw = polynomialTwo * dInversePhaseChi_dMw - dViscosity_dCriticalRatio * phaseCriticalVolume / pow(phaseMolarWeight,2);
      real64 dViscosity_dVc = dViscosity_dCriticalRatio / phaseMolarWeight;

      phaseViscosity.derivs[p][Deriv::dP] += dViscosity_dDensity * phaseDensity.derivs[p][Deriv::dP];
      phaseViscosity.derivs[p][Deriv::dT] += dViscosity_dDensity * phaseDensity.derivs[p][Deriv::dT];

      for( integer c=0; c<numComponents; ++c)
      {
        phaseViscosity.derivs[p][Deriv::dC+c] += dViscosity_dDensity * phaseDensity.derivs[p][Deriv::dC+c]
                                                 + dViscosity_dPc * m_componentCriticalPressure[c]
                                                 + dViscosity_dTc * m_componentCriticalTemperature[c]
                                                 + dViscosity_dMw * m_componentMolarWeight[c]
                                                 + dViscosity_dVc * m_componentCriticalVolume[c];
      }
    }
  }

protected:

};


class LohrenzBrayClarkViscosity : public PVTCompositionalFunctionBase
{
public:

  LohrenzBrayClarkViscosity( string const & name,
                             string_array const & componentNames,
                             array1d< real64 > const & componentMolarWeight,
                             array1d< real64 > const & componentCriticalPressure,
                             array1d< real64 > const & componentCriticalTemperature,
                             array1d< real64 > const & componentCriticalVolume,
                             array1d< real64 > const & componentAcentricFactor,
                             array1d< real64 > const & componentVolumeShift,
                             array2d< real64 > const & componentBinaryCoeff );

  virtual ~LohrenzBrayClarkViscosity() override = default;

  static string catalogName() { return "LohrenzBrayClarkViscosity"; }

  virtual string getCatalogName() const override final { return catalogName(); }

  virtual PVTCompositionalFunctionType functionType() const override
  {
    return PVTCompositionalFunctionType::VISCOSITY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = LohrenzBrayClarkViscosityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

private:

};

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_LOHRENZBRAYCLARKVISCOSITY_HPP_
