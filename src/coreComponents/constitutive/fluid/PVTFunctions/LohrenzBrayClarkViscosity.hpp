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
  template< int USD1, int USD2 >
  GEOSX_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseDensity,
                arraySlice2d< real64 const, USD2 > const & phaseComposition,
                arraySlice1d< real64,       USD1 > const & phaseViscosity ) const
  {
    GEOSX_UNUSED_VAR( phaseDensity );

    integer const numPhases = phaseComposition.size(0); 
    integer const numComponents = phaseComposition.size(1);

    // Estimate pure component properties at "near-atmospheric" pressure using
    // Stiel and Thodos [1961] correlation: https://doi.org/10.1002/aic.690070416
    // Units are converted so componentViscosity is in centipoise to match original reference

    array1d< real64 > componentViscosityAtm( numComponents ); 
    for(localIndex c=0; c<numComponents; ++c)
    {
      real64 reducedTemperature = temperature / m_componentCriticalTemperature[c];
      real64 componentChi = chiParameter( m_componentCriticalTemperature[c], m_componentCriticalPressure[c], m_componentMolarWeight[c] );

      if( m_componentMolarWeight[c] < 2.1e-3) // hydrogen correlation, Stiel & Thodos, 1961, Eq. 12
      {
        componentViscosityAtm[c] = 90.71e-5 * pow( 0.1375*temperature-1.67, 5.0/8.0 );
      }
      else if( reducedTemperature <= 1.5 ) // nonpolar gas correlation at low temp, Eq. 9
      {
        componentViscosityAtm[c] = 34e-5 * pow( reducedTemperature, 0.94 ) / componentChi;
      }
      else // nonpolar gas correlation at high temp, Eq. 10
      {
        componentViscosityAtm[c] = 17.78e-5 * pow( 4.58*reducedTemperature-1.67, 0.625 ) / componentChi;
      }
    }

    // Estimate phase viscosity (in cp) at near-atmospheric pressure using the 
    // Herning and Zipperer [1936] mixture rule.  Store in final phaseViscosity array.

    for( localIndex p=0; p<numPhases; ++p )
    {
      real64 A = 0;
      real64 B = 0; 
      for( localIndex c=0; c<numComponents; ++c)
      {
        A += phaseComposition[p][c] * sqrt(m_componentMolarWeight[c]) * componentViscosityAtm[c];
        B += phaseComposition[p][c] * sqrt(m_componentMolarWeight[c]);
      }
      phaseViscosity[p] = A/B;
    }

    // Estimate phase viscosity at given conditions using LBC [1964] correlation.
    // This is an additional term added to the atmospheric pressure estimate above.

    for( localIndex p=0; p<numPhases; ++p )
    {
      // compute phase pseudo properties via Kay's mixing rule
      real64 phaseCriticalTemperature = 0;
      real64 phaseCriticalPressure = 0;
      real64 phaseCriticalVolume = 0;
      real64 phaseMolarWeight = 0;

      for( localIndex c=0; c<numComponents; ++c)
      {
        phaseCriticalTemperature += phaseComposition[p][c] * m_componentCriticalTemperature[c];
        phaseCriticalPressure    += phaseComposition[p][c] * m_componentCriticalPressure[c];
        phaseCriticalVolume      += phaseComposition[p][c] * m_componentCriticalVolume[c];
        phaseMolarWeight         += phaseComposition[p][c] * m_componentMolarWeight[c];
      }

      // compute LBC polynomial
      real64 reducedTemperature = temperature / phaseCriticalTemperature;
      real64 reducedPressure = pressure / phaseCriticalPressure;
      real64 reducedDensity = phaseDensity[p] * phaseCriticalVolume / phaseMolarWeight;
      real64 phaseChi = chiParameter( phaseCriticalTemperature, phaseCriticalPressure, phaseMolarWeight );

      real64 polynomial = 0.1023 
                          + 0.023364*reducedDensity 
                          + 0.058533*pow(reducedDensity,2) 
                          - 0.040758*pow(reducedDensity,3) 
                          + 0.0093324*pow(reducedDensity,4);

      // compute and scale final phase viscosity
      phaseViscosity[p] += (pow(polynomial,4)-1e-4) / phaseChi;
      phaseViscosity[p] *= 1e-3; // convert centipoise to pascal.seconds
    }
  }

  /// Compute "chi" parameter from [ST 1961, LBC 1964]
  /// Using units of (K, atm, amu)
  GEOSX_HOST_DEVICE
  real64 chiParameter( real64 const criticalTemperature, 
                       real64 const criticalPressure, 
                       real64 const molarWeight ) const
  {
    real64 T = pow( criticalTemperature, 1.0/6.0 );
    real64 M = pow( 1000 * molarWeight, 0.5 ); // note: kg/mol to atomic mass units
    real64 P = pow( criticalPressure / 101325., 2.0/3.0 ); // note: pascal to atm conversion
    return T/(M*P);
  }


  /*
  template< int USD1, int USD2, int USD3, int USD4 >
  GEOSX_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseComposition,
                arraySlice1d< real64 const, USD2 > const & dPhaseComposition_dPressure,
                arraySlice1d< real64 const, USD2 > const & dPhaseComposition_dTemperature,
                arraySlice2d< real64 const, USD3 > const & dPhaseComposition_dGlobalCompFraction,
                real64 & value,
                real64 & dValue_dPressure,
                real64 & dValue_dTemperature,
                arraySlice1d< real64, USD4 > const & dValue_dGlobalCompFraction,
                bool useMass ) const;
  */

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
