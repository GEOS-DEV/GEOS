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
 * @file KValueFlashModel.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_KVALUEFLASHMODEL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_KVALUEFLASHMODEL_HPP_

#include "FunctionBase.hpp"

#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/RachfordRice.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/Utilities.hpp"

#include "functions/TableFunction.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

template< integer NUM_PHASE >
class KValueFlashParameters;

template< integer NUM_PHASE >
class KValueFlashModelUpdate final : public FunctionBaseUpdate
{
  static constexpr integer numPhases = NUM_PHASE;
  // Must be 2-phase or 3-phase
  static_assert( NUM_PHASE == 2 || NUM_PHASE == 3, "KValue flash must be 2-phase or 3-phase" );

public:
  using PhaseProp = MultiFluidVar< real64, 3, multifluid::LAYOUT_PHASE, multifluid::LAYOUT_PHASE_DC >;
  using PhaseComp = MultiFluidVar< real64, 4, multifluid::LAYOUT_PHASE_COMP, multifluid::LAYOUT_PHASE_COMP_DC >;
  using Deriv = constitutive::multifluid::DerivativeOffset;

  KValueFlashModelUpdate( integer const numComponents,
                          integer const liquidIndex,
                          integer const vapourIndex,
                          integer const aqueousIndex,
                          TableFunction const & pressureTable,
                          TableFunction const & temperatureTable,
                          arrayView1d< integer const > const & presentComponents,
                          arrayView4d< real64 const > const & kValues );

  // Mark as a 2-phase or 3-phase flash
  GEOS_HOST_DEVICE
  static constexpr integer getNumberOfPhases() { return numPhases; }

  template< int USD1, int USD2 >
  GEOS_HOST_DEVICE
  void compute( ComponentProperties::KernelWrapper const & componentProperties,
                real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & compFraction,
                arraySlice2d< real64, USD2 > const & kValues,
                PhaseProp::SliceType const phaseFraction,
                PhaseComp::SliceType const phaseCompFraction ) const;

private:
  integer const m_liquidIndex;
  integer const m_vapourIndex;
  integer const m_aquoesIndex;

  integer const m_numComponents{0};
  integer const m_numPressurePoints{0};
  integer const m_numTemperaturePoints{0};

  /// Table with pressure indices
  TableFunction::KernelWrapper m_pressureTable;

  /// Table with temperature indices
  TableFunction::KernelWrapper m_temperatureTable;

  /// Index of available components
  arrayView1d< integer const > m_presentComponents;

  /// Hypercube of k-values
  arrayView4d< real64 const > m_kValues;
};

template< integer NUM_PHASE >
class KValueFlashModel : public FunctionBase
{
public:
  KValueFlashModel( string const & name,
                    ComponentProperties const & componentProperties,
                    ModelParameters const & modelParameters,
                    arrayView1d< integer const > const phaseTypes );

  static string catalogName();

  FunctionType functionType() const override
  {
    return FunctionType::FLASH;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = KValueFlashModelUpdate< NUM_PHASE >;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  // Create parameters unique to this model
  static std::unique_ptr< ModelParameters > createParameters( std::unique_ptr< ModelParameters > parameters );

  // Determine phase ordering
  static void calculatePhaseOrdering( arrayView1d< integer const > const & phaseTypes,
                                      arrayView1d< integer > const & phaseOrder );
private:
  /// Index of present component
  array1d< integer > m_presentComponents;

  KValueFlashParameters< NUM_PHASE > const * m_parameters{};
  arrayView1d< integer const > const m_phaseTypes;
};

template< integer NUM_PHASE >
template< int USD1, int USD2 >
GEOS_HOST_DEVICE
void
KValueFlashModelUpdate< NUM_PHASE >::compute( ComponentProperties::KernelWrapper const & componentProperties,
                                              real64 const & pressure,
                                              real64 const & temperature,
                                              arraySlice1d< real64 const, USD1 > const & compFraction,
                                              arraySlice2d< real64, USD2 > const & kValues,
                                              PhaseProp::SliceType const phaseFraction,
                                              PhaseComp::SliceType const phaseCompFraction ) const
{
  GEOS_UNUSED_VAR( componentProperties );

  integer const numDofs = 2 + m_numComponents;

  // This code block is responsible for retrieving the K-equilibrium coefficients
  // ($K_i = y_i / x_i$) for all components in all secondary phases, and calculating
  // their partial derivatives with respect to pressure (P) and temperature (T).
  //
  // Since the K-values are stored in a 4D hypercube (`m_kValues`) defined on a
  // non-uniform grid of pressure and temperature points, bilinear interpolation is used.
  //
  // The overall process is:
  // 1. Use pre-built index tables (`m_pressureTable`, `m_temperatureTable`) to map the physical
  // `pressure` and `temperature` values to normalized, fractional coordinates ($p_{\alpha}$, $t_{\alpha}$). These
  // tables return:
  // - The integer part of the coordinate, which is the required lower **grid index** ($p_n, t_n$).
  // - The fractional part of the coordinate, which is the **interpolation weight** ($p_a, t_a$).
  // - The derivative of the normalized coordinate with respect to the physical
  // coordinate, required for the **chain rule** in derivative calculation ($dp_{\alpha}/dP, dt_{\alpha}/dT$).
  // 2. Loop through all secondary phases and components, read the four bounding K-values from the
  // hypercube, and apply the bilinear formula to calculate the interpolated $K_i$ and its
  // derivatives $dK_i/dP$ and $dK_i/dT$.

  real64 pa = 0.0;
  real64 pa_dp = 0.0;
  pa = m_pressureTable.compute( &pressure, &pa_dp );
  // The integer part of the normalized coordinate is the lower pressure index (pn) for interpolation.
  integer const pn = LvArray::math::min( static_cast< integer >(pa), m_numPressurePoints - 2 );
  // The fractional part is the interpolation weight (pa).
  pa -= pn;

  real64 ta = 0.0;
  real64 ta_dt = 0.0;
  ta = m_temperatureTable.compute( &temperature, &ta_dt );
  // The integer part of the normalized coordinate is the lower temperature index (tn) for interpolation.
  integer const tn = LvArray::math::min( static_cast< integer >(ta), m_numTemperaturePoints - 2 );
  // The fractional part is the interpolation weight (ta).
  ta -= tn;

  // Allocate temporary space for derivatives: 2 (for P and T) x (NUM_PHASE-1) x (m_numComponents)
  stackArray3d< real64, 2*(NUM_PHASE-1)*MultiFluidConstants::MAX_NUM_COMPONENTS > kValueDerivatives( 2, NUM_PHASE-1, m_numComponents );
  arraySlice3d< real64 > dK = kValueDerivatives.toSlice();
  for( integer ip = 0; ip < NUM_PHASE-1; ++ip )
  {
    for( integer ic = 0; ic < m_numComponents; ic++ )
    {
      // Retrieve the four bounding K-values for the 2x2 stencil at the grid indices (pn, tn)
      real64 const k00 = m_kValues( ip, ic, pn, tn );     // K at (P_n, T_n)
      real64 const k01 = m_kValues( ip, ic, pn, tn+1 );   // K at (P_n, T_{n+1})
      real64 const k10 = m_kValues( ip, ic, pn+1, tn );   // K at (P_{n+1}, T_n)
      real64 const k11 = m_kValues( ip, ic, pn+1, tn+1 ); // K at (P_{n+1}, T_{n+1})

      kValues( ip, ic ) = (1.0-pa)*(1.0-ta)*k00 + (1.0-pa)*ta*k01 + pa*(1.0-ta)*k10 + pa*ta*k11;

      // Derivative with respect to pressure (P): dK/dP = dK/d(pa) * d(pa)/dP
      dK( Deriv::dP, ip, ic ) = pa_dp*( -(1.0-ta)*k00 - ta*k01 + (1.0-ta)*k10 + ta*k11 );

      // Derivative with respect to temperature (T): dK/dT = dK/d(ta) * d(ta)/dT
      dK( Deriv::dT, ip, ic ) = ta_dt*( -(1.0-pa)*k00 + (1.0-pa)*k01 - pa*k10 + pa*k11 );
    }
  }

  // Find out which components are active
  stackArray1d< integer, MultiFluidConstants::MAX_NUM_COMPONENTS > componentIndices( m_numComponents );
  calculatePresentComponents( m_numComponents, compFraction, componentIndices );
  auto const presentComponents = componentIndices.toSliceConst();

  if constexpr (NUM_PHASE == 2)
  {
    // Check if the Rachford-Rice will actually solve
    // If all the k-values are on one side of unity then the Rachford-Rice will return a trivial solution
    integer constexpr LIQUID_PRESENT = 1;
    integer constexpr VAPOUR_PRESENT = 2;
    integer phasesPresent = 0;
    for( integer const ic : presentComponents )
    {
      real64 const k = kValues( 0, ic ) - 1.0;
      phasesPresent |= ((k < 0.0) ? LIQUID_PRESENT : VAPOUR_PRESENT);
    }

    real64 vapourFraction = 0.0;
    if( phasesPresent == (VAPOUR_PRESENT|LIQUID_PRESENT) )
    {
      vapourFraction = RachfordRice::solve( kValues[0].toSliceConst(), compFraction, presentComponents );

      // Calculate derivatives implicitly from the Rachford-Rice equation
      real64 denominator = 0.0;
      real64 pressureNumerator = 0.0;
      real64 temperatureNumerator = 0.0;
      for( integer ic = 0; ic < m_numComponents; ++ic )
      {
        real64 const k = kValues( 0, ic ) - 1.0;
        real64 const r = 1.0 / ( 1.0 + vapourFraction*k );
        pressureNumerator += compFraction[ic] * dK( Deriv::dP, 0, ic ) * r * r;
        temperatureNumerator += compFraction[ic] * dK( Deriv::dT, 0, ic ) * r * r;
        denominator += compFraction[ic] * k * k * r * r;
        phaseFraction.derivs( m_vapourIndex, Deriv::dC + ic ) = k * r;
      }

      GEOS_ERROR_IF( denominator < MultiFluidConstants::epsilon,
                     "Failed to calculate derivatives for the Rachford-Rice equation." );
      real64 const invDenominator = 1.0 / denominator;
      phaseFraction.derivs( m_vapourIndex, Deriv::dP ) = pressureNumerator * invDenominator;
      phaseFraction.derivs( m_vapourIndex, Deriv::dT ) = temperatureNumerator * invDenominator;
      for( integer ic = 0; ic < m_numComponents; ++ic )
      {
        phaseFraction.derivs( m_vapourIndex, Deriv::dC + ic ) *= invDenominator;
      }

      // Calculate phase compositions
      for( integer ic = 0; ic < m_numComponents; ++ic )
      {
        real64 const k = kValues( 0, ic ) - 1.0;
        real64 const r = 1.0 / ( 1.0 + vapourFraction*k );
        real64 const xi = compFraction[ic] * r;
        real64 const yi = xi * kValues( 0, ic );
        phaseCompFraction.value( m_liquidIndex, ic ) = xi;
        phaseCompFraction.value( m_vapourIndex, ic ) = yi;

        for( integer const idof : {Deriv::dP, Deriv::dT} )
        {
          real64 const dV = phaseFraction.derivs( m_vapourIndex, idof );
          real64 const dKi = dK( idof, 0, ic );
          real64 const dxi = r*(-vapourFraction*xi*dKi + (xi-yi)*dV);
          phaseCompFraction.derivs( m_liquidIndex, ic, idof ) = dxi;
          phaseCompFraction.derivs( m_vapourIndex, ic, idof ) = xi*dKi + dxi*kValues( 0, ic );
        }
        for( integer jc = 0; jc < m_numComponents; ++jc )
        {
          integer const idof = Deriv::dC + jc;

          real64 const dV = phaseFraction.derivs( m_vapourIndex, idof );
          real64 const dz = (jc == ic) ? 1.0 : 0.0;
          real64 const dxi = r*(dz + (xi-yi)*dV);
          phaseCompFraction.derivs( m_liquidIndex, ic, idof ) = dxi;
          phaseCompFraction.derivs( m_vapourIndex, ic, idof ) = dxi*kValues( 0, ic );
        }
      }

      if( vapourFraction < MultiFluidConstants::epsilon || 1.0-vapourFraction < MultiFluidConstants::epsilon )
      {
        vapourFraction = LvArray::math::min( LvArray::math::max( vapourFraction, 0.0 ), 1.0 );
        phaseFraction.value[m_vapourIndex] = vapourFraction;
        LvArray::forValuesInSlice( phaseFraction.derivs[m_vapourIndex], setZero );

        integer const singlePhaseIndex = vapourFraction < MultiFluidConstants::epsilon ? m_liquidIndex : m_vapourIndex;
        LvArray::forValuesInSlice( phaseCompFraction.derivs[singlePhaseIndex], setZero );
        for( integer ic = 0; ic < m_numComponents; ++ic )
        {
          phaseCompFraction.value( singlePhaseIndex, ic ) = compFraction[ic];
          phaseCompFraction.derivs( singlePhaseIndex, ic, Deriv::dC + ic ) = 1.0;
        }
      }
    }
    else
    {
      LvArray::forValuesInSlice( phaseFraction.derivs, setZero );
      LvArray::forValuesInSlice( phaseCompFraction.derivs, setZero );
      vapourFraction = (phasesPresent == LIQUID_PRESENT) ? 0.0 : 1.0;
      for( integer ic = 0; ic < m_numComponents; ++ic )
      {
        phaseCompFraction.value( m_liquidIndex, ic ) = compFraction[ic];
        phaseCompFraction.value( m_vapourIndex, ic ) = compFraction[ic];
        phaseCompFraction.derivs( m_liquidIndex, ic, Deriv::dC + ic ) = 1.0;
        phaseCompFraction.derivs( m_vapourIndex, ic, Deriv::dC + ic ) = 1.0;
      }
    }

    // Complete by calculating liquid phase fraction
    phaseFraction.value[m_vapourIndex] = vapourFraction;
    phaseFraction.value[m_liquidIndex] = 1.0 - vapourFraction;
    for( integer ic = 0; ic < numDofs; ic++ )
    {
      phaseFraction.derivs[m_liquidIndex][ic] = -phaseFraction.derivs[m_vapourIndex][ic];
    }
  }
  else
  {
    GEOS_ERROR( GEOS_FMT( "Rachford-Rice solve for {} phases not implemented.", NUM_PHASE ));
  }
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_KVALUEFLASHMODEL_HPP_
