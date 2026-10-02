/*
 * SPDX-License-Identifier: LGPL-2.1-only
 * Copyright (c) 2026 GEOS/GEOSX Contributors
 */
#include "SelfConsistentCapillaryEquilibrium.hpp"
#include "functions/FunctionManager.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace geos
{
namespace isothermalCompositionalMultiphaseBaseKernels
{

SelfConsistentCapillaryEquilibrium::SelfConsistentCapillaryEquilibrium(
  constitutive::TableCapillaryPressure & capillary,
  integer numComponents, integer maxIterations, real64 tolerance,
  real64 gravity, EvaluateFluid evaluateFluid, EvaluateRelativePermeability evaluateRelativePermeability )
  : m_numPhases( capillary.numFluidPhases() ), m_numComponents( numComponents ),
  m_maxIterations( maxIterations ), m_tolerance( tolerance ), m_gravity( gravity ),
  m_primaryPhase( -1 ),
  m_capillaryWrapper( capillary.createKernelWrapper() ), m_evaluateFluid( std::move( evaluateFluid ) ), m_evaluateRelativePermeability( std::move( evaluateRelativePermeability ) )
{
  GEOS_THROW_IF( !std::isfinite( m_tolerance ) || m_tolerance <= 0.0 || m_maxIterations <= 0 || !std::isfinite( m_gravity ),
                 "Self-consistent hydrostatic equilibrium requires finite gravity, positive finite tolerance and positive iteration count", InputError );
  GEOS_THROW_IF( m_numPhases < 2 || m_numPhases > MAX_PHASES || m_numComponents != m_numPhases,
                 "Self-consistent capillary equilibrium requires two or three fixed-composition phases/components", InputError );
  auto const order = capillary.phaseOrder();
  order.move( hostMemorySpace, false );
  using PT = constitutive::CapillaryPressureBase::PhaseType;
  integer n = 0;
  for( integer const role : { PT::WATER, PT::OIL, PT::GAS } )
  {
    if( order[role] >= 0 ) m_order[n++] = order[role];
  }
  GEOS_THROW_IF( n != m_numPhases, "Self-consistent capillary equilibrium requires gas/oil/water phase roles", InputError );
  // GEOS uses oil reference pressure whenever oil is present in the model;
  // gas is the primary reference for gas/water.
  m_primaryPhase = order[PT::OIL] >= 0 ? order[PT::OIL] : order[PT::GAS];
  auto validateTable = [&]( char const * key, bool increasing )
  {
    string const & name = capillary.getReference< string >( key );
    TableFunction const & table = FunctionManager::getInstance().getGroup< TableFunction >( name );
    auto const coordinates = table.getCoordinates();
    auto const values = table.getValues();
    GEOS_THROW_IF( coordinates.size() != 1 || values.size() < 2 || coordinates[0].size() != values.size(),
                   "Invalid self-consistent capillary table dimensions", InputError, table.getDataContext() );
    for( localIndex i = 0; i < values.size(); ++i )
    {
      GEOS_THROW_IF( !std::isfinite( coordinates[0][i] ) || coordinates[0][i] < 0.0 || coordinates[0][i] > 1.0 || !std::isfinite( values[i] ) ||
                     (i > 0 && (coordinates[0][i] <= coordinates[0][i-1] || (increasing ? values[i] < values[i-1] : values[i] > values[i-1]))),
                     "Self-consistent capillary tables require finite bounded saturations, finite monotone pressures, and strictly increasing saturation coordinates",
                     InputError, table.getDataContext() );
    }
  };
  using Keys = constitutive::TableCapillaryPressure::viewKeyStruct;
  if( m_numPhases == 3 )
  {
    validateTable( Keys::wettingIntermediateCapPresTableNameString(), false );
    validateTable( Keys::nonWettingIntermediateCapPresTableNameString(), true );
  }
  else validateTable( Keys::wettingNonWettingCapPresTableNameString(), order[PT::WATER] < 0 );
  auto const minimum = capillary.phaseMinVolumeFraction();
  minimum.move( hostMemorySpace, false );
  real64 sum = 0.0;
  for( integer ip = 0; ip < m_numPhases; ++ip )
  {
    m_capillaryMinimum[ip] = minimum[ip];
    GEOS_THROW_IF( !std::isfinite( minimum[ip] ) || minimum[ip] < 0.0 || minimum[ip] > 1.0,
                   "Invalid capillary endpoint saturation", InputError );
    sum += minimum[ip];
  }
  GEOS_THROW_IF( sum >= 1.0, "Capillary residual saturations leave no free pore volume", InputError );
}

real64 SelfConsistentCapillaryEquilibrium::roundoff( PhaseVector const & pressure ) const
{
  real64 scale = 1.0;
  for( integer ip = 0; ip < m_numPhases; ++ip ) scale = std::max( scale, std::abs( pressure[ip] ) );
  return 128.0 * std::numeric_limits< real64 >::epsilon() * scale;
}

SelfConsistentCapillaryEquilibrium::State
SelfConsistentCapillaryEquilibrium::evaluate( real64 coordinate, PhaseVector const & phasePressure, bool checkCompatibility ) const
{
  using namespace constitutive;
  StackArray< real64, 3, MAX_PHASES, cappres::LAYOUT_CAPPRES > targetPc( 1, 1, m_numPhases );
  StackArray< real64, 3, MAX_PHASES, cappres::LAYOUT_CAPPRES > actualPc( 1, 1, m_numPhases );
  StackArray< real64, 2, MAX_PHASES, compflow::LAYOUT_PHASE > saturation( 1, m_numPhases );
  StackArray< real64, 4, MAX_PHASES * MAX_PHASES, cappres::LAYOUT_CAPPRES_DS > derivatives( 1, 1, m_numPhases, m_numPhases );
  real64 sumMinimum = 0.0;
  for( integer ip = 0; ip < m_numPhases; ++ip )
  {
    GEOS_THROW_IF( !std::isfinite( phasePressure[ip] ), "Non-finite integrated hydrostatic phase pressure", InputError );
    targetPc[0][0][ip] = phasePressure[m_primaryPhase] - phasePressure[ip];
    saturation[0][ip] = m_capillaryMinimum[ip];
    sumMinimum += m_capillaryMinimum[ip];
    actualPc[0][0][ip] = 0.0;
  }
  // Table laws are separable and nonincreasing in their native Pc slots. Use
  // exact forward endpoints and a pressure-controlled safeguarded Newton step;
  // the generic inverse's fixed 1e-6 Pa endpoint snap is not a valid tolerance.
  auto invert = [&]( integer ip, integer complement, real64 low, real64 high, real64 target )
  {
    auto value = [&]( real64 trial, real64 & derivative )
    {
      saturation[0][ip] = trial;
      if( complement >= 0 ) saturation[0][complement] = 1.0 - trial;
      m_capillaryWrapper.compute( saturation[0].toSliceConst(), actualPc[0][0], derivatives[0][0] );
      real64 result = actualPc[0][0][ip];
      derivative = derivatives[0][0][ip][ip];
      if( complement >= 0 )
      {
        result -= actualPc[0][0][complement];
        derivative += derivatives[0][0][complement][complement];
      }
      GEOS_THROW_IF( !std::isfinite( result ) || !std::isfinite( derivative ) || derivative > 0.0,
                     "Invalid or nonmonotone capillary table law", InputError );
      return result;
    };
    real64 derivative = 0.0;
    real64 const lowValue = value( low, derivative );
    real64 const highValue = value( high, derivative );
    GEOS_THROW_IF( lowValue < highValue, "Nonmonotone capillary endpoints", InputError );
    if( target >= lowValue ) { value( low, derivative ); return; }
    if( target <= highValue ) { value( high, derivative ); return; }
    real64 trial = 0.5 * (low + high);
    for( integer iteration = 0; iteration < m_maxIterations; ++iteration )
    {
      real64 const actual = value( trial, derivative );
      real64 const error = actual - target;
      real64 const allowance = 0.1 * m_tolerance + 8.0 * std::numeric_limits< real64 >::epsilon() *
                               std::max( 1.0, std::max( std::abs( actual ), std::abs( target ) ) );
      if( std::abs( error ) <= allowance ) return;
      if( error > 0.0 ) low = trial; else high = trial;
      real64 const proposed = derivative < 0.0 ? trial - error / derivative : 0.5 * (low + high);
      real64 const next = proposed > low && proposed < high ? proposed : 0.5 * (low + high);
      if( next == trial ) break;
      trial = next;
    }
    GEOS_THROW( "Self-consistent hydrostatic capillary inversion did not converge to the absolute pressure tolerance", InputError );
  };
  real64 sumIndependent = 0.0;
  for( integer ip = 0; ip < m_numPhases; ++ip )
  {
    if( ip == m_primaryPhase ) continue;
    invert( ip, -1, m_capillaryMinimum[ip], 1.0 - sumMinimum + m_capillaryMinimum[ip], targetPc[0][0][ip] );
    sumIndependent += saturation[0][ip];
  }
  saturation[0][m_primaryPhase] = 1.0 - sumIndependent;
  if( m_numPhases == 3 && saturation[0][m_primaryPhase] < 0.0 )
  {
    integer const gas = m_order[2];
    integer const water = m_order[0];
    saturation[0][m_primaryPhase] = 0.0;
    invert( gas, water, m_capillaryMinimum[gas], 1.0 - m_capillaryMinimum[water],
            targetPc[0][0][gas] - targetPc[0][0][water] );
  }
  constexpr real64 saturationTolerance = 64.0 * std::numeric_limits< real64 >::epsilon();
  real64 saturationSum = 0.0;
  for( integer ip = 0; ip < m_numPhases; ++ip )
  {
    real64 & value = saturation[0][ip];
    GEOS_THROW_IF( !std::isfinite( value ) || value < -saturationTolerance || value > 1.0 + saturationTolerance,
                   "Invalid saturation in self-consistent hydrostatic equilibrium", InputError );
    // Correct only out-of-domain floating-point roundoff. In particular never
    // zero a small positive saturation, which may have substantial mobility.
    value = std::max( 0.0, std::min( 1.0, value ) );
    saturationSum += value;
  }
  GEOS_THROW_IF( std::abs( saturationSum - 1.0 ) > saturationTolerance,
                 "Invalid saturation sum in self-consistent hydrostatic equilibrium", InputError );
  for( integer ip = 0; ip < m_numPhases; ++ip ) saturation[0][ip] /= saturationSum;
  m_capillaryWrapper.compute( saturation[0].toSliceConst(), actualPc[0][0], derivatives[0][0] );
  State state;
  state.phasePressure = phasePressure;
  integer reference = -1;
  real64 sumSaturation = 0.0;
  for( integer ip = 0; ip < m_numPhases; ++ip )
  {
    state.saturation[ip] = saturation[0][ip];
    state.capillaryPressure[ip] = actualPc[0][0][ip];
    GEOS_THROW_IF( !std::isfinite( state.saturation[ip] ) || state.saturation[ip] < -saturationTolerance ||
                   state.saturation[ip] > 1.0 + saturationTolerance || !std::isfinite( state.capillaryPressure[ip] ),
                   "Invalid saturation or capillary pressure in self-consistent hydrostatic equilibrium", InputError );
    sumSaturation += state.saturation[ip];
  }
  PhaseVector relativePermeability{};
  m_evaluateRelativePermeability( state.saturation, relativePermeability );
  auto mobile = [&]( integer ip ) { return state.saturation[ip] > 0.0 && relativePermeability[ip] > 0.0; };
  for( integer ip = 0; ip < m_numPhases; ++ip )
  {
    GEOS_THROW_IF( !std::isfinite( relativePermeability[ip] ) || relativePermeability[ip] < 0.0,
                   "Invalid relative permeability in self-consistent hydrostatic equilibrium", InputError );
    // Match PhaseMobilityKernel: present phase AND actual positive kr. Neither
    // a tiny saturation nor a reported residual saturation proves immobility
    // for all interpolated three-phase laws.
    if( mobile( ip ) && (reference < 0 || state.saturation[ip] > state.saturation[reference]) ) reference = ip;
  }
  GEOS_THROW_IF( reference < 0 || std::abs( sumSaturation - 1.0 ) > saturationTolerance,
                 "No mobile reference phase or invalid saturation sum in self-consistent hydrostatic equilibrium", InputError );
  state.primaryPressure = phasePressure[reference] + state.capillaryPressure[reference];
  GEOS_THROW_IF( !std::isfinite( state.primaryPressure ) || state.primaryPressure <= 0.0,
                 "Invalid primary pressure in self-consistent hydrostatic equilibrium", InputError );
  if( checkCompatibility )
  {
    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      real64 const residual = (phasePressure[ip] - phasePressure[reference]) +
                              (state.capillaryPressure[ip] - state.capillaryPressure[reference]);
      real64 const allowance = m_tolerance + roundoff( phasePressure );
      GEOS_THROW_IF( mobile( ip ) ? std::abs( residual ) > allowance : residual > allowance,
                     "Self-consistent hydrostatic equilibrium has inconsistent mobile-phase capillary pressure or endpoint complementarity; check capillary endpoints and residual saturations",
                     InputError );
    }
  }
  m_evaluateFluid( coordinate, state.primaryPressure, state );
  for( integer ip = 0; ip < m_numPhases; ++ip )
  {
    GEOS_THROW_IF( !std::isfinite( state.massDensity[ip] ) || state.massDensity[ip] <= 0.0 ||
                   !std::isfinite( state.density[ip] ) || state.density[ip] <= 0.0,
                   "Non-finite or non-positive EOS phase density in self-consistent hydrostatic equilibrium", InputError );
    real64 sum = 0.0;
    for( integer ic = 0; ic < m_numComponents; ++ic )
    {
      real64 const value = state.composition[ip][ic];
      GEOS_THROW_IF( !std::isfinite( value ) || value < 0.0 || value > 1.0,
                     "Invalid EOS phase composition in self-consistent hydrostatic equilibrium", InputError );
      sum += value;
    }
    GEOS_THROW_IF( std::abs( sum - 1.0 ) > 1.0e-12, "Invalid EOS phase composition sum in self-consistent hydrostatic equilibrium", InputError );
  }
  return state;
}

bool SelfConsistentCapillaryEquilibrium::step( real64 from, State const & reference, real64 to, State & result ) const
{
  PhaseVector pressure = reference.phasePressure;
  real64 const potentialDifference = m_gravity * (to - from);
  for( integer iteration = 0; iteration < m_maxIterations; ++iteration )
  {
    State const trial = evaluate( to, pressure, false );
    PhaseVector updated = pressure;
    real64 error = 0.0;
    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      updated[ip] = reference.phasePressure[ip] + 0.5 * (reference.massDensity[ip] + trial.massDensity[ip]) * potentialDifference;
      error = std::max( error, std::abs( updated[ip] - pressure[ip] ) );
    }
    // Retain exactly the state at the pressure used in the EOS evaluation.
    // The discrepancy is the actual trapezoid equation residual, not a relative
    // pressure change or a saturation tolerance translated into pressure.
    if( error <= 0.1 * m_tolerance + roundoff( pressure ) ) { result = trial; return true; }
    bool accepted = false;
    for( integer line = 0; line < 20; ++line )
    {
      PhaseVector proposed = pressure;
      real64 const scale = std::ldexp( 1.0, -line );
      for( integer ip = 0; ip < m_numPhases; ++ip ) proposed[ip] += scale * (updated[ip] - pressure[ip]);
      try
      {
        evaluate( to, proposed, false );
        pressure = proposed;
        accepted = true;
        break;
      }
      catch( InputError const & )
      {
        // Keep the previously admissible EOS iterate and shorten this step.
      }
    }
    if( !accepted ) return false;
  }
  return false;
}

bool SelfConsistentCapillaryEquilibrium::march( std::vector< real64 > const & coordinates, localIndex datumIndex,
                                               PhaseVector const & datumPressures, std::vector< State > & states ) const
{
  states.resize( coordinates.size() );
  states[datumIndex] = evaluate( coordinates[datumIndex], datumPressures, false );
  for( integer const direction : { -1, 1 } )
  {
    for( localIndex i = datumIndex + direction; i >= 0 && i < static_cast< localIndex >( coordinates.size() ); i += direction )
    {
      if( !step( coordinates[i-direction], states[i-direction], coordinates[i], states[i] ) ) return false;
    }
  }
  return true;
}

std::vector< SelfConsistentCapillaryEquilibrium::State >
SelfConsistentCapillaryEquilibrium::solve( std::vector< real64 > const & coordinates,
                                         arrayView1d< real64 const > contacts,
                                         real64 datumCoordinate, real64 datumPressure ) const
{
  GEOS_THROW_IF( coordinates.empty() || !std::isfinite( datumCoordinate ) || !std::isfinite( datumPressure ) || datumPressure <= 0.0 ||
                 contacts.size() != m_numPhases - 1,
                 "Invalid datum or contact count in self-consistent hydrostatic equilibrium", InputError );
  for( size_t i = 0; i < coordinates.size(); ++i )
    GEOS_THROW_IF( !std::isfinite( coordinates[i] ) || (i > 0 && coordinates[i] <= coordinates[i-1]),
                   "Self-consistent hydrostatic coordinates must be finite and strictly increasing", InputError );
  for( localIndex i = 0; i < contacts.size(); ++i )
    GEOS_THROW_IF( !std::isfinite( contacts[i] ) || (i > 0 && contacts[i] < contacts[i-1]),
                   "Self-consistent hydrostatic contacts must be finite and ordered", InputError );
  auto indexOf = [&]( real64 value ) -> localIndex
  {
    auto const found = std::lower_bound( coordinates.begin(), coordinates.end(), value );
    GEOS_THROW_IF( found == coordinates.end() || *found != value,
                   "Missing datum/contact node in self-consistent hydrostatic table", InputError );
    return std::distance( coordinates.begin(), found );
  };
  localIndex const datumIndex = indexOf( datumCoordinate );
  std::array< localIndex, MAX_PHASES-1 > contactIndex{};
  integer datumPhase = m_order[0];
  for( integer i = 0; i < m_numPhases-1; ++i )
  {
    contactIndex[i] = indexOf( contacts[i] );
    if( datumCoordinate >= contacts[i] ) datumPhase = m_order[i+1];
  }
  std::array< integer, MAX_PHASES-1 > freePhase{};
  integer n = 0;
  for( integer ip = 0; ip < m_numPhases; ++ip ) if( ip != datumPhase ) freePhase[n++] = ip;
  PhaseVector datumPressures{};
  datumPressures.fill( datumPressure );
  std::vector< State > states;
  auto residual = [&]( std::vector< State > const & values, PhaseVector & result )
  {
    real64 norm = 0.0;
    for( integer j = 0; j < n; ++j )
    {
      auto const & p = values[contactIndex[j]].phasePressure;
      result[j] = p[m_order[j+1]] - p[m_order[j]];
      norm = std::max( norm, std::abs( result[j] ) );
    }
    return norm;
  };
  bool converged = false;
  for( integer iteration = 0; iteration < m_maxIterations; ++iteration )
  {
    GEOS_THROW_IF( !march( coordinates, datumIndex, datumPressures, states ),
                   "Self-consistent hydrostatic EOS fixed-point iteration did not converge", InputError );
    PhaseVector r{};
    real64 const norm = residual( states, r );
    if( norm <= 0.5 * m_tolerance ) { converged = true; break; }
    real64 jacobian[2][2]{};
    for( integer j = 0; j < n; ++j )
    {
      real64 const magnitude = std::max( 0.01, std::sqrt( std::numeric_limits< real64 >::epsilon() ) * std::abs( datumPressures[freePhase[j]] ) );
      bool evaluated = false;
      for( integer const sign : { 1, -1 } )
      {
        PhaseVector perturbed = datumPressures;
        real64 const delta = sign * magnitude;
        perturbed[freePhase[j]] += delta;
        std::vector< State > trial;
        try
        {
          if( !march( coordinates, datumIndex, perturbed, trial ) ) continue;
        }
        catch( InputError const & )
        {
          // At a PVT range boundary use the admissible one-sided derivative.
          continue;
        }
        PhaseVector rp{};
        residual( trial, rp );
        for( integer i = 0; i < n; ++i ) jacobian[i][j] = (rp[i] - r[i]) / delta;
        evaluated = true;
        break;
      }
      GEOS_THROW_IF( !evaluated, "No admissible EOS perturbation for self-consistent hydrostatic contact solve", InputError );
    }
    PhaseVector change{};
    real64 const determinant = n == 1 ? jacobian[0][0] : jacobian[0][0]*jacobian[1][1] - jacobian[0][1]*jacobian[1][0];
    GEOS_THROW_IF( !std::isfinite( determinant ) || std::abs( determinant ) < 1.0e-12,
                   "Singular self-consistent hydrostatic contact constraints", InputError );
    change[0] = n == 1 ? -r[0] / determinant : (-r[0]*jacobian[1][1]+r[1]*jacobian[0][1])/determinant;
    if( n == 2 ) change[1] = (-jacobian[0][0]*r[1]+jacobian[1][0]*r[0])/determinant;
    bool accepted = false;
    // Bounded backtracking never changes the datum anchor or the requested
    // pressure tolerance. Intermediate shooting states need not satisfy Pc
    // compatibility; all accepted final table/cell states are checked below.
    for( integer line = 0; line < 20; ++line )
    {
      PhaseVector proposed = datumPressures;
      real64 const scale = std::ldexp( 1.0, -line );
      for( integer j = 0; j < n; ++j ) proposed[freePhase[j]] += scale * change[j];
      std::vector< State > trial;
      try
      {
        if( !march( coordinates, datumIndex, proposed, trial ) ) continue;
      }
      catch( InputError const & )
      {
        // An inadmissible shooting trial is rejected, never clipped into an
        // invented EOS state. The already-admissible iterate stays unchanged.
        continue;
      }
      PhaseVector rt{};
      if( residual( trial, rt ) < norm ) { datumPressures = proposed; accepted = true; break; }
    }
    GEOS_THROW_IF( !accepted, "Self-consistent hydrostatic contact solve failed to reduce its absolute pressure residual", InputError );
  }
  GEOS_THROW_IF( !converged, "Self-consistent hydrostatic contact solve did not converge", InputError );
  for( localIndex i = 0; i < static_cast< localIndex >( coordinates.size() ); ++i )
  {
    states[i] = evaluate( coordinates[i], states[i].phasePressure );
  }
  return states;
}

}
}
