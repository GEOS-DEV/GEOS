/*
 * SPDX-License-Identifier: LGPL-2.1-only
 * Copyright (c) 2026 GEOS/GEOSX Contributors
 */
#ifndef GEOS_SELFCONSISTENTCAPILLARYEQUILIBRIUM_HPP
#define GEOS_SELFCONSISTENTCAPILLARYEQUILIBRIUM_HPP

#include "constitutive/capillaryPressure/TableCapillaryPressure.hpp"
#include <array>
#include <functional>
#include <vector>

namespace geos
{
namespace isothermalCompositionalMultiphaseBaseKernels
{

/**
 * Host-only boundary-value solve for fixed-phase-composition fluids and a
 * spatially uniform TableCapillaryPressure. Every constitutive evaluation uses
 * the same primary pressure that will be used by the flow solver. The datum
 * fixes its contact-ordered phase pressure, and all contacts are solved together.
 * The callback deliberately avoids instantiating fluid/capillary cross products.
 */
class SelfConsistentCapillaryEquilibrium
{
public:
  static constexpr integer MAX_PHASES = 3;
  using PhaseVector = std::array< real64, MAX_PHASES >;
  struct State
  {
    PhaseVector phasePressure{}, saturation{}, capillaryPressure{}, massDensity{}, density{};
    std::array< PhaseVector, MAX_PHASES > composition{};
    real64 primaryPressure = 0.0;
  };
  using EvaluateFluid = std::function< void( real64, real64, State & ) >;
  using EvaluateRelativePermeability = std::function< void( PhaseVector const &, PhaseVector & ) >;

  SelfConsistentCapillaryEquilibrium( constitutive::TableCapillaryPressure & capillary,
                                     integer numComponents,
                                     integer maxIterations,
                                     real64 tolerance,
                                     real64 gravity,
                                     EvaluateFluid evaluateFluid,
                                     EvaluateRelativePermeability evaluateRelativePermeability );

  std::vector< State > solve( std::vector< real64 > const & coordinates,
                             arrayView1d< real64 const > contacts,
                             real64 datumCoordinate,
                             real64 datumPressure ) const;

  State evaluate( real64 coordinate, PhaseVector const & phasePressure, bool checkCompatibility = true ) const;

private:
  bool step( real64 from, State const & reference, real64 to, State & result ) const;
  bool march( std::vector< real64 > const & coordinates, localIndex datumIndex,
              PhaseVector const & datumPressures, std::vector< State > & states ) const;
  real64 roundoff( PhaseVector const & pressure ) const;

  integer m_numPhases, m_numComponents, m_maxIterations;
  real64 m_tolerance, m_gravity;
  std::array< integer, MAX_PHASES > m_order{};
  integer m_primaryPhase;
  PhaseVector m_capillaryMinimum{};
  constitutive::TableCapillaryPressure::KernelWrapper m_capillaryWrapper;
  EvaluateFluid m_evaluateFluid;
  EvaluateRelativePermeability m_evaluateRelativePermeability;
};

}
}
#endif
