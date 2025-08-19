/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TemplatedFlowSolvers.hpp
 * @brief Template-based flow solver hierarchy for pressure and transport equations
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_TEMPLATEDFLOWSOLVERS_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_TEMPLATEDFLOWSOLVERS_HPP_

#include "FlowSolverBase.hpp"
#include "FluxTypes.hpp"
#include "EnhancedPhaseProperties.hpp"

namespace geos
{

/**
 * @brief Template enumeration for equation types
 */
enum class EquationType
{
  PRESSURE,    ///< Pressure equation
  TRANSPORT    ///< Transport equation
};

/**
 * @brief Template enumeration for discretization types
 */
enum class DiscretizationType
{
  MIMETIC          ///< Hybrid-mixed mimetic methods
};

/**
 * @brief Base template for templated flow solvers
 * @tparam NUM_PHASES Number of phases in the system
 * @tparam EQUATION_TYPE Type of equation (pressure or transport)
 * @tparam HAS_CAPILLARITY Flag for capillarity effects
 * @tparam HAS_THERMAL Flag for thermal effects
 */
template< int NUM_PHASES,
          EquationType EQUATION_TYPE,
          bool HAS_CAPILLARITY = false,
          bool HAS_THERMAL = false >
class TemplatedFlowSolver : public FlowSolverBase
{
public:
  static constexpr int numPhases() { return NUM_PHASES; }
  static constexpr EquationType equationType() { return EQUATION_TYPE; }
  static constexpr DiscretizationType discretizationType() { return DiscretizationType::MIMETIC; }
  static constexpr bool hasCapillarity() { return HAS_CAPILLARITY; }
  static constexpr bool hasThermal() { return HAS_THERMAL; }

  // Determine flux type based on capabilities
  using FluxType = typename std::conditional<
    HAS_CAPILLARITY || HAS_THERMAL,
    fluxTypes::ParabolicScalar< NUM_PHASES >,
    fluxTypes::HyperbolicScalar< NUM_PHASES >
  >::type;

  using PhaseModel = phaseProperties::EnhancedPhaseModel< NUM_PHASES >;

  /**
   * @brief Constructor
   * @param[in] name Solver name
   * @param[in] parent Parent group
   */
  TemplatedFlowSolver( string const & name, Group * const parent );

  /**
   * @brief Virtual destructor
   */
  virtual ~TemplatedFlowSolver() = default;

protected:
  /**
   * @brief Compute flux contributions using mimetic discretization
   * @param[in] subRegion Element subregion
   * @param[in] stencil Computational stencil
   */
  virtual void computeFluxContributions( ElementSubRegionBase & subRegion,
                                         localIndex const stencilSize ) = 0;

  /**
   * @brief Update phase properties
   * @param[in] subRegion Element subregion
   */
  virtual void updatePhaseProperties( ElementSubRegionBase & subRegion ) = 0;
};

/**
 * @brief Pressure equation solver template
 * @tparam NUM_PHASES Number of phases
 * @tparam HAS_CAPILLARITY Capillarity flag
 * @tparam HAS_THERMAL Thermal flag
 */
template< int NUM_PHASES,
          bool HAS_CAPILLARITY = false,
          bool HAS_THERMAL = false >
class PressureEquationSolver : public TemplatedFlowSolver< NUM_PHASES,
                                                           EquationType::PRESSURE,
                                                           HAS_CAPILLARITY,
                                                           HAS_THERMAL >
{
public:
  using Base = TemplatedFlowSolver< NUM_PHASES, EquationType::PRESSURE, HAS_CAPILLARITY, HAS_THERMAL >;
  using FluxType = typename Base::FluxType;
  using PhaseModel = typename Base::PhaseModel;

  /**
   * @brief Constructor
   * @param[in] name Solver name
   * @param[in] parent Parent group
   */
  PressureEquationSolver( string const & name, Group * const parent );

  /**
   * @brief Assemble the pressure equation system using mimetic discretization
   * @param[in] time Current time
   * @param[in] dt Time step size
   * @param[in] domain Domain partition
   * @param[in] dofManager Degree of freedom manager
   * @param[in] localMatrix Local system matrix
   * @param[in] localRhs Local right-hand side vector
   */
  virtual void assembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;

protected:
  void computeFluxContributions( ElementSubRegionBase & subRegion,
                                 localIndex const stencilSize ) override;

  void updatePhaseProperties( ElementSubRegionBase & subRegion ) override;

  /**
   * @brief Compute pressure equation residual using mimetic discretization
   * @param[in] subRegion Element subregion
   * @param[in] pressure Current pressure
   * @param[in] pressure_n Previous time step pressure
   * @param[out] residual Computed residual
   */
  void computePressureResidual( ElementSubRegionBase & subRegion,
                                arrayView1d< real64 const > const & pressure,
                                arrayView1d< real64 const > const & pressure_n,
                                arrayView1d< real64 > const & residual );
};

/**
 * @brief Transport equation solver template
 * @tparam NUM_PHASES Number of phases
 * @tparam TRANSPORT_INDEX Index of the transport equation (0 to NUM_PHASES-2)
 * @tparam HAS_CAPILLARITY Capillarity flag
 * @tparam HAS_THERMAL Thermal flag
 */
template< int NUM_PHASES,
          int TRANSPORT_INDEX,
          bool HAS_CAPILLARITY = false,
          bool HAS_THERMAL = false >
class TransportEquationSolver : public TemplatedFlowSolver< NUM_PHASES,
                                                            EquationType::TRANSPORT,
                                                            HAS_CAPILLARITY,
                                                            HAS_THERMAL >
{
  static_assert( TRANSPORT_INDEX < NUM_PHASES - 1, "Transport index must be less than NUM_PHASES - 1" );

public:
  using Base = TemplatedFlowSolver< NUM_PHASES, EquationType::TRANSPORT, HAS_CAPILLARITY, HAS_THERMAL >;
  using FluxType = typename Base::FluxType;
  using PhaseModel = typename Base::PhaseModel;

  static constexpr int transportIndex() { return TRANSPORT_INDEX; }

  /**
   * @brief Constructor
   * @param[in] name Solver name
   * @param[in] parent Parent group
   */
  TransportEquationSolver( string const & name, Group * const parent );

  /**
   * @brief Assemble the transport equation system using mimetic discretization
   * @param[in] time Current time
   * @param[in] dt Time step size
   * @param[in] domain Domain partition
   * @param[in] dofManager Degree of freedom manager
   * @param[in] localMatrix Local system matrix
   * @param[in] localRhs Local right-hand side vector
   */
  virtual void assembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;

protected:
  void computeFluxContributions( ElementSubRegionBase & subRegion,
                                 localIndex const stencilSize ) override;

  void updatePhaseProperties( ElementSubRegionBase & subRegion ) override;

  /**
   * @brief Compute transport equation residual using mimetic discretization
   * @param[in] subRegion Element subregion
   * @param[in] saturation Current saturation for this phase
   * @param[in] saturation_n Previous time step saturation
   * @param[out] residual Computed residual
   */
  void computeTransportResidual( ElementSubRegionBase & subRegion,
                                 arrayView1d< real64 const > const & saturation,
                                 arrayView1d< real64 const > const & saturation_n,
                                 arrayView1d< real64 > const & residual );
};

/**
 * @brief Specialization for single-phase flow (NUM_PHASES = 1)
 * Only pressure equation, no transport equations needed
 */
template< bool HAS_THERMAL >
class SinglePhaseFlowSolver : public PressureEquationSolver< 1, false, HAS_THERMAL >
{
public:
  using Base = PressureEquationSolver< 1, false, HAS_THERMAL >;

  SinglePhaseFlowSolver( string const & name, Group * const parent );

  // Override to implement single-phase specific logic using mimetic discretization
  void assembleSystem( real64 const time,
                       real64 const dt,
                       DomainPartition & domain,
                       DofManager const & dofManager,
                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                       arrayView1d< real64 > const & localRhs ) override;
};

/**
 * @brief Specialization for two-phase flow (NUM_PHASES = 2)
 * Pressure equation + 1 transport equation
 */
template< bool HAS_CAPILLARITY,
          bool HAS_THERMAL >
class TwoPhaseFlowSolver : public PressureEquationSolver< 2, HAS_CAPILLARITY, HAS_THERMAL >
{
public:
  using Base = PressureEquationSolver< 2, HAS_CAPILLARITY, HAS_THERMAL >;
  using TransportSolver = TransportEquationSolver< 2, 0, HAS_CAPILLARITY, HAS_THERMAL >;

  TwoPhaseFlowSolver( string const & name, Group * const parent );

  // Combined pressure and transport assembly using mimetic discretization
  void assembleSystem( real64 const time,
                       real64 const dt,
                       DomainPartition & domain,
                       DofManager const & dofManager,
                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                       arrayView1d< real64 > const & localRhs ) override;

private:
  std::unique_ptr< TransportSolver > m_transportSolver;
};

/**
 * @brief Specialization for three-phase flow (NUM_PHASES = 3)
 * Pressure equation + 2 transport equations
 */
template< bool HAS_CAPILLARITY,
          bool HAS_THERMAL >
class ThreePhaseFlowSolver : public PressureEquationSolver< 3, HAS_CAPILLARITY, HAS_THERMAL >
{
public:
  using Base = PressureEquationSolver< 3, HAS_CAPILLARITY, HAS_THERMAL >;
  using TransportSolver1 = TransportEquationSolver< 3, 0, HAS_CAPILLARITY, HAS_THERMAL >;
  using TransportSolver2 = TransportEquationSolver< 3, 1, HAS_CAPILLARITY, HAS_THERMAL >;

  ThreePhaseFlowSolver( string const & name, Group * const parent );

  // Combined pressure and transport assembly using mimetic discretization
  void assembleSystem( real64 const time,
                       real64 const dt,
                       DomainPartition & domain,
                       DofManager const & dofManager,
                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                       arrayView1d< real64 > const & localRhs ) override;

private:
  std::unique_ptr< TransportSolver1 > m_transportSolver1;
  std::unique_ptr< TransportSolver2 > m_transportSolver2;
};

// Convenient type aliases for common use cases

// Single-phase solvers
using SinglePhaseMimetic = SinglePhaseFlowSolver< false >;
using SinglePhaseThermalMimetic = SinglePhaseFlowSolver< true >;

// Two-phase solvers
using TwoPhaseMimetic = TwoPhaseFlowSolver< false, false >;
using TwoPhaseCapillaryMimetic = TwoPhaseFlowSolver< true, false >;
using TwoPhaseThermalMimetic = TwoPhaseFlowSolver< false, true >;
using TwoPhaseThermalCapillaryMimetic = TwoPhaseFlowSolver< true, true >;

// Three-phase solvers
using ThreePhaseMimetic = ThreePhaseFlowSolver< false, false >;
using ThreePhaseCapillaryMimetic = ThreePhaseFlowSolver< true, false >;
using ThreePhaseThermalMimetic = ThreePhaseFlowSolver< false, true >;
using ThreePhaseThermalCapillaryMimetic = ThreePhaseFlowSolver< true, true >;

} // namespace geos

#endif /* GEOS_PHYSICSSOLVERS_FLUIDFLOW_TEMPLATEDFLOWSOLVERS_HPP_ */
