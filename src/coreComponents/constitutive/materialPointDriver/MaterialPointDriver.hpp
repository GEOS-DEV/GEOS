/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MaterialPointDriver.hpp
 *
 * Stand-alone material-point driver utilities for exercising GEOS/MPM
 * constitutive models outside the MPM grid solver.  These classes are not
 * referenced by SolidMechanicsMPM and are built only when the optional
 * material-point driver target is enabled.
 */

#ifndef GEOS_CONSTITUTIVE_MATERIALPOINTDRIVER_MATERIALPOINTDRIVER_HPP_
#define GEOS_CONSTITUTIVE_MATERIALPOINTDRIVER_MATERIALPOINTDRIVER_HPP_

#include "common/DataTypes.hpp"
#include "constitutive/ContinuumBase.hpp"

#include <array>
#include <functional>
#include <string>
#include <vector>

namespace geos
{
namespace constitutive
{
namespace materialPointDriver
{

/** GEOS solid Voigt order: xx, yy, zz, yz, xz, xy. */
static constexpr localIndex numVoigtComponents = 6;
static constexpr localIndex numSpatialDims = 3;

using Vector3 = std::array< real64, 3 >;
using Matrix3 = std::array< std::array< real64, 3 >, 3 >;
using Voigt6 = std::array< real64, 6 >;

/**
 * Controls how the row-wise material frame is transported by the deformation.
 */
enum class MaterialFrameUpdateMode
{
  fixed,       ///< Keep the supplied reference frame fixed.
  rotation,    ///< Push the reference frame by the polar rotation.
  fiber,       ///< Push rows by F, then normalize.
  normal,      ///< Push row 0 by cof(F); rebuild tangents.
  graphite,    ///< Push row 0 by cof(F), tangent by F, rebuild row 2.
  mpmCofactor  ///< MPM non-fiber update: push all rows by cof(F), then normalize.
};

/** Energy integration mode owned by this driver. */
enum class EnergyMode
{
  off,
  stressPower,
  materialOwned
};

/** Temperature mode used while feeding constitutive dependencies. */
enum class TemperatureMode
{
  prescribed,
  isothermal,
  adiabatic,
  materialOwned
};

/** How a tensor component is controlled over one time step. */
enum class ControlMode
{
  strainIncrement,
  strainRate,
  trueStrainRate,
  stress
};

/** Per-step scalar control for one GEOS Voigt component. */
struct ComponentControl
{
  localIndex component = 0;       ///< GEOS Voigt component index.
  ControlMode mode = ControlMode::strainIncrement;
  real64 value = 0.0;             ///< Strain increment or rate, depending on mode.
  real64 target = 0.0;            ///< Stress target for stress-controlled components.
};

/** Newton controls used for mixed stress/strain loading. */
struct NonlinearSolveOptions
{
  localIndex maxIterations = 24;
  real64 absoluteTolerance = 1.0e-10;
  real64 relativeTolerance = 1.0e-8;
  real64 finiteDifferenceScale = 1.0e-7;
  real64 minimumFiniteDifference = 1.0e-12;
  real64 maximumCorrectionNorm = 0.25;
  bool useLineSearch = true;
};

/** Energy and temperature options for the driver-owned energy update. */
struct EnergyOptions
{
  EnergyMode mode = EnergyMode::stressPower;
  TemperatureMode temperatureMode = TemperatureMode::isothermal;
  real64 retentionFactor = 1.0;
  real64 heatCapacity = 1.0;
  bool outputAccumulatedStressPower = true;
};

/** Diagnostics from one attempted step. */
struct StepDiagnostics
{
  bool converged = false;
  localIndex iterations = 0;
  real64 initialResidualNorm = 0.0;
  real64 finalResidualNorm = 0.0;
  real64 energyIncrement = 0.0;
  std::string message;
};

/** Complete state carried by a single material point. */
struct MaterialPointState
{
  real64 time = 0.0;
  real64 dt = 0.0;

  Matrix3 deformationGradient{};     ///< F_n or F_{n+1} after a step.
  Matrix3 oldDeformationGradient{};  ///< Beginning-of-step F_n.
  Matrix3 fDot{};                    ///< (F_{n+1} - F_n) / dt.
  Matrix3 velocityGradient{};        ///< Step velocity gradient L.
  Matrix3 rotationBeginning{};       ///< Polar rotation of F_n.
  Matrix3 rotationEnd{};             ///< Polar rotation of F_{n+1}.

  Matrix3 referenceMaterialFrame{};  ///< User supplied row-wise GEOS/PFW material frame.
  Matrix3 materialFrame{};           ///< Current row-wise material frame.
  MaterialFrameUpdateMode materialFrameUpdate = MaterialFrameUpdateMode::fixed;

  Voigt6 strainIncrement{};          ///< sym(L) dt in GEOS Voigt engineering-shear convention.
  Voigt6 accumulatedStrain{};        ///< Accumulated small/log-style diagnostic strain.
  Voigt6 oldStress{};                ///< Stress at the start of the accepted step.
  Voigt6 stress{};                   ///< Stress at the end of the accepted step.

  real64 referenceDensity = 1.0;
  real64 density = 1.0;
  real64 referenceVolume = 1.0;
  real64 volume = 1.0;
  real64 jacobian = 1.0;
  real64 lengthScale = 1.0;
  real64 strengthScale = 1.0;

  real64 temperature = 300.0;
  real64 prescribedTemperature = 300.0;
  real64 temperatureRate = 0.0;
  real64 specificInternalEnergy = 0.0;
  real64 accumulatedStressPower = 0.0;
};

/** One load step specification. */
struct StepPlan
{
  real64 dt = 0.0;
  std::vector< ComponentControl > controls;
  NonlinearSolveOptions solveOptions;
  EnergyOptions energyOptions;
};

/**
 * Helper for building and transporting a GEOS row-wise material frame.
 */
class MaterialFrameController
{
public:
  static Matrix3 identityFrame();

  static Matrix3 fromNormal( Vector3 const & normal,
                             Vector3 tangentHint = { { 1.0, 0.0, 0.0 } } );

  static Matrix3 normalizedFrame( Matrix3 const & frame );

  static Matrix3 update( Matrix3 const & referenceFrame,
                         Matrix3 const & deformationGradient,
                         Matrix3 const & polarRotation,
                         MaterialFrameUpdateMode mode );
};

/** Driver-owned stress power and internal-energy update. */
class EnergyIntegrator
{
public:
  static real64 stressPower( Voigt6 const & stress,
                             Matrix3 const & velocityGradient );

  static real64 integrateStressPowerTrapezoid( MaterialPointState & state,
                                               Voigt6 const & stressBeginning,
                                               Voigt6 const & stressEnd,
                                               EnergyOptions const & options );
};

/**
 * Captures the common mutable arrays of a single-point GEOS continuum model.
 * This is used to make nonlinear residual evaluations side-effect free.
 */
class ConstitutivePointSnapshot
{
public:
  void capture( ContinuumBase & model );
  void restore( ContinuumBase & model ) const;

private:
  struct RealArraySnapshot
  {
    string name;
    int rank = 0;
    localIndex dims[3] = { 0, 0, 0 };
    std::vector< real64 > values;
  };

  bool m_hasStress = false;
  Voigt6 m_stress{};
  Voigt6 m_oldStress{};
  std::vector< RealArraySnapshot > m_realArrays;
};

/**
 * Local material-point driver that can be called from an executable, a unit
 * test, or a fitting harness.  It does not own, construct, or register the
 * constitutive model; callers provide a one-point GEOS ContinuumBase clone.
 */
class MaterialPointDriver
{
public:
  /** Initialize state tensors and copy common dependencies into a one-point model. */
  static void initializeOnePointModel( ContinuumBase & model,
                                       MaterialPointState & state );

  /** Apply one mixed stress/strain-control step and commit the material state. */
  static StepDiagnostics advanceOneStep( ContinuumBase & model,
                                         MaterialPointState & state,
                                         StepPlan const & plan );

  /** Copy the current driver state into common constitutive dependency wrappers. */
  static void pushDependencies( ContinuumBase & model,
                                MaterialPointState const & state );

  /** Pull common material-owned outputs back into the driver state. */
  static void pullOutputs( ContinuumBase & model,
                           MaterialPointState & state );

private:
  static void buildVelocityGradient( StepPlan const & plan,
                                     std::vector< real64 > const & unknownStrainIncrements,
                                     Matrix3 & velocityGradient );

  static void updateKinematicsForTrial( MaterialPointState & state,
                                        real64 dt );

  static void callConstitutiveUpdate( ContinuumBase & model,
                                      MaterialPointState & state,
                                      real64 dt,
                                      bool commit );

  static std::vector< localIndex > stressControlledComponents( StepPlan const & plan );

  static std::vector< real64 > residualForUnknowns( ContinuumBase & model,
                                                    MaterialPointState const & stateBeginning,
                                                    ConstitutivePointSnapshot const & modelSnapshot,
                                                    StepPlan const & plan,
                                                    std::vector< real64 > const & unknownStrainIncrements,
                                                    bool commit,
                                                    MaterialPointState * trialStateOut );
};

} // namespace materialPointDriver
} // namespace constitutive
} // namespace geos

#endif // GEOS_CONSTITUTIVE_MATERIALPOINTDRIVER_MATERIALPOINTDRIVER_HPP_
