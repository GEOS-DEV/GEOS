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
 * @file MultiFluidBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDBASE_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDBASE_HPP_

#include "common/DataLayouts.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/fluid/layouts.hpp"
#include "constitutive/fluid/MultiFluidUtils.hpp"

namespace geosx
{
namespace constitutive
{

class MultiFluidBase : public ConstitutiveBase
{
public:

  MultiFluidBase( string const & name,
                  Group * const parent );

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  // *** MultiFluid-specific interface

  /**
   * @brief Maximum supported number of fluid components (species)
   *
   * @note This puts an upper bound on memory use, allowing to optimize code better
   */
  static constexpr integer MAX_NUM_COMPONENTS = 16;

  /**
   * @brief Maximum supported number of fluid phases
   *
   * @note This puts an upper bound on memory use, allowing to optimize code better
   */
  static constexpr integer MAX_NUM_PHASES = 4;

  /**
   * @return number of fluid components (species) in the model
   */
  integer numFluidComponents() const { return LvArray::integerConversion< integer >( m_componentNames.size() ); }

  /**
   * @brief Getter for the fluid component names
   * @return an array storing the component names
   */
  arrayView1d< string const > componentNames() const { return m_componentNames; }

  /**
   * @return number of fluid phases in the model
   */
  integer numFluidPhases() const { return LvArray::integerConversion< integer >( m_phaseNames.size() ); }

  /**
   * @brief Getter for the fluid phase names
   * @return an array storing the phase names
   */
  arrayView1d< string const > phaseNames() const { return m_phaseNames; }

  /**
   * @brief Getter for the water phase index
   * @return the water phase index
   */
  virtual integer getWaterPhaseIndex() const = 0;

  /**
   * @brief Get the mass flag.
   * @return boolean value indicating whether the model is using mass-based quantities (as opposed to mole-based)
   */
  bool getMassFlag() const { return m_useMass; }

  /**
   * @brief Set the mass flag.
   * @param flag boolean value indicating whether the model should use mass-based quantities (as opposed to mole-based)
   *
   * @note This affects both input (compositions) and output quantities. The flag should be set prior to calling
   * any compute or state update methods.
   */
  void setMassFlag( bool const flag ) { m_useMass = flag; }

  arrayView3d< real64 const, multifluid::USD_PHASE > phaseFraction() const
  { return m_phaseFraction.value; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseFraction_dPressure() const
  { return m_phaseFraction.dPres; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseFraction_dTemperature() const
  { return m_phaseFraction.dTemp; }

  arrayView4d< real64 const, multifluid::USD_PHASE_DC > dPhaseFraction_dGlobalCompFraction() const
  { return m_phaseFraction.dComp; }

  arrayView3d< real64 const, multifluid::USD_PHASE > phaseDensity() const
  { return m_phaseDensity.value; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseDensity_dPressure() const
  { return m_phaseDensity.dPres; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseDensity_dTemperature() const
  { return m_phaseDensity.dTemp; }

  arrayView4d< real64 const, multifluid::USD_PHASE_DC > dPhaseDensity_dGlobalCompFraction() const
  { return m_phaseDensity.dComp; }

  arrayView3d< real64 const, multifluid::USD_PHASE > phaseMassDensity() const
  { return m_phaseMassDensity.value; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseMassDensity_dPressure() const
  { return m_phaseMassDensity.dPres; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseMassDensity_dTemperature() const
  { return m_phaseMassDensity.dTemp; }

  arrayView4d< real64 const, multifluid::USD_PHASE_DC > dPhaseMassDensity_dGlobalCompFraction() const
  { return m_phaseMassDensity.dComp; }

  arrayView3d< real64 const, multifluid::USD_PHASE > phaseViscosity() const
  { return m_phaseViscosity.value; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseViscosity_dPressure() const
  { return m_phaseViscosity.dPres; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseViscosity_dTemperature() const
  { return m_phaseViscosity.dTemp; }

  arrayView4d< real64 const, multifluid::USD_PHASE_DC > dPhaseViscosity_dGlobalCompFraction() const
  { return m_phaseViscosity.dComp; }

  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > phaseCompFraction() const
  { return m_phaseCompFraction.value; }

  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > dPhaseCompFraction_dPressure() const
  { return m_phaseCompFraction.dPres; }

  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > dPhaseCompFraction_dTemperature() const
  { return m_phaseCompFraction.dTemp; }

  arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > dPhaseCompFraction_dGlobalCompFraction() const
  { return m_phaseCompFraction.dComp; }

  arrayView2d< real64 const, multifluid::USD_FLUID > totalDensity() const
  { return m_totalDensity.value; }

  arrayView2d< real64 const, multifluid::USD_FLUID > dTotalDensity_dPressure() const
  { return m_totalDensity.dPres; }

  arrayView2d< real64 const, multifluid::USD_FLUID > dTotalDensity_dTemperature() const
  { return m_totalDensity.dTemp; }

  arrayView3d< real64 const, multifluid::USD_FLUID_DC > dTotalDensity_dGlobalCompFraction() const
  { return m_totalDensity.dComp; }

  arrayView2d< real64 const, multifluid::USD_FLUID > initialTotalMassDensity() const
  { return m_initialTotalMassDensity.toViewConst(); }

  arrayView3d< real64 const, multifluid::USD_PHASE > phaseEnthalpy() const
  { return m_phaseEnthalpy.value; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseEnthalpy_dPressure() const
  { return m_phaseEnthalpy.dPres; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseEnthalpy_dTemperature() const
  { return m_phaseEnthalpy.dTemp; }

  arrayView4d< real64 const, multifluid::USD_PHASE_DC > dPhaseEnthalpy_dGlobalCompFraction() const
  { return m_phaseEnthalpy.dComp; }

  arrayView3d< real64 const, multifluid::USD_PHASE > phaseInternalEnergy() const
  { return m_phaseInternalEnergy.value; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseInternalEnergy_dPressure() const
  { return m_phaseInternalEnergy.dPres; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseInternalEnergy_dTemperature() const
  { return m_phaseInternalEnergy.dTemp; }

  arrayView4d< real64 const, multifluid::USD_PHASE_DC > dPhaseInternalEnergy_dGlobalCompFraction() const
  { return m_phaseInternalEnergy.dComp; }

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * componentNamesString() { return "componentNames"; }
    static constexpr char const * componentMolarWeightString() { return "componentMolarWeight"; }

    static constexpr char const * phaseNamesString() { return "phaseNames"; }

    static constexpr char const * phaseFractionString() { return "phaseFraction"; } // xi_p
    static constexpr char const * dPhaseFraction_dPressureString() { return "dPhaseFraction_dPressure"; } // dXi_p/dP
    static constexpr char const * dPhaseFraction_dTemperatureString() { return "dPhaseFraction_dTemperature"; } // dXi_p/dT
    static constexpr char const * dPhaseFraction_dGlobalCompFractionString() { return "dPhaseFraction_dGlobalCompFraction"; } // dXi_p/dz

    static constexpr char const * phaseDensityString() { return "phaseDensity"; } // rho_p
    static constexpr char const * dPhaseDensity_dPressureString() { return "dPhaseDensity_dPressure"; } // dRho_p/dP
    static constexpr char const * dPhaseDensity_dTemperatureString() { return "dPhaseDensity_dTemperature"; } // dRho_p/dT
    static constexpr char const * dPhaseDensity_dGlobalCompFractionString() { return "dPhaseDensity_dGlobalCompFraction"; } // dRho_p/dz

    static constexpr char const * phaseMassDensityString() { return "phaseMassDensity"; } // rho_p
    static constexpr char const * dPhaseMassDensity_dPressureString() { return "dPhaseMassDensity_dPressure"; } // dRho_p/dP
    static constexpr char const * dPhaseMassDensity_dTemperatureString() { return "dPhaseMassDensity_dTemperature"; } // dRho_p/dT
    static constexpr char const * dPhaseMassDensity_dGlobalCompFractionString() { return "dPhaseMassDensity_dGlobalCompFraction"; } // dRho_p/dz

    static constexpr char const * phaseViscosityString() { return "phaseViscosity"; } // mu_p
    static constexpr char const * dPhaseViscosity_dPressureString() { return "dPhaseViscosity_dPressure"; } // dMu_p/dP
    static constexpr char const * dPhaseViscosity_dTemperatureString() { return "dPhaseViscosity_dTemperature"; } // dMu_p/dT
    static constexpr char const * dPhaseViscosity_dGlobalCompFractionString() { return "dPhaseViscosity_dGlobalCompFraction"; } // dMu_p/dz

    static constexpr char const * phaseCompFractionString() { return "phaseCompFraction"; } // x_cp
    static constexpr char const * dPhaseCompFraction_dPressureString() { return "dPhaseCompFraction_dPressure"; } // dx_cp/dP
    static constexpr char const * dPhaseCompFraction_dTemperatureString() { return "dPhaseCompFraction_dTemperature"; } // dx_cp/dT
    static constexpr char const * dPhaseCompFraction_dGlobalCompFractionString() { return "dPhaseCompFraction_dGlobalCompFraction"; } // dx_cp/dz

    static constexpr char const * totalDensityString() { return "totalDensity"; } // rho_t
    static constexpr char const * dTotalDensity_dPressureString() { return "dTotalDensity_dPressure"; } // dRho_t/dP
    static constexpr char const * dTotalDensity_dTemperatureString() { return "dTotalDensity_dTemperature"; } // dRho_t/dT
    static constexpr char const * dTotalDensity_dGlobalCompFractionString() { return "dTotalDensity_dGlobalCompFraction"; } // dRho_t/dz

    static constexpr char const * initialTotalMassDensityString() { return "initialTotalMassDensity"; } // rho^int_t

    static constexpr char const * phaseEnthalpyString() { return "phaseEnthalpy"; } // H_p
    static constexpr char const * dPhaseEnthalpy_dPressureString() { return "dPhaseEnthalpy_dPressure"; } // dH_p/dP
    static constexpr char const * dPhaseEnthalpy_dTemperatureString() { return "dPhaseEnthalpy_dTemperature"; } // dH_p/dT
    static constexpr char const * dPhaseEnthalpy_dGlobalCompFractionString() { return "dPhaseEnthalpy_dGlobalCompFraction"; } // dH_p/dz

    static constexpr char const * phaseInternalEnergyString() { return "phaseInternalEnergy"; } // U_p
    static constexpr char const * dPhaseInternalEnergy_dPressureString() { return "dPhaseInternalEnergy_dPressure"; } // dU_p/dP
    static constexpr char const * dPhaseInternalEnergy_dTemperatureString() { return "dPhaseInternalEnergy_dTemperature"; } // dU_p/dT
    static constexpr char const * dPhaseInternalEnergy_dGlobalCompFractionString() { return "dPhaseInternalEnergy_dGlobalCompFraction"; } // dU_p/dz

    static constexpr char const * useMassString() { return "useMass"; }
  };

protected:

  using PhaseProp = MultiFluidVar< real64, 3, multifluid::LAYOUT_PHASE, multifluid::LAYOUT_PHASE_DC >;
  using PhaseComp = MultiFluidVar< real64, 4, multifluid::LAYOUT_PHASE_COMP, multifluid::LAYOUT_PHASE_COMP_DC >;
  using FluidProp = MultiFluidVar< real64, 2, multifluid::LAYOUT_FLUID, multifluid::LAYOUT_FLUID_DC >;

  class KernelWrapper
  {
public:

    /// @cond DO_NOT_DOCUMENT
    /// We need these SMFs to avoid host-device errors with CUDA.
    KernelWrapper() = default;
    KernelWrapper( KernelWrapper const & ) = default;
    KernelWrapper & operator=( KernelWrapper const & ) = default;
    KernelWrapper & operator=( KernelWrapper && ) = default;
    /// @endcond

    /**
     * @brief Get number of elements in this wrapper.
     * @return number of elements
     */
    GEOSX_HOST_DEVICE
    localIndex numElems() const { return m_phaseFraction.value.size( 0 ); }

    /**
     * @brief Get number of gauss points per element.
     * @return number of gauss points per element
     */
    GEOSX_HOST_DEVICE
    localIndex numGauss() const { return m_phaseFraction.value.size( 1 ); }

    /**
     * @brief Get number of fluid components.
     * @return number of components
     */
    GEOSX_HOST_DEVICE
    integer numComponents() const { return LvArray::integerConversion< integer >( m_componentMolarWeight.size() ); }

    /**
     * @brief Get number of fluid phases.
     * @return number of phases
     */
    GEOSX_HOST_DEVICE
    integer numPhases() const { return LvArray::integerConversion< integer >( m_phaseFraction.value.size( 2 ) ); }

protected:

    KernelWrapper( arrayView1d< real64 const > componentMolarWeight,
                   bool const useMass,
                   PhaseProp::ViewType phaseFraction,
                   PhaseProp::ViewType phaseDensity,
                   PhaseProp::ViewType phaseMassDensity,
                   PhaseProp::ViewType phaseViscosity,
                   PhaseProp::ViewType phaseEnthalpy,
                   PhaseProp::ViewType phaseInternalEnergy,
                   PhaseComp::ViewType phaseCompFraction,
                   FluidProp::ViewType totalDensity )
      : m_componentMolarWeight( std::move( componentMolarWeight ) ),
      m_useMass( useMass ),
      m_phaseFraction( std::move( phaseFraction ) ),
      m_phaseDensity( std::move( phaseDensity ) ),
      m_phaseMassDensity( std::move( phaseMassDensity ) ),
      m_phaseViscosity( std::move( phaseViscosity ) ),
      m_phaseEnthalpy( std::move( phaseEnthalpy ) ),
      m_phaseInternalEnergy( std::move( phaseInternalEnergy ) ),
      m_phaseCompFraction( std::move( phaseCompFraction ) ),
      m_totalDensity( std::move( totalDensity ) )
    { }

    arrayView1d< real64 const > m_componentMolarWeight;

    bool m_useMass;

    PhaseProp::ViewType m_phaseFraction;
    PhaseProp::ViewType m_phaseDensity;
    PhaseProp::ViewType m_phaseMassDensity;
    PhaseProp::ViewType m_phaseViscosity;
    PhaseProp::ViewType m_phaseEnthalpy;
    PhaseProp::ViewType m_phaseInternalEnergy;
    PhaseComp::ViewType m_phaseCompFraction;
    FluidProp::ViewType m_totalDensity;

private:

    GEOSX_HOST_DEVICE
    virtual void compute( real64 const pressure,
                          real64 const temperature,
                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseEnthalpy,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseInternalEnergy,
                          arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & phaseCompFraction,
                          real64 & totalDensity ) const = 0;

    GEOSX_HOST_DEVICE
    virtual void compute( real64 const pressure,
                          real64 const temperature,
                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                          PhaseProp::SliceType const phaseFraction,
                          PhaseProp::SliceType const phaseDensity,
                          PhaseProp::SliceType const phaseMassDensity,
                          PhaseProp::SliceType const phaseViscosity,
                          PhaseProp::SliceType const phaseEnthalpy,
                          PhaseProp::SliceType const phaseInternalEnergy,
                          PhaseComp::SliceType const phaseCompFraction,
                          FluidProp::SliceType const totalDensity ) const = 0;

    GEOSX_HOST_DEVICE
    virtual void update( localIndex const k,
                         localIndex const q,
                         real64 const pressure,
                         real64 const temperature,
                         arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const = 0;
  };

private:

  /**
   * @brief Function called internally to resize member arrays
   * @param size primary dimension (e.g. number of cells)
   * @param numPts secondary dimension (e.g. number of gauss points per cell)
   */
  void resizeFields( localIndex const size, localIndex const numPts );

  /**
   * @brief Called internally to set array dim labels.
   */
  void setLabels();

protected:

  virtual void postProcessInput() override;

  // flag indicating whether input/output component fractions are treated as mass fractions
  int m_useMass;

  // general fluid composition information

  array1d< string > m_componentNames;
  array1d< real64 > m_componentMolarWeight;
  array1d< string > m_phaseNames;

  // constitutive data

  PhaseProp m_phaseFraction;
  PhaseProp m_phaseDensity;
  PhaseProp m_phaseMassDensity;
  PhaseProp m_phaseViscosity;
  PhaseProp m_phaseEnthalpy;
  PhaseProp m_phaseInternalEnergy;
  PhaseComp m_phaseCompFraction;
  FluidProp m_totalDensity;

  // initial data (used to compute the body force in the poromechanics solver)

  array2d< real64, multifluid::LAYOUT_FLUID > m_initialTotalMassDensity;

};

} //namespace constitutive

} //namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDBASE_HPP_
