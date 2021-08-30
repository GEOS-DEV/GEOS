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
 * @file MultiFluidBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDBASE_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDBASE_HPP_

#include "common/DataLayouts.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/fluid/layouts.hpp"

#include "RAJA/RAJA.hpp"

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
   * @param ic component index
   * @return name of ic-th fluid component
   */
  arrayView1d< string const > componentNames() const { return m_componentNames; }

  /**
   * @return number of fluid phases in the model
   */
  integer numFluidPhases() const { return LvArray::integerConversion< integer >( m_phaseNames.size() ); }

  /**
   * @param ip phase index
   * @return name of ip-th fluid phase
   */
  arrayView1d< string const > phaseNames() const { return m_phaseNames; }

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
  { return m_phaseFraction; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseFraction_dPressure() const
  { return m_dPhaseFraction_dPressure; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseFraction_dTemperature() const
  { return m_dPhaseFraction_dTemperature; }

  arrayView4d< real64 const, multifluid::USD_PHASE_DC > dPhaseFraction_dGlobalCompFraction() const
  { return m_dPhaseFraction_dGlobalCompFraction; }

  arrayView3d< real64 const, multifluid::USD_PHASE > phaseDensity() const
  { return m_phaseDensity; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseDensity_dPressure() const
  { return m_dPhaseDensity_dPressure; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseDensity_dTemperature() const
  { return m_dPhaseDensity_dTemperature; }

  arrayView4d< real64 const, multifluid::USD_PHASE_DC > dPhaseDensity_dGlobalCompFraction() const
  { return m_dPhaseDensity_dGlobalCompFraction; }

  arrayView3d< real64 const, multifluid::USD_PHASE > phaseMassDensity() const
  { return m_phaseMassDensity; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseMassDensity_dPressure() const
  { return m_dPhaseMassDensity_dPressure; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseMassDensity_dTemperature() const
  { return m_dPhaseMassDensity_dTemperature; }

  arrayView4d< real64 const, multifluid::USD_PHASE_DC > dPhaseMassDensity_dGlobalCompFraction() const
  { return m_dPhaseMassDensity_dGlobalCompFraction; }

  arrayView3d< real64 const, multifluid::USD_PHASE > phaseViscosity() const
  { return m_phaseViscosity; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseViscosity_dPressure() const
  { return m_dPhaseViscosity_dPressure; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseViscosity_dTemperature() const
  { return m_dPhaseViscosity_dTemperature; }

  arrayView4d< real64 const, multifluid::USD_PHASE_DC > dPhaseViscosity_dGlobalCompFraction() const
  { return m_dPhaseViscosity_dGlobalCompFraction; }

  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > phaseCompFraction() const
  { return m_phaseCompFraction; }

  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > dPhaseCompFraction_dPressure() const
  { return m_dPhaseCompFraction_dPressure; }

  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > dPhaseCompFraction_dTemperature() const
  { return m_dPhaseCompFraction_dTemperature; }

  arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > dPhaseCompFraction_dGlobalCompFraction() const
  { return m_dPhaseCompFraction_dGlobalCompFraction; }

  arrayView2d< real64 const, multifluid::USD_FLUID > totalDensity() const
  { return m_totalDensity; }

  arrayView2d< real64 const, multifluid::USD_FLUID > dTotalDensity_dPressure() const
  { return m_dTotalDensity_dPressure; }

  arrayView2d< real64 const, multifluid::USD_FLUID > dTotalDensity_dTemperature() const
  { return m_dTotalDensity_dTemperature; }

  arrayView3d< real64 const, multifluid::USD_FLUID_DC > dTotalDensity_dGlobalCompFraction() const
  { return m_dTotalDensity_dGlobalCompFraction; }

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

    static constexpr char const * useMassString() { return "useMass"; }
  };

protected:

  class KernelWrapper
  {
public:

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

    /**
     * @brief Struct holding views into fluid data, used to simplify parameter passing in kernel wrapper constructors.
     * @tparam NDIM number of dimensions
     * @tparam USD unit-stride-dim of primary property and derivatives
     * @tparam USD_DC unit-stride-dim of compositional derivatives
     */
    template< int NDIM, int USD, int USD_DC >
    struct PropViews
    {
      ArrayView< real64, NDIM, USD > value;                      ///< View into property values
      ArrayView< real64, NDIM, USD > dPressure;                  ///< View into property pressure derivatives
      ArrayView< real64, NDIM, USD > dTemperature;               ///< View into property temperature derivatives
      ArrayView< real64, NDIM + 1, USD_DC > dGlobalCompFraction; ///< View into property compositional derivatives
    };

    /// Alias for phase property views struct
    using PhasePropViews = PropViews< 3, multifluid::USD_PHASE, multifluid::USD_PHASE_DC >;

    /// Alias for phase composition views struct
    using PhaseCompViews = PropViews< 4, multifluid::USD_PHASE_COMP, multifluid::USD_PHASE_COMP_DC >;

    /// Alias for fluid property views struct
    using FluidPropViews = PropViews< 2, multifluid::USD_FLUID, multifluid::USD_FLUID_DC >;

protected:

    KernelWrapper( arrayView1d< real64 const > const & componentMolarWeight,
                   bool const useMass,
                   PhasePropViews const & phaseFraction,
                   PhasePropViews const & phaseDensity,
                   PhasePropViews const & phaseMassDensity,
                   PhasePropViews const & phaseViscosity,
                   PhaseCompViews const & phaseCompFraction,
                   FluidPropViews const & totalDensity )
      : m_componentMolarWeight( componentMolarWeight ),
      m_useMass( useMass ),
      m_phaseFraction( phaseFraction ),
      m_phaseDensity( phaseDensity ),
      m_phaseMassDensity( phaseMassDensity ),
      m_phaseViscosity( phaseViscosity ),
      m_phaseCompFraction( phaseCompFraction ),
      m_totalDensity( totalDensity )
    { }

    arrayView1d< real64 const > m_componentMolarWeight;

    bool m_useMass;

    PhasePropViews m_phaseFraction;
    PhasePropViews m_phaseDensity;
    PhasePropViews m_phaseMassDensity;
    PhasePropViews m_phaseViscosity;
    PhaseCompViews m_phaseCompFraction;
    FluidPropViews m_totalDensity;

private:

    GEOSX_HOST_DEVICE
    virtual void compute( real64 const pressure,
                          real64 const temperature,
                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
                          arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & phaseCompFraction,
                          real64 & totalDensity ) const = 0;

    GEOSX_HOST_DEVICE
    virtual void compute( real64 const pressure,
                          real64 const temperature,
                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseFraction_dPressure,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseFraction_dTemperature,
                          arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseFraction_dGlobalCompFraction,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseDensity_dPressure,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseDensity_dTemperature,
                          arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseDensity_dGlobalCompFraction,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseMassDensity_dPressure,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseMassDensity_dTemperature,
                          arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseMassDensity_dGlobalCompFraction,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseViscosity_dPressure,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseViscosity_dTemperature,
                          arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseViscosity_dGlobalCompFraction,
                          arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & phaseCompFraction,
                          arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & dPhaseCompFraction_dPressure,
                          arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & dPhaseCompFraction_dTemperature,
                          arraySlice3d< real64, multifluid::USD_PHASE_COMP_DC-2 > const & dPhaseCompFraction_dGlobalCompFraction,
                          real64 & totalDensity,
                          real64 & dTotalDensity_dPressure,
                          real64 & dTotalDensity_dTemperature,
                          arraySlice1d< real64, multifluid::USD_FLUID_DC - 2 > const & dTotalDensity_dGlobalCompFraction ) const = 0;

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

  array3d< real64, multifluid::LAYOUT_PHASE > m_phaseFraction;
  array3d< real64, multifluid::LAYOUT_PHASE > m_dPhaseFraction_dPressure;
  array3d< real64, multifluid::LAYOUT_PHASE > m_dPhaseFraction_dTemperature;
  array4d< real64, multifluid::LAYOUT_PHASE_DC > m_dPhaseFraction_dGlobalCompFraction;

  array3d< real64, multifluid::LAYOUT_PHASE > m_phaseDensity;
  array3d< real64, multifluid::LAYOUT_PHASE > m_dPhaseDensity_dPressure;
  array3d< real64, multifluid::LAYOUT_PHASE > m_dPhaseDensity_dTemperature;
  array4d< real64, multifluid::LAYOUT_PHASE_DC > m_dPhaseDensity_dGlobalCompFraction;

  array3d< real64, multifluid::LAYOUT_PHASE > m_phaseMassDensity;
  array3d< real64, multifluid::LAYOUT_PHASE > m_dPhaseMassDensity_dPressure;
  array3d< real64, multifluid::LAYOUT_PHASE > m_dPhaseMassDensity_dTemperature;
  array4d< real64, multifluid::LAYOUT_PHASE_DC > m_dPhaseMassDensity_dGlobalCompFraction;

  array3d< real64, multifluid::LAYOUT_PHASE > m_phaseViscosity;
  array3d< real64, multifluid::LAYOUT_PHASE > m_dPhaseViscosity_dPressure;
  array3d< real64, multifluid::LAYOUT_PHASE > m_dPhaseViscosity_dTemperature;
  array4d< real64, multifluid::LAYOUT_PHASE_DC > m_dPhaseViscosity_dGlobalCompFraction;

  array4d< real64, multifluid::LAYOUT_PHASE_COMP > m_phaseCompFraction;
  array4d< real64, multifluid::LAYOUT_PHASE_COMP > m_dPhaseCompFraction_dPressure;
  array4d< real64, multifluid::LAYOUT_PHASE_COMP > m_dPhaseCompFraction_dTemperature;
  array5d< real64, multifluid::LAYOUT_PHASE_COMP_DC > m_dPhaseCompFraction_dGlobalCompFraction;

  array2d< real64, multifluid::LAYOUT_FLUID > m_totalDensity;
  array2d< real64, multifluid::LAYOUT_FLUID > m_dTotalDensity_dPressure;
  array2d< real64, multifluid::LAYOUT_FLUID > m_dTotalDensity_dTemperature;
  array3d< real64, multifluid::LAYOUT_FLUID_DC > m_dTotalDensity_dGlobalCompFraction;

};

} //namespace constitutive

} //namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDBASE_HPP_
