/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MultiFluidBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDBASE_HPP
#define GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDBASE_HPP

#include "constitutive/ConstitutiveBase.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{

namespace constitutive
{

namespace detail
{

template< typename T, int DIM >
struct ArraySlice_helper
{
  using type = ArraySlice< T, DIM >;
};

// an array slice of DIM=0 decays to a reference to scalar
template< typename T >
struct ArraySlice_helper< T, 0 >
{
  using type = T &;
};

// an array1 slice of DIM=1 uses specialization (possibly a raw pointer)
template< typename T >
struct ArraySlice_helper< T, 1 >
{
  using type = arraySlice1d< T >;
};

}

template< int DIM >
using real_ArraySlice = typename detail::ArraySlice_helper< real64, DIM >::type;

template< int DIM >
using real_array_const_slice = typename detail::ArraySlice_helper< real64 const, DIM >::type;

// helper struct to represent a var and its derivatives
template< int DIM >
struct CompositionalVarContainer
{
  real_ArraySlice< DIM >   value; // variable value
  real_ArraySlice< DIM >   dPres; // derivative w.r.t. pressure
  real_ArraySlice< DIM >   dTemp; // derivative w.r.t. temperature
  real_ArraySlice< DIM+1 > dComp; // derivative w.r.t. composition
};

template< int DIM >
struct CompositionalVarConstContainer
{
  real_array_const_slice< DIM >   value; // variable value
  real_array_const_slice< DIM >   dPres; // derivative w.r.t. pressure
  real_array_const_slice< DIM >   dTemp; // derivative w.r.t. temperature
  real_array_const_slice< DIM+1 > dComp; // derivative w.r.t. composition
};

class MultiFluidBase : public ConstitutiveBase
{
public:

  MultiFluidBase( std::string const & name, Group * const parent );

  virtual ~MultiFluidBase() override;

  virtual void DeliverClone( string const & name,
                             Group * const parent,
                             std::unique_ptr< ConstitutiveBase > & clone ) const override;

  virtual void AllocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  // *** MultiFluid-specific interface

  /**
   * @brief Maximum supported number of fluid components (species)
   *
   * @note This puts an upper bound on memory use, allowing to optimize code better
   */
  static constexpr localIndex MAX_NUM_COMPONENTS = 32;

  /**
   * @brief Maximum supported number of fluid phases
   *
   * @note This puts an upper bound on memory use, allowing to optimize code better
   */
  static constexpr localIndex MAX_NUM_PHASES = 4;

  /**
   * @brief Perform a single point constitutive update.
   * @param[in] pressure target pressure value
   * @param[in] temperature target temperature value
   * @param[in] composition target fluid composition array
   * @param[in] k first constitutive index (e.g. elem index)
   * @param[in] q second constitutive index (e.g. quadrature index)
   *
   * @note This function should generally not be called from a kernel, use BatchUpdate instead
   */
  virtual void PointUpdate( real64 const & pressure,
                            real64 const & temperature,
                            arraySlice1d< real64 const > const & composition,
                            localIndex const k,
                            localIndex const q ) = 0;

  /**
   * @brief Perform a batch constitutive update (all points).
   * @param[in] pressure array containing target pressure values
   * @param[in] temperature array containing target temperature values
   * @param[in] composition 2D array containing target fluid composition values
   */
  virtual void BatchUpdate( arrayView1d< real64 const > const & pressure,
                            arrayView1d< real64 const > const & temperature,
                            arrayView2d< real64 const > const & composition ) = 0;

  /**
   * @brief Compute constitutive values at a single point.
   * @param[in]  pressure target pressure value
   * @param[in]  temperature target temperature value
   * @param[in]  composition target fluid composition array
   * @param[out] phaseFraction phase fractions
   * @param[out] dPhaseFraction_dPressure derivatives of phase fractions w.r.t. pressure
   * @param[out] dPhaseFraction_dTemperature derivatives of phase fractions w.r.t. temperature
   * @param[out] dPhaseFraction_dGlobalCompFraction derivatives of phase fractions w.r.t. composition
   * @param[out] phaseDensity phase densitites
   * @param[out] dPhaseDensity_dPressure derivatives of phase densitites w.r.t. pressure
   * @param[out] dPhaseDensity_dTemperature derivatives of phase densitites w.r.t. temperature
   * @param[out] dPhaseDensity_dGlobalCompFraction derivatives of phase densitites w.r.t. composition
   * @param[out] phaseViscosity phase viscosities
   * @param[out] dPhaseViscosity_dPressure derivatives of phase viscosities w.r.t. pressure
   * @param[out] dPhaseViscosity_dTemperature derivatives of phase viscosities w.r.t. temperature
   * @param[out] dPhaseViscosity_dGlobalCompFraction derivatives of phase viscosities w.r.t. composition
   * @param[out] phaseCompFraction phase compositions
   * @param[out] dPhaseCompFraction_dPressure derivatives of phase compositions w.r.t. pressure
   * @param[out] dPhaseCompFraction_dTemperature derivatives of phase compositions w.r.t. temperature
   * @param[out] dPhaseCompFraction_dGlobalCompFraction derivatives of phase compositions w.r.t. composition
   * @param[out] totalDensity total fluid mixture density
   * @param[out] dTotalDensity_dPressure derivatives of total density w.r.t. pressure
   * @param[out] dTotalDensity_dTemperature derivatives of total density w.r.t. temperature
   * @param[out] dTotalDensity_dGlobalCompFraction derivatives of total density w.r.t. composition
   *
   * @note This function should only be called in extremely rare cases, when constitutive state
   * needs to be evaluated at a point where constitutive model does not have storage allocated.
   * It should not be called from kernels since it is virtual.
   */
  virtual void Compute( real64 const & pressure,
                        real64 const & temperature,
                        arraySlice1d< real64 const > const & composition,
                        arraySlice1d< real64 > const & phaseFraction,
                        arraySlice1d< real64 > const & dPhaseFraction_dPressure,
                        arraySlice1d< real64 > const & dPhaseFraction_dTemperature,
                        arraySlice2d< real64 > const & dPhaseFraction_dGlobalCompFraction,
                        arraySlice1d< real64 > const & phaseDensity,
                        arraySlice1d< real64 > const & dPhaseDensity_dPressure,
                        arraySlice1d< real64 > const & dPhaseDensity_dTemperature,
                        arraySlice2d< real64 > const & dPhaseDensity_dGlobalCompFraction,
                        arraySlice1d< real64 > const & phaseViscosity,
                        arraySlice1d< real64 > const & dPhaseViscosity_dPressure,
                        arraySlice1d< real64 > const & dPhaseViscosity_dTemperature,
                        arraySlice2d< real64 > const & dPhaseViscosity_dGlobalCompFraction,
                        arraySlice2d< real64 > const & phaseCompFraction,
                        arraySlice2d< real64 > const & dPhaseCompFraction_dPressure,
                        arraySlice2d< real64 > const & dPhaseCompFraction_dTemperature,
                        arraySlice3d< real64 > const & dPhaseCompFraction_dGlobalCompFraction,
                        real64 & totalDensity,
                        real64 & dTotalDensity_dPressure,
                        real64 & dTotalDensity_dTemperature,
                        arraySlice1d< real64 > const & dTotalDensity_dGlobalCompFraction ) const = 0;

  /**
   * @return number of fluid components (species) in the model
   */
  localIndex numFluidComponents() const { return m_componentNames.size(); }

  /**
   * @param ic component index
   * @return name of ic-th fluid component
   */
  arrayView1d< string const > const & componentNames() const { return m_componentNames; }

  /**
   * @return number of fluid phases in the model
   */
  localIndex numFluidPhases() const { return m_phaseNames.size(); }

  /**
   * @param ip phase index
   * @return name of ip-th fluid phase
   */
  arrayView1d< string const > const & phaseNames() const { return m_phaseNames; }

  /**
   * @brief Get the mass flag.
   * @return boolean value indicating whether the model is using mass-based quantities (as opposed to mole-based)
   */
  bool getMassFlag() const;

  /**
   * @brief Set the mass flag.
   * @param flag boolean value indicating whether the model should use mass-based quantities (as opposed to mole-based)
   *
   * @note This affects both input (compositions) and output quantities. The flag should be set prior to calling
   * any compute or state update methods.
   */
  void setMassFlag( bool flag );

  arrayView3d< real64 const > const & phaseFraction() const { return m_phaseFraction; }
  arrayView3d< real64 const > const & dPhaseFraction_dPressure() const { return m_dPhaseFraction_dPressure; }
  arrayView3d< real64 const > const & dPhaseFraction_dTemperature() const { return m_dPhaseFraction_dTemperature; }
  arrayView4d< real64 const > const & dPhaseFraction_dGlobalCompFraction() const { return m_dPhaseFraction_dGlobalCompFraction; }

  arrayView3d< real64 const > const & phaseDensity() const { return m_phaseDensity; }
  arrayView3d< real64 const > const & dPhaseDensity_dPressure() const { return m_dPhaseDensity_dPressure; }
  arrayView3d< real64 const > const & dPhaseDensity_dTemperature() const { return m_dPhaseDensity_dTemperature; }
  arrayView4d< real64 const > const & dPhaseDensity_dGlobalCompFraction() const { return m_dPhaseDensity_dGlobalCompFraction; }

  arrayView3d< real64 const > const & phaseViscosity() const { return m_phaseViscosity; }
  arrayView3d< real64 const > const & dPhaseViscosity_dPressure() const { return m_dPhaseViscosity_dPressure; }
  arrayView3d< real64 const > const & dPhaseViscosity_dTemperature() const { return m_dPhaseViscosity_dTemperature; }
  arrayView4d< real64 const > const & dPhaseViscosity_dGlobalCompFraction() const { return m_dPhaseViscosity_dGlobalCompFraction; }

  arrayView4d< real64 const > const & phaseCompFraction() const { return m_phaseCompFraction; }
  arrayView4d< real64 const > const & dPhaseCompFraction_dPressure() const { return m_dPhaseCompFraction_dPressure; }
  arrayView4d< real64 const > const & dPhaseCompFraction_dTemperature() const { return m_dPhaseCompFraction_dTemperature; }
  arrayView5d< real64 const > const & dPhaseCompFraction_dGlobalCompFraction() const { return m_dPhaseCompFraction_dGlobalCompFraction; }

  arrayView2d< real64 const > const & totalDensity() const { return m_totalDensity; }
  arrayView2d< real64 const > const & dTotalDensity_dPressure() const { return m_dTotalDensity_dPressure; }
  arrayView2d< real64 const > const & dTotalDensity_dTemperature() const { return m_dTotalDensity_dTemperature; }
  arrayView3d< real64 const > const & dTotalDensity_dGlobalCompFraction() const { return m_dTotalDensity_dGlobalCompFraction; }

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr auto componentNamesString       = "componentNames";
    static constexpr auto componentMolarWeightString = "componentMolarWeight";

    static constexpr auto phaseNamesString     = "phaseNames";

    static constexpr auto phaseFractionString                            = "phaseFraction";                          // xi_p
    static constexpr auto dPhaseFraction_dPressureString                 = "dPhaseFraction_dPressure";               // dXi_p/dP
    static constexpr auto dPhaseFraction_dTemperatureString              = "dPhaseFraction_dTemperature";            // dXi_p/dT
    static constexpr auto dPhaseFraction_dGlobalCompFractionString       = "dPhaseFraction_dGlobalCompFraction";     // dXi_p/dz

    static constexpr auto phaseDensityString                             = "phaseDensity";                           // rho_p
    static constexpr auto dPhaseDensity_dPressureString                  = "dPhaseDensity_dPressure";                // dRho_p/dP
    static constexpr auto dPhaseDensity_dTemperatureString               = "dPhaseDensity_dTemperature";             // dRho_p/dT
    static constexpr auto dPhaseDensity_dGlobalCompFractionString        = "dPhaseDensity_dGlobalCompFraction";      // dRho_p/dz

    static constexpr auto phaseViscosityString                           = "phaseViscosity";                         // mu_p
    static constexpr auto dPhaseViscosity_dPressureString                = "dPhaseViscosity_dPressure";              // dMu_p/dP
    static constexpr auto dPhaseViscosity_dTemperatureString             = "dPhaseViscosity_dTemperature";           // dMu_p/dT
    static constexpr auto dPhaseViscosity_dGlobalCompFractionString      = "dPhaseViscosity_dGlobalCompFraction";    // dMu_p/dz

    static constexpr auto phaseCompFractionString                        = "phaseCompFraction";                      // x_cp
    static constexpr auto dPhaseCompFraction_dPressureString             = "dPhaseCompFraction_dPressure";           // dx_cp/dP
    static constexpr auto dPhaseCompFraction_dTemperatureString          = "dPhaseCompFraction_dTemperature";        // dx_cp/dT
    static constexpr auto dPhaseCompFraction_dGlobalCompFractionString   = "dPhaseCompFraction_dGlobalCompFraction"; // dx_cp/dz

    static constexpr auto totalDensityString                             = "totalDensity";                           // rho_t
    static constexpr auto dTotalDensity_dPressureString                  = "dTotalDensity_dPressure";                // dRho_t/dP
    static constexpr auto dTotalDensity_dTemperatureString               = "dTotalDensity_dTemperature";             // dRho_t/dT
    static constexpr auto dTotalDensity_dGlobalCompFractionString        = "dTotalDensity_dGlobalCompFraction";      // dRho_t/dz
  } viewKeysMultiFluidBase;

protected:

  virtual void PostProcessInput() override;

  /**
   * @brief Function to batch process constitutive updates via a kernel launch.
   * @tparam LEAFCLASS The derived class that provides the functions for usein the kernel
   * @tparam ARGS Parameter pack for arbitrary number of arbitrary types for the function parameter list
   * @param pressure array containing the pressure values
   * @param temperature array containing the temperature values
   * @param composition array containing the fluid composition
   * @param args arbitrary number of arbitrary types that are passed to the kernel
   */
  template< typename LEAFCLASS, typename POLICY=serialPolicy, typename ... ARGS >
  void BatchUpdateKernel( arrayView1d< real64 const > const & pressure,
                          arrayView1d< real64 const > const & temperature,
                          arrayView2d< real64 const > const & composition,
                          ARGS && ... args );

  /**
   * @brief Function called internally to resize member arrays
   * @param size primary dimension (e.g. number of cells)
   * @param numPts secondary dimension (e.g. number of gauss points per cell)
   */
  void ResizeFields( localIndex const size, localIndex const numPts );

  // flag indicating whether input/output component fractions are treated as mass fractions
  bool m_useMass;

  // general fluid composition information
  string_array m_componentNames;
  array1d< real64 > m_componentMolarWeight;

  string_array m_phaseNames;

  array3d< real64 > m_phaseFraction;
  array3d< real64 > m_dPhaseFraction_dPressure;
  array3d< real64 > m_dPhaseFraction_dTemperature;
  array4d< real64 > m_dPhaseFraction_dGlobalCompFraction;

  array3d< real64 > m_phaseDensity;
  array3d< real64 > m_dPhaseDensity_dPressure;
  array3d< real64 > m_dPhaseDensity_dTemperature;
  array4d< real64 > m_dPhaseDensity_dGlobalCompFraction;

  array3d< real64 > m_phaseViscosity;
  array3d< real64 > m_dPhaseViscosity_dPressure;
  array3d< real64 > m_dPhaseViscosity_dTemperature;
  array4d< real64 > m_dPhaseViscosity_dGlobalCompFraction;

  array4d< real64 > m_phaseCompFraction;
  array4d< real64 > m_dPhaseCompFraction_dPressure;
  array4d< real64 > m_dPhaseCompFraction_dTemperature;
  array5d< real64 > m_dPhaseCompFraction_dGlobalCompFraction;

  array2d< real64 > m_totalDensity;
  array2d< real64 > m_dTotalDensity_dPressure;
  array2d< real64 > m_dTotalDensity_dTemperature;
  array3d< real64 > m_dTotalDensity_dGlobalCompFraction;

};

template< typename LEAFCLASS, typename POLICY, typename ... ARGS >
void MultiFluidBase::BatchUpdateKernel( arrayView1d< real64 const > const & pressure,
                                        arrayView1d< real64 const > const & temperature,
                                        arrayView2d< real64 const > const & composition,
                                        ARGS && ... args )
{
  localIndex const numElem = m_phaseDensity.size( 0 );
  localIndex const numQ    = m_phaseDensity.size( 1 );

  localIndex const NC = numFluidComponents();
  localIndex const NP = numFluidPhases();
  bool const useMass = m_useMass;

  arrayView1d< string const > const & phaseNames = m_phaseNames;
  arrayView1d< real64 const > const & componentMolarWeight = m_componentMolarWeight;

  arrayView3d< real64 > const & phaseFraction = m_phaseFraction;
  arrayView3d< real64 > const & dPhaseFraction_dPressure = m_dPhaseFraction_dPressure;
  arrayView3d< real64 > const & dPhaseFraction_dTemperature = m_dPhaseFraction_dTemperature;
  arrayView4d< real64 > const & dPhaseFraction_dGlobalCompFraction = m_dPhaseFraction_dGlobalCompFraction;

  arrayView3d< real64 > const & phaseDensity = m_phaseDensity;
  arrayView3d< real64 > const & dPhaseDensity_dPressure = m_dPhaseDensity_dPressure;
  arrayView3d< real64 > const & dPhaseDensity_dTemperature = m_dPhaseDensity_dTemperature;
  arrayView4d< real64 > const & dPhaseDensity_dGlobalCompFraction = m_dPhaseDensity_dGlobalCompFraction;

  arrayView3d< real64 > const & phaseViscosity = m_phaseViscosity;
  arrayView3d< real64 > const & dPhaseViscosity_dPressure = m_dPhaseViscosity_dPressure;
  arrayView3d< real64 > const & dPhaseViscosity_dTemperature = m_dPhaseViscosity_dTemperature;
  arrayView4d< real64 > const & dPhaseViscosity_dGlobalCompFraction = m_dPhaseViscosity_dGlobalCompFraction;

  arrayView4d< real64 > const & phaseCompFraction = m_phaseCompFraction;
  arrayView4d< real64 > const & dPhaseCompFraction_dPressure = m_dPhaseCompFraction_dPressure;
  arrayView4d< real64 > const & dPhaseCompFraction_dTemperature = m_dPhaseCompFraction_dTemperature;
  arrayView5d< real64 > const & dPhaseCompFraction_dGlobalCompFraction = m_dPhaseCompFraction_dGlobalCompFraction;

  arrayView2d< real64 > const & totalDensity = m_totalDensity;
  arrayView2d< real64 > const & dTotalDensity_dPressure = m_dTotalDensity_dPressure;
  arrayView2d< real64 > const & dTotalDensity_dTemperature = m_dTotalDensity_dTemperature;
  arrayView3d< real64 > const & dTotalDensity_dGlobalCompFraction = m_dTotalDensity_dGlobalCompFraction;

  forAll< POLICY >( numElem, [=] ( localIndex const k )
  {
    for( localIndex q = 0; q < numQ; ++q )
    {
      LEAFCLASS::Compute( NC, NP, useMass,
                          phaseNames,
                          componentMolarWeight,
                          pressure[k],
                          temperature[k],
                          composition[k],
                          phaseFraction[k][q],
                          dPhaseFraction_dPressure[k][q],
                          dPhaseFraction_dTemperature[k][q],
                          dPhaseFraction_dGlobalCompFraction[k][q],
                          phaseDensity[k][q],
                          dPhaseDensity_dPressure[k][q],
                          dPhaseDensity_dTemperature[k][q],
                          dPhaseDensity_dGlobalCompFraction[k][q],
                          phaseViscosity[k][q],
                          dPhaseViscosity_dPressure[k][q],
                          dPhaseViscosity_dTemperature[k][q],
                          dPhaseViscosity_dGlobalCompFraction[k][q],
                          phaseCompFraction[k][q],
                          dPhaseCompFraction_dPressure[k][q],
                          dPhaseCompFraction_dTemperature[k][q],
                          dPhaseCompFraction_dGlobalCompFraction[k][q],
                          totalDensity[k][q],
                          dTotalDensity_dPressure[k][q],
                          dTotalDensity_dTemperature[k][q],
                          dTotalDensity_dGlobalCompFraction[k][q],
                          args ... );
    }
  } );
}

} //namespace constitutive

} //namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDBASE_HPP
