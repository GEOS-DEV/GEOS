/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
  * @file MultiFluidBase.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_MULTIFLUIDBASE_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_MULTIFLUIDBASE_HPP

#include "constitutive/ConstitutiveBase.hpp"

namespace geosx
{

namespace constitutive
{

namespace detail
{

template<typename T, int DIM>
struct array_slice_helper
{
  using type = array_slice<T, DIM>;
};

// an array slice of DIM=0 decays to a reference to scalar
template<typename T>
struct array_slice_helper<T, 0>
{
  using type = T &;
};

// an array1 slice of DIM=1 uses specialization (possibly a raw pointer)
template<typename T>
struct array_slice_helper<T, 1>
{
  using type = arraySlice1d<T>;
};

}

template<int DIM>
using real_array_slice = typename detail::array_slice_helper<real64, DIM>::type;

template<int DIM>
using real_array_const_slice = typename detail::array_slice_helper<real64 const, DIM>::type;

// helper struct to represent a var and its derivatives
template<int DIM>
struct CompositionalVarContainer
{
  real_array_slice<DIM>   value; // variable value
  real_array_slice<DIM>   dPres; // derivative w.r.t. pressure
  real_array_slice<DIM>   dTemp; // derivative w.r.t. temperature
  real_array_slice<DIM+1> dComp; // derivative w.r.t. composition
};

template<int DIM>
struct CompositionalVarConstContainer
{
  real_array_const_slice<DIM>   value; // variable value
  real_array_const_slice<DIM>   dPres; // derivative w.r.t. pressure
  real_array_const_slice<DIM>   dTemp; // derivative w.r.t. temperature
  real_array_const_slice<DIM+1> dComp; // derivative w.r.t. composition
};

class MultiFluidBase : public ConstitutiveBase
{
public:

  MultiFluidBase( std::string const & name, ManagedGroup * const parent );

  virtual ~MultiFluidBase() override;

  virtual void AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
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
                            arraySlice1d<real64 const> const & composition,
                            localIndex const k,
                            localIndex const q ) = 0;

  /**
   * @brief Perform a batch constitutive update (all points).
   * @param[in] pressure array containing target pressure values
   * @param[in] temperature array containing target temperature values
   * @param[in] composition 2D array containing target fluid composition values
   */
  virtual void BatchUpdate( arrayView1d<real64 const> const & pressure,
                            arrayView1d<real64 const> const & temperature,
                            arrayView2d<real64 const> const & composition ) = 0;

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
                        arraySlice1d<real64 const> const & composition,
                        arraySlice1d<real64> const & phaseFraction,
                        arraySlice1d<real64> const & dPhaseFraction_dPressure,
                        arraySlice1d<real64> const & dPhaseFraction_dTemperature,
                        arraySlice2d<real64> const & dPhaseFraction_dGlobalCompFraction,
                        arraySlice1d<real64> const & phaseDensity,
                        arraySlice1d<real64> const & dPhaseDensity_dPressure,
                        arraySlice1d<real64> const & dPhaseDensity_dTemperature,
                        arraySlice2d<real64> const & dPhaseDensity_dGlobalCompFraction,
                        arraySlice1d<real64> const & phaseViscosity,
                        arraySlice1d<real64> const & dPhaseViscosity_dPressure,
                        arraySlice1d<real64> const & dPhaseViscosity_dTemperature,
                        arraySlice2d<real64> const & dPhaseViscosity_dGlobalCompFraction,
                        arraySlice2d<real64> const & phaseCompFraction,
                        arraySlice2d<real64> const & dPhaseCompFraction_dPressure,
                        arraySlice2d<real64> const & dPhaseCompFraction_dTemperature,
                        arraySlice3d<real64> const & dPhaseCompFraction_dGlobalCompFraction,
                        real64 & totalDensity,
                        real64 & dTotalDensity_dPressure,
                        real64 & dTotalDensity_dTemperature,
                        arraySlice1d<real64> const & dTotalDensity_dGlobalCompFraction ) const = 0;

  /**
   * @return number of fluid components (species) in the model
   */
  localIndex numFluidComponents() const;

  /**
   * @param ic component index
   * @return name of ic-th fluid component
   */
  string const & componentName( localIndex ic ) const;

  /**
   * @return number of fluid phases in the model
   */
  localIndex numFluidPhases() const;

  /**
   * @param ip phase index
   * @return name of ip-th fluid phase
   */
  string const & phaseName( localIndex ip ) const;

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

    using ViewKey = dataRepository::ViewKey;

    ViewKey componentNames       = { componentNamesString };
    ViewKey componentMolarWeight = { componentMolarWeightString };

    ViewKey phaseNames     = { phaseNamesString };

    ViewKey phaseFraction                            = { phaseFractionString };                            // xi_p
    ViewKey dPhaseFraction_dPressure                 = { dPhaseFraction_dPressureString };                 // dXi_p/dP
    ViewKey dPhaseFraction_dTemperature              = { dPhaseFraction_dTemperatureString };              // dXi_p/dT
    ViewKey dPhaseFraction_dGlobalCompFraction       = { dPhaseFraction_dGlobalCompFractionString };       // dXi_p/dz

    ViewKey phaseDensity                             = { phaseDensityString };                             // rho_p
    ViewKey dPhaseDensity_dPressure                  = { dPhaseDensity_dPressureString };                  // dRho_p/dP
    ViewKey dPhaseDensity_dTemperature               = { dPhaseDensity_dTemperatureString };               // dRho_p/dT
    ViewKey dPhaseDensity_dGlobalCompFraction        = { dPhaseDensity_dGlobalCompFractionString };        // dRho_p/dz

    ViewKey phaseViscosity                           = { phaseViscosityString };                           // rho_p
    ViewKey dPhaseViscosity_dPressure                = { dPhaseViscosity_dPressureString };                // dRho_p/dP
    ViewKey dPhaseViscosity_dTemperature             = { dPhaseViscosity_dTemperatureString };             // dRho_p/dT
    ViewKey dPhaseViscosity_dGlobalCompFraction      = { dPhaseViscosity_dGlobalCompFractionString };      // dRho_p/dz

    ViewKey phaseCompFraction                        = { phaseCompFractionString };                        // x_cp
    ViewKey dPhaseCompFraction_dPressure             = { dPhaseCompFraction_dPressureString };             // dx_cp/dP
    ViewKey dPhaseCompFraction_dTemperature          = { dPhaseCompFraction_dTemperatureString };          // dx_cp/dT
    ViewKey dPhaseCompFraction_dGlobalCompFraction   = { dPhaseCompFraction_dGlobalCompFractionString };   // dx_cp/dz

    ViewKey totalDensity                             = { totalDensityString };                             // rho_t
    ViewKey dTotalDensity_dPressure                  = { dTotalDensity_dPressureString };                  // dRho_t/dP
    ViewKey dTotalDensity_dTemperature               = { dTotalDensity_dTemperatureString };               // dRho_t/dT
    ViewKey dTotalDensity_dGlobalCompFraction        = { dTotalDensity_dGlobalCompFractionString };        // dRho_t/dz

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
  template< typename LEAFCLASS, typename POLICY=elemPolicy, typename ... ARGS >
  void BatchUpdateKernel( arrayView1d<real64 const> const & pressure,
                          arrayView1d<real64 const> const & temperature,
                          arrayView2d<real64 const> const & composition,
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
  string_array    m_componentNames;
  array1d<real64> m_componentMolarWeight;

  string_array    m_phaseNames;

  array3d<real64> m_phaseFraction;
  array3d<real64> m_dPhaseFraction_dPressure;
  array3d<real64> m_dPhaseFraction_dTemperature;
  array4d<real64> m_dPhaseFraction_dGlobalCompFraction;

  array3d<real64> m_phaseDensity;
  array3d<real64> m_dPhaseDensity_dPressure;
  array3d<real64> m_dPhaseDensity_dTemperature;
  array4d<real64> m_dPhaseDensity_dGlobalCompFraction;

  array3d<real64> m_phaseViscosity;
  array3d<real64> m_dPhaseViscosity_dPressure;
  array3d<real64> m_dPhaseViscosity_dTemperature;
  array4d<real64> m_dPhaseViscosity_dGlobalCompFraction;

  array4d<real64> m_phaseCompFraction;
  array4d<real64> m_dPhaseCompFraction_dPressure;
  array4d<real64> m_dPhaseCompFraction_dTemperature;
  array5d<real64> m_dPhaseCompFraction_dGlobalCompFraction;

  array2d<real64> m_totalDensity;
  array2d<real64> m_dTotalDensity_dPressure;
  array2d<real64> m_dTotalDensity_dTemperature;
  array3d<real64> m_dTotalDensity_dGlobalCompFraction;

};

template<typename LEAFCLASS, typename POLICY, typename ... ARGS>
void MultiFluidBase::BatchUpdateKernel( arrayView1d<real64 const> const & pressure,
                                        arrayView1d<real64 const> const & temperature,
                                        arrayView2d<real64 const> const & composition,
                                        ARGS && ... args )
{
  localIndex const numElem = m_phaseDensity.size(0);
  localIndex const numQ    = m_phaseDensity.size(1);

  localIndex const NC = numFluidComponents();
  localIndex const NP = numFluidPhases();
  bool const useMass = m_useMass;

  arrayView1d<string const> const & phaseNames = m_phaseNames;
  arrayView1d<real64 const> const & componentMolarWeight = m_componentMolarWeight;

  arrayView3d<real64> const & phaseFraction = m_phaseFraction;
  arrayView3d<real64> const & dPhaseFraction_dPressure = m_dPhaseFraction_dPressure;
  arrayView3d<real64> const & dPhaseFraction_dTemperature = m_dPhaseFraction_dTemperature;
  arrayView4d<real64> const & dPhaseFraction_dGlobalCompFraction = m_dPhaseFraction_dGlobalCompFraction;

  arrayView3d<real64> const & phaseDensity = m_phaseDensity;
  arrayView3d<real64> const & dPhaseDensity_dPressure = m_dPhaseDensity_dPressure;
  arrayView3d<real64> const & dPhaseDensity_dTemperature = m_dPhaseDensity_dTemperature;
  arrayView4d<real64> const & dPhaseDensity_dGlobalCompFraction = m_dPhaseDensity_dGlobalCompFraction;

  arrayView3d<real64> const & phaseViscosity = m_phaseViscosity;
  arrayView3d<real64> const & dPhaseViscosity_dPressure = m_dPhaseViscosity_dPressure;
  arrayView3d<real64> const & dPhaseViscosity_dTemperature = m_dPhaseViscosity_dTemperature;
  arrayView4d<real64> const & dPhaseViscosity_dGlobalCompFraction = m_dPhaseViscosity_dGlobalCompFraction;

  arrayView4d<real64> const & phaseCompFraction = m_phaseCompFraction;
  arrayView4d<real64> const & dPhaseCompFraction_dPressure = m_dPhaseCompFraction_dPressure;
  arrayView4d<real64> const & dPhaseCompFraction_dTemperature = m_dPhaseCompFraction_dTemperature;
  arrayView5d<real64> const & dPhaseCompFraction_dGlobalCompFraction = m_dPhaseCompFraction_dGlobalCompFraction;

  arrayView2d<real64> const & totalDensity = m_totalDensity;
  arrayView2d<real64> const & dTotalDensity_dPressure = m_dTotalDensity_dPressure;
  arrayView2d<real64> const & dTotalDensity_dTemperature = m_dTotalDensity_dTemperature;
  arrayView3d<real64> const & dTotalDensity_dGlobalCompFraction = m_dTotalDensity_dGlobalCompFraction;

  forall_in_range<POLICY>( 0, numElem, GEOSX_LAMBDA ( localIndex const k )
  {
    for (localIndex q = 0; q < numQ; ++q)
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
                          args... );
    }
  });
}

} //namespace constitutive

} //namespace geosx

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_MULTIFLUIDBASE_HPP
