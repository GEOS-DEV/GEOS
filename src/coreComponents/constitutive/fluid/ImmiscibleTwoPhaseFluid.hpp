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
  * @file ImmiscibleTwoPhaseFluid.hpp
  */

#ifndef GEOSX_CONSTITUTIVE_FLUID_IMMISCIBLETWOPHASEFLUID_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_IMMISCIBLETWOPHASEFLUID_HPP_

#include "constitutive/fluid/MultiFluidBase.hpp"

#include "constitutive/ExponentialRelation.hpp"

namespace geosx
{

  
namespace constitutive
{

class ImmiscibleTwoPhaseFluid : public MultiFluidBase
{
public:

  static constexpr localIndex NUM_PHASES     = 2;
  static constexpr localIndex NUM_COMPONENTS = 2;
  
  ImmiscibleTwoPhaseFluid( std::string const & name, Group * const parent );

  virtual ~ImmiscibleTwoPhaseFluid() override;

  virtual void DeliverClone( string const & name,
                             Group * const parent,
                             std::unique_ptr<ConstitutiveBase> & clone ) const override;

  static std::string CatalogName() { return "ImmiscibleTwoPhaseFluid"; }

  virtual string GetCatalogName() override { return CatalogName(); }

  virtual void PointUpdate( real64 const & pressure,
                            real64 const & temperature,
                            arraySlice1d<real64 const> const & composition,
                            localIndex const k,
                            localIndex const q ) override;

  virtual void BatchUpdate( arrayView1d<real64 const> const & pressure,
                            arrayView1d<real64 const> const & temperature,
                            arrayView2d<real64 const> const & composition ) override;

  static void Compute( localIndex const NC,
                       localIndex const NP,
                       bool const useMass,
                       arrayView1d<string const> const & phaseNames,
                       arrayView1d<real64 const> const & componentMolarWeight,
                       real64 const & pressure,
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
                       arraySlice1d<real64> const & dTotalDensity_dGlobalCompFraction,
                       real64 const & referencePressure,
                       arrayView1d<real64 const> const referencePhaseDensity,
                       arrayView1d<real64 const> const referencePhaseViscosity,
                       arrayView1d<real64 const> const phaseCompressibility,
                       arrayView1d<real64 const> const phaseViscosibility );

  virtual void Compute( real64 const & pressure, real64 const & temperature,
                        arraySlice1d<double const> const & composition,
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
                        real64 & totalDensity, real64 & dTotalDensity_dPressure,
                        real64 & dTotalDensity_dTemperature,
                        arraySlice1d<real64> const & dTotalDensity_dGlobalCompFraction) const override;

  
  /**
   * @brief Compute kernel for the partial constitutive update (single property)
   * @tparam EAT exponential approximation type
   * @param[in]  pressure
   * @param[out] property
   * @param[out] dProperty_dPressure
   * @param[in]  propertyRelation property exponential relation
   */
  template<ExponentApproximationType EAT>
  inline static void ComputePhaseProperty( real64 const & pressure,
                                           real64 & property,
                                           real64 & dProperty_dPressure,
                                           ExponentialRelation<real64, EAT> const & propertyRelation );


  struct viewKeyStruct : MultiFluidBase::viewKeyStruct
  {
    static constexpr auto phaseCompressibilityString    = "phaseCompressibility";
    static constexpr auto phaseViscosibilityString      = "phaseViscosibility";
    
    static constexpr auto referencePhaseDensityString   = "referencePhaseDensity";
    static constexpr auto referencePhaseViscosityString = "referencePhaseViscosity";
    static constexpr auto referencePressureString       = "referencePressure";
    
    dataRepository::ViewKey phaseCompressibility    = { phaseCompressibilityString  };
    dataRepository::ViewKey phaseViscosibility      = { phaseViscosibilityString    };
    
    dataRepository::ViewKey referencePhaseDensity   = { referencePhaseDensityString   };
    dataRepository::ViewKey referencePhaseViscosity = { referencePhaseViscosityString };
    dataRepository::ViewKey referencePressure       = { referencePressureString       };
    
  } viewKeysImmiscibleTwoPhaseFluid;

  
protected:
  
  virtual void PostProcessInput() override;

  virtual void InitializePostSubGroups( Group * const group ) override;

private:

  /// default phase densities
  array1d<real64> m_defaultPhaseDensity;
  /// default phase viscosities
  array1d<real64> m_defaultPhaseViscosity;

  /// phase compressibilities
  array1d<real64> m_phaseCompressibility;
  /// phase viscosity exponential coefficient
  array1d<real64> m_phaseViscosibility;

  /// reference density parameter
  array1d<real64> m_referencePhaseDensity;
  /// reference viscosity parameter
  array1d<real64> m_referencePhaseViscosity;

  /// reference pressure parameter
  real64 m_referencePressure;

};

template<ExponentApproximationType EAT>
inline void ImmiscibleTwoPhaseFluid::ComputePhaseProperty( real64 const & pressure,
                                                           real64 & phaseProperty, real64 & dPhaseProperty_dPressure,
                                                           ExponentialRelation<real64, EAT> const & propertyRelation )
{
  propertyRelation.Compute( pressure, phaseProperty, dPhaseProperty_dPressure );
}

  
} //namespace constitutive

} //namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_IMMISCIBLETWOPHASEFLUID_HPP_
