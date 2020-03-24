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
 * @file CompressibleSinglePhaseFluid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_COMPRESSIBLESINGLEPHASEFLUID_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_COMPRESSIBLESINGLEPHASEFLUID_HPP_

#include "constitutive/fluid/SingleFluidBase.hpp"

#include "constitutive/ExponentialRelation.hpp"

namespace geosx
{

namespace constitutive
{

class CompressibleSinglePhaseFluid : public SingleFluidBase
{
public:

  CompressibleSinglePhaseFluid( std::string const & name, Group * const parent );

  virtual ~CompressibleSinglePhaseFluid() override;

  // *** ConstitutiveBase interface

  virtual void DeliverClone( string const & name,
                             Group * const parent,
                             std::unique_ptr< ConstitutiveBase > & clone ) const override;

  static std::string CatalogName() { return "CompressibleSinglePhaseFluid"; }

  virtual string GetCatalogName() override { return CatalogName(); }

  virtual void AllocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  // *** SingleFluidBase interface

  virtual void PointUpdate( real64 const & pressure, localIndex const k, localIndex const q ) override;

  virtual void BatchUpdate( arrayView1d< real64 const > const & pressure ) override;

  virtual void Compute( real64 const & pressure,
                        real64 & density,
                        real64 & dDensity_dPressure,
                        real64 & viscosity,
                        real64 & dViscosity_dPressure ) const override;

  // *** Compute kernels

  /**
   * @brief Compute kernel for the partial constitutive update (single property)
   * @tparam EAT exponential approximation type
   * @param[in]  pressure
   * @param[out] property
   * @param[out] dProperty_dPressure
   * @param[in]  propertyRelation property exponential relation
   */
  template< ExponentApproximationType EAT >
  inline static void Compute( real64 const & pressure,
                              real64 & property,
                              real64 & dProperty_dPressure,
                              ExponentialRelation< real64, EAT > const & propertyRelation );

  /**
   * @brief Compute kernel for the full constitutive update
   * @tparam DENS_EAT density exponential approximation type
   * @tparam VISC_EAT viscosity exponential appeoximation type
   * @param[in]  pressure target pressure
   * @param[out] density fluid density
   * @param[out] dDensity_dPressure fluid density derivative w.r.t. pressure
   * @param[out] viscosity fluid viscosity
   * @param[out] dViscosity_dPressure fluid viscosity derivative w.r.t. pressure
   * @param[in]  densityRelation density exponential relation
   * @param[in]  viscosityRelation viscosity exponential relation
   */
  template< ExponentApproximationType DENS_EAT, ExponentApproximationType VISC_EAT >
  inline static void Compute( real64 const & pressure,
                              real64 & density,
                              real64 & dDensity_dPressure,
                              real64 & viscosity,
                              real64 & dViscosity_dPressure,
                              ExponentialRelation< real64, DENS_EAT > const & densityRelation,
                              ExponentialRelation< real64, VISC_EAT > const & viscosityRelation );

  // *** Data repository keys

  struct viewKeyStruct : public SingleFluidBase::viewKeyStruct
  {
    static constexpr auto compressibilityString    = "compressibility";
    static constexpr auto viscosibilityString      = "viscosibility";
    static constexpr auto referencePressureString  = "referencePressure";
    static constexpr auto referenceDensityString   = "referenceDensity";
    static constexpr auto referenceViscosityString = "referenceViscosity";
    static constexpr auto densityModelString       = "densityModel";
    static constexpr auto viscosityModelString     = "viscosityModel";

    dataRepository::ViewKey compressibility    = { compressibilityString    };
    dataRepository::ViewKey viscosibility      = { viscosibilityString      };
    dataRepository::ViewKey referencePressure  = { referencePressureString  };
    dataRepository::ViewKey referenceDensity   = { referenceDensityString   };
    dataRepository::ViewKey referenceViscosity = { referenceViscosityString };
    dataRepository::ViewKey densityModel       = { densityModelString       };
    dataRepository::ViewKey viscosityModel     = { viscosityModelString     };

  } viewKeysCompressibleSinglePhaseFluid;

protected:
  virtual void PostProcessInput() override;

private:

  /// scalar fluid bulk modulus parameter
  real64 m_compressibility;

  /// scalar fluid viscosity exponential coefficient
  real64 m_viscosibility;

  /// reference pressure parameter
  real64 m_referencePressure;

  /// reference density parameter
  real64 m_referenceDensity;

  /// reference viscosity parameter
  real64 m_referenceViscosity;

  /// input string for type of density model (linear, quadratic, exponential)
  string m_densityModelString;

  /// input string for type of viscosity model (linear, quadratic, exponential)
  string m_viscosityModelString;

  /// type of density model (linear, quadratic, exponential)
  ExponentApproximationType m_densityModelType;

  /// type of viscosity model (linear, quadratic, exponential)
  ExponentApproximationType m_viscosityModelType;
};


template< ExponentApproximationType EAT >
inline void CompressibleSinglePhaseFluid::Compute( real64 const & pressure,
                                                   real64 & property, real64 & dProperty_dPressure,
                                                   ExponentialRelation< real64, EAT > const & propertyRelation )
{
  propertyRelation.Compute( pressure, property, dProperty_dPressure );
}

template< ExponentApproximationType DENS_EAT, ExponentApproximationType VISC_EAT >
inline void CompressibleSinglePhaseFluid::Compute( real64 const & pressure,
                                                   real64 & density, real64 & dDensity_dPressure,
                                                   real64 & viscosity, real64 & dViscosity_dPressure,
                                                   ExponentialRelation< real64, DENS_EAT > const & densityRelation,
                                                   ExponentialRelation< real64, VISC_EAT > const & viscosityRelation )
{
  Compute( pressure, density, dDensity_dPressure, densityRelation );
  Compute( pressure, viscosity, dViscosity_dPressure, viscosityRelation );
}

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_FLUID_COMPRESSIBLESINGLEPHASEFLUID_HPP_ */
