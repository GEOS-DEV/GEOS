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
  * @file CompressibleSinglePhaseFluid.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_COMPRESSIBLESINGLEPHASEFLUID_HPP_
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_COMPRESSIBLESINGLEPHASEFLUID_HPP_

#include "constitutive/Fluid/SingleFluidBase.hpp"

#include "constitutive/ExponentialRelation.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const compressibleSinglePhaseFluid = "CompressibleSinglePhaseFluid";
}
}

namespace constitutive
{

class CompressibleSinglePhaseFluid : public SingleFluidBase
{
public:

  CompressibleSinglePhaseFluid( std::string const & name, ManagedGroup * const parent );

  virtual ~CompressibleSinglePhaseFluid() override;

  // *** ConstitutiveBase interface

  std::unique_ptr<ConstitutiveBase> DeliverClone( string const & name,
                                                  ManagedGroup * const parent ) const override;

  static std::string CatalogName() { return dataRepository::keys::compressibleSinglePhaseFluid; }

  virtual string GetCatalogName() override { return CatalogName(); }

  virtual void AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  // *** SingleFluidBase interface

  virtual void PointUpdate( real64 const & pressure, localIndex const k, localIndex const q ) override;

  virtual void BatchUpdate( arrayView1d<real64 const> const & pressure ) override;

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
  template<ExponentApproximationType EAT>
  inline static void Compute( real64 const & pressure,
                              real64 & property,
                              real64 & dProperty_dPressure,
                              ExponentialRelation<real64, EAT> const & propertyRelation );

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
  template<ExponentApproximationType DENS_EAT, ExponentApproximationType VISC_EAT>
  inline static void Compute( real64 const & pressure,
                              real64 & density,
                              real64 & dDensity_dPressure,
                              real64 & viscosity,
                              real64 & dViscosity_dPressure,
                              ExponentialRelation<real64, DENS_EAT> const & densityRelation,
                              ExponentialRelation<real64, VISC_EAT> const & viscosityRelation );

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


template<ExponentApproximationType EAT>
inline void CompressibleSinglePhaseFluid::Compute( real64 const & pressure,
                                                   real64 & property, real64 & dProperty_dPressure,
                                                   ExponentialRelation<real64, EAT> const & propertyRelation )
{
  propertyRelation.Compute( pressure, property, dProperty_dPressure );
}

template<ExponentApproximationType DENS_EAT, ExponentApproximationType VISC_EAT>
inline void CompressibleSinglePhaseFluid::Compute( real64 const & pressure,
                                                   real64 & density, real64 & dDensity_dPressure,
                                                   real64 & viscosity, real64 & dViscosity_dPressure,
                                                   ExponentialRelation<real64, DENS_EAT> const & densityRelation,
                                                   ExponentialRelation<real64, VISC_EAT> const & viscosityRelation )
{
  Compute( pressure, density, dDensity_dPressure, densityRelation );
  Compute( pressure, viscosity, dViscosity_dPressure, viscosityRelation );
}

} /* namespace constitutive */

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_COMPRESSIBLESINGLEPHASEFLUID_HPP_ */
