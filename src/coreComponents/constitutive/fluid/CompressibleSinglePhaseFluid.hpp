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

/**
 * @brief Lightweight compute-only class for the model for device capture.
 * @tparam DENS_EAT type of density exponent approximation
 * @tparam VISC_EAT type of viscosity exponent approximation
 */
template< ExponentApproximationType DENS_EAT, ExponentApproximationType VISC_EAT >
class CompressibleSinglePhaseCompute
{
public:

  using DensRelationType = ExponentialRelation< real64, DENS_EAT >;
  using ViscRelationType = ExponentialRelation< real64, VISC_EAT >;

  CompressibleSinglePhaseCompute( DensRelationType const & densRelation,
                                  ViscRelationType const & viscRelation )
  : m_densRelation( densRelation ),
    m_viscRelation( viscRelation )
  {}

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void Compute( real64 const pres, real64 & dens, real64 & visc ) const
  {
    m_densRelation.Compute( pres, dens );
    m_viscRelation.Compute( pres, visc );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void Compute( real64 const pres, real64 & dens, real64 & dDens_dPres, real64 & visc, real64 & dVisc_dPres ) const
  {
    m_densRelation.Compute( pres, dens, dDens_dPres );
    m_viscRelation.Compute( pres, visc, dVisc_dPres );
  }

private:

  DensRelationType m_densRelation;
  ViscRelationType m_viscRelation;
};

/**
 * @brief Update class for the model for device capture.
 * @tparam DENS_EAT type of density exponent approximation
 * @tparam VISC_EAT type of viscosity exponent approximation
 */
template< ExponentApproximationType DENS_EAT, ExponentApproximationType VISC_EAT >
class CompressibleSinglePhaseUpdate
{
public:

  using ComputeType = CompressibleSinglePhaseCompute< DENS_EAT, VISC_EAT >;

  CompressibleSinglePhaseUpdate( typename ComputeType::DensRelationType const & densRelation,
                                 typename ComputeType::ViscRelationType const & viscRelation,
                                 arrayView2d< real64 > const & density,
                                 arrayView2d< real64 > const & dDens_dPres,
                                 arrayView2d< real64 > const & viscosity,
                                 arrayView2d< real64 > const & dVisc_dPres )
    : m_compute ( densRelation, viscRelation ),
    m_density( density ),
    m_dDens_dPres( dDens_dPres ),
    m_viscosity( viscosity ),
    m_dVisc_dPres( dVisc_dPres )
  {}

  /// Default copy constructor
  CompressibleSinglePhaseUpdate( CompressibleSinglePhaseUpdate const & ) = default;

  /// Default move constructor
  CompressibleSinglePhaseUpdate( CompressibleSinglePhaseUpdate && ) = default;

  /// Deleted default constructor
  CompressibleSinglePhaseUpdate() = delete;

  /// Deleted copy assignment operator
  CompressibleSinglePhaseUpdate & operator=( CompressibleSinglePhaseUpdate const & ) = delete;

  /// Deleted move assignment operator
  CompressibleSinglePhaseUpdate & operator=( CompressibleSinglePhaseUpdate && ) =  delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void Update( real64 const pres, localIndex const k, localIndex const q ) const
  {
    m_compute.Compute( pres, m_density[k][q], m_dDens_dPres[k][q], m_viscosity[k][q], m_dVisc_dPres[k][q] );
  }

private:

  ComputeType m_compute;

  arrayView2d< real64 > m_density;
  arrayView2d< real64 > m_dDens_dPres;

  arrayView2d< real64 > m_viscosity;
  arrayView2d< real64 > m_dVisc_dPres;
};

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

  virtual void PointUpdate( real64 const pressure, localIndex const k, localIndex const q ) override;

  virtual void BatchUpdate( arrayView1d< real64 const > const & pressure,
                            arrayView1d< real64 const > const & deltaPressure ) override;

  virtual void Compute( real64 const pressure,
                        real64 const deltaPressure,
                        real64 & density,
                        real64 & dDensity_dPressure,
                        real64 & viscosity,
                        real64 & dViscosity_dPressure ) const override;

  // *** Compute kernels

  /// Type of kernel wrapper for in-kernel compute (TODO: support multiple EAT, not just linear)
  using ComputeWrapper = CompressibleSinglePhaseCompute< ExponentApproximationType::Linear, ExponentApproximationType::Linear >;

  /**
   * @brief Create a compute-only kernel wrapper.
   * @return the wrapper
   */
  ComputeWrapper createComputeWrapper() const
  {
    ComputeWrapper::DensRelationType densRelation( m_referencePressure, m_referenceDensity, m_compressibility );
    ComputeWrapper::ViscRelationType viscRelation( m_referencePressure, m_referenceViscosity, m_viscosibility );
    return ComputeWrapper( densRelation, viscRelation );
  }

  /// Type of kernel wrapper for in-kernel update (TODO: support multiple EAT, not just linear)
  using UpdateWrapper = CompressibleSinglePhaseUpdate< ExponentApproximationType::Linear, ExponentApproximationType::Linear >;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  UpdateWrapper createUpdateWrapper()
  {
    UpdateWrapper::ComputeType::DensRelationType densRelation( m_referencePressure, m_referenceDensity, m_compressibility );
    UpdateWrapper::ComputeType::ViscRelationType viscRelation( m_referencePressure, m_referenceViscosity, m_viscosibility );
    return UpdateWrapper( densRelation, viscRelation,
                          m_density.toView(), m_dDensity_dPressure.toView(),
                          m_viscosity.toView(), m_dViscosity_dPressure.toView() );
  }

  /**
   * @brief Compute kernel for the partial constitutive update (single property)
   * @tparam EAT exponential approximation type
   * @param[in]  pressure
   * @param[out] property
   * @param[out] dProperty_dPressure
   * @param[in]  propertyRelation property exponential relation
   */
  template< ExponentApproximationType EAT >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void Compute( real64 const pressure,
                       real64 const deltaPressure,
                       real64 & property,
                       real64 & dProperty_dPressure,
                       ExponentialRelation< real64, EAT > const & propertyRelation )
  {
    propertyRelation.Compute( pressure + deltaPressure, property, dProperty_dPressure );
  }

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
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void Compute( real64 const pressure,
                       real64 const deltaPressure,
                       real64 & density,
                       real64 & dDensity_dPressure,
                       real64 & viscosity,
                       real64 & dViscosity_dPressure,
                       ExponentialRelation< real64, DENS_EAT > const & densityRelation,
                       ExponentialRelation< real64, VISC_EAT > const & viscosityRelation )
  {
    Compute( pressure + deltaPressure, density, dDensity_dPressure, densityRelation );
    Compute( pressure + deltaPressure, viscosity, dViscosity_dPressure, viscosityRelation );
  }

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

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_FLUID_COMPRESSIBLESINGLEPHASEFLUID_HPP_ */
