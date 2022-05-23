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
 * @brief Update class for the model suitable for lambda capture.
 * @tparam DENS_PRES_EAT type of density exponent approximation for the pressure part
 * @tparam DENS_TEMP_EAT type of density exponent approximation for the temperature part
 * @tparam VISC_EAT type of viscosity exponent approximation
 * @tparam INTENERGY_EAT type of internal energy exponent approximation 
 */
template< ExponentApproximationType DENS_PRES_EAT, ExponentApproximationType DENS_TEMP_EAT, 
          ExponentApproximationType VISC_EAT, ExponentApproximationType INTENERGY_EAT >
class CompressibleSinglePhaseUpdate final : public SingleFluidBaseUpdate
{
public:

  using DensPresRelationType  = ExponentialRelation< real64, DENS_PRES_EAT >;
  using DensTempRelationType  = ExponentialRelation< real64, DENS_TEMP_EAT >; 
  using ViscRelationType      = ExponentialRelation< real64, VISC_EAT >;
  using IntEnergyRelationType = ExponentialRelation< real64, INTENERGY_EAT >; 

  CompressibleSinglePhaseUpdate( DensPresRelationType const & densPresRelation, 
                                 DensTempRelationType const & densTempRelation,
                                 ViscRelationType const & viscRelation, 
                                 IntEnergyRelationType const & intEnergyRelation, 
                                 arrayView2d< real64 > const & density,
                                 arrayView2d< real64 > const & dDens_dPres,
                                 arrayView2d< real64 > const & dDens_dTemp, 
                                 arrayView2d< real64 > const & viscosity,
                                 arrayView2d< real64 > const & dVisc_dPres, 
                                 arrayView2d< real64 > const & dVisc_dTemp,
                                 arrayView2d< real64 > const & internalEnergy, 
                                 arrayView2d< real64 > const & dIntEnergy_dPres, 
                                 arrayView2d< real64 > const & dIntEnergy_dTemp,
                                 arrayView2d< real64 > const & enthalpy, 
                                 arrayView2d< real64 > const & dEnthalpy_dPres, 
                                 arrayView2d< real64 > const & dEnthalpy_dTemp,
                                 real64 const & refIntEnergy, 
                                 integer const isThermal )
    : SingleFluidBaseUpdate( density,
                             dDens_dPres,
                             dDens_dTemp, 
                             viscosity,
                             dVisc_dPres,
                             dVisc_dTemp, 
                             internalEnergy, 
                             dIntEnergy_dPres,
                             dIntEnergy_dTemp,
                             enthalpy, 
                             dEnthalpy_dPres, 
                             dEnthalpy_dTemp ),
    m_densPresRelation( densPresRelation ), 
    m_densTempRelation( densTempRelation ), 
    m_viscRelation( viscRelation ), 
    m_intEnergyRelation( intEnergyRelation ), 
    m_refIntEnergy( refIntEnergy ), 
    m_isThermal( isThermal )
  {}

  /// Default copy constructor
  CompressibleSinglePhaseUpdate( CompressibleSinglePhaseUpdate const & ) = default;

  /// Default move constructor
  CompressibleSinglePhaseUpdate( CompressibleSinglePhaseUpdate && ) = default;

  /// Deleted copy assignment operator
  CompressibleSinglePhaseUpdate & operator=( CompressibleSinglePhaseUpdate const & ) = delete;

  /// Deleted move assignment operator
  CompressibleSinglePhaseUpdate & operator=( CompressibleSinglePhaseUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void compute( real64 const pressure,
                        real64 & density,
                        real64 & viscosity) const override
  {
    m_densPresRelation.compute( pressure, density );
    m_viscRelation.compute( pressure, viscosity );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void compute( real64 const pressure,
                        real64 & density,
                        real64 & dDensity_dPressure,
                        real64 & viscosity,
                        real64 & dViscosity_dPressure ) const override
  {
    m_densPresRelation.compute( pressure, density, dDensity_dPressure );
    m_viscRelation.compute( pressure, viscosity, dViscosity_dPressure );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void compute( real64 const pressure,
                        real64 const temperature, 
                        real64 & density,
                        real64 & viscosity ) const override
  {
    m_viscRelation.compute( pressure, viscosity );

    if ( m_isThermal )
    {
      real64 density_pressurePart, density_temperaturePart; 

      m_densPresRelation.compute( pressure, density_pressurePart ); 
      m_densTempRelation.compute( temperature, density_temperaturePart ); 

      density = density_pressurePart * density_temperaturePart; 
    }
    else
    {
      m_densPresRelation.compute( pressure, density ); 
    }
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void compute( real64 const pressure,
                        real64 const temperature, 
                        real64 & density,
                        real64 & dDensity_dPressure,
                        real64 & dDensity_dTemperature, 
                        real64 & viscosity,
                        real64 & dViscosity_dPressure, 
                        real64 & dViscosity_dTemperature, 
                        real64 & internalEnergy, 
                        real64 & dInternalEnergy_dPressure, 
                        real64 & dInternalEnergy_dTemperature,
                        real64 & enthalpy, 
                        real64 & dEnthalpy_dPressure,
                        real64 & dEnthalpy_dTemperature ) const override
  {
    m_viscRelation.compute( pressure, viscosity, dViscosity_dPressure );
    dViscosity_dTemperature = 0.0; 

    if ( m_isThermal )
    {
      real64 density_pressurePart, density_temperaturePart; 
      real64 density_pressurePart_deriv, density_temperaturePart_deriv; 

      m_densPresRelation.compute( pressure, density_pressurePart, density_pressurePart_deriv ); 
      m_densTempRelation.compute( temperature, density_temperaturePart, density_temperaturePart_deriv ); 

      density = density_pressurePart * density_temperaturePart; 

      dDensity_dPressure = density_pressurePart_deriv * density_temperaturePart; 
      dDensity_dTemperature = density_temperaturePart_deriv * density_pressurePart; 

      /// Compute the internal energy (only sensitive to temperature)
      m_intEnergyRelation.compute( temperature, internalEnergy, dInternalEnergy_dTemperature ); 
      dInternalEnergy_dPressure = 0.0; 

      enthalpy = internalEnergy - m_refIntEnergy; 
      dEnthalpy_dPressure = 0.0; 
      dEnthalpy_dTemperature = dInternalEnergy_dTemperature; 
    }
    else
    {
      m_densPresRelation.compute( pressure, density, dDensity_dPressure ); 
      dDensity_dTemperature = 0.0; 
    }

  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const pressure ) const override
  {
    compute( pressure,
             m_density[k][q],
             m_dDens_dPres[k][q],
             m_viscosity[k][q],
             m_dVisc_dPres[k][q] );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const pressure,
                       real64 const temperature ) const override
  {
    compute( pressure,
             temperature, 
             m_density[k][q],
             m_dDens_dPres[k][q],
             m_dDens_dTemp[k][q], 
             m_viscosity[k][q],
             m_dVisc_dPres[k][q],
             m_dVisc_dTemp[k][q], 
             m_internalEnergy[k][q],
             m_dIntEnergy_dPres[k][q],
             m_dIntEnergy_dTemp[k][q],
             m_enthalpy[k][q], 
             m_dEnthalpy_dPres[k][q],
             m_dEnthalpy_dTemp[k][q] );
  }

private:

  /// Relationship between the fluid density and pressure 
  DensPresRelationType m_densPresRelation; 

  /// Relationship between the fluid density and temperature  
  DensTempRelationType m_densTempRelation; 

  /// Relationship between the fluid viscosity and pressure 
  ViscRelationType m_viscRelation;

  /// Relationship between the fluid internal energy and temperature 
  IntEnergyRelationType m_intEnergyRelation; 

  /// Reference internal energy of the fluid
  real64 const m_refIntEnergy; 

  /// Flag to determine whether it is a nonisothermal fluid
  integer m_isThermal; 
};

class CompressibleSinglePhaseFluid : public SingleFluidBase
{
public:

  CompressibleSinglePhaseFluid( string const & name, Group * const parent );

  virtual ~CompressibleSinglePhaseFluid() override;

  static string catalogName() { return "CompressibleSinglePhaseFluid"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual bool isThermal() const override; 

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /// Type of kernel wrapper for in-kernel update (TODO: support multiple EAT, not just linear)
  using KernelWrapper = CompressibleSinglePhaseUpdate< ExponentApproximationType::Full, ExponentApproximationType::Full,
                                                       ExponentApproximationType::Full, ExponentApproximationType::Linear >;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  struct viewKeyStruct
  {
    static constexpr char const * isThermalString() { return "isThermal"; }
    static constexpr char const * defaultDensityString() { return "defaultDensity"; }
    static constexpr char const * defaultViscosityString() { return "defaultViscosity"; }
    static constexpr char const * compressibilityString() { return "compressibility"; }
    static constexpr char const * thermalExpansionCoeffString() { return "thermalExpansionCoeff"; }
    static constexpr char const * viscosibilityString() { return "viscosibility"; }
    static constexpr char const * volumetricHeatCapacityString() { return "volumetricHeatCapacity"; }
    static constexpr char const * referencePressureString() { return "referencePressure"; }
    static constexpr char const * referenceTemperatureString() { return "referenceTemperature"; }
    static constexpr char const * referenceDensityString() { return "referenceDensity"; }
    static constexpr char const * referenceViscosityString() { return "referenceViscosity"; }
    static constexpr char const * referenceInternalEnergyString() { return "referenceInternalEnergy"; }
    static constexpr char const * densityPressureModelTypeString() { return "densityPressureModelType"; }
    static constexpr char const * densityTemperatureModelTypeString() { return "densityTemperatureModelType"; }
    static constexpr char const * viscosityModelTypeString() { return "viscosityModelType"; }
    static constexpr char const * internalEnergyModelTypeString() { return "internalEnergyModelType"; }
  };

  real64 defaultDensity() const override final { return m_defaultDensity; }
  real64 defaultViscosity() const override final { return m_defaultViscosity; }

protected:

  virtual void postProcessInput() override;

private:

  /// flag to determine if it is a nonisothermal fluid
  integer m_isThermal; 

  /// default density value
  real64 m_defaultDensity;

  /// default viscosity value
  real64 m_defaultViscosity;

  /// scalar fluid bulk modulus parameter
  real64 m_compressibility;

  /// scalar fluid thermal expansion coefficient
  real64 m_thermalExpansionCoeff; 

  /// scalar fluid viscosity exponential coefficient
  real64 m_viscosibility;

  /// scalar fluid volumetric heat capacity coefficient 
  real64 m_volumetricHeatCapacity;

  /// reference pressure parameter
  real64 m_referencePressure;

  /// reference temperature parameter
  real64 m_referenceTemperature;

  /// reference density parameter
  real64 m_referenceDensity;

  /// reference viscosity parameter
  real64 m_referenceViscosity;

  /// reference internal energy parameter
  real64 m_referenceInternalEnergy; 

  /// type of density model in terms of pressure 
  ExponentApproximationType m_densityPressureModelType;

  /// type of density model in terms of temperature 
  ExponentApproximationType m_densityTemperatureModelType;

  /// type of viscosity model 
  ExponentApproximationType m_viscosityModelType;

  /// type of internal energy model 
  ExponentApproximationType m_internalEnergyModelType;
};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_FLUID_COMPRESSIBLESINGLEPHASEFLUID_HPP_ */
