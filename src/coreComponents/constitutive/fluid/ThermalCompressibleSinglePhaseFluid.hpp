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
 * @file ThermalCompressibleSinglePhaseFluid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_THERMALCOMPRESSIBLESINGLEPHASEFLUID_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_THERMALCOMPRESSIBLESINGLEPHASEFLUID_HPP_

#include "constitutive/fluid/SingleFluidBase.hpp"
#include "constitutive/fluid/CompressibleSinglePhaseFluid.hpp"

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
class ThermalCompressibleSinglePhaseUpdate : public SingleFluidBaseUpdate
{
public:

  using DensPresRelationType  = ExponentialRelation< real64, DENS_PRES_EAT >;
  using DensTempRelationType  = ExponentialRelation< real64, DENS_TEMP_EAT >;
  using ViscRelationType      = ExponentialRelation< real64, VISC_EAT >;
  using IntEnergyRelationType = ExponentialRelation< real64, INTENERGY_EAT >;

  ThermalCompressibleSinglePhaseUpdate( DensPresRelationType const & densPresRelation,
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
                                        real64 const & refIntEnergy )
    : SingleFluidBaseUpdate( density,
                             dDens_dPres,
                             viscosity,
                             dVisc_dPres ),
    m_dDens_dTemp( dDens_dTemp ),
    m_dVisc_dTemp( dVisc_dTemp ),
    m_internalEnergy( internalEnergy ),
    m_dIntEnergy_dPres( dIntEnergy_dPres ),
    m_dIntEnergy_dTemp( dIntEnergy_dTemp ),
    m_enthalpy( enthalpy ),
    m_dEnthalpy_dPres( dEnthalpy_dPres ),
    m_dEnthalpy_dTemp( dEnthalpy_dTemp ),
    m_densPresRelation( densPresRelation ),
    m_densTempRelation( densTempRelation ),
    m_viscRelation( viscRelation ),
    m_intEnergyRelation( intEnergyRelation ),
    m_refIntEnergy( refIntEnergy )
  {}

  /// Default copy constructor
  ThermalCompressibleSinglePhaseUpdate( ThermalCompressibleSinglePhaseUpdate const & ) = default;

  /// Default move constructor
  ThermalCompressibleSinglePhaseUpdate( ThermalCompressibleSinglePhaseUpdate && ) = default;

  /// Deleted copy assignment operator
  ThermalCompressibleSinglePhaseUpdate & operator=( ThermalCompressibleSinglePhaseUpdate const & ) = delete;

  /// Deleted move assignment operator
  ThermalCompressibleSinglePhaseUpdate & operator=( ThermalCompressibleSinglePhaseUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void compute( real64 const pressure,
                        real64 & density,
                        real64 & viscosity ) const override
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

  /// Derivative of density w.r.t. temperature
  arrayView2d< real64 > m_dDens_dTemp;

  /// Derivative of viscosity w.r.t. temperature
  arrayView2d< real64 > m_dVisc_dTemp;

  /// Fluid internal energy
  arrayView2d< real64 > m_internalEnergy;

  /// Derivative of internal energy w.r.t. pressure
  arrayView2d< real64 > m_dIntEnergy_dPres;

  /// Derivative of internal energy w.r.t. temperature
  arrayView2d< real64 > m_dIntEnergy_dTemp;

  /// Fluid enthalpy
  arrayView2d< real64 > m_enthalpy;

  /// Derivative of enthalpy w.r.t. pressure
  arrayView2d< real64 > m_dEnthalpy_dPres;

  /// Derivative of enthalpy w.r.t. temperature
  arrayView2d< real64 > m_dEnthalpy_dTemp;

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

};

class ThermalCompressibleSinglePhaseFluid : public CompressibleSinglePhaseFluid
{
public:

  ThermalCompressibleSinglePhaseFluid( string const & name, Group * const parent );

  virtual ~ThermalCompressibleSinglePhaseFluid() override;

  static string catalogName() { return "ThermalCompressibleSinglePhaseFluid"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /// Type of kernel wrapper for in-kernel update (TODO: support multiple EAT, not just linear)
  using KernelWrapper = ThermalCompressibleSinglePhaseUpdate< ExponentApproximationType::Full, ExponentApproximationType::Full,
                                                              ExponentApproximationType::Linear, ExponentApproximationType::Linear >;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  struct viewKeyStruct : public CompressibleSinglePhaseFluid::viewKeyStruct
  {
    static constexpr char const * thermalExpansionCoeffString() { return "thermalExpansionCoeff"; }
    static constexpr char const * volumetricHeatCapacityString() { return "volumetricHeatCapacity"; }
    static constexpr char const * referenceTemperatureString() { return "referenceTemperature"; }
    static constexpr char const * referenceInternalEnergyString() { return "referenceInternalEnergy"; }
    static constexpr char const * densityPressureModelTypeString() { return "densityPressureModelType"; }
    static constexpr char const * densityTemperatureModelTypeString() { return "densityTemperatureModelType"; }
    static constexpr char const * internalEnergyModelTypeString() { return "internalEnergyModelType"; }
  };

protected:

  virtual void postProcessInput() override;

private:

  /// scalar fluid thermal expansion coefficient
  real64 m_thermalExpansionCoeff;

  /// scalar fluid volumetric heat capacity coefficient
  real64 m_volumetricHeatCapacity;

  /// reference temperature parameter
  real64 m_referenceTemperature;

  /// reference internal energy parameter
  real64 m_referenceInternalEnergy;

  /// type of density model in terms of pressure
  ExponentApproximationType m_densityPressureModelType;

  /// type of density model in terms of temperature
  ExponentApproximationType m_densityTemperatureModelType;

  /// type of internal energy model
  ExponentApproximationType m_internalEnergyModelType;
};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_FLUID_THERMALCOMPRESSIBLESINGLEPHASEFLUID_HPP_ */
