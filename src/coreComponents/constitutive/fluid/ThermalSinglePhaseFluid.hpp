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
 * @file ThermalSinglePhaseFluid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_THERMALSINGLEPHASEFLUID_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_THERMALSINGLEPHASEFLUID_HPP_

#include "constitutive/fluid/SingleFluidBase.hpp"

namespace geosx
{

namespace constitutive
{

/**
 * @brief Update class for the model suitable for lambda capture.
 */
class ThermalSinglePhaseUpdate final : public SingleFluidBaseUpdate
{
public:

  ThermalSinglePhaseUpdate( real64 const & compressibility, 
                            real64 const & thermalExpansionCoeff, 
                            real64 const & viscosibility, 
                            real64 const & volumetricHeatCapacity, 
                            real64 const & referencePressure, 
                            real64 const & referenceTemperature, 
                            real64 const & referenceDensity, 
                            real64 const & referenceViscosity, 
                            real64 const & referenceInternalEnergy, 
                            arrayView2d< real64 > const & density,
                            arrayView2d< real64 > const & dDens_dPres,
                            arrayView2d< real64 > const & dDens_dTemp, 
                            arrayView2d< real64 > const & viscosity,
                            arrayView2d< real64 > const & dVisc_dPres, 
                            arrayView2d< real64 > const & dVisc_dTemp,
                            arrayView2d< real64 > const & internalEnergy, 
                            arrayView2d< real64 > const & dIntEnergy_dPres, 
                            arrayView2d< real64 > const & dIntEnergy_dTemp )
    : SingleFluidBaseUpdate( density,
                             dDens_dPres,
                             dDens_dTemp, 
                             viscosity,
                             dVisc_dPres,
                             dVisc_dTemp, 
                             internalEnergy, 
                             dIntEnergy_dPres,
                             dIntEnergy_dTemp ),
    m_compressibility( compressibility ),
    m_thermalExpansionCoeff( thermalExpansionCoeff ),
    m_viscosibility( viscosibility ),  
    m_volumetricHeatCapacity( volumetricHeatCapacity ), 
    m_referencePressure( referencePressure ), 
    m_referenceTemperature( referenceTemperature ),
    m_referenceDensity( referenceDensity ), 
    m_referenceViscosity( referenceViscosity ),  
    m_referenceInternalEnergy( referenceInternalEnergy )
  {}

  /// Default copy constructor
  ThermalSinglePhaseUpdate( ThermalSinglePhaseUpdate const & ) = default;

  /// Default move constructor
  ThermalSinglePhaseUpdate( ThermalSinglePhaseUpdate && ) = default;

  /// Deleted copy assignment operator
  ThermalSinglePhaseUpdate & operator=( ThermalSinglePhaseUpdate const & ) = delete;

  /// Deleted move assignment operator
  ThermalSinglePhaseUpdate & operator=( ThermalSinglePhaseUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void compute( real64 const pressure,
                        real64 & density,
                        real64 & viscosity ) const override
  {
    density = m_referenceDensity * exp( m_compressibility * (pressure - m_referencePressure) ); 
    viscosity = m_referenceViscosity * exp( m_viscosibility * (pressure - m_referencePressure) ); 
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void compute( real64 const pressure,
                        real64 & density,
                        real64 & dDensity_dPressure,
                        real64 & viscosity,
                        real64 & dViscosity_dPressure ) const override
  {
    compute( pressure, density, viscosity );

    dDensity_dPressure = density * m_compressibility; 
    dViscosity_dPressure = viscosity * m_viscosibility;  
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
                        real64 & dInternalEnergy_dTemperature ) const override
  {
    /// Compute the density and viscosity 
    density = m_referenceDensity * exp( m_compressibility * (pressure - m_referencePressure)
                                        + m_thermalExpansionCoeff * (temperature - m_referenceTemperature) ); 

    dDensity_dPressure = density * m_compressibility; 
    dDensity_dTemperature = density * m_thermalExpansionCoeff; 

    viscosity = m_referenceViscosity * exp( m_viscosibility * (pressure - m_referencePressure) ); 

    dViscosity_dPressure = viscosity * m_viscosibility; 
    dViscosity_dTemperature = 0.0; 

    /// Compute the internal energy (only sensitive to temperature)
    internalEnergy = m_referenceInternalEnergy + m_volumetricHeatCapacity * ( temperature - m_referenceTemperature ); 
    dInternalEnergy_dPressure = 0.0; 
    dInternalEnergy_dTemperature = m_volumetricHeatCapacity; 
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
             m_dIntEnergy_dTemp[k][q] );
  }

private:
  /// Fluid compressibility 
  real64 m_compressibility; 

  /// Fluid thermal expansion coefficient 
  real64 m_thermalExpansionCoeff; 

  /// Fluid viscosibility 
  real64 m_viscosibility; 

  /// Fluid volumetric heat capacity
  real64 m_volumetricHeatCapacity;

  /// reference pressure parameter
  real64 m_referencePressure;

  /// reference temperature parameter
  real64 m_referenceTemperature;

  /// reference density parameter
  real64 m_referenceDensity;

  /// reference viscosity parameter
  real64 m_referenceViscosity;

  /// Fluid reference internal energy
  real64 m_referenceInternalEnergy;

};

class ThermalSinglePhaseFluid : public SingleFluidBase
{
public:

  ThermalSinglePhaseFluid( string const & name, Group * const parent );

  virtual ~ThermalSinglePhaseFluid() override;

  static string catalogName() { return "ThermalSinglePhaseFluid"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /// Type of kernel wrapper for in-kernel update (TODO: support multiple EAT, not just linear)
  using KernelWrapper = ThermalSinglePhaseUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  struct viewKeyStruct
  {
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
  };

  real64 defaultDensity() const override final { return m_defaultDensity; }
  real64 defaultViscosity() const override final { return m_defaultViscosity; }

protected:

  virtual void postProcessInput() override;

private:

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
};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_FLUID_COMPRESSIBLESINGLEPHASEFLUID_HPP_ */
