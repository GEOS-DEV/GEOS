/* SPDX-License-Identifier: LGPL-2.1-only */
#ifndef GEOS_CONSTITUTIVE_SOLID_THERMALLAW_HPP_
#define GEOS_CONSTITUTIVE_SOLID_THERMALLAW_HPP_

#include "constitutive/solid/ThermomechanicalLaw.hpp"
#include <cstddef>
#include <cmath>
#include <stdexcept>

namespace geos
{
namespace constitutive
{

/** Caloric and Fourier transport policies. Energies are specific, not volumetric.
 * These policies are independent of GEOS storage and of the mechanical kernel.
 */
struct ConstantThermalLaw
{
  static char const * prefix() { return "Thermal"; }
  double capacity = 1.0;
  double conductivity = 0.0;
  double referenceTemperature = 293.15;
  double capacitySlope = 0.0;
  double conductivitySlope = 0.0;

  double specificHeat( double const ) const { return capacity; }
  double thermalConductivity( double const ) const { return conductivity; }
  double energy( double const temperature ) const
  { return capacity * ( temperature - referenceTemperature ); }
  double temperature( double const specificEnergy ) const
  { return referenceTemperature + specificEnergy / capacity; }
};

/** c(T)=c0+c1(T-Tr), k(T)=k0+k1(T-Tr). Invert only on the c(T)>0 branch. */
struct LinearThermalLaw : ConstantThermalLaw
{
  static char const * prefix() { return "LinearThermal"; }
  double specificHeat( double const t ) const
  { return capacity + capacitySlope * ( t - referenceTemperature ); }
  double thermalConductivity( double const t ) const
  { return conductivity + conductivitySlope * ( t - referenceTemperature ); }
  double energy( double const t ) const
  {
    double const d = t - referenceTemperature;
    return d * ( capacity + 0.5 * capacitySlope * d );
  }
  double temperature( double const u ) const
  {
    if( !( std::abs( capacitySlope ) > 0.0 ) ) return ConstantThermalLaw::temperature( u );
    double const discriminant = capacity * capacity + 2.0 * capacitySlope * u;
    if( !( discriminant > 0.0 ) || !std::isfinite( discriminant ) )
      throw std::domain_error( "Thermal energy is outside the positive-capacity branch" );
    return referenceTemperature + 2.0 * u / ( capacity + std::sqrt( discriminant ) );
  }
};

/** Host-side capability queried once per particle subregion, never in a device loop. */
class ThermalMaterial
{
public:
  virtual ~ThermalMaterial() = default;
  virtual double thermalSpecificHeat( double temperature ) const = 0;
  virtual double thermalConductivity( double temperature ) const = 0;
  virtual double thermalEnergy( double temperature ) const = 0;
  virtual double thermalTemperature( double specificEnergy ) const = 0;
  virtual bool solverThermalExpansion() const = 0;
  virtual bool supportsThermomechanicalEnergy() const = 0;
  virtual thermomechanics::Matrix thermomechanicalCompliance( std::size_t particle ) const = 0;
  virtual double plasticHeatFraction() const = 0;
  virtual void prepareThermomechanicalTrial( std::size_t particle, double temperature ) = 0;

  // Temperature is recovered from TOTAL energy with the mechanical state held
  // fixed, never by interpreting elastic or retained plastic energy as heat.
  double internalEnergy( double t, thermomechanics::EnergyState const & state ) const
  { return thermalEnergy( t )+state.stored( t ); }
  double stateHeatCapacity( double t, thermomechanics::EnergyState const & state ) const
  { return thermalSpecificHeat( t )+state.thermoelasticCoefficient; }
  double temperatureFromInternalEnergy( double u, double initial,
                                        thermomechanics::EnergyState const & state ) const
  {
    return thermomechanics::solveTemperature( initial,
      1.0e-11*std::max( 1.0, std::abs( thermalSpecificHeat( initial )*std::max( 1.0, initial ) ) ), 100,
      [&]( double t )
      {
        if( !( thermalSpecificHeat( t ) > 0.0 ) || !( stateHeatCapacity( t, state ) > 0.0 ) )
          throw std::domain_error( "Nonpositive heat capacity of the frozen mechanical state" );
        return internalEnergy( t, state )-u;
      } );
  }
};

} // namespace constitutive
} // namespace geos
#endif
