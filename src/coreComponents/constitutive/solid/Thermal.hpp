/* SPDX-License-Identifier: LGPL-2.1-only */
#ifndef GEOS_CONSTITUTIVE_SOLID_THERMAL_HPP_
#define GEOS_CONSTITUTIVE_SOLID_THERMAL_HPP_

#include "constitutive/solid/ThermalLaw.hpp"
#include "constitutive/solid/ElasticIsotropic.hpp"
#include "constitutive/solid/PerfectlyPlastic.hpp"
#include <type_traits>

namespace geos
{
namespace constitutive
{

/** Add a thermal policy to a solid without changing its mechanical kernel type.
 * ConstitutivePassThruMPM deliberately dispatches to BASE for mechanical updates.
 * Thermal transport uses the independent ThermalMaterial capability. This also
 * preserves createDerivedKernelUpdates forwarding for nested mechanical models.
 */
template< typename BASE, typename LAW = ConstantThermalLaw >
class Thermal : public BASE, public ThermalMaterial
{
public:
  using KernelWrapper = typename BASE::KernelWrapper;

  Thermal( string const & name, dataRepository::Group * const parent ) : BASE( name, parent )
  {
    using dataRepository::InputFlags;
    this->registerWrapper( "plasticHeatFraction", &m_plasticHeatFraction ).
      setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( 0.9 ).
      setDescription( "Fraction of plastic work converted to heat; the remainder is retained internal energy." );
    this->registerWrapper( "thermalYieldSoftening", &m_yieldSoftening ).
      setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( 0.0 ).
      setDescription( "J2 yield multiplier exp[-softening*(T-Tr)]; units 1/K. Supported for PerfectlyPlastic." );
    this->registerWrapper( "thermalReferenceYield", &m_referenceYield ).
      setSizedFromParent( 1 ).setInputFlag( InputFlags::FALSE ).setApplyDefaultValue( -1.0 ).
      setPlotLevel( dataRepository::PlotLevel::NOPLOT ).setRestartFlags( dataRepository::RestartFlags::WRITE_AND_READ );
    this->registerWrapper( "thermalSpecificHeat", &m_law.capacity ).
      setInputFlag( InputFlags::REQUIRED ).setDescription( "Specific heat at the reference temperature [energy/(mass K)]" );
    this->registerWrapper( "thermalConductivity", &m_law.conductivity ).
      setInputFlag( InputFlags::REQUIRED ).setDescription( "Isotropic conductivity at the reference temperature [power/(length K)]" );
    this->registerWrapper( "thermalReferenceTemperature", &m_law.referenceTemperature ).
      setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( 293.15 );
    this->registerWrapper( "thermalHeatCapacitySlope", &m_law.capacitySlope ).
      setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( 0.0 );
    this->registerWrapper( "thermalConductivitySlope", &m_law.conductivitySlope ).
      setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( 0.0 );
    this->registerWrapper( "solverThermalExpansion", &m_solverExpansion ).
      setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( 0 ).
      setDescription( "Apply isotropic thermal strain using defaultDrainedLinearTEC. Leave off for materials owning thermal expansion." );
  }

  static string catalogName() { return string( LAW::prefix() ) + BASE::catalogName(); }
  string getCatalogName() const override { return catalogName(); }
  string getMechanicalCatalogName() const override { return BASE::catalogName(); }
  LAW createThermalKernelUpdates() const { return m_law; }
  double thermalSpecificHeat( double const t ) const override { return m_law.specificHeat( t ); }
  double thermalConductivity( double const t ) const override { return m_law.thermalConductivity( t ); }
  double thermalEnergy( double const t ) const override { return m_law.energy( t ); }
  double thermalTemperature( double const u ) const override { return m_law.temperature( u ); }
  bool solverThermalExpansion() const override { return m_solverExpansion == 1; }
  double plasticHeatFraction() const override { return m_plasticHeatFraction; }

  static constexpr bool energySupported =
    std::is_same< BASE, ElasticIsotropic >::value ||
    std::is_same< BASE, PerfectlyPlastic >::value;
  bool supportsThermomechanicalEnergy() const override { return energySupported; }

  thermomechanics::Matrix thermomechanicalCompliance( std::size_t const p ) const override
  {
    if constexpr( energySupported )
    {
      // Query elasticity, never the algorithmic elastoplastic tangent.
      real64 c[6][6] = {};
      auto const wrapper = BASE::createKernelUpdates();
      wrapper.getElasticStiffness( static_cast< localIndex >( p ), 0, c );
      thermomechanics::Matrix matrix{};
      for( integer i = 0; i < 6; ++i )
        for( integer j = 0; j < 6; ++j ) matrix[i][j] = c[i][j];
      return thermomechanics::compliance( matrix );
    }
    throw std::domain_error( "Mechanical model needs a thermomechanical stored-energy adapter" );
  }

  void prepareThermomechanicalTrial( std::size_t const p, double const t ) override
  {
    static_cast< void >( p ); static_cast< void >( t );
    if constexpr( std::is_same< BASE, PerfectlyPlastic >::value )
    {
      char const * const key = "yieldStress";
      auto & yield = this->template getReference< array1d< real64 > >( key );
      GEOS_ERROR_IF( p >= static_cast< std::size_t >( m_referenceYield.size() ), "Missing thermal reference-yield state" );
      if( m_referenceYield[p] < 0.0 ) m_referenceYield[p] = yield[p];
      double const scale = std::exp( -m_yieldSoftening*( t-m_law.referenceTemperature ) );
      GEOS_ERROR_IF( !std::isfinite( scale ) || !( scale > 0.0 ), "Inadmissible thermal yield multiplier; reduce timestep or softening." );
      yield[p] = m_referenceYield[p]*scale;
    }
  }

protected:
  void postInputInitialization() override
  {
    BASE::postInputInitialization();
    GEOS_ERROR_IF( !std::isfinite( m_law.capacity ) || !( m_law.capacity > 0.0 ) ||
                   !std::isfinite( m_law.conductivity ) || m_law.conductivity < 0.0 ||
                   !std::isfinite( m_law.referenceTemperature ) || m_law.referenceTemperature < 0.0 ||
                   !std::isfinite( m_law.capacitySlope ) || !std::isfinite( m_law.conductivitySlope ),
                   this->getDataContext() << ": invalid thermal properties" );
    if constexpr( std::is_same< LAW, ConstantThermalLaw >::value )
    {
      GEOS_ERROR_IF( std::abs( m_law.capacitySlope ) > 0.0 || std::abs( m_law.conductivitySlope ) > 0.0,
                     "Use a LinearThermal material for temperature-dependent thermal properties." );
    }
    GEOS_ERROR_IF( !std::isfinite( m_plasticHeatFraction ) || m_plasticHeatFraction < 0.0 || m_plasticHeatFraction > 1.0 ||
                   !std::isfinite( m_yieldSoftening ) || m_yieldSoftening < 0.0,
                   "plasticHeatFraction must be in [0,1]; thermalYieldSoftening must be finite and nonnegative." );
    GEOS_ERROR_IF( m_yieldSoftening > 0.0 &&
                   !( std::is_same< BASE, PerfectlyPlastic >::value ),
                   "thermalYieldSoftening requires a PerfectlyPlastic base." );
    GEOS_ERROR_IF( m_solverExpansion != 0 && m_solverExpansion != 1,
                   "solverThermalExpansion must be 0 or 1." );
    GEOS_ERROR_IF( m_solverExpansion == 1 &&
                   ( BASE::catalogName() == "Graphite" || BASE::catalogName() == "Geomechanics" ),
                   "This mechanical model already owns thermal expansion; leave solverThermalExpansion=0." );
  }

private:
  LAW m_law;
  integer m_solverExpansion = 0;
  real64 m_plasticHeatFraction = 0.9;
  real64 m_yieldSoftening = 0.0;
  array1d< real64 > m_referenceYield;
};

using ThermalElasticIsotropic = Thermal< ElasticIsotropic >;

} // namespace constitutive
} // namespace geos
#endif
