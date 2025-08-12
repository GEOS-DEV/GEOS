#include "InvariantImmiscibleFluid.hpp"
#include "common/format/StringUtilities.hpp"
#include "constitutive/ConstitutiveBase.hpp"

namespace geos
{
namespace constitutive
{

InvariantImmiscibleFluid::InvariantImmiscibleFluid( string const & name, Group * parent )
  : MultiFluidBase( name, parent )
{
  // Override input flags for mandatory options
  
  registerWrapper( viewKeyStruct::componentNamesString(), &m_componentNames )
    .setInputFlag( dataRepository::InputFlags::REQUIRED )
    .setDescription( "List of fluid components (e.g. CH4, H2O, C5H12)" );
  
  registerWrapper( viewKeyStruct::phaseNamesString(), &m_phaseNames )
    .setInputFlag( dataRepository::InputFlags::REQUIRED )
    .setDescription( "List of fluid phases (e.g. gas, water, oil)" );

  registerWrapper( viewKeyStruct::componentMolarWeightString(), &m_componentMolarWeight )
    .setInputFlag( dataRepository::InputFlags::REQUIRED )
    .setDescription( "Molar weights of components" );

  // Densities: constant phase densities
  registerWrapper( "densities", &m_densities )
    .setInputFlag( dataRepository::InputFlags::REQUIRED )
    .setDescription( "Constant phase densities" );

  // Viscosities: constant phase viscosities
  registerWrapper( "viscosities", &m_viscosities )
    .setInputFlag( dataRepository::InputFlags::REQUIRED )
    .setDescription( "Constant phase viscosities" );
}

integer InvariantImmiscibleFluid::getWaterPhaseIndex() const
{
  // find index of water phase in user ordering
  for( size_t i=0; i< m_phaseNames.size(); ++i )
  {
    if( stringutilities::toLower( m_phaseNames[i] ) == "water" )
      return static_cast< integer >(i);
  }
  return -1;
}

void InvariantImmiscibleFluid::postInputInitialization()
{
  // check base inputs
  MultiFluidBase::postInputInitialization();

  integer numPhase = numFluidPhases();
  // check densities and viscosities size
  GEOS_THROW_IF_NE_MSG( m_densities.size(), numPhase,
                        GEOS_FMT( "%s: 'Densities' must have %d values", getFullName(), numPhase ), InputError );
  GEOS_THROW_IF_NE_MSG( m_viscosities.size(), numPhase,
                        GEOS_FMT( "%s: 'Viscosities' must have %d values", getFullName(), numPhase ), InputError );
}

GEOS_HOST_DEVICE
InvariantImmiscibleFluid::KernelWrapper InvariantImmiscibleFluid::createKernelWrapper() const
{
  return KernelWrapper(
    m_densities.toViewConst(),
    m_viscosities.toViewConst(),
    m_componentMolarWeight.toViewConst(),
    m_useMass,
    const_cast< MultiFluidBase::PhaseProp & >(m_phaseFraction).toView(),
    const_cast< MultiFluidBase::PhaseProp & >(m_phaseDensity).toView(),
    const_cast< MultiFluidBase::PhaseProp & >(m_phaseMassDensity).toView(),
    const_cast< MultiFluidBase::PhaseProp & >(m_phaseViscosity).toView(),
    const_cast< MultiFluidBase::PhaseProp & >(m_phaseEnthalpy).toView(),
    const_cast< MultiFluidBase::PhaseProp & >(m_phaseInternalEnergy).toView(),
    const_cast< MultiFluidBase::PhaseComp & >(m_phaseCompFraction).toView(),
    const_cast< MultiFluidBase::FluidProp & >(m_totalDensity).toView()
    );
}


GEOS_HOST_DEVICE
void InvariantImmiscibleFluid::update( localIndex const k,
                                       localIndex const q,
                                       real64 const pressure,
                                       real64 const temperature,
                                       arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const
{
  // Create a kernel wrapper for this fluid model (works on both host and device)
  KernelWrapper kernelWrapper = createKernelWrapper();

  // Call the kernel wrapper's update method
  kernelWrapper.update( k, q, pressure, temperature, composition );
}

void InvariantImmiscibleFluid::checkTablesParameters( real64 pressure, real64 temperature ) const
{
  // No PVT tables for this model, so nothing to check
  GEOS_UNUSED_VAR( pressure );
  GEOS_UNUSED_VAR( temperature );
}

GEOS_HOST_DEVICE
void InvariantImmiscibleFluid::KernelWrapper::compute( real64 const pressure,
                                                       real64 const temperature,
                                                       arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                                                       PhaseProp::SliceType const phaseFraction,
                                                       PhaseProp::SliceType const phaseDensity,
                                                       PhaseProp::SliceType const phaseMassDensity,
                                                       PhaseProp::SliceType const phaseViscosity,
                                                       PhaseProp::SliceType const phaseEnthalpy,
                                                       PhaseProp::SliceType const phaseInternalEnergy,
                                                       PhaseComp::SliceType const phaseCompFraction,
                                                       FluidProp::SliceType const totalDensity ) const
{
  // Note: pressure and temperature are marked as unused as this model uses constant properties
  GEOS_UNUSED_VAR( pressure, temperature );

  using Deriv = constitutive::multifluid::DerivativeOffset;
  
  integer nPhase = phaseDensity.value.size();
  integer nComp = phaseCompFraction.value.size( 1 );

  for( integer ip=0; ip<nPhase; ++ip )
  {
    // Phase fractions are equal to the composition for each phase
    phaseFraction.value[ip] = composition[ip];
    phaseFraction.derivs[ip][Deriv::dP] = 0.0;
    for( integer ic = 0; ic < nComp; ++ic )
    {
      phaseFraction.derivs[ip][Deriv::dC+ic] = (ip == ic) ? 1.0 : 0.0;
    }
      
    // densities and viscosities constant
    real64 const mult = m_useMass ? 1.0 : 1.0 / m_componentMolarWeight[ip];
    phaseDensity.value[ip] = m_densities[ip] * mult;
    phaseMassDensity.value[ip] = m_densities[ip];
    phaseViscosity.value[ip] = m_viscosities[ip];
    // derivatives
    for( integer ic = 0; ic < nComp; ++ic )
    {
      phaseDensity.derivs[ip][Deriv::dC+ic] = 0.0;
      phaseMassDensity.derivs[ip][Deriv::dC+ic] = 0.0;
      phaseViscosity.derivs[ip][Deriv::dC+ic] = 0.0;
    }

    // zero enthalpy/internal energy
    phaseEnthalpy.value[ip] = 0.0;
    phaseInternalEnergy.value[ip] = 0.0;
    for( integer d=0; d<phaseEnthalpy.derivs[ip].size( 0 ); ++d )
    {
      phaseEnthalpy.derivs[ip][d] = 0.0;
      phaseInternalEnergy.derivs[ip][d] = 0.0;
    }

    // phase composition identity: each phase pure component
    for( integer ic=0; ic<nComp; ++ic )
    {
      phaseCompFraction.value[ip][ic] = (ic==ip ? 1.0 : 0.0);
      for( integer d=0; d<phaseCompFraction.derivs[ip][ic].size( 0 ); ++d )
        phaseCompFraction.derivs[ip][ic][d] = 0.0;
    }
  }

  // total density: sum over phases of phaseMassDensity * phaseFraction
  computeTotalDensity( phaseFraction,
                       phaseDensity,
                       totalDensity );
}

GEOS_HOST_DEVICE
void InvariantImmiscibleFluid::KernelWrapper::update( localIndex const k,
                                                      localIndex const q,
                                                      real64 const pressure,
                                                      real64 const temperature,
                                                      arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const
{
  // Create slice views of the member variables at cell k, quadrature point q
  PhaseProp::SliceType phaseFractionSlice{m_phaseFraction.value[k][q], m_phaseFraction.derivs[k][q]};
  PhaseProp::SliceType phaseDensitySlice{m_phaseDensity.value[k][q], m_phaseDensity.derivs[k][q]};
  PhaseProp::SliceType phaseMassDensitySlice{m_phaseMassDensity.value[k][q], m_phaseMassDensity.derivs[k][q]};
  PhaseProp::SliceType phaseViscositySlice{m_phaseViscosity.value[k][q], m_phaseViscosity.derivs[k][q]};
  PhaseProp::SliceType phaseEnthalpySlice{m_phaseEnthalpy.value[k][q], m_phaseEnthalpy.derivs[k][q]};
  PhaseProp::SliceType phaseInternalEnergySlice{m_phaseInternalEnergy.value[k][q], m_phaseInternalEnergy.derivs[k][q]};
  PhaseComp::SliceType phaseCompFractionSlice{m_phaseCompFraction.value[k][q], m_phaseCompFraction.derivs[k][q]};
  FluidProp::SliceType totalDensitySlice{m_totalDensity.value[k][q], m_totalDensity.derivs[k][q]};

  // call compute with slice views
  compute( pressure, temperature, composition,
           phaseFractionSlice,
           phaseDensitySlice,
           phaseMassDensitySlice,
           phaseViscositySlice,
           phaseEnthalpySlice,
           phaseInternalEnergySlice,
           phaseCompFractionSlice,
           totalDensitySlice );
}


// Register this fluid model
REGISTER_CATALOG_ENTRY( ConstitutiveBase, InvariantImmiscibleFluid, string const &, geos::dataRepository::Group * const )

} // namespace constitutive
} // namespace geos
