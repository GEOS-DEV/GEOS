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
#ifndef GEOS_USE_DEVICE_CODE
  // Create a kernel wrapper for this fluid model
  KernelWrapper kernelWrapper = createKernelWrapper();

  // Call the kernel wrapper's update method
  kernelWrapper.update( k, q, pressure, temperature, composition );
#else
  // For device code path, create slice views directly
  PhaseProp::SliceType phaseFractionSlice{m_phaseFraction.value[k][q], m_phaseFraction.derivs[k][q]};
  PhaseProp::SliceType phaseDensitySlice{m_phaseDensity.value[k][q], m_phaseDensity.derivs[k][q]};
  PhaseProp::SliceType phaseMassDensitySlice{m_phaseMassDensity.value[k][q], m_phaseMassDensity.derivs[k][q]};
  PhaseProp::SliceType phaseViscositySlice{m_phaseViscosity.value[k][q], m_phaseViscosity.derivs[k][q]};
  PhaseProp::SliceType phaseEnthalpySlice{m_phaseEnthalpy.value[k][q], m_phaseEnthalpy.derivs[k][q]};
  PhaseProp::SliceType phaseInternalEnergySlice{m_phaseInternalEnergy.value[k][q], m_phaseInternalEnergy.derivs[k][q]};
  PhaseComp::SliceType phaseCompFractionSlice{m_phaseCompFraction.value[k][q], m_phaseCompFraction.derivs[k][q]};
  FluidProp::SliceType totalDensitySlice{m_totalDensity.value[k][q], m_totalDensity.derivs[k][q]};

  // Directly call compute with slice views
  compute( pressure, temperature, composition,
           phaseFractionSlice,
           phaseDensitySlice,
           phaseMassDensitySlice,
           phaseViscositySlice,
           phaseEnthalpySlice,
           phaseInternalEnergySlice,
           phaseCompFractionSlice,
           totalDensitySlice );
#endif
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
  GEOS_UNUSED_VAR( pressure, temperature );

  integer nPhase = phaseDensity.value.size();
  integer nComp = phaseCompFraction.value.size( 1 );
  
  // Create local storage vector for s with nPhase size
  real64 s[3];
  // Create container for derivatives of s with respect to composition
  real64 ds_dcomp[3][3];
  
  if (nPhase == 2) {
          // Precompute common terms
          real64 denom = m_densities[0] - composition[0] * m_densities[0] + composition[0] * m_densities[1];
          real64 denom_squared = denom * denom;

          // Loop without conditionals
          for (integer ip = 0; ip < 2; ++ip) {
              real64 sign = (ip == 0) ? 1.0 : -1.0;
              real64 numerator = (ip == 0) ? composition[0] * m_densities[1] : (-m_densities[0] + composition[0] * m_densities[0]);
              s[ip] = sign * numerator / denom;
              
              // Calculate derivatives with respect to composition[0]
              // For ip=0: d(s[0])/d(comp[0]) = d(comp[0]*m_densities[1]/denom)/d(comp[0])
              // For ip=1: d(s[1])/d(comp[0]) = d((-m_densities[0]+comp[0]*m_densities[0])/denom)/d(comp[0])
              real64 dnumerator_dcomp0 = (ip == 0) ? m_densities[1] : m_densities[0];
              real64 ddenom_dcomp0 = -m_densities[0] + m_densities[1];
              
              // Quotient rule: d(num/denom)/dx = (denom*dnum/dx - num*ddenom/dx)/denom^2
              ds_dcomp[ip][0] = sign * (denom * dnumerator_dcomp0 - numerator * ddenom_dcomp0) / denom_squared;
              
              // For 2-phase system, we have only one composition variable (comp[0])
              // Zero out other derivatives
              for (integer ic = 1; ic < nComp; ++ic) {
                  ds_dcomp[ip][ic] = 0.0;
              }
          }

    } else if (nPhase == 3) {
        // Precompute common terms
        real64 denom = composition[2]*m_densities[0]*m_densities[1] + composition[1]*m_densities[0]*m_densities[2] + m_densities[1]*m_densities[2]
                     - composition[1]*m_densities[1]*m_densities[2] - composition[2]*m_densities[1]*m_densities[2];
        real64 denom_squared = denom * denom;

        // Numerators for each ip
        real64 numerators[3];
        numerators[0] = -( -m_densities[1]*m_densities[2] + composition[1]*m_densities[1]*m_densities[2] + composition[2]*m_densities[1]*m_densities[2] );
        numerators[1] = composition[1]*m_densities[0]*m_densities[2];
        numerators[2] = composition[2]*m_densities[0]*m_densities[1];

        // Calculate derivatives of denominator with respect to each composition
        real64 ddenom_dcomp[3];
        ddenom_dcomp[0] = 0.0; // Derivative with respect to comp[0] - not present in denom
        ddenom_dcomp[1] = m_densities[0]*m_densities[2] - m_densities[1]*m_densities[2]; // d(denom)/d(comp[1])
        ddenom_dcomp[2] = m_densities[0]*m_densities[1] - m_densities[1]*m_densities[2]; // d(denom)/d(comp[2])

        // Loop without conditionals
        for (integer ip = 0; ip < 3; ++ip) {
            real64 sign = (ip == 0) ? -1.0 : -1.0; // all negatives
            s[ip] = sign * numerators[ip] / denom;
            
            // Calculate derivatives of numerator with respect to compositions
            real64 dnumerator_dcomp[3];
            
            if (ip == 0) {
                dnumerator_dcomp[0] = 0.0; // No comp[0] in numerator[0]
                dnumerator_dcomp[1] = -m_densities[1]*m_densities[2]; // d(numerator[0])/d(comp[1])
                dnumerator_dcomp[2] = -m_densities[1]*m_densities[2]; // d(numerator[0])/d(comp[2])
            } else if (ip == 1) {
                dnumerator_dcomp[0] = 0.0; // No comp[0] in numerator[1]
                dnumerator_dcomp[1] = m_densities[0]*m_densities[2]; // d(numerator[1])/d(comp[1])
                dnumerator_dcomp[2] = 0.0; // No comp[2] in numerator[1]
            } else { // ip == 2
                dnumerator_dcomp[0] = 0.0; // No comp[0] in numerator[2]
                dnumerator_dcomp[1] = 0.0; // No comp[1] in numerator[2]
                dnumerator_dcomp[2] = m_densities[0]*m_densities[1]; // d(numerator[2])/d(comp[2])
            }
            
            // Apply quotient rule for each composition variable
            for (integer ic = 0; ic < nComp; ++ic) {
                // Quotient rule: d(num/denom)/dx = (denom*dnum/dx - num*ddenom/dx)/denom^2
                ds_dcomp[ip][ic] = sign * (denom * dnumerator_dcomp[ic] - numerators[ip] * ddenom_dcomp[ic]) / denom_squared;
            }
        }

    } else {
        // Use device-compatible error handling
        #ifndef GEOS_USE_DEVICE_CODE
        std::cerr << "Unsupported number of phases: " << nPhase << std::endl;
        #endif
        
        // Initialize ds_dcomp to zeros for unsupported cases
        for (integer ip = 0; ip < 3; ++ip) {
            for (integer ic = 0; ic < 3; ++ic) {
                ds_dcomp[ip][ic] = 0.0;
            }
        }
    }
  

  for( integer ip=0; ip<nPhase; ++ip )
  {
    phaseFraction.value[ip] = s[ip];
    // Set derivatives of phase fraction with respect to composition
    for( integer d=0; d<phaseFraction.derivs[ip].size( 0 ); ++d )
      phaseFraction.derivs[ip][d] = 0.0; //ds_dcomp[ip][d];

    // densities and viscosities constant
    phaseDensity.value[ip] = m_densities[ip];
    phaseMassDensity.value[ip] = m_densities[ip];
    phaseViscosity.value[ip] = m_viscosities[ip];
    // derivatives
    for( integer d=0; d<phaseDensity.derivs[ip].size( 0 ); ++d )
    {
      phaseDensity.derivs[ip][d] = 0.0;
      phaseMassDensity.derivs[ip][d] = 0.0;
      phaseViscosity.derivs[ip][d] = 0.0;
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
