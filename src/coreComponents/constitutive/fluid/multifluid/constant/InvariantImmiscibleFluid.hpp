#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_INVARIANTIMMISCIBLEFLUID_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_INVARIANTIMMISCIBLEFLUID_HPP_

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief Fluid model with constant phase properties (density, viscosity) independent of P, T.
 *        Formation volume factors are assumed = 1.
 */
class InvariantImmiscibleFluid : public MultiFluidBase
{
public:
  InvariantImmiscibleFluid( string const & name, Group * parent );

  static string catalogName() { return "InvariantImmiscibleFluid"; }
  virtual string getCatalogName() const override { return catalogName(); }

  static constexpr bool isThermalType() { return false; }

  virtual integer getWaterPhaseIndex() const override;

  void checkTablesParameters( real64 pressure, real64 temperature ) const override;

  /**
   * @brief Kernel wrapper class for InvariantImmiscibleFluid
   *        This kernel can be called on the GPU
   */
  class KernelWrapper final : public MultiFluidBase::KernelWrapper
  {
public:

    /**
     * @brief Compute fluid properties for an invariant immiscible fluid model
     *
     * This function calculates the properties of multiple immiscible fluid phases:
     * - Phase fractions directly from composition (1:1 mapping)
     * - Phase densities from constant input values (adjusted for molar/mass basis)
     * - Phase viscosities from constant input values
     * - Phase enthalpy and internal energy (set to zero in this simple model)
     * - Phase component fractions (identity mapping - each phase is a pure component)
     * - Total fluid density as weighted sum of phase densities
     *
     *
     * @param pressure        The pressure at which to evaluate fluid properties (unused in this model)
     * @param temperature     The temperature at which to evaluate fluid properties (unused in this model)
     * @param composition     Array of composition values (directly maps to phase fractions)
     * @param phaseFraction   Output view for phase fractions and their derivatives
     * @param phaseDensity    Output view for phase densities and their derivatives
     * @param phaseMassDensity Output view for phase mass densities and their derivatives
     * @param phaseViscosity  Output view for phase viscosities and their derivatives
     * @param phaseEnthalpy   Output view for phase enthalpies and their derivatives
     * @param phaseInternalEnergy Output view for phase internal energies and their derivatives
     * @param phaseCompFraction Output view for phase component fractions and their derivatives
     * @param totalDensity    Output view for total fluid density and its derivatives
     */
    GEOS_HOST_DEVICE
    void compute( real64 const pressure,
                  real64 const temperature,
                  arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                  MultiFluidBase::PhaseProp::SliceType const phaseFraction,
                  MultiFluidBase::PhaseProp::SliceType const phaseDensity,
                  MultiFluidBase::PhaseProp::SliceType const phaseMassDensity,
                  MultiFluidBase::PhaseProp::SliceType const phaseViscosity,
                  MultiFluidBase::PhaseProp::SliceType const phaseEnthalpy,
                  MultiFluidBase::PhaseProp::SliceType const phaseInternalEnergy,
                  MultiFluidBase::PhaseComp::SliceType const phaseCompFraction,
                  MultiFluidBase::FluidProp::SliceType const totalDensity ) const override;

    GEOS_HOST_DEVICE
    void update( localIndex const k,
                 localIndex const q,
                 real64 const pressure,
                 real64 const temperature,
                 arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const override;

private:
    friend class InvariantImmiscibleFluid;

    arrayView1d< real64 const > m_densities;
    arrayView1d< real64 const > m_viscosities;

    /**
     * @brief Constructor for the class doing in-kernel fluid property updates
     * @param[in] densities constant phase densities
     * @param[in] viscosities constant phase viscosities
     * @param[in] componentMolarWeight component molecular weights
     * @param[in] useMass flag to decide whether we return mass or molar densities
     * @param[in] phaseFraction phase fractions (+ derivatives) in the cell
     * @param[in] phaseDensity phase mass/molar densities (+ derivatives) in the cell
     * @param[in] phaseMassDensity phase mass densities (+ derivatives) in the cell
     * @param[in] phaseViscosity phase viscosities (+ derivatives) in the cell
     * @param[in] phaseEnthalpy phase enthalpies (+ derivatives) in the cell
     * @param[in] phaseInternalEnergy phase internal energies (+ derivatives) in the cell
     * @param[in] phaseCompFraction phase component fractions (+ derivatives) in the cell
     * @param[in] totalDensity total density in the cell
     */
    KernelWrapper( arrayView1d< real64 const > densities,
                   arrayView1d< real64 const > viscosities,
                   arrayView1d< real64 const > componentMolarWeight,
                   bool const useMass,
                   MultiFluidBase::PhaseProp::ViewType const phaseFraction,
                   MultiFluidBase::PhaseProp::ViewType const phaseDensity,
                   MultiFluidBase::PhaseProp::ViewType const phaseMassDensity,
                   MultiFluidBase::PhaseProp::ViewType const phaseViscosity,
                   MultiFluidBase::PhaseProp::ViewType const phaseEnthalpy,
                   MultiFluidBase::PhaseProp::ViewType const phaseInternalEnergy,
                   MultiFluidBase::PhaseComp::ViewType const phaseCompFraction,
                   MultiFluidBase::FluidProp::ViewType const totalDensity )
      : MultiFluidBase::KernelWrapper( std::move( componentMolarWeight ),
                                       useMass,
                                       std::move( phaseFraction ),
                                       std::move( phaseDensity ),
                                       std::move( phaseMassDensity ),
                                       std::move( phaseViscosity ),
                                       std::move( phaseEnthalpy ),
                                       std::move( phaseInternalEnergy ),
                                       std::move( phaseCompFraction ),
                                       std::move( totalDensity )),
      m_densities( std::move( densities ) ),
      m_viscosities( std::move( viscosities ) )
    {}

  };

  /**
   * @brief Create a KernelWrapper for this fluid model
   * @return A KernelWrapper instance for this fluid model
   */
  KernelWrapper createKernelWrapper() const;

protected:
  virtual void postInputInitialization() override;

private:
  array1d< real64 > m_densities;
  array1d< real64 > m_viscosities;
};


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
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

  integer const nPhase = numPhases();
  integer const nComp = numComponents();

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
GEOS_FORCE_INLINE
void InvariantImmiscibleFluid::KernelWrapper::update( localIndex const k,
                                                      localIndex const q,
                                                      real64 const pressure,
                                                      real64 const temperature,
                                                      arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const
{
  // Compute with slice views at cell k, quadrature point q
  compute( pressure, temperature, composition,
           m_phaseFraction( k, q ),
           m_phaseDensity( k, q ),
           m_phaseMassDensity( k, q ),
           m_phaseViscosity( k, q ),
           m_phaseEnthalpy( k, q ),
           m_phaseInternalEnergy( k, q ),
           m_phaseCompFraction( k, q ),
           m_totalDensity( k, q ) );
}

} // namespace constitutive
} // namespace geos

#endif // GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_INVARIANTIMMISCIBLEFLUID_HPP_
