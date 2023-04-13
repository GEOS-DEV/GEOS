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
 * @file ReactiveBrineFluid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_REACTIVEBRINEFLUID_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_REACTIVEBRINEFLUID_HPP_

#include "codingUtilities/EnumStrings.hpp"
#include "constitutive/fluid/ReactiveMultiFluid.hpp"
#include "constitutive/fluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/PhaseModel.hpp"
#include "constitutive/fluid/PVTFunctions/BrineEnthalpy.hpp"
#include "constitutive/fluid/PVTFunctions/NoOpPVTFunction.hpp"
#include "constitutive/fluid/PVTFunctions/WaterDensity.hpp"
#include "constitutive/fluid/PVTFunctions/PhillipsBrineViscosity.hpp"
#include "PVTFunctions/PureWaterProperties.hpp"



#include <memory>

namespace geos
{

namespace constitutive
{

template< typename PHASE >
class ReactiveBrineFluid : public ReactiveMultiFluid
{
public:

  using exec_policy = parallelDevicePolicy<>;

  ReactiveBrineFluid( string const & name,
                      Group * const parent );

  virtual std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                Group * const parent ) const override;

  static string catalogName();

  virtual string getCatalogName() const override { return catalogName(); }

  virtual bool isThermal() const override final;

  /**
   * @brief Kernel wrapper class for ReactiveBrineFluid.
   */
  class KernelWrapper final : public ReactiveMultiFluid::KernelWrapper
  {
public:

    GEOSX_HOST_DEVICE
    virtual void compute( real64 const pressure,
                          real64 const temperature,
                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseEnthalpy,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseInternalEnergy,
                          arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & phaseCompFraction,
                          real64 & totalDensity ) const override;

    GEOSX_HOST_DEVICE
    virtual void compute( real64 const pressure,
                          real64 const temperature,
                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                          PhaseProp::SliceType const phaseFraction,
                          PhaseProp::SliceType const phaseDensity,
                          PhaseProp::SliceType const phaseMassDensity,
                          PhaseProp::SliceType const phaseViscosity,
                          PhaseProp::SliceType const phaseEnthalpy,
                          PhaseProp::SliceType const phaseInternalEnergy,
                          PhaseComp::SliceType const phaseCompFraction,
                          FluidProp::SliceType const totalDensity ) const override;

    GEOSX_HOST_DEVICE
    virtual void update( localIndex const k,
                         localIndex const q,
                         real64 const pressure,
                         real64 const temperature,
                         arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const override;

    virtual void updateChemistry( localIndex const k,
                                  localIndex const q,
                                  real64 const pressure,
                                  real64 const temperature,
                                  arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const override;

private:

    friend class ReactiveBrineFluid;

    KernelWrapper( PHASE const & phase,
                   arrayView1d< real64 const > componentMolarWeight,
                   bool const useMass,
                   bool const isThermal,
                   PhaseProp::ViewType phaseFraction,
                   PhaseProp::ViewType phaseDensity,
                   PhaseProp::ViewType phaseMassDensity,
                   PhaseProp::ViewType phaseViscosity,
                   PhaseProp::ViewType phaseEnthalpy,
                   PhaseProp::ViewType phaseInternalEnergy,
                   PhaseComp::ViewType phaseCompFraction,
                   FluidProp::ViewType totalDensity,
                   integer const numPrimarySpecies,
                   chemicalReactions::EquilibriumReactions const & equilibriumReactions,
                   chemicalReactions::KineticReactions const & kineticReactions,
                   arrayView2d< real64, compflow::USD_COMP > const & primarySpeciesConcentration,
                   arrayView2d< real64, compflow::USD_COMP > const & secondarySpeciesConcentration,
                   arrayView2d< real64, compflow::USD_COMP > const & primarySpeciesTotalConcentration,
                   arrayView2d< real64, compflow::USD_COMP > const & kineticReactionRates );


    /// Flag to specify whether the model is thermal or not
    bool m_isThermal;

    /// Brine constitutive kernel wrappers
    typename PHASE::KernelWrapper m_phase;

  };

  virtual integer getWaterPhaseIndex() const override final;

  /**
   * @brief Names of the submodels for input
   */
  enum class SubModelInputNames : integer
  {
    DENSITY,         ///< the keyword for the density model
    VISCOSITY,       ///< the keyword for the viscosity model
    ENTHALPY         ///< the keyword for the enthalpy model
  };

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  struct viewKeyStruct : ReactiveMultiFluid::viewKeyStruct
  {
    static constexpr char const * phasePVTParaFilesString() { return "phasePVTParaFiles"; }
  };

protected:

  virtual void postProcessInput() override;

private:

  void createPVTModels();

  /// Names of the files defining the viscosity and density models
  path_array m_phasePVTParaFiles;

  /// Brine constitutive models
  std::unique_ptr< PHASE > m_phase;

};

// these aliases are useful in constitutive dispatch
using ReactiveBrine =
  ReactiveBrineFluid< PhaseModel< PVTProps::WaterDensity, PVTProps::PhillipsBrineViscosity, PVTProps::NoOpPVTFunction > >;
using ReactiveBrineThermal =
  ReactiveBrineFluid< PhaseModel< PVTProps::WaterDensity, PVTProps::PhillipsBrineViscosity, PVTProps::BrineEnthalpy > >;

template< typename PHASE >
GEOSX_HOST_DEVICE
inline void
ReactiveBrineFluid< PHASE >::KernelWrapper::
  compute( real64 pressure,
           real64 temperature,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseEnthalpy,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseInternalEnergy,
           arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & phaseCompFraction,
           real64 & totalDensity ) const
{
  integer constexpr numComp = chemicalReactions::ReactionsBase::maxNumPrimarySpecies;
  integer constexpr numPhase = 1;

  // 1. We perform a sort of single phase flash
  stackArray1d< real64, numComp > compMoleFrac( composition.size() );
  phaseFraction[0] = 1.0; // it's a single phase system
  for( integer ic = 0; ic < composition.size(); ++ic )
  {
    compMoleFrac[ic] = composition[ic];
    phaseCompFraction[0][ic] = composition[ic];
  }

  // 2. Compute phase fractions and phase component fractions
  real64 const temperatureInCelsius = temperature - 273.15;

  // 3. Compute phase density and phase viscoisty
  m_phase.density.compute( pressure,
                           temperatureInCelsius,
                           phaseCompFraction[0].toSliceConst(),
                           phaseDensity[0],
                           m_useMass );
  m_phase.viscosity.compute( pressure,
                             temperatureInCelsius,
                             phaseCompFraction[0].toSliceConst(),
                             phaseViscosity[0],
                             m_useMass );

  // for now, we have to compute the phase mass density here
  m_phase.density.compute( pressure,
                           temperatureInCelsius,
                           phaseCompFraction[0].toSliceConst(),
                           phaseMassDensity[0],
                           true );

  // 4. Compute enthalpy and internal energy
  if( m_isThermal )
  {

    m_phase.enthalpy.compute( pressure,
                              temperatureInCelsius,
                              phaseCompFraction[0].toSliceConst(),
                              phaseEnthalpy[0],
                              m_useMass );

    computeInternalEnergy< numComp, numPhase >( pressure,
                                                phaseFraction,
                                                phaseMassDensity,
                                                phaseEnthalpy,
                                                phaseInternalEnergy );

  }

  // 6. Compute total fluid mass/molar density
  computeTotalDensity< numComp, numPhase >( phaseFraction,
                                            phaseDensity,
                                            totalDensity );
}

template< typename PHASE >
GEOSX_HOST_DEVICE
inline void
ReactiveBrineFluid< PHASE >::KernelWrapper::
  compute( real64 const pressure,
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
  integer constexpr numComp = chemicalReactions::ReactionsBase::maxNumPrimarySpecies;

  // 1. We perform a sort of single phase flash
  stackArray1d< real64, numComp > compMoleFrac( composition.size() );

  phaseFraction.value[0] = 1.0; // it's a single phase system
  for( integer ic = 0; ic < composition.size(); ++ic )
  {
    compMoleFrac[ic] = composition[ic];
    phaseCompFraction.value[0][ic] = composition[ic];
  }

  real64 const temperatureInCelsius = temperature - 273.15;

  // 2. Compute phase densities and phase viscosities
  m_phase.density.compute( pressure,
                           temperatureInCelsius,
                           phaseCompFraction.value[0].toSliceConst(), phaseCompFraction.derivs[0].toSliceConst(),
                           phaseDensity.value[0], phaseDensity.derivs[0],
                           m_useMass );

  m_phase.viscosity.compute( pressure,
                             temperatureInCelsius,
                             phaseCompFraction.value[0].toSliceConst(), phaseCompFraction.derivs[0].toSliceConst(),
                             phaseViscosity.value[0], phaseViscosity.derivs[0],
                             m_useMass );


  // for now, we have to compute the phase mass density here
  m_phase.density.compute( pressure,
                           temperatureInCelsius,
                           phaseCompFraction.value[0].toSliceConst(), phaseCompFraction.derivs[0].toSliceConst(),
                           phaseMassDensity.value[0], phaseMassDensity.derivs[0],
                           true );

  // 3. Compute enthalpy and internal energy
  if( m_isThermal )
  {
    m_phase.enthalpy.compute( pressure,
                              temperatureInCelsius,
                              phaseCompFraction.value[0].toSliceConst(), phaseCompFraction.derivs[0].toSliceConst(),
                              phaseEnthalpy.value[0], phaseEnthalpy.derivs[0],
                              m_useMass );

    computeInternalEnergy( pressure,
                           phaseFraction,
                           phaseMassDensity,
                           phaseEnthalpy,
                           phaseInternalEnergy );
  }

  // 6. Compute total fluid mass/molar density and derivatives
  computeTotalDensity( phaseFraction,
                       phaseDensity,
                       totalDensity );
}

template< typename PHASE >
GEOSX_HOST_DEVICE inline void
ReactiveBrineFluid< PHASE >::KernelWrapper::
  update( localIndex const k,
          localIndex const q,
          real64 const pressure,
          real64 const temperature,
          arraySlice1d< geos::real64 const, compflow::USD_COMP - 1 > const & composition ) const
{
  compute( pressure,
           temperature,
           composition,
           m_phaseFraction( k, q ),
           m_phaseDensity( k, q ),
           m_phaseMassDensity( k, q ),
           m_phaseViscosity( k, q ),
           m_phaseEnthalpy( k, q ),
           m_phaseInternalEnergy( k, q ),
           m_phaseCompFraction( k, q ),
           m_totalDensity( k, q ) );
}

template< typename PHASE >
inline void
ReactiveBrineFluid< PHASE >::KernelWrapper::
  updateChemistry( localIndex const k,
                   localIndex const q,
                   real64 const pressure,
                   real64 const temperature,
                   arraySlice1d< geos::real64 const, compflow::USD_COMP - 1 > const & composition ) const

{
  real64 const totalMolecularWeight = PVTProps::PureWaterProperties::MOLECULAR_WEIGHT;

  convertMoleFractionToMolarity( m_totalDensity( k, q ).value,
                                 totalMolecularWeight,
                                 composition,
                                 m_primarySpeciesTotalConcentration[k] );

  computeChemistry( pressure,
                    temperature,
                    m_primarySpeciesTotalConcentration[k],
                    m_primarySpeciesConcentration[k],
                    m_secondarySpeciesConcentration[k],
                    m_kineticReactionRates[k] );
}


} // namespace constitutive

} // namespace geos

#endif //GEOSX_CONSTITUTIVE_FLUID_REACTIVEBRINEFLUID_HPP_
