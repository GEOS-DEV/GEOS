/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MultiFluidPVTPackageWrapper.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDPVTPACKAGEWRAPPER_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDPVTPACKAGEWRAPPER_HPP_

#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/fluid/MultiFluidUtils.hpp"

#include "pvt/pvt.hpp"

#include <memory>

namespace geosx
{

namespace constitutive
{

/**
 * @brief Kernel wrapper class for MultiFluidPVTPackage.
 * @note Not thread-safe, do not use with any parallel launch policy.
 */
class MultiFluidPVTPackageWrapperUpdate final : public MultiFluidBaseUpdate
{
public:

  MultiFluidPVTPackageWrapperUpdate( pvt::MultiphaseSystem & fluid,
                                     arrayView1d< pvt::PHASE_TYPE > const & phaseTypes,
                                     arrayView1d< real64 const > const & componentMolarWeight,
                                     bool useMass,
                                     arrayView3d< real64 > const & phaseFraction,
                                     arrayView3d< real64 > const & dPhaseFraction_dPressure,
                                     arrayView3d< real64 > const & dPhaseFraction_dTemperature,
                                     arrayView4d< real64 > const & dPhaseFraction_dGlobalCompFraction,
                                     arrayView3d< real64 > const & phaseDensity,
                                     arrayView3d< real64 > const & dPhaseDensity_dPressure,
                                     arrayView3d< real64 > const & dPhaseDensity_dTemperature,
                                     arrayView4d< real64 > const & dPhaseDensity_dGlobalCompFraction,
                                     arrayView3d< real64 > const & phaseMassDensity,
                                     arrayView3d< real64 > const & dPhaseMassDensity_dPressure,
                                     arrayView3d< real64 > const & dPhaseMassDensity_dTemperature,
                                     arrayView4d< real64 > const & dPhaseMassDensity_dGlobalCompFraction,
                                     arrayView3d< real64 > const & phaseViscosity,
                                     arrayView3d< real64 > const & dPhaseViscosity_dPressure,
                                     arrayView3d< real64 > const & dPhaseViscosity_dTemperature,
                                     arrayView4d< real64 > const & dPhaseViscosity_dGlobalCompFraction,
                                     arrayView4d< real64 > const & phaseCompFraction,
                                     arrayView4d< real64 > const & dPhaseCompFraction_dPressure,
                                     arrayView4d< real64 > const & dPhaseCompFraction_dTemperature,
                                     arrayView5d< real64 > const & dPhaseCompFraction_dGlobalCompFraction,
                                     arrayView2d< real64 > const & totalDensity,
                                     arrayView2d< real64 > const & dTotalDensity_dPressure,
                                     arrayView2d< real64 > const & dTotalDensity_dTemperature,
                                     arrayView3d< real64 > const & dTotalDensity_dGlobalCompFraction )
    : MultiFluidBaseUpdate( componentMolarWeight,
                            useMass,
                            phaseFraction,
                            dPhaseFraction_dPressure,
                            dPhaseFraction_dTemperature,
                            dPhaseFraction_dGlobalCompFraction,
                            phaseDensity,
                            dPhaseDensity_dPressure,
                            dPhaseDensity_dTemperature,
                            dPhaseDensity_dGlobalCompFraction,
                            phaseMassDensity,
                            dPhaseMassDensity_dPressure,
                            dPhaseMassDensity_dTemperature,
                            dPhaseMassDensity_dGlobalCompFraction,
                            phaseViscosity,
                            dPhaseViscosity_dPressure,
                            dPhaseViscosity_dTemperature,
                            dPhaseViscosity_dGlobalCompFraction,
                            phaseCompFraction,
                            dPhaseCompFraction_dPressure,
                            dPhaseCompFraction_dTemperature,
                            dPhaseCompFraction_dGlobalCompFraction,
                            totalDensity,
                            dTotalDensity_dPressure,
                            dTotalDensity_dTemperature,
                            dTotalDensity_dGlobalCompFraction ),
    m_fluid( fluid ),
    m_phaseTypes( phaseTypes )
  {}

  /// Default copy constructor
  MultiFluidPVTPackageWrapperUpdate( MultiFluidPVTPackageWrapperUpdate const & ) = default;

  /// Default move constructor
  MultiFluidPVTPackageWrapperUpdate( MultiFluidPVTPackageWrapperUpdate && ) = default;

  /// Deleted copy assignment operator
  MultiFluidPVTPackageWrapperUpdate & operator=( MultiFluidPVTPackageWrapperUpdate const & ) = delete;

  /// Deleted move assignment operator
  MultiFluidPVTPackageWrapperUpdate & operator=( MultiFluidPVTPackageWrapperUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void compute( real64 const pressure,
                        real64 const temperature,
                        arraySlice1d< real64 const > const & composition,
                        arraySlice1d< real64 > const & phaseFraction,
                        arraySlice1d< real64 > const & phaseDensity,
                        arraySlice1d< real64 > const & phaseMassDensity,
                        arraySlice1d< real64 > const & phaseViscosity,
                        arraySlice2d< real64 > const & phaseCompFraction,
                        real64 & totalDensity ) const override;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void compute( real64 const pressure,
                        real64 const temperature,
                        arraySlice1d< real64 const > const & composition,
                        arraySlice1d< real64 > const & phaseFraction,
                        arraySlice1d< real64 > const & dPhaseFraction_dPressure,
                        arraySlice1d< real64 > const & dPhaseFraction_dTemperature,
                        arraySlice2d< real64 > const & dPhaseFraction_dGlobalCompFraction,
                        arraySlice1d< real64 > const & phaseDensity,
                        arraySlice1d< real64 > const & dPhaseDensity_dPressure,
                        arraySlice1d< real64 > const & dPhaseDensity_dTemperature,
                        arraySlice2d< real64 > const & dPhaseDensity_dGlobalCompFraction,
                        arraySlice1d< real64 > const & phaseMassDensity,
                        arraySlice1d< real64 > const & dPhaseMassDensity_dPressure,
                        arraySlice1d< real64 > const & dPhaseMassDensity_dTemperature,
                        arraySlice2d< real64 > const & dPhaseMassDensity_dGlobalCompFraction,
                        arraySlice1d< real64 > const & phaseViscosity,
                        arraySlice1d< real64 > const & dPhaseViscosity_dPressure,
                        arraySlice1d< real64 > const & dPhaseViscosity_dTemperature,
                        arraySlice2d< real64 > const & dPhaseViscosity_dGlobalCompFraction,
                        arraySlice2d< real64 > const & phaseCompFraction,
                        arraySlice2d< real64 > const & dPhaseCompFraction_dPressure,
                        arraySlice2d< real64 > const & dPhaseCompFraction_dTemperature,
                        arraySlice3d< real64 > const & dPhaseCompFraction_dGlobalCompFraction,
                        real64 & totalDensity,
                        real64 & dTotalDensity_dPressure,
                        real64 & dTotalDensity_dTemperature,
                        arraySlice1d< real64 > const & dTotalDensity_dGlobalCompFraction ) const override;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const pressure,
                       real64 const temperature,
                       arraySlice1d< real64 const > const & composition ) const override
  {
    compute( pressure,
             temperature,
             composition,
             m_phaseFraction[k][q],
             m_dPhaseFraction_dPressure[k][q],
             m_dPhaseFraction_dTemperature[k][q],
             m_dPhaseFraction_dGlobalCompFraction[k][q],
             m_phaseDensity[k][q],
             m_dPhaseDensity_dPressure[k][q],
             m_dPhaseDensity_dTemperature[k][q],
             m_dPhaseDensity_dGlobalCompFraction[k][q],
             m_phaseMassDensity[k][q],
             m_dPhaseMassDensity_dPressure[k][q],
             m_dPhaseMassDensity_dTemperature[k][q],
             m_dPhaseMassDensity_dGlobalCompFraction[k][q],
             m_phaseViscosity[k][q],
             m_dPhaseViscosity_dPressure[k][q],
             m_dPhaseViscosity_dTemperature[k][q],
             m_dPhaseViscosity_dGlobalCompFraction[k][q],
             m_phaseCompFraction[k][q],
             m_dPhaseCompFraction_dPressure[k][q],
             m_dPhaseCompFraction_dTemperature[k][q],
             m_dPhaseCompFraction_dGlobalCompFraction[k][q],
             m_totalDensity[k][q],
             m_dTotalDensity_dPressure[k][q],
             m_dTotalDensity_dTemperature[k][q],
             m_dTotalDensity_dGlobalCompFraction[k][q] );
  }

private:

  pvt::MultiphaseSystem & m_fluid;

  arrayView1d< pvt::PHASE_TYPE > m_phaseTypes;

};

class MultiFluidPVTPackageWrapper : public MultiFluidBase
{
public:

  MultiFluidPVTPackageWrapper( string const & name, Group * const parent );

  virtual ~MultiFluidPVTPackageWrapper() override;

  virtual std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                Group * const parent ) const override;

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = MultiFluidPVTPackageWrapperUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper()
  {
    return KernelWrapper( *m_fluid,
                          m_phaseTypes,
                          m_componentMolarWeight,
                          m_useMass,
                          m_phaseFraction,
                          m_dPhaseFraction_dPressure,
                          m_dPhaseFraction_dTemperature,
                          m_dPhaseFraction_dGlobalCompFraction,
                          m_phaseDensity,
                          m_dPhaseDensity_dPressure,
                          m_dPhaseDensity_dTemperature,
                          m_dPhaseDensity_dGlobalCompFraction,
                          m_phaseMassDensity,
                          m_dPhaseMassDensity_dPressure,
                          m_dPhaseMassDensity_dTemperature,
                          m_dPhaseMassDensity_dGlobalCompFraction,
                          m_phaseViscosity,
                          m_dPhaseViscosity_dPressure,
                          m_dPhaseViscosity_dTemperature,
                          m_dPhaseViscosity_dGlobalCompFraction,
                          m_phaseCompFraction,
                          m_dPhaseCompFraction_dPressure,
                          m_dPhaseCompFraction_dTemperature,
                          m_dPhaseCompFraction_dGlobalCompFraction,
                          m_totalDensity,
                          m_dTotalDensity_dPressure,
                          m_dTotalDensity_dTemperature,
                          m_dTotalDensity_dGlobalCompFraction );
  }

protected:

  virtual void postProcessInput() override;

  virtual void initializePostSubGroups() override;

  /// function that populates m_fluid ptr; to be overriden by derived classes
  virtual void createFluid() = 0;

  /// PVTPackage fluid object
  std::unique_ptr< pvt::MultiphaseSystem > m_fluid;

  /// PVTPackage phase labels
  array1d< pvt::PHASE_TYPE > m_phaseTypes;
};

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void MultiFluidPVTPackageWrapperUpdate::compute( real64 const pressure,
                                                 real64 const temperature,
                                                 arraySlice1d< real64 const > const & composition,
                                                 arraySlice1d< real64 > const & phaseFrac,
                                                 arraySlice1d< real64 > const & phaseDens,
                                                 arraySlice1d< real64 > const & phaseMassDens,
                                                 arraySlice1d< real64 > const & phaseVisc,
                                                 arraySlice2d< real64 > const & phaseCompFrac,
                                                 real64 & totalDens ) const
{
#if defined(__CUDA_ARCH__)
  GEOSX_ERROR( "This function cannot be used on GPU" );
#else
  localIndex const NC = m_componentMolarWeight.size();
  localIndex const NP = m_phaseTypes.size();

  // 1. Convert input mass fractions to mole fractions and keep derivatives
  std::vector< double > compMoleFrac( NC );

  if( m_useMass )
  {
    real64 totalMolality = 0.0;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      compMoleFrac[ic] = composition[ic] / m_componentMolarWeight[ic]; // this is molality (units of mole/mass)
      totalMolality += compMoleFrac[ic];
    }

    real64 const totalMolalityInv = 1.0 / totalMolality;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      compMoleFrac[ic] *= totalMolalityInv;
    }
  }
  else
  {
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      compMoleFrac[ic] = composition[ic];
    }
  }

  // 2. Trigger PVTPackage compute and get back phase split
  m_fluid.Update( pressure, temperature, compMoleFrac );

  GEOSX_WARNING_IF( !m_fluid.hasSucceeded(),
                    "Phase equilibrium calculations not converged" );

  pvt::MultiphaseSystemProperties const & props = m_fluid.getMultiphaseSystemProperties();

  // 3. Extract phase split and phase properties from PVTPackage
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    pvt::PHASE_TYPE const & phaseType = m_phaseTypes[ip];

    auto const & frac = props.getPhaseMoleFraction( phaseType );
    auto const & comp = props.getMoleComposition( phaseType );
    auto const & dens = m_useMass ? props.getMassDensity( phaseType ) : props.getMoleDensity( phaseType );
    auto const & massDens = props.getMassDensity( phaseType );

    phaseFrac[ip] = frac.value;
    phaseDens[ip] = dens.value;
    phaseMassDens[ip] = massDens.value;
    phaseVisc[ip] = 0.001;   // TODO
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      phaseCompFrac[ip][jc] = comp.value[jc];
    }
  }

  // 4. if mass variables used instead of molar, perform the conversion
  if( m_useMass )
  {
    // 4.1. Convert phase fractions (requires two passes)
    real64 totalMass{};

    // 4.1.1. Compute mass of each phase and total mass (on a 1-mole basis)
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      auto const & phaseMW = props.getMolecularWeight( m_phaseTypes[ip] );
      phaseFrac[ip] *= phaseMW.value;
      totalMass += phaseFrac[ip];
    }

    // 4.1.2. Normalize to get mass fractions
    real64 const totalMassInv = 1.0 / totalMass;
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      phaseFrac[ip] *= totalMassInv;
    }

    // 4.2. Convert phase compositions
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      auto const & phaseMW = props.getMolecularWeight( m_phaseTypes[ip] );
      real64 const phaseMWInv = 1.0 / phaseMW.value;

      for( localIndex ic = 0; ic < NC; ++ic )
      {
        phaseCompFrac[ip][ic] = phaseCompFrac[ip][ic] * m_componentMolarWeight[ic] * phaseMWInv;
      }
    }
  }

  // 5. Compute total fluid mass/molar density
  {
    totalDens = 0.0;

    // 5.1. Sum mass/molar fraction/density ratio over all phases to get the inverse of density
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      totalDens += phaseFrac[ip] / phaseDens[ip];
    }

    // 5.2. Invert the previous quantity to get actual density
    totalDens = 1.0 / totalDens;
  }
#endif
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void MultiFluidPVTPackageWrapperUpdate::compute( real64 pressure,
                                                 real64 temperature,
                                                 arraySlice1d< real64 const > const & composition,
                                                 arraySlice1d< real64 > const & phaseFraction,
                                                 arraySlice1d< real64 > const & dPhaseFraction_dPressure,
                                                 arraySlice1d< real64 > const & dPhaseFraction_dTemperature,
                                                 arraySlice2d< real64 > const & dPhaseFraction_dGlobalCompFraction,
                                                 arraySlice1d< real64 > const & phaseDensity,
                                                 arraySlice1d< real64 > const & dPhaseDensity_dPressure,
                                                 arraySlice1d< real64 > const & dPhaseDensity_dTemperature,
                                                 arraySlice2d< real64 > const & dPhaseDensity_dGlobalCompFraction,
                                                 arraySlice1d< real64 > const & phaseMassDensity,
                                                 arraySlice1d< real64 > const & dPhaseMassDensity_dPressure,
                                                 arraySlice1d< real64 > const & dPhaseMassDensity_dTemperature,
                                                 arraySlice2d< real64 > const & dPhaseMassDensity_dGlobalCompFraction,
                                                 arraySlice1d< real64 > const & phaseViscosity,
                                                 arraySlice1d< real64 > const & dPhaseViscosity_dPressure,
                                                 arraySlice1d< real64 > const & dPhaseViscosity_dTemperature,
                                                 arraySlice2d< real64 > const & dPhaseViscosity_dGlobalCompFraction,
                                                 arraySlice2d< real64 > const & phaseCompFraction,
                                                 arraySlice2d< real64 > const & dPhaseCompFraction_dPressure,
                                                 arraySlice2d< real64 > const & dPhaseCompFraction_dTemperature,
                                                 arraySlice3d< real64 > const & dPhaseCompFraction_dGlobalCompFraction,
                                                 real64 & totalDensity,
                                                 real64 & dTotalDensity_dPressure,
                                                 real64 & dTotalDensity_dTemperature,
                                                 arraySlice1d< real64 > const & dTotalDensity_dGlobalCompFraction ) const
{
#if defined(__CUDA_ARCH__)
  GEOSX_ERROR( "This function cannot be used on GPU" );
#else
// 0. make shortcut structs to avoid long names (TODO maybe remove)
  CompositionalVarContainer< 1 > phaseFrac{
    phaseFraction,
    dPhaseFraction_dPressure,
    dPhaseFraction_dTemperature,
    dPhaseFraction_dGlobalCompFraction
  };

  CompositionalVarContainer< 1 > phaseDens{
    phaseDensity,
    dPhaseDensity_dPressure,
    dPhaseDensity_dTemperature,
    dPhaseDensity_dGlobalCompFraction
  };

  CompositionalVarContainer< 1 > phaseMassDens {
    phaseMassDensity,
    dPhaseMassDensity_dPressure,
    dPhaseMassDensity_dTemperature,
    dPhaseMassDensity_dGlobalCompFraction
  };

  CompositionalVarContainer< 1 > phaseVisc {
    phaseViscosity,
    dPhaseViscosity_dPressure,
    dPhaseViscosity_dTemperature,
    dPhaseViscosity_dGlobalCompFraction
  };

  CompositionalVarContainer< 2 > phaseCompFrac{
    phaseCompFraction,
    dPhaseCompFraction_dPressure,
    dPhaseCompFraction_dTemperature,
    dPhaseCompFraction_dGlobalCompFraction
  };

  CompositionalVarContainer< 0 > totalDens{
    totalDensity,
    dTotalDensity_dPressure,
    dTotalDensity_dTemperature,
    dTotalDensity_dGlobalCompFraction
  };

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex const NC = numComponents();
  localIndex const NP = numPhases();

  // 1. Convert input mass fractions to mole fractions and keep derivatives

  std::vector< double > compMoleFrac( NC );
  stackArray2d< real64, maxNumComp * maxNumComp > dCompMoleFrac_dCompMassFrac( NC, NC );

  if( m_useMass )
  {
    dCompMoleFrac_dCompMassFrac.resize( NC, NC );
    dCompMoleFrac_dCompMassFrac.setValues< serialPolicy >( 0.0 );

    real64 totalMolality = 0.0;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      real64 const mwInv = 1.0 / m_componentMolarWeight[ic];
      compMoleFrac[ic] = composition[ic] * mwInv; // this is molality (units of mole/mass)
      dCompMoleFrac_dCompMassFrac[ic][ic] = mwInv;
      totalMolality += compMoleFrac[ic];
    }

    real64 const totalMolalityInv = 1.0 / totalMolality;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      compMoleFrac[ic] *= totalMolalityInv;

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dCompMoleFrac_dCompMassFrac[ic][jc] -= compMoleFrac[ic] / m_componentMolarWeight[jc];
        dCompMoleFrac_dCompMassFrac[ic][jc] *= totalMolalityInv;
      }
    }
  }
  else
  {
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      compMoleFrac[ic] = composition[ic];
    }
  }

  // 2. Trigger PVTPackage compute and get back phase split
  m_fluid.Update( pressure, temperature, compMoleFrac );

  GEOSX_WARNING_IF( !m_fluid.hasSucceeded(),
                    "Phase equilibrium calculations not converged" );

  pvt::MultiphaseSystemProperties const & props = m_fluid.getMultiphaseSystemProperties();

  // 3. Extract phase split, phase properties and derivatives from PVTPackage
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    pvt::PHASE_TYPE const & phaseType = m_phaseTypes[ip];

    auto const & frac = props.getPhaseMoleFraction( phaseType );
    auto const & comp = props.getMoleComposition( phaseType );
    auto const & dens = m_useMass ? props.getMassDensity( phaseType ) : props.getMoleDensity( phaseType );
    auto const & massDens = props.getMassDensity( phaseType );

    phaseFrac.value[ip] = frac.value;
    phaseFrac.dPres[ip] = frac.dP;
    phaseFrac.dTemp[ip] = frac.dT;

    phaseDens.value[ip] = dens.value;
    phaseDens.dPres[ip] = dens.dP;
    phaseDens.dTemp[ip] = dens.dT;

    phaseMassDens.value[ip] = massDens.value;
    phaseMassDens.dPres[ip] = massDens.dP;
    phaseMassDens.dTemp[ip] = massDens.dT;

    // TODO
    phaseVisc.value[ip] = 0.001;
    phaseVisc.dPres[ip] = 0.0;
    phaseVisc.dTemp[ip] = 0.0;

    for( localIndex jc = 0; jc < NC; ++jc )
    {
      phaseFrac.dComp[ip][jc]     = frac.dz[jc];
      phaseDens.dComp[ip][jc]     = dens.dz[jc];
      phaseMassDens.dComp[ip][ip] = massDens.dz[jc];
      phaseVisc.dComp[ip][jc]     = 0.0; // TODO

      phaseCompFrac.value[ip][jc] = comp.value[jc];
      phaseCompFrac.dPres[ip][jc] = comp.dP[jc];
      phaseCompFrac.dTemp[ip][jc] = comp.dT[jc];

      for( localIndex ic = 0; ic < NC; ++ic )
      {
        phaseCompFrac.dComp[ip][ic][jc] = comp.dz[ic][jc];
      }
    }
  }

  // 4. if mass variables used instead of molar, perform the conversion
  if( m_useMass )
  {
    // 4.1. Convert phase fractions (requires two passes)
    real64 totalMass{};
    real64 dTotalMass_dP{};
    real64 dTotalMass_dT{};
    real64 dTotalMass_dC[maxNumComp]{};

    // 4.1.1. Compute mass of each phase and total mass (on a 1-mole basis)
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      auto const & phaseMW = props.getMolecularWeight( m_phaseTypes[ip] );
      real64 const nu = phaseFrac.value[ip];

      phaseFrac.value[ip] *= phaseMW.value;
      phaseFrac.dPres[ip] = phaseFrac.dPres[ip] * phaseMW.value + nu * phaseMW.dP;
      phaseFrac.dTemp[ip] = phaseFrac.dTemp[ip] * phaseMW.value + nu * phaseMW.dT;

      totalMass += phaseFrac.value[ip];
      dTotalMass_dP += phaseFrac.dPres[ip];
      dTotalMass_dT += phaseFrac.dTemp[ip];

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        phaseFrac.dComp[ip][jc] = phaseFrac.dComp[ip][jc] * phaseMW.value + nu * phaseMW.dz[jc];
        dTotalMass_dC[jc] += phaseFrac.dComp[ip][jc];
      }
    }

    // 4.1.2. Normalize to get mass fractions
    real64 const totalMassInv = 1.0 / totalMass;
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      phaseFrac.value[ip] *= totalMassInv;
      phaseFrac.dPres[ip] = ( phaseFrac.dPres[ip] - phaseFrac.value[ip] * dTotalMass_dP ) * totalMassInv;
      phaseFrac.dTemp[ip] = ( phaseFrac.dTemp[ip] - phaseFrac.value[ip] * dTotalMass_dT ) * totalMassInv;

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        phaseFrac.dComp[ip][jc] = ( phaseFrac.dComp[ip][jc] - phaseFrac.value[ip] * dTotalMass_dC[jc] ) * totalMassInv;
      }
    }

    // 4.2. Convert phase compositions
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      auto const & phaseMW = props.getMolecularWeight( m_phaseTypes[ip] );
      real64 const phaseMWInv = 1.0 / phaseMW.value;

      for( localIndex ic = 0; ic < NC; ++ic )
      {

        real64 const compMW = m_componentMolarWeight[ic];

        phaseCompFrac.value[ip][ic] = phaseCompFrac.value[ip][ic] * compMW * phaseMWInv;
        phaseCompFrac.dPres[ip][ic] =
          ( phaseCompFrac.dPres[ip][ic] * compMW - phaseCompFrac.value[ip][ic] * phaseMW.dP ) * phaseMWInv;
        phaseCompFrac.dTemp[ip][ic] =
          ( phaseCompFrac.dTemp[ip][ic] * compMW - phaseCompFrac.value[ip][ic] * phaseMW.dT ) * phaseMWInv;

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          phaseCompFrac.dComp[ip][ic][jc] =
            ( phaseCompFrac.dComp[ip][ic][jc] * compMW - phaseCompFrac.value[ip][ic] * phaseMW.dz[jc] ) * phaseMWInv;
        }
      }
    }

    // 4.3. Update derivatives w.r.t. mole fractions to derivatives w.r.t mass fractions
    array1d< real64 > work( NC );
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      applyChainRuleInPlace( NC, dCompMoleFrac_dCompMassFrac, phaseFrac.dComp[ip], work );
      applyChainRuleInPlace( NC, dCompMoleFrac_dCompMassFrac, phaseDens.dComp[ip], work );

      for( localIndex ic = 0; ic < NC; ++ic )
      {
        applyChainRuleInPlace( NC, dCompMoleFrac_dCompMassFrac, phaseCompFrac.dComp[ip][ic], work );
      }
    }
  }

  // 5. Compute total fluid mass/molar density and derivatives
  {
    totalDens.value = 0.0;
    totalDens.dPres = 0.0;
    totalDens.dTemp = 0.0;
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      totalDens.dComp[jc] = 0.0;
    }

    // 5.1. Sum mass/molar fraction/density ratio over all phases to get the inverse of density
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      real64 const densInv = 1.0 / phaseDens.value[ip];
      real64 const value = phaseFrac.value[ip] * densInv;

      totalDens.value += value;
      totalDens.dPres += ( phaseFrac.dPres[ip] - value * phaseDens.dPres[ip] ) * densInv;
      totalDens.dTemp += ( phaseFrac.dTemp[ip] - value * phaseDens.dTemp[ip] ) * densInv;

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        totalDens.dComp[jc] += ( phaseFrac.dComp[ip][jc] - value * phaseDens.dComp[ip][jc] ) * densInv;
      }
    }

    // 5.2. Invert the previous quantity to get actual density
    totalDens.value = 1.0 / totalDens.value;
    real64 const minusDens2 = -totalDens.value * totalDens.value;
    totalDens.dPres *= minusDens2;
    totalDens.dTemp *= minusDens2;

    for( localIndex jc = 0; jc < NC; ++jc )
    {
      totalDens.dComp[jc] *= minusDens2;
    }
  }
#endif
}

} //namespace constitutive

} //namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDPVTPACKAGEWRAPPER_HPP_
