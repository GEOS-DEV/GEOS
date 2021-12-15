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
 * @file BlackOilFluid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_BLACKOILFLUID_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_BLACKOILFLUID_HPP_

#include "constitutive/fluid/BlackOilFluidBase.hpp"
#include "constitutive/fluid/PVTOData.hpp"
#include "math/interpolation/Interpolation.hpp"

namespace geosx
{

namespace constitutive
{

class BlackOilFluid : public BlackOilFluidBase
{
public:

  static constexpr real64 minForPhasePresence = 1e-10;

  static constexpr integer NC_BO = 3;
  static constexpr integer NP_BO = 3;
  static constexpr integer HNC_BO = NC_BO - 1;

  BlackOilFluid( string const & name, Group * const parent );

  static string catalogName() { return "BlackOilFluid"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /**
   * @brief Kernel wrapper class for BlackOilFluid
   *        This kernel can be called on the GPU
   */
  class KernelWrapper final : public BlackOilFluidBase::KernelWrapper
  {
public:

    GEOSX_HOST_DEVICE
    virtual void compute( real64 pressure,
                          real64 temperature,
                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
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
                          PhaseComp::SliceType const phaseCompFraction,
                          FluidProp::SliceType const totalDensity ) const override;

    GEOSX_HOST_DEVICE
    virtual void update( localIndex k,
                         localIndex q,
                         real64 pressure,
                         real64 temperature,
                         arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const override;

private:

    friend class BlackOilFluid;

    KernelWrapper( PVTOData const & PVTO,
                   arrayView1d< integer const > phaseTypes,
                   arrayView1d< integer const > phaseOrder,
                   arrayView1d< integer const > hydrocarbonPhaseOrder,
                   arrayView1d< real64 const > surfacePhaseMassDensity,
                   arrayView1d< TableFunction::KernelWrapper const > formationVolFactorTables,
                   arrayView1d< TableFunction::KernelWrapper const > viscosityTables,
                   BlackOilFluidBase::WaterParams const waterParams,
                   arrayView1d< real64 const > componentMolarWeight,
                   bool const useMass,
                   PhaseProp::ViewType phaseFraction,
                   PhaseProp::ViewType phaseDensity,
                   PhaseProp::ViewType phaseMassDensity,
                   PhaseProp::ViewType phaseViscosity,
                   PhaseComp::ViewType phaseCompFraction,
                   FluidProp::ViewType totalDensity );

    GEOSX_HOST_DEVICE
    void computeDensitiesViscosities( real64 const pressure,
                                      bool const needDerivs,
                                      real64 const composition[],
                                      arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & phaseFrac,
                                      PhaseProp::SliceType const & phaseDens,
                                      PhaseProp::SliceType const & phaseMassDens,
                                      PhaseProp::SliceType const & phaseVisc,
                                      real64 phaseMW[NP_BO],
                                      real64 dphaseMW_dPres[],
                                      real64 dphaseMW_dGlobalCompFrac[] ) const;

    GEOSX_HOST_DEVICE
    void computeEquilibrium( real64 const pressure,
                             bool const needDerivs,
                             real64 const composition[],
                             PhaseProp::SliceType const & phaseFraction,
                             PhaseComp::SliceType const & phaseCompFraction ) const;

    // TODO: move this function elsewhere so we can use it in other models
    GEOSX_HOST_DEVICE
    void convertToMolar( arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                         real64 compMoleFrac[],
                         real64 dCompMoleFrac_dCompMassFrac[NC_BO][NC_BO] ) const;

    GEOSX_HOST_DEVICE
    void convertToMass( real64 const dCompMoleFrac_dCompMassFrac[NC_BO][NC_BO],
                        real64 const phaseMW[NP_BO],
                        real64 const dPhaseMW_dPressure[NP_BO],
                        real64 const dPhaseMW_dGlobalCompFraction[NP_BO *NC_BO],
                        PhaseProp::SliceType const & phaseFrac,
                        PhaseComp::SliceType const & phaseCompFrac,
                        arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseDens_dGlobalCompFrac,
                        arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseVisc_dGlobalCompFrac ) const;

    GEOSX_HOST_DEVICE
    void computeRs( real64 const presBub,
                    real64 & Rs,
                    real64 & dRs_dPres ) const;

    GEOSX_HOST_DEVICE
    void computeSaturatedBoViscosity( real64 const Rs,
                                      real64 const dRs_dPres,
                                      real64 & Bo,
                                      real64 & dBo_dPres,
                                      real64 & visc,
                                      real64 & dVisc_dPres )  const;

    GEOSX_HOST_DEVICE
    void computeUndersaturatedBoViscosity( bool const needDerivs,
                                           integer const numHydrocarbonComp,
                                           real64 const pres,
                                           real64 const Rs,
                                           real64 const dRs_dComp[],
                                           real64 & Bo,
                                           real64 & dBo_dPres,
                                           real64 dBo_dComp[],
                                           real64 & visc,
                                           real64 & dVisc_dPres,
                                           real64 dvisc_dComp[] ) const;

    GEOSX_HOST_DEVICE
    void computeUndersaturatedBoViscosity( real64 const Rs,
                                           real64 const pres,
                                           real64 & Bo,
                                           real64 & visc ) const;

    GEOSX_HOST_DEVICE
    void computeMassMoleDensity( bool const needDerivs,
                                 bool const useMass,
                                 integer const numHydrocarbonComp,
                                 real64 const Rs,
                                 real64 const dRs_dPres,
                                 real64 const dRs_dComp[],
                                 real64 const Bo,
                                 real64 const dBo_dPres,
                                 real64 const dBo_dComp[],
                                 real64 & dens,
                                 real64 & dDens_dPres,
                                 arraySlice1d< real64, multifluid::USD_PHASE_DC - 3 > const & dDens_dC ) const;

    /// Data needed to update the oil phase properties
    PVTOData::KernelWrapper m_PVTOView;
  };

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

protected:

  virtual void postProcessInput() override;

private:

  virtual void readInputDataFromTableFunctions() override;

  virtual void readInputDataFromPVTFiles() override;

  /**
   * @brief Read all the PVT table provided by the user in Eclipse format
   * @param[in] oilTable the oil table data read from file
   * @param[in] oilSurfaceMassDensity the oil phase surface mass density
   * @param[in] oilSurfaceMolecularWeight the oil phase surface molecular weight
   * @param[in] gasSurfaceMassDensity the oil phase surface mass density
   * @param[in] gasSurfaceMolecularWeight the oil phase surface molecular weight
   */
  void fillPVTOData( array1d< array1d< real64 > > const & oilTable,
                     real64 oilSurfaceMassDensity,
                     real64 oilSurfaceMolecularWeight,
                     real64 gasSurfaceMassDensity,
                     real64 gasSurfaceMolecularWeight );

  /**
   * @brief Form the tables:
   *    - m_undersaturatedPressure
   *    - m_undersaturatedBO
   *    - m_undersaturatedVisc
   *  using the input data
   */
  void createUndersaturatedProperties();

  /**
   * @brief Extrapolate the tables:
   *    - m_undersaturatedPressure
   *    - m_undersaturatedBO
   *    - m_undersaturatedVisc
   *  up to the value of m_maxRelativePressure (if the data is missing for a given branch)
   */
  void extendUndersaturatedProperties();

  /**
   * @brief Refine the undersaturated tables and copy the data into the final 2D arrays that will be used in the kernel
   */
  void refineUndersaturatedTables( integer numRefinedPresPoints );

  /**
   * @brief Check the monotonicity of the PVTO table values
   */
  void checkTableConsistency() const;

  /// The data needed to compute the oil phase properties
  PVTOData m_PVTO;

};

GEOSX_HOST_DEVICE
inline void
BlackOilFluid::KernelWrapper::
  compute( real64 pressure,
           real64 temperature,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
           arraySlice2d< real64, multifluid::USD_PHASE_COMP - 2 > const & phaseCompFraction,
           real64 & totalDens ) const
{
  GEOSX_UNUSED_VAR( pressure );
  GEOSX_UNUSED_VAR( temperature );
  GEOSX_UNUSED_VAR( composition );
  GEOSX_UNUSED_VAR( phaseFraction );
  GEOSX_UNUSED_VAR( phaseDensity );
  GEOSX_UNUSED_VAR( phaseMassDensity );
  GEOSX_UNUSED_VAR( phaseViscosity );
  GEOSX_UNUSED_VAR( phaseCompFraction );
  GEOSX_UNUSED_VAR( totalDens );

  GEOSX_ERROR( "BlackOilFluid: this compute function is not implemented" );
}

GEOSX_HOST_DEVICE
inline void
BlackOilFluid::KernelWrapper::
  compute( real64 pressure,
           real64 temperature,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
           PhaseProp::SliceType const phaseFraction,
           PhaseProp::SliceType const phaseDensity,
           PhaseProp::SliceType const phaseMassDensity,
           PhaseProp::SliceType const phaseViscosity,
           PhaseComp::SliceType const phaseCompFraction,
           FluidProp::SliceType const totalDensity ) const
{
  GEOSX_UNUSED_VAR( temperature );

  real64 compMoleFrac[NC_BO]{};
  real64 dCompMoleFrac_dCompMassFrac[NC_BO][NC_BO]{};
  real64 phaseMW[NP_BO]{};
  real64 dPhaseMW_dPressure[NP_BO]{};
  real64 dPhaseMW_dGlobalCompFraction[NP_BO*NC_BO]{};

  // 1. Convert to mass if necessary

  if( m_useMass )
  {
    convertToMolar( composition, compMoleFrac, dCompMoleFrac_dCompMassFrac );
  }
  else
  {
    for( integer ic = 0; ic < NC_BO; ++ic )
    {
      compMoleFrac[ic] = composition[ic];
    }
  }

  // 2. Compute phase fractions and phase component fractions

  computeEquilibrium( pressure,
                      true,
                      compMoleFrac,
                      phaseFraction,
                      phaseCompFraction );

  // 3. Compute phase densities and viscosities

  computeDensitiesViscosities( pressure,
                               true,
                               compMoleFrac,
                               phaseFraction.value,
                               phaseDensity,
                               phaseMassDensity,
                               phaseViscosity,
                               phaseMW,
                               dPhaseMW_dPressure,
                               dPhaseMW_dGlobalCompFraction );

  // 4. If mass variables used instead of molar, perform the conversion
  if( m_useMass )
  {
    convertToMass( dCompMoleFrac_dCompMassFrac,
                   phaseMW,
                   dPhaseMW_dPressure,
                   dPhaseMW_dGlobalCompFraction,
                   phaseFraction,
                   phaseCompFraction,
                   phaseDensity.dComp,
                   phaseViscosity.dComp );
  }

  // 5. Compute total fluid mass/molar density and derivatives
  totalDensity.value = 0.0;
  totalDensity.dPres = 0.0;
  for( integer ic = 0; ic < NC_BO; ++ic )
  {
    totalDensity.dComp[ic] = 0.0;
  }

  for( integer ip = 0; ip < NP_BO; ++ip )
  {
    if( phaseFraction.value[ip] <= 0. )
    {
      continue;
    }
    real64 const densInv = 1.0 / phaseDensity.value[ip];
    real64 const value = phaseFraction.value[ip] * densInv;

    totalDensity.value += value;
    totalDensity.dPres += ( phaseFraction.dPres[ip] - value * phaseDensity.dPres[ip] ) * densInv;
    for( integer ic = 0; ic < NC_BO; ++ic )
    {
      totalDensity.dComp[ic] += ( phaseFraction.dComp[ip][ic] - value * phaseDensity.dComp[ip][ic] ) * densInv;
    }
  }

  totalDensity.value = 1.0 / totalDensity.value;
  real64 const minusDens2 = -totalDensity.value * totalDensity.value;
  totalDensity.dPres *= minusDens2;
  for( integer ic = 0; ic < NC_BO; ++ic )
  {
    totalDensity.dComp[ic] *= minusDens2;
  }
}

GEOSX_HOST_DEVICE
inline void
BlackOilFluid::KernelWrapper::
  computeEquilibrium( real64 const pressure,
                      bool const needDerivs,
                      real64 const composition[],
                      PhaseProp::SliceType const & phaseFraction,
                      PhaseComp::SliceType const & phaseCompFraction ) const
{
  using PT = BlackOilFluid::PhaseType;

  integer const ipOil   = m_phaseOrder[PT::OIL];
  integer const ipGas   = m_phaseOrder[PT::GAS];
  integer const ipWater = m_phaseOrder[PT::WATER];

  integer const icOil   = ipOil;
  integer const icGas   = ipGas;
  integer const icWater = ipWater;

  real64 const zo = composition[icOil];
  real64 const zg = composition[icGas];
  real64 const zw = composition[icWater];

  // 1. Make everything zero first

  for( integer ip = 0; ip < NP_BO; ++ip )
  {
    phaseFraction.value[ip] = 0.;
    if( needDerivs )
    {
      phaseFraction.dPres[ip] = 0.0;
    }
    for( integer ic = 0; ic < NC_BO; ++ic )
    {
      phaseCompFraction.value[ip][ic] = 0.0;
      if( needDerivs )
      {
        phaseFraction.dComp[ip][ic] = 0.;
        phaseCompFraction.dPres[ip][ic] = 0.0;
        for( integer jc = 0; jc < NC_BO; ++jc )
        {
          phaseCompFraction.dComp[ip][ic][jc] = 0.0;
        }
      }
    }
  }

  // 2. Check feed first, and if only water is present (e.g., water inj), then skip

  if( zw >= 1. - minForPhasePresence )
  {
    phaseFraction.value[ipWater] = zw;
    if( needDerivs )
    {
      phaseFraction.dComp[ipWater][icWater] = 1.;
    }
    phaseCompFraction.value[ipWater][icWater] = 1.0;
    return;
  }

  // 3. Compute Rs (saturated) as a function of pressure

  // oil
  real64 const & oilSurfaceMoleDensity = m_PVTOView.m_surfaceMoleDensity[PT::OIL];
  real64 RsSat = 0.0;
  real64 dRsSat_dP = 0.0;
  computeRs( pressure, RsSat, dRsSat_dP );

  // gas
  real64 const & gasSurfaceMoleDensity = m_PVTOView.m_surfaceMoleDensity[PT::GAS];

  if( RsSat < minForPhasePresence )
  {
    RsSat = minForPhasePresence;
  }

  // 4. Use the saturated Rs and density to compute the gas phase fraction
  //    Note : we assume the oil component cannot enter the gas phase

  real64 const Kg = ( oilSurfaceMoleDensity + gasSurfaceMoleDensity * RsSat ) / ( RsSat * gasSurfaceMoleDensity );
  real64 const dKg_dP = -oilSurfaceMoleDensity / gasSurfaceMoleDensity * dRsSat_dP / ( RsSat * RsSat );
  real64 const gasPhaseFraction = zo / ( 1. - Kg ) + zg;
  real64 const dGasPhaseFraction_dP = zo / ( ( 1. - Kg ) * ( 1. - Kg ) ) *  dKg_dP;
  real64 const dGasPhaseFraction_dzo = 1. / ( 1. - Kg );
  real64 const dGasPhaseFraction_dzg = 1.;

  // 4. Update phase fraction and phase component fractions

  // 4.1 The gas phase is present
  if( ( gasPhaseFraction > 0 ) /*&& ( gasPhaseFraction < 1 )*/ )
  {

    // phase fractions
    phaseFraction.value[ipOil] = 1. - gasPhaseFraction - zw;
    phaseFraction.value[ipGas] = gasPhaseFraction;
    phaseFraction.value[ipWater] =  zw;

    if( needDerivs )
    {
      phaseFraction.dPres[ipOil] = -dGasPhaseFraction_dP;
      phaseFraction.dPres[ipGas] = dGasPhaseFraction_dP;
      phaseFraction.dComp[ipOil][icOil] = -dGasPhaseFraction_dzo;
      phaseFraction.dComp[ipOil][icGas] = -dGasPhaseFraction_dzg;
      phaseFraction.dComp[ipOil][icWater] = -1.;
      phaseFraction.dComp[ipGas][icOil] = dGasPhaseFraction_dzo;
      phaseFraction.dComp[ipGas][icGas] = dGasPhaseFraction_dzg;
      phaseFraction.dComp[ipWater][icWater] = 1.;
    }

    // oil
    real64 const tmp = ( oilSurfaceMoleDensity + gasSurfaceMoleDensity * RsSat );
    real64 const tmpOil = oilSurfaceMoleDensity / tmp;
    real64 const dTmpOil_dP = -oilSurfaceMoleDensity * gasSurfaceMoleDensity * dRsSat_dP / ( tmp * tmp );
    phaseCompFraction.value[ipOil][icOil] = tmpOil;
    phaseCompFraction.value[ipOil][icGas] = 1. - tmpOil;
    phaseCompFraction.value[ipOil][icWater] = 0.;

    if( needDerivs )
    {
      phaseCompFraction.dPres[ipOil][icOil] = dTmpOil_dP;
      phaseCompFraction.dPres[ipOil][icGas] = -dTmpOil_dP;
    }

    // gas
    real64 const tmpGas = gasSurfaceMoleDensity / ( gasSurfaceMoleDensity );
    phaseCompFraction.value[ipGas][icOil] = 1. - tmpGas;
    phaseCompFraction.value[ipGas][icGas] = tmpGas;
    phaseCompFraction.value[ipGas][icWater] = 0.;

    // water
    phaseCompFraction.value[ipWater][icOil] = 0;
    phaseCompFraction.value[ipWater][icGas] = 0;
    phaseCompFraction.value[ipWater][icWater] = 1.;

  }
  // 4.2 The gas phase is absent
  else
  {

    // phase fractions
    phaseFraction.value[ipOil] = 1. - zw;
    phaseFraction.value[ipGas] = 0.;
    phaseFraction.value[ipWater] =  zw;

    // oil
    phaseCompFraction.value[ipOil][icOil] = zo;
    phaseCompFraction.value[ipOil][icGas] = zg;
    phaseCompFraction.value[ipOil][icWater] = 0.;

    // gas
    phaseCompFraction.value[ipWater][icOil] = 0;
    phaseCompFraction.value[ipWater][icGas] = 0;
    phaseCompFraction.value[ipWater][icWater] = 1.;

    if( needDerivs )
    {
      phaseFraction.dComp[ipOil][icWater] = -1.;
      phaseFraction.dComp[ipWater][icWater] = 1.;
      phaseCompFraction.dComp[ipOil][icOil][icOil] = 1.;
      phaseCompFraction.dComp[ipOil][icGas][icGas] = 1.;
    }
  }

}

GEOSX_HOST_DEVICE
inline void
BlackOilFluid::KernelWrapper::
  computeDensitiesViscosities( real64 const pressure,
                               bool const needDerivs,
                               real64 const composition[],
                               arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & phaseFrac,
                               PhaseProp::SliceType const & phaseDens,
                               PhaseProp::SliceType const & phaseMassDens,
                               PhaseProp::SliceType const & phaseVisc,
                               real64 phaseMW[NP_BO],
                               real64 dPhaseMW_dPres[],
                               real64 dPhaseMW_dGlobalCompFrac[] ) const
{
  using PT = BlackOilFluid::PhaseType;

  integer const ipOil = m_phaseOrder[PT::OIL];
  integer const ipGas = m_phaseOrder[PT::GAS];
  integer const ipWater = m_phaseOrder[PT::WATER];

  integer const icOil = ipOil;
  integer const icGas = ipGas;

  bool const isWater = (ipWater >= 0 && phaseFrac[ipWater] > 0);
  bool const isGas = (ipGas >= 0 && phaseFrac[ipGas] > 0);
  bool const isOil = (ipOil >= 0 && phaseFrac[ipOil] > 0);

  for( integer ip = 0; ip < NP_BO; ++ip )
  {
    phaseMassDens.value[ip] = 0.;
    phaseDens.value[ip]     = 0.;
    phaseVisc.value[ip]     = 0.;
    phaseDens.dPres[ip]     = 0.;
    phaseMassDens.dPres[ip] = 0.;
    phaseVisc.dPres[ip]     = 0.;
    for( integer ic = 0; ic < NC_BO; ++ic )
    {
      phaseDens.dComp[ip][ic] = 0.;
      phaseMassDens.dComp[ip][ic] = 0.;
      phaseVisc.dComp[ip][ic]     = 0.;
    }
  }

  // 1. Gas phase: look up in the formation vol factor tables

  if( isGas )
  {
    // interpolate in the table to get the phase formation vol factor and its derivative wrt pressure
    real64 derivative = 0.0;
    real64 const fvf = m_formationVolFactorTables[0].compute( &pressure, &derivative );

    // we are ready to update the densities
    real64 const fvfInv = 1.0 / fvf;

    phaseMassDens.value[ipGas] = m_surfacePhaseMassDensity[ipGas] * fvfInv;
    real64 const mult = m_useMass ? 1.0 : 1.0 / m_componentMolarWeight[ipGas];
    phaseDens.value[ipGas] = phaseMassDens.value[ipGas] * mult;
    phaseMW[ipGas] = m_componentMolarWeight[ipGas];

    if( needDerivs )
    {
      phaseMassDens.dPres[ipGas] = -derivative * phaseMassDens.value[ipGas] * fvfInv;
      phaseDens.dPres[ipGas] = phaseMassDens.dPres[ipGas] * mult;
      dPhaseMW_dPres[ipGas] = 0.;
      for( integer ic = 0; ic < NC_BO; ic++ )
      {
        phaseMassDens.dComp[ipGas][ic] = 0.;
        phaseDens.dComp[ipGas][ic]     = 0.;
        dPhaseMW_dGlobalCompFrac[ipGas*NC_BO+ic]  = 0.;
      }
    }

    phaseVisc.value[ipGas] = m_viscosityTables[0].compute( &pressure, &(phaseVisc.dPres)[ipGas] );
  }

  // 2. Water phase: use the constant formation volume factor and compressibility provided by the user

  if( isWater )
  {
    // if water is present
    real64 const expCompDeltaPres = std::exp( -m_waterParams.compressibility * ( pressure - m_waterParams.referencePressure ) );
    real64 const dExpCompDeltaPres_dPres = -m_waterParams.compressibility * expCompDeltaPres;
    real64 const denom = m_waterParams.formationVolFactor * expCompDeltaPres;
    real64 const dDenom_dPres = m_waterParams.formationVolFactor * dExpCompDeltaPres_dPres;
    real64 const denomInv = 1.0 / denom;
    phaseMassDens.value[ipWater] = m_surfacePhaseMassDensity[ipWater] * denomInv;
    real64 const mult = m_useMass ? 1.0 : 1.0 / m_componentMolarWeight[ipWater];
    phaseDens.value[ipWater] = phaseMassDens.value[ipWater] * mult;
    phaseVisc.value[ipWater] = m_waterParams.viscosity;
    phaseMW[ipWater] = m_componentMolarWeight[ipWater];

    if( needDerivs )
    {
      phaseMassDens.dPres[ipWater] = -dDenom_dPres * phaseMassDens.value[ipWater] * denomInv;
      phaseDens.dPres[ipWater] = phaseMassDens.dPres[ipWater] * mult;
      dPhaseMW_dPres[ipWater] = 0.;
      for( integer ic = 0; ic < NC_BO; ic++ )
      {
        phaseMassDens.dComp[ipWater][ic] = 0.;
        phaseDens.dComp[ipWater][ic] = 0.0;
        dPhaseMW_dGlobalCompFrac[ipWater*NC_BO+ic] = 0.;
      }
    }

  }

  // 3. Oil phase: make the distinction between saturated and unsaturated conditions

  if( isOil )
  {

    real64 Rs = 0.0;
    real64 dRs_dC[2]{};
    real64 dRs_dP = 0.0;
    real64 Bo = 0.0;
    real64 dBo_dP = 0.0;
    real64 dBo_dC[2]{};
    real64 visc = 0.0;
    real64 dVisc_dP = 0.0;
    real64 dVisc_dC[2]{};

    // saturated conditions
    if( isGas )
    {

      // compute Rs as a function of pressure
      computeRs( pressure, Rs, dRs_dP );

      // compute saturated properties (Bo, viscosity) as a function of Rs
      computeSaturatedBoViscosity( Rs, dRs_dP, Bo, dBo_dP, visc, dVisc_dP );

    }
    // unsaturated conditions
    else
    {

      // compute Rs as a function of composition
      real64 const densRatio = m_PVTOView.m_surfaceMoleDensity[PT::OIL] / m_PVTOView.m_surfaceMoleDensity[PT::GAS];
      Rs = densRatio * composition[icGas] / composition[icOil];
      dRs_dC[icOil] = -densRatio * composition[icGas] / (composition[icOil] * composition[icOil]);
      dRs_dC[icGas] = densRatio  / composition[icOil];

      if( isZero( composition[icOil] ) )
      {
	std::cout << "composition = " << composition << std::endl;
      }
      
      // compute undersaturated properties (Bo, viscosity) by two-step interpolation in undersaturated tables
      // this part returns numerical derivatives
      computeUndersaturatedBoViscosity( needDerivs, HNC_BO, pressure, Rs, dRs_dC, Bo, dBo_dP, dBo_dC,
                                        visc, dVisc_dP, dVisc_dC );

    }

    // compute densities
    computeMassMoleDensity( needDerivs, true, HNC_BO, Rs, dRs_dP, dRs_dC, Bo, dBo_dP, dBo_dC,
                            phaseMassDens.value[ipOil], phaseMassDens.dPres[ipOil],
                            phaseMassDens.dComp[ipOil] );
    computeMassMoleDensity( needDerivs, false, HNC_BO, Rs, dRs_dP, dRs_dC, Bo, dBo_dP, dBo_dC,
                            phaseDens.value[ipOil], phaseDens.dPres[ipOil],
                            phaseDens.dComp[ipOil] );

    phaseMW[ipOil] = phaseMassDens.value[ipOil] / phaseDens.value[ipOil];
    real64 const tmp = 1. / ( phaseDens.value[ipOil] * phaseDens.value[ipOil] );
    dPhaseMW_dPres[ipOil] = ( phaseMassDens.dPres[ipOil] * phaseDens.value[ipOil] - phaseDens.dPres[ipOil] * phaseMassDens.value[ipOil] ) * tmp;
    for( integer ic = 0; ic < NC_BO; ++ic )
    {
      dPhaseMW_dGlobalCompFrac[ipOil*NC_BO+ic] = ( phaseMassDens.dComp[ipOil][ic] * phaseDens.value[ipOil]
                                                   - phaseDens.dComp[ipOil][ic] * phaseMassDens.value[ipOil] ) * tmp;
    }

    if( m_useMass )
    {
      phaseDens.value[ipOil] = phaseMassDens.value[ipOil];
      if( needDerivs )
      {
        phaseDens.dPres[ipOil] = phaseMassDens.dPres[ipOil];
        for( integer i = 0; i < phaseMassDens.dComp[ipOil].size(); ++i )
        {
          phaseDens.dComp[ipOil][i] = phaseMassDens.dComp[ipOil][i];
        }
      }
    }

    // copy viscosity into the final array
    // TODO: skip this step
    phaseVisc.value[ipOil] = visc;
    if( needDerivs )
    {
      phaseVisc.dPres[ipOil] = dVisc_dP;
      for( integer i = 0; i < HNC_BO; ++i )
      {
        phaseVisc.dComp[ipOil][i] = dVisc_dC[i];
      }
    }
  }
}

GEOSX_HOST_DEVICE
inline void
BlackOilFluid::KernelWrapper::
  convertToMolar( arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                  real64 compMoleFrac[],
                  real64 dCompMoleFrac_dCompMassFrac[NC_BO][NC_BO] ) const
{
  for( integer ic = 0; ic < NC_BO; ++ic )
  {
    for( integer jc = 0; jc < NC_BO; ++jc )
    {
      dCompMoleFrac_dCompMassFrac[ic][jc] = 0.0;
    }
  }

  real64 totalMolality = 0.0;
  for( integer ic = 0; ic < NC_BO; ++ic )
  {
    real64 const mwInv = 1.0 / m_componentMolarWeight[ic];
    compMoleFrac[ic] = composition[ic] * mwInv;
    dCompMoleFrac_dCompMassFrac[ic][ic] = mwInv;
    totalMolality += compMoleFrac[ic];
  }

  real64 const totalMolalityInv = 1.0 / totalMolality;
  for( integer ic = 0; ic < NC_BO; ++ic )
  {
    compMoleFrac[ic] *= totalMolalityInv;

    for( integer jc = 0; jc < NC_BO; ++jc )
    {
      dCompMoleFrac_dCompMassFrac[ic][jc] -= compMoleFrac[ic] / m_componentMolarWeight[jc];
      dCompMoleFrac_dCompMassFrac[ic][jc] *= totalMolalityInv;
    }
  }

}


GEOSX_HOST_DEVICE
inline void
BlackOilFluid::KernelWrapper::
  convertToMass( real64 const dCompMoleFrac_dCompMassFrac[NC_BO][NC_BO],
                 real64 const phaseMW[NP_BO],
                 real64 const dPhaseMW_dPres[NP_BO],
                 real64 const dPhaseMW_dGlobalCompFrac[NP_BO * NC_BO],
                 PhaseProp::SliceType const & phaseFrac,
                 PhaseComp::SliceType const & phaseCompFrac,
                 arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseDens_dGlobalCompFrac,
                 arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseVisc_dGlobalCompFrac ) const
{
  // 1. Convert phase fractions (requires two passes)

  real64 totalMass{};
  real64 dTotalMass_dP{};
  real64 dTotalMass_dC[NC_BO]{};

  // 1.1. Compute mass of each phase and total mass (on a 1-mole basis)

  for( integer ip = 0; ip < NP_BO; ++ip )
  {
    bool const phaseExists = (phaseFrac.value[ip] > 0);
    if( !phaseExists )
    {
      continue;
    }

    real64 const nu = phaseFrac.value[ip];

    phaseFrac.value[ip] *= phaseMW[ip];
    phaseFrac.dPres[ip] = phaseFrac.dPres[ip] * phaseMW[ip] + nu * dPhaseMW_dPres[ip];
    totalMass += phaseFrac.value[ip];
    dTotalMass_dP += phaseFrac.dPres[ip];

    for( integer jc = 0; jc < NC_BO; ++jc )
    {
      phaseFrac.dComp[ip][jc] = phaseFrac.dComp[ip][jc] * phaseMW[ip] + nu * dPhaseMW_dGlobalCompFrac[jc + ip*NC_BO];
      dTotalMass_dC[jc] += phaseFrac.dComp[ip][jc];
    }
  }

  // 1.2. Normalize to get mass fractions

  real64 const totalMassInv = 1.0 / totalMass;
  for( integer ip = 0; ip < NC_BO; ++ip )
  {
    bool const phaseExists = (phaseFrac.value[ip] > 0);
    if( !phaseExists )
    {
      continue;
    }

    phaseFrac.value[ip] *= totalMassInv;
    phaseFrac.dPres[ip] = ( phaseFrac.dPres[ip] - phaseFrac.value[ip] * dTotalMass_dP ) * totalMassInv;

    for( integer jc = 0; jc < NC_BO; ++jc )
    {
      phaseFrac.dComp[ip][jc] = ( phaseFrac.dComp[ip][jc] - phaseFrac.value[ip] * dTotalMass_dC[jc] ) * totalMassInv;
    }
  }

  // 2. Convert phase compositions

  for( integer ip = 0; ip < NC_BO; ++ip )
  {
    bool const phaseExists = (phaseFrac.value[ip] > 0);
    if( !phaseExists )
    {
      continue;
    }

    real64 const phaseMWInv = 1.0 / phaseMW[ip];

    for( integer ic = 0; ic < NC_BO; ++ic )
    {

      real64 const compMW = m_componentMolarWeight[ic];

      phaseCompFrac.value[ip][ic] = phaseCompFrac.value[ip][ic] * compMW * phaseMWInv;
      phaseCompFrac.dPres[ip][ic] = ( phaseCompFrac.dPres[ip][ic] * compMW - phaseCompFrac.value[ip][ic] * dPhaseMW_dPres[ip] ) * phaseMWInv;

      for( integer jc = 0; jc < NC_BO; ++jc )
      {
        phaseCompFrac.dComp[ip][ic][jc] =
          ( phaseCompFrac.dComp[ip][ic][jc] * compMW - phaseCompFrac.value[ip][ic] * dPhaseMW_dGlobalCompFrac[jc + ip*NC_BO] ) * phaseMWInv;
      }
    }
  }

  // 3. Update derivatives w.r.t. mole fractions to derivatives w.r.t mass fractions

  real64 work[NC_BO];
  for( integer ip = 0; ip < NC_BO; ++ip )
  {
    bool const phaseExists = (phaseFrac.value[ip] > 0);
    if( !phaseExists )
    {
      continue;
    }

    applyChainRuleInPlace( NC_BO, dCompMoleFrac_dCompMassFrac, phaseFrac.dComp[ip], work );
    applyChainRuleInPlace( NC_BO, dCompMoleFrac_dCompMassFrac, dPhaseDens_dGlobalCompFrac[ip], work );
    applyChainRuleInPlace( NC_BO, dCompMoleFrac_dCompMassFrac, dPhaseVisc_dGlobalCompFrac[ip], work );
    for( integer ic = 0; ic < NC_BO; ++ic )
    {
      applyChainRuleInPlace( NC_BO, dCompMoleFrac_dCompMassFrac, phaseCompFrac.dComp[ip][ic], work );
    }
  }
}

GEOSX_HOST_DEVICE
inline void
BlackOilFluid::KernelWrapper::
  computeRs( real64 const presBub,
             real64 & Rs,
             real64 & dRs_dPres ) const
{

  integer const idx = LvArray::sortedArrayManipulation::find( m_PVTOView.m_bubblePressure.begin(),
                                                              m_PVTOView.m_bubblePressure.size(),
                                                              presBub );
  integer const iUp  = LvArray::math::min( LvArray::math::max( idx, 1 ), LvArray::integerConversion< integer >( m_PVTOView.m_bubblePressure.size()-1 ) );
  integer const iLow = iUp-1;
  interpolation::linearInterpolation( presBub - m_PVTOView.m_bubblePressure[iLow], m_PVTOView.m_bubblePressure[iUp] - presBub,
                                      m_PVTOView.m_Rs[iLow], m_PVTOView.m_Rs[iUp],
                                      Rs, dRs_dPres );
}

GEOSX_HOST_DEVICE
inline void
BlackOilFluid::KernelWrapper::
  computeSaturatedBoViscosity( real64 const Rs,
                               real64 const dRs_dPres,
                               real64 & Bo,
                               real64 & dBo_dPres,
                               real64 & visc,
                               real64 & dVisc_dPres ) const
{
  arrayView1d< real64 const > const & RsVec = m_PVTOView.m_Rs;
  arrayView1d< real64 const > const & BoVec = m_PVTOView.m_saturatedBo;
  arrayView1d< real64 const > const & viscVec = m_PVTOView.m_saturatedViscosity;

  integer const idx = LvArray::sortedArrayManipulation::find( RsVec.begin(),
                                                              RsVec.size(),
                                                              Rs );
  integer const iUp  = LvArray::math::min( LvArray::math::max( idx, 1 ), LvArray::integerConversion< integer >( RsVec.size()-1 ) );
  integer const iLow = iUp-1;

  interpolation::linearInterpolation( Rs - RsVec[iLow], RsVec[iUp] - Rs,
                                      BoVec[iLow], BoVec[iUp],
                                      Bo, dBo_dPres );
  interpolation::linearInterpolation( Rs - RsVec[iLow], RsVec[iUp] - Rs,
                                      viscVec[iLow], viscVec[iUp],
                                      visc, dVisc_dPres );

  // chain rule
  dBo_dPres *= dRs_dPres;
  dVisc_dPres *= dRs_dPres;
}

GEOSX_HOST_DEVICE
inline void
BlackOilFluid::KernelWrapper::
  computeUndersaturatedBoViscosity( bool const needDerivs,
                                    integer const numHydrocarbonComp,
                                    real64 const P,
                                    real64 const Rs,
                                    real64 const dRs_dComp[],
                                    real64 & Bo,
                                    real64 & dBo_dPres,
                                    real64 dBo_dComp[],
                                    real64 & visc,
                                    real64 & dVisc_dPres,
                                    real64 dVisc_dComp[] ) const
{
  computeUndersaturatedBoViscosity( Rs, P, Bo, visc );

  if( needDerivs )
  {
    // numerical derivatives

    // 1. dPres
    real64 const eps = 1e-6;
    real64 const inv_eps = 1.0 / eps;
    real64 const P_eps = P + eps;
    real64 Bo_eps   = 0.0;
    real64 visc_eps = 0.0;
    computeUndersaturatedBoViscosity( Rs, P_eps, Bo_eps, visc_eps );
    dBo_dPres = (Bo_eps - Bo) * inv_eps;
    dVisc_dPres = (visc_eps - visc) * inv_eps;

    // 2. dRs
    real64 const Rs_eps = Rs + eps;
    computeUndersaturatedBoViscosity( Rs_eps, P, Bo_eps, visc_eps );
    real64 const dBo_dRs =   (Bo_eps - Bo) * inv_eps;
    real64 const dVisc_dRs = (visc_eps - visc) * inv_eps;

    // 3. chainrule to dComp
    for( integer i = 0; i < numHydrocarbonComp; ++i )
    {
      dBo_dComp[i] = dBo_dRs * dRs_dComp[i];
      dVisc_dComp[i] = dVisc_dRs * dRs_dComp[i];
    }
  }

}

GEOSX_HOST_DEVICE
inline void
BlackOilFluid::KernelWrapper::
  computeUndersaturatedBoViscosity( real64 const Rs,
                                    real64 const pres,
                                    real64 & Bo,
                                    real64 & visc ) const
{
  // Step 1: interpolate for presBub
  arrayView1d< real64 const > const RsVec = m_PVTOView.m_Rs;
  integer idx = LvArray::sortedArrayManipulation::find( RsVec.begin(),
                                                        RsVec.size(),
                                                        Rs );
  integer const iUp  = LvArray::math::min( LvArray::math::max( idx, 1 ), LvArray::integerConversion< integer >( RsVec.size()-1 ) );
  integer const iLow = iUp-1;

  real64 const presBub = interpolation::linearInterpolation( Rs - m_PVTOView.m_Rs[iLow], m_PVTOView.m_Rs[iUp] - Rs,
                                                             m_PVTOView.m_bubblePressure[iLow], m_PVTOView.m_bubblePressure[iUp] );
  real64 const deltaPres = pres - presBub;

  // Step 2: get indices in undersatured pressure table
  real64 const deltaRsUp = LvArray::math::abs( m_PVTOView.m_Rs[iUp] - Rs );
  real64 const deltaRsLow = LvArray::math::abs( Rs - m_PVTOView.m_Rs[iLow] );
  arraySlice1d< real64 const > const & presUp = m_PVTOView.m_undersaturatedPressure2d[iUp];
  arraySlice1d< real64 const > const & presLow = m_PVTOView.m_undersaturatedPressure2d[iLow];

  idx = LvArray::sortedArrayManipulation::find( presUp.begin(),
                                                presUp.size(),
                                                deltaPres );
  integer const iUpP  = LvArray::math::min( LvArray::math::max( idx, 1 ), LvArray::integerConversion< integer >( presUp.size()-1 ) );
  integer const iLowP = iUpP-1;

  // Step 3: interpolate for Bo
  arraySlice1d< real64 const > const & BoUp = m_PVTOView.m_undersaturatedBo2d[iUp];
  arraySlice1d< real64 const > const & BoLow = m_PVTOView.m_undersaturatedBo2d[iLow];

  real64 const BoInterpLow = interpolation::linearInterpolation( deltaPres-presLow[iLowP], presLow[iUpP]-deltaPres,
                                                                 BoLow[iLowP], BoLow[iUpP] );
  real64 const BoInterpUp = interpolation::linearInterpolation( deltaPres-presUp[iLowP], presUp[iUpP]-deltaPres,
                                                                BoUp[iLowP], BoUp[iUpP] );
  Bo = interpolation::linearInterpolation( deltaRsLow, deltaRsUp, BoInterpLow, BoInterpUp );

  // Step 4: interpolate for viscosity
  arraySlice1d< real64 const > const & viscUp = m_PVTOView.m_undersaturatedViscosity2d[iUp];
  arraySlice1d< real64 const > const & viscLow = m_PVTOView.m_undersaturatedViscosity2d[iLow];

  real64 const viscInterpLow = interpolation::linearInterpolation( deltaPres-presLow[iLowP], presLow[iUpP]-deltaPres,
                                                                   viscLow[iLowP], viscLow[iUpP] );
  real64 const viscInterpUp = interpolation::linearInterpolation( deltaPres-presUp[iLowP], presUp[iUpP]-deltaPres,
                                                                  viscUp[iLowP], viscUp[iUpP] );
  visc = interpolation::linearInterpolation( deltaRsLow, deltaRsUp, viscInterpLow, viscInterpUp );
}

GEOSX_HOST_DEVICE
inline void
BlackOilFluid::KernelWrapper::
  computeMassMoleDensity( bool const needDerivs,
                          bool const useMass,
                          integer const numHydrocarbonComp,
                          real64 const Rs,
                          real64 const dRs_dPres,
                          real64 const dRs_dComp[],
                          real64 const Bo,
                          real64 const dBo_dPres,
                          real64 const dBo_dComp[],
                          real64 & dens,
                          real64 & dDens_dPres,
                          arraySlice1d< real64, multifluid::USD_PHASE_DC - 3 > const & dDens_dComp ) const
{
  using PT = BlackOilFluid::PhaseType;

  real64 const oilDens = (useMass)? m_PVTOView.m_surfaceMassDensity[PT::OIL]:
                         m_PVTOView.m_surfaceMoleDensity[PT::OIL];
  real64 const gasDens = (useMass)? m_PVTOView.m_surfaceMassDensity[PT::GAS]:
                         m_PVTOView.m_surfaceMoleDensity[PT::GAS];

  real64 const Binv = 1. / Bo;
  real64 const tmp = ( oilDens + gasDens * Rs );
  dens =  Binv * tmp;
  if( needDerivs )
  {
    dDens_dPres = Binv * Binv * (Bo * gasDens * dRs_dPres - tmp * dBo_dPres);
    for( integer i = 0; i < numHydrocarbonComp; ++i )
    {
      dDens_dComp[i] = Binv * Binv * (Bo * gasDens * dRs_dComp[i] - tmp * dBo_dComp[i]);
    }
  }
}

GEOSX_HOST_DEVICE
inline void
BlackOilFluid::KernelWrapper::
  update( localIndex const k,
          localIndex const q,
          real64 const pressure,
          real64 const temperature,
          arraySlice1d< geosx::real64 const, compflow::USD_COMP - 1 > const & composition ) const
{
  compute( pressure,
           temperature,
           composition,
           m_phaseFraction( k, q ),
           m_phaseDensity( k, q ),
           m_phaseMassDensity( k, q ),
           m_phaseViscosity( k, q ),
           m_phaseCompFraction( k, q ),
           m_totalDensity( k, q ) );
}

} // namespace constitutive

} // namespace geosx

#endif // GEOSX_CONSTITUTIVE_FLUID_BLACKOILFLUID_HPP_
