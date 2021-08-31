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
 * @file BlackOilFluid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_BLACKOILFLUID_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_BLACKOILFLUID_HPP_

#include "constitutive/fluid/BlackOilFluidBase.hpp"
#include "constitutive/fluid/PVTOData.hpp"
#include "functions/TableFunction.hpp"
//#include "codingUtilities/Utilities.hpp"
#include "math/interpolation/Interpolation.hpp"

namespace geosx
{

namespace constitutive
{

/**
 * @brief Kernel wrapper class for BlackOilFluid
 *        This kernel can be called on the GPU
 */
class BlackOilFluidUpdate final : public MultiFluidBaseUpdate
{
public:

  static constexpr real64 minForPhasePresence = 1e-10;

  static constexpr localIndex NC_BO = 3;
  static constexpr localIndex NP_BO = 3;
  static constexpr localIndex HNC_BO = NC_BO - 1;

  BlackOilFluidUpdate( arrayView1d< integer const > const & phaseTypes,
                       arrayView1d< integer const > const & phaseOrder,
                       arrayView1d< integer const > const & hydrocarbonPhaseOrder,
                       arrayView1d< real64 const > const & surfacePhaseMassDensity,
                       arrayView1d< TableFunction::KernelWrapper const > const & formationVolFactorTables,
                       arrayView1d< TableFunction::KernelWrapper const > const & viscosityTables,
                       real64 const waterRefPressure,
                       real64 const waterFormationVolFactor,
                       real64 const waterCompressibility,
                       real64 const waterViscosity,
                       arrayView1d< real64 const > const & componentMolarWeight,
                       bool useMass,
                       arrayView3d< real64, multifluid::USD_PHASE > const & phaseFraction,
                       arrayView3d< real64, multifluid::USD_PHASE > const & dPhaseFraction_dPressure,
                       arrayView3d< real64, multifluid::USD_PHASE > const & dPhaseFraction_dTemperature,
                       arrayView4d< real64, multifluid::USD_PHASE_DC > const & dPhaseFraction_dGlobalCompFraction,
                       arrayView3d< real64, multifluid::USD_PHASE > const & phaseDensity,
                       arrayView3d< real64, multifluid::USD_PHASE > const & dPhaseDensity_dPressure,
                       arrayView3d< real64, multifluid::USD_PHASE > const & dPhaseDensity_dTemperature,
                       arrayView4d< real64, multifluid::USD_PHASE_DC > const & dPhaseDensity_dGlobalCompFraction,
                       arrayView3d< real64, multifluid::USD_PHASE > const & phaseMassDensity,
                       arrayView3d< real64, multifluid::USD_PHASE > const & dPhaseMassDensity_dPressure,
                       arrayView3d< real64, multifluid::USD_PHASE > const & dPhaseMassDensity_dTemperature,
                       arrayView4d< real64, multifluid::USD_PHASE_DC > const & dPhaseMassDensity_dGlobalCompFraction,
                       arrayView3d< real64, multifluid::USD_PHASE > const & phaseViscosity,
                       arrayView3d< real64, multifluid::USD_PHASE > const & dPhaseViscosity_dPressure,
                       arrayView3d< real64, multifluid::USD_PHASE > const & dPhaseViscosity_dTemperature,
                       arrayView4d< real64, multifluid::USD_PHASE_DC > const & dPhaseViscosity_dGlobalCompFraction,
                       arrayView4d< real64, multifluid::USD_PHASE_COMP > const & phaseCompFraction,
                       arrayView4d< real64, multifluid::USD_PHASE_COMP > const & dPhaseCompFraction_dPressure,
                       arrayView4d< real64, multifluid::USD_PHASE_COMP > const & dPhaseCompFraction_dTemperature,
                       arrayView5d< real64, multifluid::USD_PHASE_COMP_DC > const & dPhaseCompFraction_dGlobalCompFraction,
                       arrayView2d< real64, multifluid::USD_FLUID > const & totalDensity,
                       arrayView2d< real64, multifluid::USD_FLUID > const & dTotalDensity_dPressure,
                       arrayView2d< real64, multifluid::USD_FLUID > const & dTotalDensity_dTemperature,
                       arrayView3d< real64, multifluid::USD_FLUID_DC > const & dTotalDensity_dGlobalCompFraction,
                       PVTOData::KernelWrapper const PVTO )
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
    m_phaseTypes( phaseTypes ),
    m_phaseOrder( phaseOrder ),
    m_hydrocarbonPhaseOrder( hydrocarbonPhaseOrder ),
    m_surfacePhaseMassDensity( surfacePhaseMassDensity ),
    m_formationVolFactorTables( formationVolFactorTables ),
    m_viscosityTables( viscosityTables ),
    m_waterRefPressure( waterRefPressure ),
    m_waterFormationVolFactor( waterFormationVolFactor ),
    m_waterCompressibility( waterCompressibility ),
    m_waterViscosity( waterViscosity ),
    m_PVTOView( PVTO )
  {}

  GEOSX_HOST_DEVICE
  virtual void compute( real64 const pressure,
                        real64 const temperature,
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
                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseFraction_dPressure,
                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseFraction_dTemperature,
                        arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseFraction_dGlobalCompFraction,
                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseDensity_dPressure,
                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseDensity_dTemperature,
                        arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseDensity_dGlobalCompFraction,
                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseMassDensity_dPressure,
                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseMassDensity_dTemperature,
                        arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseMassDensity_dGlobalCompFraction,
                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseViscosity_dPressure,
                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseViscosity_dTemperature,
                        arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseViscosity_dGlobalCompFraction,
                        arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & phaseCompFraction,
                        arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & dPhaseCompFraction_dPressure,
                        arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & dPhaseCompFraction_dTemperature,
                        arraySlice3d< real64, multifluid::USD_PHASE_COMP_DC-2 > const & dPhaseCompFraction_dGlobalCompFraction,
                        real64 & totalDensity,
                        real64 & dTotalDensity_dPressure,
                        real64 & dTotalDensity_dTemperature,
                        arraySlice1d< real64, multifluid::USD_FLUID_DC - 2 > const & dTotalDensity_dGlobalCompFraction ) const override;

  GEOSX_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const pressure,
                       real64 const temperature,
                       arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const override
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

  GEOSX_HOST_DEVICE
  void computeDensitiesViscosities( real64 pressure,
                                    bool needDerivs,
                                    real64 const composition[],
                                    arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & phaseFrac,
                                    arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDens,
                                    arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseDens_dPres,
                                    arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseDens_dGlobalCompFrac,
                                    arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDens,
                                    arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseMassDens_dPres,
                                    arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseMassDens_dGlobalCompFrac,
                                    arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseVisc,
                                    arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseVisc_dPres,
                                    arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseVisc_dGlobalCompFrac,
                                    real64 phaseMW[NP_BO],
                                    real64 dphaseMW_dPres[],
                                    real64 dphaseMW_dGlobalCompFrac[] ) const;

  GEOSX_HOST_DEVICE
  void computeEquilibrium( real64 pressure,
                           bool needDerivs,
                           real64 const composition[],
                           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
                           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseFraction_dPressure,
                           arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseFraction_dGlobalCompFraction,
                           arraySlice2d< real64, multifluid::USD_PHASE_COMP - 2 > const & phaseCompFraction,
                           arraySlice2d< real64, multifluid::USD_PHASE_COMP - 2 > const & dPhaseCompFraction_dPressure,
                           arraySlice3d< real64, multifluid::USD_PHASE_COMP_DC - 2 > const & dPhaseCompFraction_dGlobalCompFraction ) const;

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
                      arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
                      arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseFraction_dPressure,
                      arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseFraction_dGlobalCompFraction,
                      arraySlice2d< real64, multifluid::USD_PHASE_COMP - 2 > const & phaseCompFraction,
                      arraySlice2d< real64, multifluid::USD_PHASE_COMP - 2 > const & dPhaseCompFraction_dPressure,
                      arraySlice3d< real64, multifluid::USD_PHASE_COMP_DC - 2 > const & dPhaseCompFraction_dGlobalCompFraction,
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
                                         localIndex const numHydrocarbonComp,
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
                               localIndex const numHydrocarbonComp,
                               real64 const Rs,
                               real64 const dRs_dPres,
                               real64 const dRs_dComp[],
                               real64 const Bo,
                               real64 const dBo_dPres,
                               real64 const dBo_dComp[],
                               real64 & dens,
                               real64 & dDens_dPres,
                               arraySlice1d< real64, multifluid::USD_PHASE_DC - 3 > const & dDens_dC ) const;

  /// Phase ordering info
  arrayView1d< integer const > m_phaseTypes;
  arrayView1d< integer const > m_phaseOrder;
  arrayView1d< integer const > m_hydrocarbonPhaseOrder;

  /// Surface mass density for each phase
  arrayView1d< real64 const > m_surfacePhaseMassDensity;

  /// Table kernel wrappers to interpolate in the oil and gas (B vs p) tables
  arrayView1d< TableFunction::KernelWrapper const > m_formationVolFactorTables;

  /// Table kernel wrappers to interpolate in the oil and gas (\mu vs p) tables
  arrayView1d< TableFunction::KernelWrapper const > m_viscosityTables;

  /// Water reference pressure
  real64 m_waterRefPressure;

  /// Water formation volume factor
  real64 m_waterFormationVolFactor;

  /// Water compressibility
  real64 m_waterCompressibility;

  /// Water viscosity
  real64 m_waterViscosity;

  /// Data needed to update the oil phase properties
  PVTOData::KernelWrapper m_PVTOView;
};


class BlackOilFluid : public BlackOilFluidBase
{
public:

  BlackOilFluid( string const & name, Group * const parent );

  virtual ~BlackOilFluid() override = default;

  virtual std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                Group * const parent ) const override;

  static string catalogName() { return "BlackOilFluid"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = BlackOilFluidUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper()
  {
    return KernelWrapper( m_phaseTypes,
                          m_phaseOrder,
                          m_hydrocarbonPhaseOrder,
                          m_surfacePhaseMassDensity,
                          m_formationVolFactorTables,
                          m_viscosityTables,
                          m_waterRefPressure,
                          m_waterFormationVolFactor,
                          m_waterCompressibility,
                          m_waterViscosity,
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
                          m_dTotalDensity_dGlobalCompFraction,
                          m_PVTO.createKernelWrapper());
  }

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
                     real64 const oilSurfaceMassDensity,
                     real64 const oilSurfaceMolecularWeight,
                     real64 const gasSurfaceMassDensity,
                     real64 const gasSurfaceMolecularWeight );

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
  void refineUndersaturatedTables( localIndex const numRefinedPresPoints );
  void refineTableAndCopyOld( localIndex const nLevels );

  /**
   * @brief Check the monotonicity of the PVTO table values
   */
  void checkTableConsistency() const;

  /// The data needed to compute the oil phase properties
  PVTOData m_PVTO;

};

GEOSX_HOST_DEVICE
inline void
BlackOilFluidUpdate::compute( real64 pressure,
                              real64 temperature,
                              arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                              arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
                              arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
                              arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
                              arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
                              arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & phaseCompFraction,
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
BlackOilFluidUpdate::compute( real64 pressure,
                              real64 temperature,
                              arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                              arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
                              arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseFraction_dPressure,
                              arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseFraction_dTemperature,
                              arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseFraction_dGlobalCompFraction,
                              arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
                              arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseDensity_dPressure,
                              arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseDensity_dTemperature,
                              arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseDensity_dGlobalCompFraction,
                              arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
                              arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseMassDensity_dPressure,
                              arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseMassDensity_dTemperature,
                              arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseMassDensity_dGlobalCompFraction,
                              arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
                              arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseViscosity_dPressure,
                              arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseViscosity_dTemperature,
                              arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseViscosity_dGlobalCompFraction,
                              arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & phaseCompFraction,
                              arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & dPhaseCompFraction_dPressure,
                              arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & dPhaseCompFraction_dTemperature,
                              arraySlice3d< real64, multifluid::USD_PHASE_COMP_DC-2 > const & dPhaseCompFraction_dGlobalCompFraction,
                              real64 & totalDensity,
                              real64 & dTotalDensity_dPressure,
                              real64 & dTotalDensity_dTemperature,
                              arraySlice1d< real64, multifluid::USD_FLUID_DC - 2 > const & dTotalDensity_dGlobalCompFraction ) const
{
  GEOSX_UNUSED_VAR( temperature );
  GEOSX_UNUSED_VAR( dPhaseFraction_dTemperature );
  GEOSX_UNUSED_VAR( dPhaseDensity_dTemperature );
  GEOSX_UNUSED_VAR( dPhaseMassDensity_dTemperature );
  GEOSX_UNUSED_VAR( dPhaseViscosity_dTemperature );
  GEOSX_UNUSED_VAR( dPhaseCompFraction_dTemperature );
  GEOSX_UNUSED_VAR( dTotalDensity_dTemperature );

  real64 compMoleFrac[NC_BO];
  real64 dCompMoleFrac_dCompMassFrac[NC_BO][NC_BO];
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
    for( localIndex ic = 0; ic < NC_BO; ++ic )
    {
      compMoleFrac[ic] = composition[ic];
    }
  }

  // 2. Compute phase fractions and phase component fractions

  computeEquilibrium( pressure,
                      true,
                      compMoleFrac,
                      phaseFraction,
                      dPhaseFraction_dPressure,
                      dPhaseFraction_dGlobalCompFraction,
                      phaseCompFraction,
                      dPhaseCompFraction_dPressure,
                      dPhaseCompFraction_dGlobalCompFraction );

  // 3. Compute phase densities and viscosities

  computeDensitiesViscosities( pressure,
                               true,
                               compMoleFrac,
                               phaseFraction,
                               phaseDensity,
                               dPhaseDensity_dPressure,
                               dPhaseDensity_dGlobalCompFraction,
                               phaseMassDensity,
                               dPhaseMassDensity_dPressure,
                               dPhaseMassDensity_dGlobalCompFraction,
                               phaseViscosity,
                               dPhaseViscosity_dPressure,
                               dPhaseViscosity_dGlobalCompFraction,
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
                   dPhaseFraction_dPressure,
                   dPhaseFraction_dGlobalCompFraction,
                   phaseCompFraction,
                   dPhaseCompFraction_dPressure,
                   dPhaseCompFraction_dGlobalCompFraction,
                   dPhaseDensity_dGlobalCompFraction,
                   dPhaseViscosity_dGlobalCompFraction );
  }

  // 5. Compute total fluid mass/molar density and derivatives
  totalDensity = 0.0;
  dTotalDensity_dPressure = 0.0;
  for( localIndex ic = 0; ic < NC_BO; ++ic )
  {
    dTotalDensity_dGlobalCompFraction[ic] = 0.0;
  }

  for( localIndex ip = 0; ip < NP_BO; ++ip )
  {
    if( phaseFraction[ip] <= 0. )
    {
      continue;
    }
    real64 const densInv = 1.0 / phaseDensity[ip];
    real64 const value = phaseFraction[ip] * densInv;

    totalDensity += value;
    dTotalDensity_dPressure += ( dPhaseFraction_dPressure[ip] - value * dPhaseDensity_dPressure[ip] ) * densInv;
    for( localIndex ic = 0; ic < NC_BO; ++ic )
    {
      dTotalDensity_dGlobalCompFraction[ic] += ( dPhaseFraction_dGlobalCompFraction[ip][ic]
                                                 - value * dPhaseDensity_dGlobalCompFraction[ip][ic] ) * densInv;
    }
  }

  totalDensity = 1.0 / totalDensity;
  real64 const minusDens2 = -totalDensity * totalDensity;
  dTotalDensity_dPressure *= minusDens2;
  for( localIndex ic = 0; ic < NC_BO; ++ic )
  {
    dTotalDensity_dGlobalCompFraction[ic] *= minusDens2;
  }
}

GEOSX_HOST_DEVICE
inline void
BlackOilFluidUpdate::computeEquilibrium( real64 pressure,
                                         bool needDerivs,
                                         real64 const composition[],
                                         arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
                                         arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseFraction_dPressure,
                                         arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseFraction_dGlobalCompFraction,
                                         arraySlice2d< real64, multifluid::USD_PHASE_COMP - 2 > const & phaseCompFraction,
                                         arraySlice2d< real64, multifluid::USD_PHASE_COMP - 2 > const & dPhaseCompFraction_dPressure,
                                         arraySlice3d< real64, multifluid::USD_PHASE_COMP_DC - 2 > const & dPhaseCompFraction_dGlobalCompFraction ) const
{
  using PT = BlackOilFluid::PhaseType;

  localIndex const ipOil   = m_phaseOrder[PT::OIL];
  localIndex const ipGas   = m_phaseOrder[PT::GAS];
  localIndex const ipWater = m_phaseOrder[PT::WATER];

  localIndex const icOil   = ipOil;
  localIndex const icGas   = ipGas;
  localIndex const icWater = ipWater;

  real64 const zo = composition[icOil];
  real64 const zg = composition[icGas];
  real64 const zw = composition[icWater];

  // 1. Make everything zero first

  for( localIndex ip = 0; ip < NP_BO; ++ip )
  {
    phaseFraction[ip] = 0.;
    if( needDerivs )
    {
      dPhaseFraction_dPressure[ip] = 0.0;
    }
    for( localIndex ic = 0; ic < NC_BO; ++ic )
    {
      phaseCompFraction[ip][ic] = 0.0;
      if( needDerivs )
      {
        dPhaseFraction_dGlobalCompFraction[ip][ic] = 0.;
        dPhaseCompFraction_dPressure[ip][ic] = 0.0;
        for( localIndex jc = 0; jc < NC_BO; ++jc )
        {
          dPhaseCompFraction_dGlobalCompFraction[ip][ic][jc] = 0.0;
        }
      }
    }
  }

  // 2. Check feed first, and if only water is present (e.g., water inj), then skip

  if( zw >= 1. - minForPhasePresence )
  {
    phaseFraction[ipWater] = zw;
    if( needDerivs )
    {
      dPhaseFraction_dGlobalCompFraction[ipWater][icWater] = 1.;
    }
    phaseCompFraction[ipWater][icWater] = 1.0;
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
  if( ( gasPhaseFraction > 0 ) && ( gasPhaseFraction < 1 ) )
  {

    // phase fractions
    phaseFraction[ipOil] = 1. - gasPhaseFraction - zw;
    phaseFraction[ipGas] = gasPhaseFraction;
    phaseFraction[ipWater] =  zw;

    if( needDerivs )
    {
      dPhaseFraction_dPressure[ipOil] = -dGasPhaseFraction_dP;
      dPhaseFraction_dPressure[ipGas] = dGasPhaseFraction_dP;
      dPhaseFraction_dGlobalCompFraction[ipOil][icOil] = -dGasPhaseFraction_dzo;
      dPhaseFraction_dGlobalCompFraction[ipOil][icGas] = -dGasPhaseFraction_dzg;
      dPhaseFraction_dGlobalCompFraction[ipOil][icWater] = -1.;
      dPhaseFraction_dGlobalCompFraction[ipGas][icOil] = dGasPhaseFraction_dzo;
      dPhaseFraction_dGlobalCompFraction[ipGas][icGas] = dGasPhaseFraction_dzg;
      dPhaseFraction_dGlobalCompFraction[ipWater][icWater] = 1.;
    }

    // oil
    real64 const tmp = ( oilSurfaceMoleDensity + gasSurfaceMoleDensity * RsSat );
    real64 const tmpOil = oilSurfaceMoleDensity / tmp;
    real64 const dTmpOil_dP = -oilSurfaceMoleDensity * gasSurfaceMoleDensity * dRsSat_dP / ( tmp * tmp );
    phaseCompFraction[ipOil][icOil] = tmpOil;
    phaseCompFraction[ipOil][icGas] = 1. - tmpOil;
    phaseCompFraction[ipOil][icWater] = 0.;

    if( needDerivs )
    {
      dPhaseCompFraction_dPressure[ipOil][icOil] = dTmpOil_dP;
      dPhaseCompFraction_dPressure[ipOil][icGas] = -dTmpOil_dP;
    }

    // gas
    real64 const tmpGas = gasSurfaceMoleDensity / ( gasSurfaceMoleDensity );
    phaseCompFraction[ipGas][icOil] = 1. - tmpGas;
    phaseCompFraction[ipGas][icGas] = tmpGas;
    phaseCompFraction[ipGas][icWater] = 0.;

    // water
    phaseCompFraction[ipWater][icOil] = 0;
    phaseCompFraction[ipWater][icGas] = 0;
    phaseCompFraction[ipWater][icWater] = 1.;

  }
  // 4.2 The gas phase is absent
  else
  {

    // phase fractions
    phaseFraction[ipOil] = 1. - zw;
    phaseFraction[ipGas] = 0.;
    phaseFraction[ipWater] =  zw;

    // oil
    phaseCompFraction[ipOil][icOil] = zo;
    phaseCompFraction[ipOil][icGas] = zg;
    phaseCompFraction[ipOil][icWater] = 0.;

    // gas
    phaseCompFraction[ipWater][icOil] = 0;
    phaseCompFraction[ipWater][icGas] = 0;
    phaseCompFraction[ipWater][icWater] = 1.;

    if( needDerivs )
    {
      dPhaseFraction_dGlobalCompFraction[ipOil][icWater] = -1.;
      dPhaseFraction_dGlobalCompFraction[ipWater][icWater] = 1.;
      dPhaseCompFraction_dGlobalCompFraction[ipOil][icOil][icOil] = 1.;
      dPhaseCompFraction_dGlobalCompFraction[ipOil][icGas][icGas] = 1.;
    }
  }

}

GEOSX_HOST_DEVICE
inline void
BlackOilFluidUpdate::computeDensitiesViscosities( real64 pressure,
                                                  bool needDerivs,
                                                  real64 const composition[],
                                                  arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & phaseFrac,
                                                  arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDens,
                                                  arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseDens_dPres,
                                                  arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseDens_dGlobalCompFrac,
                                                  arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDens,
                                                  arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseMassDens_dPres,
                                                  arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseMassDens_dGlobalCompFrac,
                                                  arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseVisc,
                                                  arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseVisc_dPres,
                                                  arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseVisc_dGlobalCompFrac,
                                                  real64 phaseMW[NP_BO],
                                                  real64 dPhaseMW_dPres[],
                                                  real64 dPhaseMW_dGlobalCompFrac[] ) const
{
  using PT = BlackOilFluid::PhaseType;

  real64 fvf = 0.0;
  real64 derivative = 0.0;

  localIndex const ipOil = m_phaseOrder[PT::OIL];
  localIndex const ipGas = m_phaseOrder[PT::GAS];
  localIndex const ipWater = m_phaseOrder[PT::WATER];

  localIndex const icOil = ipOil;
  localIndex const icGas = ipGas;

  bool const isWater = (ipWater >= 0 && phaseFrac[ipWater] > 0);
  bool const isGas = (ipGas >= 0 && phaseFrac[ipGas] > 0);
  bool const isOil = (ipOil >= 0 && phaseFrac[ipOil] > 0);

  for( localIndex ip = 0; ip < NP_BO; ++ip )
  {
    phaseMassDens[ip] = 0.;
    phaseDens[ip]     = 0.;
    phaseVisc[ip]     = 0.;
    dPhaseDens_dPres[ip]     = 0.;
    dPhaseMassDens_dPres[ip] = 0.;
    dPhaseVisc_dPres[ip]     = 0.;
    for( localIndex ic = 0; ic < NC_BO; ++ic )
    {
      dPhaseMassDens_dGlobalCompFrac[ip][ic] = 0.;
      dPhaseMassDens_dGlobalCompFrac[ip][ic] = 0.;
      dPhaseVisc_dGlobalCompFrac[ip][ic]     = 0.;
    }
  }

  // 1. Gas phase: look up in the formation vol factor tables

  if( isGas )
  {
    // interpolate in the table to get the phase formation vol factor and its derivative wrt pressure
    m_formationVolFactorTables[0].compute( &pressure, fvf, &derivative );

    // we are ready to update the densities
    real64 const fvfInv = 1.0 / fvf;

    phaseMassDens[ipGas] = m_surfacePhaseMassDensity[ipGas] * fvfInv;
    real64 const mult = m_useMass ? 1.0 : 1.0 / m_componentMolarWeight[ipGas];
    phaseDens[ipGas] = phaseMassDens[ipGas] * mult;
    phaseMW[ipGas] = m_componentMolarWeight[ipGas];

    if( needDerivs )
    {
      dPhaseMassDens_dPres[ipGas] = -derivative * phaseMassDens[ipGas] * fvfInv;
      dPhaseDens_dPres[ipGas] = dPhaseMassDens_dPres[ipGas] * mult;
      dPhaseMW_dPres[ipGas] = 0.;
      for( localIndex ic = 0; ic < NC_BO; ic++ )
      {
        dPhaseMassDens_dGlobalCompFrac[ipGas][ic] = 0.;
        dPhaseDens_dGlobalCompFrac[ipGas][ic]     = 0.;
        dPhaseMW_dGlobalCompFrac[ipGas*NC_BO+ic]  = 0.;
      }
    }

    m_viscosityTables[0].compute( &pressure, phaseVisc[ipGas], &(dPhaseVisc_dPres)[ipGas] );
  }

  // 2. Water phase: use the constant formation volume factor and compressibility provided by the user

  if( isWater )
  {
    // if water is present
    real64 const expCompDeltaPres = std::exp( -m_waterCompressibility * ( pressure - m_waterRefPressure ) );
    real64 const dExpCompDeltaPres_dPres = -m_waterCompressibility * expCompDeltaPres;
    real64 const denom = m_waterFormationVolFactor * expCompDeltaPres;
    real64 const dDenom_dPres = m_waterFormationVolFactor * dExpCompDeltaPres_dPres;
    real64 const denomInv = 1.0 / denom;
    phaseMassDens[ipWater] = m_surfacePhaseMassDensity[ipWater] * denomInv;
    real64 const mult = m_useMass ? 1.0 : 1.0 / m_componentMolarWeight[ipWater];
    phaseDens[ipWater] = phaseMassDens[ipWater] * mult;
    phaseVisc[ipWater] = m_waterViscosity;
    phaseMW[ipWater] = m_componentMolarWeight[ipWater];

    if( needDerivs )
    {
      dPhaseMassDens_dPres[ipWater] = -dDenom_dPres * phaseMassDens[ipWater] * denomInv;
      dPhaseDens_dPres[ipWater] = dPhaseMassDens_dPres[ipWater] * mult;
      dPhaseMW_dPres[ipWater] = 0.;
      for( localIndex ic = 0; ic < NC_BO; ic++ )
      {
        dPhaseMassDens_dGlobalCompFrac[ipWater][ic] = 0.;
        dPhaseDens_dGlobalCompFrac[ipWater][ic] = 0.0;
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

      // compute undersaturated properties (Bo, viscosity) by two-step interpolation in undersaturated tables
      // this part returns numerical derivatives
      computeUndersaturatedBoViscosity( needDerivs, HNC_BO, pressure, Rs, dRs_dC, Bo, dBo_dP, dBo_dC,
                                        visc, dVisc_dP, dVisc_dC );

    }

    // compute densities
    computeMassMoleDensity( needDerivs, true, HNC_BO, Rs, dRs_dP, dRs_dC, Bo, dBo_dP, dBo_dC,
                            phaseMassDens[ipOil], dPhaseMassDens_dPres[ipOil],
                            dPhaseMassDens_dGlobalCompFrac[ipOil] );
    computeMassMoleDensity( needDerivs, false, HNC_BO, Rs, dRs_dP, dRs_dC, Bo, dBo_dP, dBo_dC,
                            phaseDens[ipOil], dPhaseDens_dPres[ipOil],
                            dPhaseDens_dGlobalCompFrac[ipOil] );

    phaseMW[ipOil] = phaseMassDens[ipOil] / phaseDens[ipOil];
    real64 const tmp = 1. / ( phaseDens[ipOil] * phaseDens[ipOil] );
    dPhaseMW_dPres[ipOil] = ( dPhaseMassDens_dPres[ipOil] * phaseDens[ipOil] - dPhaseDens_dPres[ipOil] * phaseMassDens[ipOil] ) * tmp;
    for( localIndex ic = 0; ic < NC_BO; ++ic )
    {
      dPhaseMW_dGlobalCompFrac[ipOil*NC_BO+ic] = ( dPhaseMassDens_dGlobalCompFrac[ipOil][ic] * phaseDens[ipOil]
                                                   - dPhaseDens_dGlobalCompFrac[ipOil][ic] * phaseMassDens[ipOil] ) * tmp;
    }

    if( m_useMass )
    {
      phaseDens[ipOil] = phaseMassDens[ipOil];
      if( needDerivs )
      {
        dPhaseDens_dPres[ipOil] = dPhaseMassDens_dPres[ipOil];
        for( localIndex i = 0; i < dPhaseMassDens_dGlobalCompFrac[ipOil].size(); ++i )
        {
          dPhaseDens_dGlobalCompFrac[ipOil][i] = dPhaseMassDens_dGlobalCompFrac[ipOil][i];
        }
      }
    }

    // copy viscosity into the final array
    // TODO: skip this step
    phaseVisc[ipOil] = visc;
    if( needDerivs )
    {
      dPhaseVisc_dPres[ipOil] = dVisc_dP;
      for( localIndex i = 0; i < HNC_BO; ++i )
      {
        dPhaseVisc_dGlobalCompFrac[ipOil][i] = dVisc_dC[i];
      }
    }
  }
}

GEOSX_HOST_DEVICE
inline void
BlackOilFluidUpdate::convertToMolar( arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                                     real64 compMoleFrac[],
                                     real64 dCompMoleFrac_dCompMassFrac[NC_BO][NC_BO] ) const
{
  for( localIndex ic = 0; ic < NC_BO; ++ic )
  {
    for( localIndex jc = 0; jc < NC_BO; ++jc )
    {
      dCompMoleFrac_dCompMassFrac[ic][jc] = 0.0;
    }
  }

  real64 totalMolality = 0.0;
  for( localIndex ic = 0; ic < NC_BO; ++ic )
  {
    real64 const mwInv = 1.0 / m_componentMolarWeight[ic];
    compMoleFrac[ic] = composition[ic] * mwInv;
    dCompMoleFrac_dCompMassFrac[ic][ic] = mwInv;
    totalMolality += compMoleFrac[ic];
  }

  real64 const totalMolalityInv = 1.0 / totalMolality;
  for( localIndex ic = 0; ic < NC_BO; ++ic )
  {
    compMoleFrac[ic] *= totalMolalityInv;

    for( localIndex jc = 0; jc < NC_BO; ++jc )
    {
      dCompMoleFrac_dCompMassFrac[ic][jc] -= compMoleFrac[ic] / m_componentMolarWeight[jc];
      dCompMoleFrac_dCompMassFrac[ic][jc] *= totalMolalityInv;
    }
  }

}


GEOSX_HOST_DEVICE
inline void
BlackOilFluidUpdate::convertToMass( real64 const dCompMoleFrac_dCompMassFrac[NC_BO][NC_BO],
                                    real64 const phaseMW[NP_BO],
                                    real64 const dPhaseMW_dPres[NP_BO],
                                    real64 const dPhaseMW_dGlobalCompFrac[NP_BO *NC_BO],
                                    arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFrac,
                                    arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseFrac_dPres,
                                    arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseFrac_dGlobalCompFrac,
                                    arraySlice2d< real64, multifluid::USD_PHASE_COMP - 2 > const & phaseCompFrac,
                                    arraySlice2d< real64, multifluid::USD_PHASE_COMP - 2 > const & dPhaseCompFrac_dPres,
                                    arraySlice3d< real64, multifluid::USD_PHASE_COMP_DC - 2 > const & dPhaseCompFrac_dGlobalCompFrac,
                                    arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseDens_dGlobalCompFrac,
                                    arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseVisc_dGlobalCompFrac ) const
{
  // 1. Convert phase fractions (requires two passes)

  real64 totalMass{};
  real64 dTotalMass_dP{};
  real64 dTotalMass_dC[NC_BO]{};

  // 1.1. Compute mass of each phase and total mass (on a 1-mole basis)

  for( localIndex ip = 0; ip < NP_BO; ++ip )
  {
    bool const phaseExists = (phaseFrac[ip] > 0);
    if( !phaseExists )
    {
      continue;
    }

    real64 const nu = phaseFrac[ip];

    phaseFrac[ip] *= phaseMW[ip];
    dPhaseFrac_dPres[ip] = dPhaseFrac_dPres[ip] * phaseMW[ip] + nu * dPhaseMW_dPres[ip];
    totalMass += phaseFrac[ip];
    dTotalMass_dP += dPhaseFrac_dPres[ip];

    for( localIndex jc = 0; jc < NC_BO; ++jc )
    {
      dPhaseFrac_dGlobalCompFrac[ip][jc] = dPhaseFrac_dGlobalCompFrac[ip][jc] * phaseMW[ip] + nu * dPhaseMW_dGlobalCompFrac[jc + ip*NC_BO];
      dTotalMass_dC[jc] += dPhaseFrac_dGlobalCompFrac[ip][jc];
    }
  }

  // 1.2. Normalize to get mass fractions

  real64 const totalMassInv = 1.0 / totalMass;
  for( localIndex ip = 0; ip < NC_BO; ++ip )
  {
    bool const phaseExists = (phaseFrac[ip] > 0);
    if( !phaseExists )
    {
      continue;
    }

    phaseFrac[ip] *= totalMassInv;
    dPhaseFrac_dPres[ip] = ( dPhaseFrac_dPres[ip] - phaseFrac[ip] * dTotalMass_dP ) * totalMassInv;

    for( localIndex jc = 0; jc < NC_BO; ++jc )
    {
      dPhaseFrac_dGlobalCompFrac[ip][jc] = ( dPhaseFrac_dGlobalCompFrac[ip][jc] - phaseFrac[ip] * dTotalMass_dC[jc] ) * totalMassInv;
    }
  }

  // 2. Convert phase compositions

  for( localIndex ip = 0; ip < NC_BO; ++ip )
  {
    bool const phaseExists = (phaseFrac[ip] > 0);
    if( !phaseExists )
    {
      continue;
    }

    real64 const phaseMWInv = 1.0 / phaseMW[ip];

    for( localIndex ic = 0; ic < NC_BO; ++ic )
    {

      real64 const compMW = m_componentMolarWeight[ic];

      phaseCompFrac[ip][ic] = phaseCompFrac[ip][ic] * compMW * phaseMWInv;
      dPhaseCompFrac_dPres[ip][ic] =
        ( dPhaseCompFrac_dPres[ip][ic] * compMW - phaseCompFrac[ip][ic] * dPhaseMW_dPres[ip] ) * phaseMWInv;

      for( localIndex jc = 0; jc < NC_BO; ++jc )
      {
        dPhaseCompFrac_dGlobalCompFrac[ip][ic][jc] =
          ( dPhaseCompFrac_dGlobalCompFrac[ip][ic][jc] * compMW - phaseCompFrac[ip][ic] * dPhaseMW_dGlobalCompFrac[jc + ip*NC_BO] ) * phaseMWInv;
      }
    }
  }

  // 3. Update derivatives w.r.t. mole fractions to derivatives w.r.t mass fractions

  real64 work[NC_BO];
  for( localIndex ip = 0; ip < NC_BO; ++ip )
  {
    bool const phaseExists = (phaseFrac[ip] > 0);
    if( !phaseExists )
    {
      continue;
    }

    applyChainRuleInPlace( NC_BO, dCompMoleFrac_dCompMassFrac, dPhaseFrac_dGlobalCompFrac[ip], work );
    applyChainRuleInPlace( NC_BO, dCompMoleFrac_dCompMassFrac, dPhaseDens_dGlobalCompFrac[ip], work );
    applyChainRuleInPlace( NC_BO, dCompMoleFrac_dCompMassFrac, dPhaseVisc_dGlobalCompFrac[ip], work );
    for( localIndex ic = 0; ic < NC_BO; ++ic )
    {
      applyChainRuleInPlace( NC_BO, dCompMoleFrac_dCompMassFrac, dPhaseCompFrac_dGlobalCompFrac[ip][ic], work );
    }
  }
}

GEOSX_HOST_DEVICE
inline void
BlackOilFluidUpdate::computeRs( real64 const presBub,
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
BlackOilFluidUpdate::computeSaturatedBoViscosity( real64 const Rs,
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
BlackOilFluidUpdate::computeUndersaturatedBoViscosity( bool const needDerivs,
                                                       localIndex const numHydrocarbonComp,
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
    for( localIndex i = 0; i < numHydrocarbonComp; ++i )
    {
      dBo_dComp[i] = dBo_dRs * dRs_dComp[i];
      dVisc_dComp[i] = dVisc_dRs * dRs_dComp[i];
    }
  }

}

GEOSX_HOST_DEVICE
inline void
BlackOilFluidUpdate::computeUndersaturatedBoViscosity( real64 const Rs,
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
BlackOilFluidUpdate::computeMassMoleDensity( bool const needDerivs,
                                             bool const useMass,
                                             localIndex const numHydrocarbonComp,
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
    for( localIndex i = 0; i < numHydrocarbonComp; ++i )
    {
      dDens_dComp[i] = Binv * Binv * (Bo * gasDens * dRs_dComp[i] - tmp * dBo_dComp[i]);
    }
  }
}

} // namespace constitutive

} // namespace geosx

#endif // GEOSX_CONSTITUTIVE_FLUID_BLACKOILFLUID_HPP_
