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

#include <memory>

namespace PVTPackage
{
class MultiphaseSystem;
enum class PHASE_TYPE : int;
}

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

  MultiFluidPVTPackageWrapperUpdate( PVTPackage::MultiphaseSystem & fluid,
                                     arrayView1d< PVTPackage::PHASE_TYPE > const & phaseTypes,
                                     arrayView1d< real64 const > const & componentMolarWeight,
                                     bool useMass,
                                     arrayView3d< real64 > const & phaseFraction,
                                     arrayView3d< real64 > const & dPhaseFractionDPressure,
                                     arrayView3d< real64 > const & dPhaseFractionDTemperature,
                                     arrayView4d< real64 > const & dPhaseFractionDGlobalCompFraction,
                                     arrayView3d< real64 > const & phaseDensity,
                                     arrayView3d< real64 > const & dPhaseDensityDPressure,
                                     arrayView3d< real64 > const & dPhaseDensityDTemperature,
                                     arrayView4d< real64 > const & dPhaseDensityDGlobalCompFraction,
                                     arrayView3d< real64 > const & phaseViscosity,
                                     arrayView3d< real64 > const & dPhaseViscosityDPressure,
                                     arrayView3d< real64 > const & dPhaseViscosityDTemperature,
                                     arrayView4d< real64 > const & dPhaseViscosityDGlobalCompFraction,
                                     arrayView4d< real64 > const & phaseCompFraction,
                                     arrayView4d< real64 > const & dPhaseCompFractionDPressure,
                                     arrayView4d< real64 > const & dPhaseCompFractionDTemperature,
                                     arrayView5d< real64 > const & dPhaseCompFractionDGlobalCompFraction,
                                     arrayView2d< real64 > const & totalDensity,
                                     arrayView2d< real64 > const & dTotalDensityDPressure,
                                     arrayView2d< real64 > const & dTotalDensityDTemperature,
                                     arrayView3d< real64 > const & dTotalDensityDGlobalCompFraction )
    : MultiFluidBaseUpdate( componentMolarWeight,
                            useMass,
                            phaseFraction,
                            dPhaseFractionDPressure,
                            dPhaseFractionDTemperature,
                            dPhaseFractionDGlobalCompFraction,
                            phaseDensity,
                            dPhaseDensityDPressure,
                            dPhaseDensityDTemperature,
                            dPhaseDensityDGlobalCompFraction,
                            phaseViscosity,
                            dPhaseViscosityDPressure,
                            dPhaseViscosityDTemperature,
                            dPhaseViscosityDGlobalCompFraction,
                            phaseCompFraction,
                            dPhaseCompFractionDPressure,
                            dPhaseCompFractionDTemperature,
                            dPhaseCompFractionDGlobalCompFraction,
                            totalDensity,
                            dTotalDensityDPressure,
                            dTotalDensityDTemperature,
                            dTotalDensityDGlobalCompFraction ),
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

  virtual void Compute( real64 const pressure,
                        real64 const temperature,
                        arraySlice1d< real64 const > const & composition,
                        arraySlice1d< real64 > const & phaseFraction,
                        arraySlice1d< real64 > const & phaseDensity,
                        arraySlice1d< real64 > const & phaseViscosity,
                        arraySlice2d< real64 > const & phaseCompFraction,
                        real64 & totalDensity ) const override;

  virtual void Compute( real64 const pressure,
                        real64 const temperature,
                        arraySlice1d< real64 const > const & composition,
                        arraySlice1d< real64 > const & phaseFraction,
                        arraySlice1d< real64 > const & dPhaseFractionDPressure,
                        arraySlice1d< real64 > const & dPhaseFractionDTemperature,
                        arraySlice2d< real64 > const & dPhaseFractionDGlobalCompFraction,
                        arraySlice1d< real64 > const & phaseDensity,
                        arraySlice1d< real64 > const & dPhaseDensityDPressure,
                        arraySlice1d< real64 > const & dPhaseDensityDTemperature,
                        arraySlice2d< real64 > const & dPhaseDensityDGlobalCompFraction,
                        arraySlice1d< real64 > const & phaseViscosity,
                        arraySlice1d< real64 > const & dPhaseViscosityDPressure,
                        arraySlice1d< real64 > const & dPhaseViscosityDTemperature,
                        arraySlice2d< real64 > const & dPhaseViscosityDGlobalCompFraction,
                        arraySlice2d< real64 > const & phaseCompFraction,
                        arraySlice2d< real64 > const & dPhaseCompFractionDPressure,
                        arraySlice2d< real64 > const & dPhaseCompFractionDTemperature,
                        arraySlice3d< real64 > const & dPhaseCompFractionDGlobalCompFraction,
                        real64 & totalDensity,
                        real64 & dTotalDensityDPressure,
                        real64 & dTotalDensityDTemperature,
                        arraySlice1d< real64 > const & dTotalDensityDGlobalCompFraction ) const override;

  GEOSX_FORCE_INLINE
  virtual void Update( localIndex const k,
                       localIndex const q,
                       real64 const pressure,
                       real64 const temperature,
                       arraySlice1d< real64 const > const & composition ) const override
  {
    Compute( pressure,
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

  PVTPackage::MultiphaseSystem & m_fluid;

  arrayView1d< PVTPackage::PHASE_TYPE > m_phaseTypes;

};

class MultiFluidPVTPackageWrapper : public MultiFluidBase
{
public:

  MultiFluidPVTPackageWrapper( std::string const & name, Group * const parent );

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

  virtual void PostProcessInput() override;

  virtual void InitializePostSubGroups( Group * const group ) override;

  /// function that populates m_fluid ptr; to be overriden by derived classes
  virtual void createFluid() = 0;

  /// PVTPackage fluid object
  std::unique_ptr< PVTPackage::MultiphaseSystem > m_fluid;

  /// PVTPackage phase labels
  array1d< PVTPackage::PHASE_TYPE > m_phaseTypes;
};

} //namespace constitutive

} //namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDPVTPACKAGEWRAPPER_HPP_
