/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
  * @file MultiFluidPVTPackageWrapper.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_MULTIFLUIDPVTPACKAGEWRAPPER_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_MULTIFLUIDPVTPACKAGEWRAPPER_HPP

#include "constitutive/Fluid/MultiFluidBase.hpp"

// PVTPackage includes
#include "MultiphaseSystem/PVTEnums.hpp"

#include <memory>

namespace PVTPackage
{
class MultiphaseSystem;
}

namespace geosx
{

namespace constitutive
{

class MultiFluidPVTPackageWrapper : public MultiFluidBase
{
public:

  MultiFluidPVTPackageWrapper( std::string const & name, ManagedGroup * const parent );

  virtual ~MultiFluidPVTPackageWrapper() override;


  virtual void PointUpdate( real64 const & pressure,
                            real64 const & temperature,
                            arraySlice1d<real64 const> const & composition,
                            localIndex const k,
                            localIndex const q ) override;

  virtual void BatchUpdate( arrayView1d<real64 const> const & pressure,
                            arrayView1d<real64 const> const & temperature,
                            arrayView2d<real64 const> const & composition ) override;

  static void Compute( localIndex const NC, localIndex const NP, bool const useMass,
                       arrayView1d<string const> const & phaseNames,
                       arrayView1d<real64 const> const & componentMolarWeight,
                       real64 const & pressure,
                       real64 const & temperature,
                       arraySlice1d<real64 const> const & composition,
                       arraySlice1d<real64> const & phaseFraction,
                       arraySlice1d<real64> const & dPhaseFraction_dPressure,
                       arraySlice1d<real64> const & dPhaseFraction_dTemperature,
                       arraySlice2d<real64> const & dPhaseFraction_dGlobalCompFraction,
                       arraySlice1d<real64> const & phaseDensity,
                       arraySlice1d<real64> const & dPhaseDensity_dPressure,
                       arraySlice1d<real64> const & dPhaseDensity_dTemperature,
                       arraySlice2d<real64> const & dPhaseDensity_dGlobalCompFraction,
                       arraySlice1d<real64> const & phaseViscosity,
                       arraySlice1d<real64> const & dPhaseViscosity_dPressure,
                       arraySlice1d<real64> const & dPhaseViscosity_dTemperature,
                       arraySlice2d<real64> const & dPhaseViscosity_dGlobalCompFraction,
                       arraySlice2d<real64> const & phaseCompFraction,
                       arraySlice2d<real64> const & dPhaseCompFraction_dPressure,
                       arraySlice2d<real64> const & dPhaseCompFraction_dTemperature,
                       arraySlice3d<real64> const & dPhaseCompFraction_dGlobalCompFraction,
                       real64 & totalDensity,
                       real64 & dTotalDensity_dPressure,
                       real64 & dTotalDensity_dTemperature,
                       arraySlice1d<real64> const & dTotalDensity_dGlobalCompFraction,
                       PVTPackage::MultiphaseSystem * const fluid );

  virtual void Compute( real64 const & pressure,
                        real64 const & temperature,
                        arraySlice1d<real64 const> const & composition,
                        arraySlice1d<real64> const & phaseFraction,
                        arraySlice1d<real64> const & dPhaseFraction_dPressure,
                        arraySlice1d<real64> const & dPhaseFraction_dTemperature,
                        arraySlice2d<real64> const & dPhaseFraction_dGlobalCompFraction,
                        arraySlice1d<real64> const & phaseDensity,
                        arraySlice1d<real64> const & dPhaseDensity_dPressure,
                        arraySlice1d<real64> const & dPhaseDensity_dTemperature,
                        arraySlice2d<real64> const & dPhaseDensity_dGlobalCompFraction,
                        arraySlice1d<real64> const & phaseViscosity,
                        arraySlice1d<real64> const & dPhaseViscosity_dPressure,
                        arraySlice1d<real64> const & dPhaseViscosity_dTemperature,
                        arraySlice2d<real64> const & dPhaseViscosity_dGlobalCompFraction,
                        arraySlice2d<real64> const & phaseCompFraction,
                        arraySlice2d<real64> const & dPhaseCompFraction_dPressure,
                        arraySlice2d<real64> const & dPhaseCompFraction_dTemperature,
                        arraySlice3d<real64> const & dPhaseCompFraction_dGlobalCompFraction,
                        real64 & totalDensity,
                        real64 & dTotalDensity_dPressure,
                        real64 & dTotalDensity_dTemperature,
                        arraySlice1d<real64> const & dTotalDensity_dGlobalCompFraction ) const override;

protected:
  virtual void PostProcessInput() override;

  virtual void InitializePostSubGroups( ManagedGroup * const group ) override;

  // function that populates m_fluid ptr; to be overriden by derived classes
  virtual void createFluid() = 0;

  // PVTPackage phase type labels
  array1d<PVTPackage::PHASE_TYPE> m_pvtPackagePhaseTypes;

  // PVTPackage fluid object
  std::unique_ptr<PVTPackage::MultiphaseSystem> m_fluid;

};

} //namespace constitutive

} //namespace geosx

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_MULTIFLUIDPVTPACKAGEWRAPPER_HPP
