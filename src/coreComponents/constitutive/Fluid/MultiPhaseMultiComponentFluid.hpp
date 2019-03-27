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
  * @file MultiPhaseMultiComponentFluid.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_MULTIPHASEMULTICOMPONENTFLUID_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_MULTIPHASEMULTICOMPONENTFLUID_HPP

#include "constitutive/Fluid/MultiFluidBase.hpp"

#include <memory>

#include "PVTFunctions/FlashModelBase.hpp"
#include "PVTFunctions/PVTFunctionBase.hpp"


namespace geosx
{

  /*
using namespace PVTProps;
using namespace Utility;
  */
  
namespace dataRepository
{
namespace keys
{
string const multiPhaseMultiComponentFluid = "MultiPhaseMultiComponentFluid";
}
}

namespace constitutive
{

class MultiPhaseMultiComponentFluid : public MultiFluidBase
{
public:

  MultiPhaseMultiComponentFluid( std::string const & name, ManagedGroup * const parent );

  virtual ~MultiPhaseMultiComponentFluid() override;

  std::unique_ptr<ConstitutiveBase> DeliverClone( string const & name,
                                                  ManagedGroup * const parent ) const override;

  static std::string CatalogName() { return dataRepository::keys::multiPhaseMultiComponentFluid; }

  virtual string GetCatalogName() override { return CatalogName(); }


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
                       const PVTProps::array1dT<PVTProps::PVTFunction>& phaseDensityFuns,
                       const PVTProps::array1dT<PVTProps::PVTFunction>& phaseViscosityFuns,
                       const PVTProps::FlashModel flashModel);

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

  struct viewKeyStruct : MultiFluidBase::viewKeyStruct
  {
    static constexpr auto flashModelParaFileString = "flashModelParaFile";
    static constexpr auto phasePVTParaFilesString = "phasePVTParaFiles";
    
    using ViewKey = dataRepository::ViewKey;

    ViewKey flashModelParaFile = { flashModelParaFileString };
    ViewKey phasePVTParaFiles       = { phasePVTParaFilesString };

  } viewKeysMultiPhaseMultiComponentFluid;

  
protected:
  virtual void PostProcessInput() override;

  virtual void InitializePostSubGroups( ManagedGroup * const group ) override;

private:
  
  void CreatePVTModels();

  // phase PVT parameter filenames
  string_array m_phasePVTParaFiles;

  string m_flashModelParaFile;

  // number of entries corrosponds to number of phases
  PVTProps::array1dT<PVTProps::PVTFunction> m_phaseDensityFuns;
  PVTProps::array1dT<PVTProps::PVTFunction> m_phaseViscosityFuns;  

  PVTProps::FlashModel m_flashModel;
  
};

} //namespace constitutive

} //namespace geosx

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_MULTIPHASEMULTICOMPONENTFLUID_HPP
