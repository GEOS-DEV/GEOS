/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
  * @file MultiPhaseMultiComponentFluid.hpp
  */

#ifndef GEOSX_CONSTITUTIVE_FLUID_MULTIPHASEMULTICOMPONENTFLUID_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_MULTIPHASEMULTICOMPONENTFLUID_HPP_

#include "constitutive/fluid/MultiFluidBase.hpp"

#include <memory>

#include "PVTFunctions/FlashModelBase.hpp"
#include "PVTFunctions/PVTFunctionBase.hpp"


namespace geosx
{

// using namespace PVTProps;
// using namespace Utility;
  
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

  MultiPhaseMultiComponentFluid( std::string const & name, Group * const parent );

  virtual ~MultiPhaseMultiComponentFluid() override;

  virtual void DeliverClone( string const & name,
                             Group * const parent,
                             std::unique_ptr<ConstitutiveBase> & clone ) const override;

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
                       array1d<std::shared_ptr<PVTProps::PVTFunction>> const & phaseDensityFuns,
                       array1d<std::shared_ptr<PVTProps::PVTFunction>> const & phaseViscosityFuns,
                       std::shared_ptr<PVTProps::FlashModel> const & flashModel);

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

  virtual void InitializePostSubGroups( Group * const group ) override;

private:
  
  void CreatePVTModels();

  // phase PVT parameter filenames
  path_array m_phasePVTParaFiles;

  string m_flashModelParaFile;

  // number of entries corrosponds to number of phases
  array1d<std::shared_ptr<PVTProps::PVTFunction> > m_phaseDensityFuns;
  array1d<std::shared_ptr<PVTProps::PVTFunction> > m_phaseViscosityFuns;

  std::shared_ptr<PVTProps::FlashModel> m_flashModel;
  
};

} //namespace constitutive

} //namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_MULTIPHASEMULTICOMPONENTFLUID_HPP_
