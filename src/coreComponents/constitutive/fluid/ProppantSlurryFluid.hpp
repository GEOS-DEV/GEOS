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
 * @file ProppantSlurryFluid.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_PROPPANTSLURRYFLUID_HPP_
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_PROPPANTSLURRYFLUID_HPP_

#include "constitutive/fluid/SlurryFluidBase.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const proppantSlurryFluid = "ProppantSlurryFluid";
}
}

namespace constitutive
{

class ProppantSlurryFluid : public SlurryFluidBase
{
public:

  ProppantSlurryFluid( std::string const & name, Group * const parent );

  virtual ~ProppantSlurryFluid() override;

  // *** ConstitutiveBase interface

  virtual void DeliverClone( string const & name,
                             Group * const parent,
                             std::unique_ptr< ConstitutiveBase > & clone ) const override;

  static std::string CatalogName() { return dataRepository::keys::proppantSlurryFluid; }

  virtual string GetCatalogName() override { return CatalogName(); }

  virtual void AllocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  // *** ProppantSlurryFluid interface

  virtual void PointUpdate( real64 const & pressure, real64 const & proppantConcentration, arraySlice1d< real64 const > const & Componentconcentration,
                            real64 const & shearRate, integer const & isProppantBoundary, localIndex const k, localIndex const q ) override;

  virtual void PointUpdateFluidProperty( real64 const & pressure, arraySlice1d< real64 const > const & componentConcentration, real64 const & shearRate,
                                         localIndex const k, localIndex const q ) override;

  virtual void PointUpdateComponentDensity( real64 const & pressure, arraySlice1d< real64 const > const & componentConcentration, localIndex const k,
                                            localIndex const q ) override;

  virtual void BatchUpdate( arrayView1d< real64 const > const & pressure, arrayView1d< real64 const > const & proppantConcentration,
                            arrayView2d< real64 const > const & componentConcentration, arrayView1d< real64 const > const & shearRate ) override;

  void ComputeFluidDensity( real64 const & pressure,
                            arraySlice1d< real64 const > const & componentConcentration,
                            arraySlice1d< real64 > const & componentDensity,
                            arraySlice1d< real64 > const & dComponentDensity_dPressure,
                            arraySlice2d< real64 > const & dComponentDensity_dComponentConcentration,
                            real64 & fluidDensity,
                            real64 & dFluidDensity_dPressure,
                            arraySlice1d< real64 > const & dFluidDensity_dComponentConcentration ) const;

  void ComputeComponentDensity( real64 const & pressure,
                                arraySlice1d< real64 const > const & componentConcentration,
                                arraySlice1d< real64 > const & componentDensity,
                                arraySlice1d< real64 > const & dComponentDensity_dPressure,
                                arraySlice2d< real64 > const & dComponentDensity_dComponentConcentration ) const;


  void ComputeFluidViscosity( arraySlice1d< real64 const > const & componentDensity,
                              arraySlice1d< real64 const > const & dComponentDensity_dPressure,
                              arraySlice2d< real64 const > const & dComponentDensity_dComponentConcentration,
                              real64 const & fluidDensity,
                              real64 const & dFluidDensity_dPressure,
                              arraySlice1d< real64 const > const & dFluidDensity_dComponentConcentration,
                              real64 & fluidViscosity,
                              real64 & dFluidViscosity_dPressure,
                              arraySlice1d< real64 > const & dFluidViscosity_dComponentConcentration ) const;


  void Compute( real64 const & proppantConcentration,
                real64 const & fluidDensity,
                real64 const & dFluidDensity_dPressure,
                arraySlice1d< real64 const > const & dFluidDensity_dComponentConcentration,
                real64 const & fluidViscosity,
                real64 const & dFluidViscosity_dPressure,
                arraySlice1d< real64 const > const & dFluidViscosity_dComponentConcentration,
                integer const & isProppantBoundary,
                real64 & density,
                real64 & dDensity_dPressure,
                real64 & dDensity_dProppantConcentration,
                arraySlice1d< real64 > const & dDensity_dComponentConcentration,
                real64 & viscosity,
                real64 & dViscosity_dPressure,
                real64 & dViscosity_dProppantConcentration,
                arraySlice1d< real64 > const & dViscosity_dComponentConcentration ) const;

  // *** Data repository keys

  struct viewKeyStruct : public SlurryFluidBase::viewKeyStruct
  {
    static constexpr auto compressibilityString    = "compressibility";
    static constexpr auto referencePressureString  = "referencePressure";
    static constexpr auto referenceDensityString   = "referenceDensity";

    static constexpr auto referenceProppantDensityString   = "referenceProppantDensity";

    static constexpr auto maxProppantConcentrationString   = "maxProppantConcentration";
    static constexpr auto referenceViscosityString = "referenceViscosity";

    dataRepository::ViewKey compressibility    = { compressibilityString    };
    dataRepository::ViewKey referencePressure  = { referencePressureString  };
    dataRepository::ViewKey referenceDensity   = { referenceDensityString   };
    dataRepository::ViewKey referenceViscosity = { referenceViscosityString };
    dataRepository::ViewKey maxProppantConcentration = { maxProppantConcentrationString };
    dataRepository::ViewKey referenceProppantDensity = { referenceProppantDensityString };

  } viewKeysProppantSlurryFluid;

protected:
  virtual void PostProcessInput() override;

private:

  real64 m_compressibility;

  real64 m_referenceProppantDensity;

  real64 m_referencePressure;

  real64 m_referenceDensity;

  real64 m_referenceViscosity;

  real64 m_maxProppantConcentration;

};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_PROPPANTSLURRYFLUID_HPP_ */
