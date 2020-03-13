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
 * @file ParticleFluid.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_PARTICLEFLUID_HPP_
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_PARTICLEFLUID_HPP_

#include "constitutive/Fluid/ParticleFluidBase.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const particleFluid = "ParticleFluid";
}
}

namespace constitutive
{

class ParticleFluid : public ParticleFluidBase
{
public:

  enum class ParticleSettlingModel
  {
    Stokes,
    Intermediate
  };

  static ParticleSettlingModel stringToParticleSettlingModel( string const & str );

  ParticleFluid( std::string const & name, Group * const parent );

  virtual ~ParticleFluid() override;

  // *** ConstitutiveBase interface

  virtual void DeliverClone( string const & name,
                             Group * const parent,
                             std::unique_ptr< ConstitutiveBase > & clone ) const override;

  static std::string CatalogName() { return dataRepository::keys::particleFluid; }

  virtual string GetCatalogName() override { return CatalogName(); }

  virtual void AllocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  // *** ParticleFluid interface

  virtual void PointUpdate( localIndex const NC, real64 const & proppantConcentration, arraySlice1d< real64 const > const & componentConcentration,
                            arraySlice1d< real64 const > const & nIndex, arraySlice1d< real64 const > const & KIndex, real64 const & fluidDensity,
                            real64 const & dFluidDensity_dPressure, arraySlice1d< real64 const > const & dFluidDensity_dComponentConcentration,
                            localIndex const k ) override;

  virtual void BatchUpdate( arrayView1d< real64 const > const & concentration ) override;

  virtual void PointUpdateMob( real64 const & concentration, real64 const & aperture, localIndex const k ) override;

  virtual void BatchUpdateMob( arrayView1d< real64 const > const & concentration, arrayView1d< real64 const > const & aperture ) override;

  // *** Data repository keys

  struct viewKeyStruct : public ParticleFluidBase::viewKeyStruct
  {

    static constexpr auto fluidViscosityString    = "fluidViscosity";
    static constexpr auto proppantDiameterString    = "proppantDiameter";
    static constexpr auto proppantDensityString    = "proppantDensity";
    static constexpr auto hinderedSettlingCoefficientString    = "hinderedSettlingCoefficient";
    static constexpr auto collisionAlphaString    = "collisionAlpha";
    static constexpr auto slipConcentrationString    = "slipConcentration";
    static constexpr auto collisionBetaString    = "collisionBeta";
    static constexpr auto bridgingFactorString    = "bridgingFactor";
    static constexpr auto sphericityString    = "sphericity";

    static constexpr auto particleSettlingModelString    = "particleSettlingModel";

    dataRepository::ViewKey fluidViscosity    = { fluidViscosityString    };
    dataRepository::ViewKey proppantDiameter    = { proppantDiameterString };
    dataRepository::ViewKey proppantDensity   = { proppantDensityString };
    dataRepository::ViewKey hinderedSettlingCoefficient  = { hinderedSettlingCoefficientString };
    dataRepository::ViewKey collisionAlpha   = { collisionAlphaString };
    dataRepository::ViewKey slipConcentration = { slipConcentrationString };
    dataRepository::ViewKey collisionBeta   = { collisionBetaString };

    dataRepository::ViewKey bridgingFactor   = { bridgingFactorString };
    dataRepository::ViewKey sphericity   = { sphericityString };

    dataRepository::ViewKey particleSettlingModel   = { particleSettlingModelString };

  } viewKeysParticleFluid;

protected:

  virtual void PostProcessInput() override;

private:

  void Compute( localIndex const NC,
                real64 const & proppantConcentration,
                arraySlice1d< real64 const > const & componentConcentration,
                arraySlice1d< real64 const > const & nIndex,
                arraySlice1d< real64 const > const & KIndex,
                real64 const & fluidDensity,
                real64 const & dFluidDensity_dPressure,
                arraySlice1d< real64 const > const & dFluidDensity_dComponentConcentration,
                real64 & settlingFactor,
                real64 & dSettlingFactor_dPressure,
                real64 & dSettlingFactor_dProppantConcentration,
                arraySlice1d< real64 > const & dSettlingFactor_dComponentConcentration,
                real64 & collisionFactor,
                real64 & dCollisionFactor_dProppantConcentration ) const;

  void ComputeMob( real64 const & concentration,
                   real64 const & aperture,
                   integer & isProppantMobile,
                   real64 & proppantPackPermeability ) const;


  string m_particleSettlingModelString;

  ParticleSettlingModel m_particleSettlingModel;

  real64 m_proppantDensity;

  real64 m_fluidViscosity;

  real64 m_proppantDiameter;

  real64 m_hinderedSettlingCoefficient;

  real64 m_collisionAlpha;

  real64 m_slipConcentration;

  real64 m_collisionBeta;

  real64 m_bridgingFactor;

  real64 m_sphericity;

  real64 m_packPermeabilityCoef;

  real64 m_bridgingAperture;

};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_PARTICLEFLUID_HPP_ */
