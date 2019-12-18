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
  * @file ParticleFluid.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_PARTICLEFLUID_HPP_
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_PARTICLEFLUID_HPP_

#include "constitutive/fluid/ParticleFluidBase.hpp"

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
    Intermediate,
    Turbulence
  };

  static ParticleSettlingModel stringToParticleSettlingModel( string const & str );
  
  ParticleFluid( std::string const & name, Group * const parent );

  virtual ~ParticleFluid() override;

  // *** ConstitutiveBase interface

  virtual void DeliverClone( string const & name,
                             Group * const parent,
                             std::unique_ptr<ConstitutiveBase> & clone ) const override;

  static std::string CatalogName() { return dataRepository::keys::particleFluid; }

  virtual string GetCatalogName() override { return CatalogName(); }

  virtual void AllocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  // *** ParticleFluid interface

  virtual void PointUpdate(localIndex const NC, real64 const & proppantConcentration, arraySlice1d<real64 const> const & componentConcentration, arraySlice1d<real64 const> const & nIndex, arraySlice1d<real64 const> const & KIndex, real64 const &fluidDensity, real64 const &dFluidDensity_dPressure, arraySlice1d<real64 const> const &dFluidDensity_dComponentConcentration, localIndex const k) override; 

  virtual void PointUpdate(localIndex const NC, real64 const & proppantConcentration, real64 const &fluidDensity, real64 const &dFluidDensity_dPressure, arraySlice1d<real64 const> const &dFluidDensity_dComponentConcentration, real64 const &fluidViscosity, real64 const &dFluidViscosity_dPressure, arraySlice1d<real64 const> const &dFluidViscosity_dComponentConcentration, localIndex const k) override;  
  
  virtual void BatchUpdate( arrayView1d<real64 const> const & concentration) override;

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

    dataRepository::ViewKey sphericity   = { sphericityString };

    dataRepository::ViewKey particleSettlingModel   = { particleSettlingModelString };            

  } viewKeysParticleFluid;

protected:

  virtual void PostProcessInput() override;

private:

  void Compute( localIndex const NC,
                real64 const & proppantConcentration,
                real64 const & fluidDensity,
                real64 const & dFluidDensity_dPressure,
                arraySlice1d<real64 const> const & dFluidDensity_dComponentConcentration,
                real64 const & fluidViscosity,
                real64 const & dFluidViscosity_dPressure,
                arraySlice1d<real64 const> const & dFluidViscosity_dComponentConcentration,
                real64 & settlingFactor,
                real64 & dSettlingFactor_dPressure,
                real64 & dSettlingFactor_dProppantConcentration,
                arraySlice1d<real64> const & dSettlingFactor_dComponentConcentration,
                real64 & collisionFactor,
                real64 & dCollisionFactor_dProppantConcentration ) const;

  string m_particleSettlingModelString;

  ParticleSettlingModel m_particleSettlingModel;
  
  real64 m_proppantDensity;  

  real64 m_fluidViscosity;
  
  real64 m_proppantDiameter;  

  real64 m_hinderedSettlingCoefficient;

  real64 m_collisionAlpha;

  real64 m_slipConcentration;

  real64 m_collisionBeta;

  real64 m_sphericity;

  real64 m_packPermeabilityCoef;

};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_PARTICLEFLUID_HPP_ */
