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

  ParticleFluid( std::string const & name, ManagedGroup * const parent );

  virtual ~ParticleFluid() override;

  // *** ConstitutiveBase interface

  virtual void DeliverClone( string const & name,
                             ManagedGroup * const parent,
                             std::unique_ptr<ConstitutiveBase> & clone ) const override;

  static std::string CatalogName() { return dataRepository::keys::particleFluid; }

  virtual string GetCatalogName() override { return CatalogName(); }

  virtual void AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  // *** ParticleFluid interface

  virtual void PointUpdate(real64 const & concentration, localIndex const k) override;

  virtual void BatchUpdate( arrayView1d<real64 const> const & concentration) override {}

  virtual void PointUpdateMob(real64 const & concentration, real64 const & aperture, localIndex const k) override;

  virtual void BatchUpdateMob( arrayView1d<real64 const> const & concentration, arrayView1d<real64 const> const & aperture) override {}
  
  // *** Data repository keys

  struct viewKeyStruct : public ParticleFluidBase::viewKeyStruct
  {

    static constexpr auto fluidViscosityString    = "fluidViscosity";
    static constexpr auto proppantDiameterString    = "proppantDiameter";    
    static constexpr auto fluidDensityString    = "fluidDensity";
    static constexpr auto proppantDensityString    = "proppantDensity";
    static constexpr auto hinderedSettlingCoefficientString    = "hinderedSettlingCoefficient";
    static constexpr auto collisionAlphaString    = "collisionAlpha";
    static constexpr auto slipConcentrationString    = "slipConcentration";    
    static constexpr auto collisionBetaString    = "collisionBeta";
    static constexpr auto bridgingFactorString    = "bridgingFactor";
    static constexpr auto sphericityString    = "sphericity";        

    dataRepository::ViewKey fluidViscosity    = { fluidViscosityString    };
    dataRepository::ViewKey proppantDiameter    = { proppantDiameterString };
    dataRepository::ViewKey fluidDensity    = { fluidDensityString    };
    dataRepository::ViewKey proppantDensity   = { proppantDensityString };
    dataRepository::ViewKey hinderedSettlingCoefficient  = { hinderedSettlingCoefficientString };
    dataRepository::ViewKey collisionAlpha   = { collisionAlphaString };
    dataRepository::ViewKey slipConcentration = { slipConcentrationString };
    dataRepository::ViewKey collisionBeta   = { collisionBetaString };

    dataRepository::ViewKey bridgingFactor   = { bridgingFactorString };
    dataRepository::ViewKey sphericity   = { sphericityString };        

  } viewKeysParticleFluid;

protected:

  virtual void PostProcessInput() override;

private:

  void Compute( real64 const & concentration,
		real64 & settlingFactor,
		real64 & dSettlingFactor_dConc,
		real64 & collisionFactor,
		real64 & dCollisionFactor_dConc ) const;

  void ComputeMob( real64 const & concentration,
		   real64 const & aperture,		
		   bool & isProppantMobile,
		   real64 & proppantPackPermeability ) const;
  
  
  real64 m_fluidDensity;
  
  real64 m_proppantDensity;  

  real64 m_fluidViscosity;
  
  real64 m_proppantDiameter;  

  real64 m_singleParticleSettlingVelocity;

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
