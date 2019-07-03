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
  * @file ParticleFluidBase.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_PARTICLEFLUIDBASE_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_PARTICLEFLUIDBASE_HPP

#include "constitutive/ConstitutiveBase.hpp"

namespace geosx
{

namespace constitutive
{

class ParticleFluidBase : public ConstitutiveBase
{
public:

  ParticleFluidBase( std::string const & name, ManagedGroup * const parent );

  virtual ~ParticleFluidBase() override;

  // *** ConstitutiveBase interface

  virtual void DeliverClone( string const & name,
                             ManagedGroup * const parent,
                             std::unique_ptr<ConstitutiveBase> & clone ) const override = 0;

  virtual void AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  // *** ParticleFluidBase-specific interface

  virtual void PointUpdate(real64 const & concentration, localIndex const k) = 0;

  virtual void BatchUpdate( arrayView1d<real64 const> const & concentration ) = 0;

  virtual void PointUpdateMob(real64 const & concentration, real64 const &apeture, localIndex const k) = 0;

  virtual void BatchUpdateMob( arrayView1d<real64 const> const & concentration, arrayView1d<real64 const> const &aperture) = 0;
  
 
  // *** Data repository keys

  struct viewKeyStruct
  {

    static constexpr auto settlingFactorString    = "settlingFactor";    
    static constexpr auto dSettlingFactor_dConcString  = "dSettlingFactor_dConc";

    static constexpr auto collisionFactorString    = "collisionFactor";    
    static constexpr auto dCollisionFactor_dConcString  = "dCollisionFactor_dConc";        

    static constexpr auto maxProppantConcentrationString    = "maxProppantConcentration";

    static constexpr auto isProppantMobileString    = "isProppantMobile";
    
    static constexpr auto proppantPackPermeabilityString    = "proppantPackPermeability";
    
    using ViewKey = dataRepository::ViewKey;

    ViewKey settlingFactor  = { settlingFactorString };
    ViewKey dSettlingFactor_dConc = { dSettlingFactor_dConcString };

    ViewKey collisionFactor  = { collisionFactorString };
    ViewKey collisionFactor_dConc = { dCollisionFactor_dConcString };

    ViewKey maxProppantConcentration   = { maxProppantConcentrationString };

    ViewKey isProppantMobile   = { isProppantMobileString };
    ViewKey proppantPackPermeability   = { proppantPackPermeabilityString };        

  } viewKeysParticleFluidBase;

protected:

  virtual void PostProcessInput() override;

  array1d<real64> m_settlingFactor;
  array1d<real64> m_dSettlingFactor_dConc;

  array1d<real64> m_collisionFactor;
  array1d<real64> m_dCollisionFactor_dConc;   

  array1d<bool> m_isProppantMobile;
  array1d<real64> m_proppantPackPermeability;  
  
  real64 m_maxProppantConcentration;    
  
};

} //namespace constitutive

} //namespace geosx

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_PARTICLEFLUIDBASE_HPP
