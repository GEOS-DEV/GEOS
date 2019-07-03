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
  * @file ProppantSlurryFluid.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_PROPPANTSLURRYFLUID_HPP_
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_PROPPANTSLURRYFLUID_HPP_

#include "constitutive/Fluid/SlurryFluidBase.hpp"

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

  ProppantSlurryFluid( std::string const & name, ManagedGroup * const parent );

  virtual ~ProppantSlurryFluid() override;

  // *** ConstitutiveBase interface

  virtual void DeliverClone( string const & name,
                             ManagedGroup * const parent,
                             std::unique_ptr<ConstitutiveBase> & clone ) const override;

  static std::string CatalogName() { return dataRepository::keys::proppantSlurryFluid; }

  virtual string GetCatalogName() override { return CatalogName(); }

  virtual void AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  // *** ProppantSlurryFluid interface

  virtual void PointUpdate( real64 const & pressure, real64 const & concentration, localIndex const k, localIndex const q ) override;

  virtual void BatchUpdate( arrayView1d<real64 const> const & pressure, arrayView1d<real64 const> const & concentration) override {}

  virtual void Compute( real64 const & pressure,
			real64 const & concentration,			
                        real64 & density,
                        real64 & fluidDensity,			
                        real64 & dDensity_dPressure,
                        real64 & dDensity_dConcentration,
                        real64 & dFluidDensity_dPressure,
                        real64 & viscosity,
                        real64 & dViscosity_dPressure,			
                        real64 & dViscosity_dConcentration ) const override;

  // *** Data repository keys

  struct viewKeyStruct : public SlurryFluidBase::viewKeyStruct
  {
    static constexpr auto compressibilityString    = "compressibility";
    static constexpr auto referencePressureString  = "referencePressure";
    static constexpr auto referenceDensityString   = "referenceDensity";
    static constexpr auto referenceProppantDensityString   = "referenceProppantDensity";
    static constexpr auto referenceFluidDensityString   = "referenceFluidDensity";    
    static constexpr auto maxProppantVolumeFractionString   = "maxProppantVolumeFraction";    
    static constexpr auto referenceViscosityString = "referenceViscosity";
    static constexpr auto referenceProppantVolumeFractionString   = "referenceProppantVolumeFraction";        

    dataRepository::ViewKey compressibility    = { compressibilityString    };
    dataRepository::ViewKey referencePressure  = { referencePressureString  };
    dataRepository::ViewKey referenceDensity   = { referenceDensityString   };
    dataRepository::ViewKey referenceViscosity = { referenceViscosityString };
    dataRepository::ViewKey maxProppantVolumeFraction = { maxProppantVolumeFractionString };    
    dataRepository::ViewKey referenceProppantDensity = { referenceProppantDensityString };
    dataRepository::ViewKey referenceFluidDensity = { referenceFluidDensityString };
    dataRepository::ViewKey referenceProppantVolumeFraction = { referenceProppantVolumeFractionString };        
    
  } viewKeysProppantSlurryFluid;

protected:
  virtual void PostProcessInput() override;

private:

  real64 m_compressibility;

  real64 m_referenceProppantDensity;

  real64 m_referenceFluidDensity;  

  real64 m_referencePressure;

  real64 m_referenceDensity;

  real64 m_referenceViscosity;

  real64 m_maxProppantVolumeFraction;

  real64 m_referenceProppantVolumeFraction;    

};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_PROPPANTSLURRYFLUID_HPP_ */
