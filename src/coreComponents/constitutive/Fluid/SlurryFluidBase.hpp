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
  * @file SlurryFluidBase.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_SLURRYFLUIDBASE_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_SLURRYFLUIDBASE_HPP

#include "constitutive/ConstitutiveBase.hpp"

namespace geosx
{

namespace constitutive
{

class SlurryFluidBase : public ConstitutiveBase
{
public:

  SlurryFluidBase( std::string const & name, ManagedGroup * const parent );

  virtual ~SlurryFluidBase() override;

  // *** ConstitutiveBase interface

  virtual void DeliverClone( string const & name,
                             ManagedGroup * const parent,
                             std::unique_ptr<ConstitutiveBase> & clone ) const override = 0;

  virtual void AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  // *** SlurryFluidBase-specific interface

  virtual void PointUpdate( real64 const & pressure, real64 const & concentration, localIndex const k, localIndex const q ) = 0;

  virtual void BatchUpdate( arrayView1d<real64 const> const & pressure,  arrayView1d<real64 const> const & concentration ) = 0;

  virtual void Compute( real64 const & pressure,
			real64 const & concentration,			
                        real64 & density,
                        real64 & fluidDensity,			
                        real64 & dDensity_dPressure,
                        real64 & dDensity_dConcentration,
			real64 & dFluidDensity_dPressure,
                        real64 & viscosity,
                        real64 & dViscosity_dPressure,
                        real64 & dViscosity_dConcentration ) const = 0;
  

  array2d<real64> const & density() const { return m_density; }
  array2d<real64>       & density()       { return m_density; }

  array2d<real64> const & fluidDensity() const { return m_fluidDensity; }
  array2d<real64>       & fluidDensity()       { return m_fluidDensity; }  

  array2d<real64> const & dDensity_dPressure() const { return m_dDens_dPres; }
  array2d<real64>       & dDensity_dPressure()       { return m_dDens_dPres; }

  array2d<real64> const & dDensity_dConcentration() const { return m_dDens_dConc; }
  array2d<real64>       & dDensity_dConcentration()       { return m_dDens_dConc; }  

  array2d<real64> const & viscosity() const { return m_viscosity; }
  array2d<real64>       & viscosity()       { return m_viscosity; }

  array2d<real64> const & dViscosity_dPressure() const { return m_dVisc_dPres; }
  array2d<real64>       & dViscosity_dPressure()       { return m_dVisc_dPres; }

  array2d<real64> const & dViscosity_dConcentration() const { return m_dVisc_dConc; }
  array2d<real64>       & dViscosity_dConcentration()       { return m_dVisc_dConc; }  

  // *** Data repository keys

  struct viewKeyStruct
  {
    static constexpr auto defaultDensityString      = "defaultDensity";
    static constexpr auto densityString      = "density";
    static constexpr auto fluidDensityString      = "FluidDensity";    
    static constexpr auto dDens_dPresString  = "dDensity_dPres";
    static constexpr auto dDens_dConcString  = "dDensity_dConc";
    static constexpr auto dFluidDens_dPresString  = "dFluidDensity_dPres";    

    static constexpr auto defaultViscosityString    = "defaultViscosity";
    static constexpr auto viscosityString    = "viscosity";
    static constexpr auto dVisc_dPresString  = "dVisc_dPres";
    static constexpr auto dVisc_dConcString  = "dVisc_dConc";    

    using ViewKey = dataRepository::ViewKey;

    ViewKey density     = { densityString };
    ViewKey fluidDensity = { fluidDensityString };    
    ViewKey dDens_dPres = { dDens_dPresString };
    ViewKey dDens_dConc = { dDens_dConcString };
    ViewKey dFluidDens_dPres = { dFluidDens_dPresString };    

    ViewKey viscosity   = { viscosityString };
    ViewKey dVisc_dPres = { dVisc_dPresString };
    ViewKey dVisc_dConc = { dVisc_dConcString };    

  } viewKeysSlurryFluidBase;

protected:

  virtual void PostProcessInput() override;

  real64 m_defaultDensity;
  real64 m_defaultViscosity;

  array2d<real64> m_density;
  array2d<real64> m_fluidDensity;  
  array2d<real64> m_dDens_dPres;
  array2d<real64> m_dDens_dConc;  
  array2d<real64> m_dFluidDens_dPres;  
  
  array2d<real64> m_viscosity;
  array2d<real64> m_dVisc_dPres;
  array2d<real64> m_dVisc_dConc;  

};

} //namespace constitutive

} //namespace geosx

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_SLURRYFLUIDBASE_HPP
