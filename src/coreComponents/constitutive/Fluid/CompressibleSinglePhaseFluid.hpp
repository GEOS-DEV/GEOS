/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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
  * @file CompressibleSinglePhaseFluid.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_COMPRESSIBLESINGLEPHASEFLUID_HPP_
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_COMPRESSIBLESINGLEPHASEFLUID_HPP_

#include "constitutive/Fluid/SingleFluidBase.hpp"

#include "constitutive/ExponentialRelation.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const compressibleSinglePhaseFluid = "CompressibleSinglePhaseFluid";
}
}

namespace constitutive
{

class CompressibleSinglePhaseFluid : public SingleFluidBase
{
public:

  CompressibleSinglePhaseFluid( std::string const & name, ManagedGroup * const parent );

  virtual ~CompressibleSinglePhaseFluid() override;

  std::unique_ptr<ConstitutiveBase> DeliverClone( string const & name,
                                                  ManagedGroup * const parent ) const override;

  static std::string CatalogName() { return dataRepository::keys::compressibleSinglePhaseFluid; }

  virtual string GetCatalogName() override { return CatalogName(); }

  virtual void ProcessInputFile_PostProcess() override;

  // *** SingleFluid-specific interface

  virtual void PointUpdate( real64 const & pressure, localIndex const k, localIndex const q ) override;

  virtual void BatchUpdate( arrayView1d<real64 const> const & pressure ) override;

  inline static void Compute( real64 const & pressure,
                              real64 & density,
                              real64 & dDensity_dPressure,
                              real64 & viscosity,
                              real64 & dViscosity_dPressure,
                              ExponentialRelation<real64> densityRelation,
                              ExponentialRelation<real64> viscosityRelation );

  virtual void Compute( real64 const & pressure,
                        real64 & density,
                        real64 & dDensity_dPressure,
                        real64 & viscosity,
                        real64 & dViscosity_dPressure ) override;

  struct viewKeyStruct : public SingleFluidBase::viewKeyStruct
  {
    static constexpr auto compressibilityString    = "compressibility";
    static constexpr auto viscosibilityString      = "viscosibility";
    static constexpr auto referencePressureString  = "referencePressure";
    static constexpr auto referenceDensityString   = "referenceDensity";
    static constexpr auto referenceViscosityString = "referenceViscosity";

    dataRepository::ViewKey compressibility    = { compressibilityString    };
    dataRepository::ViewKey viscosibility      = { viscosibilityString      };
    dataRepository::ViewKey referencePressure  = { referencePressureString  };
    dataRepository::ViewKey referenceDensity   = { referenceDensityString   };
    dataRepository::ViewKey referenceViscosity = { referenceViscosityString };

  } viewKeysCompressibleSinglePhaseFluid;


private:

  /// scalar fluid bulk modulus parameter
  real64 m_compressibility;

  /// scalar fluid viscosity exponential coefficient
  real64 m_viscosibility;

  /// reference pressure parameter
  real64 m_referencePressure;

  /// reference density parameter
  real64 m_referenceDensity;

  /// reference viscosity parameter
  real64 m_referenceViscosity;

  ExponentialRelation<real64> m_densityRelation;
  ExponentialRelation<real64> m_viscosityRelation;
};


inline
void CompressibleSinglePhaseFluid::Compute( real64 const & pressure,
                                            real64 & density, real64 & dDensity_dPressure,
                                            real64 & viscosity, real64 & dViscosity_dPressure,
                                            ExponentialRelation<real64> densityRelation,
                                            ExponentialRelation<real64> viscosityRelation )
{
  densityRelation.Compute( pressure, density, dDensity_dPressure );
  viscosityRelation.Compute( pressure, viscosity, dViscosity_dPressure );
}

} /* namespace constitutive */

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_COMPRESSIBLESINGLEPHASEFLUID_HPP_ */
