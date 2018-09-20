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

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_COMPRESSIBLESINGLEPHASEFLUID_HPP_
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_COMPRESSIBLESINGLEPHASEFLUID_HPP_

#include "constitutive/ConstitutiveBase.hpp"

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

class CompressibleSinglePhaseFluid : public ConstitutiveBase
{
public:
  CompressibleSinglePhaseFluid( std::string const & name, ManagedGroup * const parent );

  virtual ~CompressibleSinglePhaseFluid() override;

  std::unique_ptr<ConstitutiveBase> DeliverClone( string const & name,
                                                  ManagedGroup * const parent ) const override;

  virtual void AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;


  static std::string CatalogName() { return dataRepository::keys::compressibleSinglePhaseFluid; }

  virtual string GetCatalogName() override { return CatalogName(); }

  virtual void StateUpdate( dataRepository::ManagedGroup const * const input,
                            dataRepository::ManagedGroup const * const parameters,
                            dataRepository::ManagedGroup * const stateVariables,
                            integer const systemAssembleFlag ) const override final {}

  virtual void FluidDensityCompute( real64 const & pres,
                                    localIndex const i,
                                    real64 & dens,
                                    real64 & dDens_dPres ) const override final;

  virtual void FluidViscosityCompute( real64 const & pres,
                                      localIndex const i,
                                      real64 & visc,
                                      real64 & dVisc_dPres ) const override final;

  virtual void PressureUpdatePoint( real64 const & pres,
                                    localIndex const k,
                                    localIndex const q ) override final;

  virtual void FillDocumentationNode() override;

  virtual void ReadXML_PostProcess() override;

  virtual void FinalInitialization( ManagedGroup * const parent ) override final;

  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    dataRepository::ViewKey compressibility    = { "compressibility"    };
    dataRepository::ViewKey viscosibility      = { "viscosibility"      };
    dataRepository::ViewKey referencePressure  = { "referencePressure"  };
    dataRepository::ViewKey referenceDensity   = { "referenceDensity"   };
    dataRepository::ViewKey referenceViscosity = { "referenceViscosity" };
  } viewKeys;

  array2d<real64> const & density() const { return m_density; }
  array2d<real64>       & density()       { return m_density; }

  array2d<real64> const & dPressure_dDensity() const { return m_dDensity_dPressure; }
  array2d<real64>       & dPressure_dDensity()       { return m_dDensity_dPressure; }

  array2d<real64> const & viscosity() const { return m_viscosity; }
  array2d<real64>       & viscosity()       { return m_viscosity; }

  array2d<real64> const & dViscosity_dDensity() const { return m_dViscosity_dPressure; }
  array2d<real64>       & dViscosity_dDensity()       { return m_dViscosity_dPressure; }


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

  array2d<real64> m_density;
  array2d<real64> m_dDensity_dPressure;

  array2d<real64> m_viscosity;
  array2d<real64> m_dViscosity_dPressure;

  ExponentialRelation<localIndex, real64> m_densityRelation;
  ExponentialRelation<localIndex, real64> m_viscosityRelation;
};


inline void CompressibleSinglePhaseFluid::PressureUpdatePoint(real64 const & pres,
                                                              localIndex const k,
                                                              localIndex const q)
{
  m_densityRelation.Compute( pres, m_density[k][q], m_dDensity_dPressure[k][q] );
  m_viscosityRelation.Compute( pres, m_viscosity[k][q], m_dViscosity_dPressure[k][q] );
}

} /* namespace constitutive */

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_COMPRESSIBLESINGLEPHASEFLUID_HPP_ */
