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

#ifndef LINEAREOS_HPP_
#define LINEAREOS_HPP_

#include "constitutive/ConstitutiveBase.hpp"

#include "constitutive/ExponentialRelation.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const linearEOS = "LinearEOS";
}
}

namespace constitutive
{

class LinearEOS : public ConstitutiveBase
{
public:
  LinearEOS( std::string const & name, ManagedGroup * const parent );

  virtual ~LinearEOS() override;

  static std::string CatalogName() { return dataRepository::keys::linearEOS; }


  virtual void SetParamStatePointers( void *& ) override final {}

  virtual void StateUpdate( dataRepository::ManagedGroup const * const input,
                            dataRepository::ManagedGroup const * const parameters,
                            dataRepository::ManagedGroup * const stateVariables,
                            integer const systemAssembleFlag ) const override final {}

  virtual UpdateFunctionPointer GetStateUpdateFunctionPointer() override final;

  virtual void FluidPressureUpdate(real64 const &dens,
                                   localIndex const i,
                                   real64 &pres,
                                   real64 &dPres_dDens) override final;

  virtual void FluidDensityUpdate(real64 const &pres,
                                  localIndex const i,
                                  real64 &dens,
                                  real64 &dDens_dPres) override final;

  virtual void FluidViscosityUpdate(real64 const &pres,
                                    localIndex const i,
                                    real64 &visc,
                                    real64 &dVisc_dPres) override final;

  virtual void SimplePorosityUpdate(real64 const &pres,
                                    real64 const &poro_ref,
                                    localIndex const i,
                                    real64 &poro,
                                    real64 &dPoro_dPres) override final;

  virtual void FillDocumentationNode() override;

  virtual void ReadXML_PostProcess() override;

  virtual void FinalInitialization(ManagedGroup * const parent) override final;

  void GetStiffness( realT c[6][6]) const override;

  struct ViewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    dataRepository::ViewKey fluidBulkModulus   = { "fluidBulkModulus"   };
    dataRepository::ViewKey solidBulkModulus   = { "solidBulkModulus"   };
    dataRepository::ViewKey fluidViscosibility = { "fluidViscosibility" };

    dataRepository::ViewKey referencePressure  = { "referencePressure"  };
    dataRepository::ViewKey referenceDensity   = { "referenceDensity"   };
    dataRepository::ViewKey referenceViscosity = { "referenceViscosity" };
  } viewKeys;


private:

  /// scalar fluid bulk modulus parameter
  real64 m_fluidBulkModulus;

  /// scalar fluid bulk modulus parameter
  real64 m_solidBulkModulus;

  /// scalar fluid viscosity exponential coefficient
  real64 m_fluidViscosibility;

  /// reference pressure parameter for EOS relation
  real64 m_referencePressure;

  /// reference density parameter for EOS relation
  real64 m_referenceDensity;

  /// reference viscosity parameter for EOS relation
  real64 m_referenceViscosity;

  ExponentialRelation<localIndex, real64> m_densityRelation;
  ExponentialRelation<localIndex, real64> m_viscosityRelation;
  ExponentialRelation<localIndex, real64> m_porosityRelation;
};



}

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_HYPOELASTICLINEAR_HPP_ */
