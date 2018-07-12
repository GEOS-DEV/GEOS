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

  virtual ~LinearEOS();

  static std::string CatalogName() { return dataRepository::keys::linearEOS; }


  virtual void SetParamStatePointers( void *& ) override final;

  virtual void StateUpdate( dataRepository::ManagedGroup const * const input,
                            dataRepository::ManagedGroup const * const parameters,
                            dataRepository::ManagedGroup * const stateVariables,
                            integer const systemAssembleFlag ) const override;

  virtual UpdateFunctionPointer GetStateUpdateFunctionPointer() override final;

  R2SymTensor  StateUpdatePoint( R2SymTensor const & D,
                                 R2Tensor const & Rot,
                                 localIndex const i,
                                 integer const systemAssembleFlag ) override;

  void EquationOfStatePressureUpdate( real64 const & dRho,
                         localIndex const i,
                         real64 & P,
                         real64 & dPdRho ) override;

  void EquationOfStateDensityUpdate( real64 const & dP,
                                     localIndex const i,
                                     real64 & dRho,
                                     real64 & dRho_dP ) override;

  virtual void FillDocumentationNode() override;

  virtual void ReadXML_PostProcess() override;

  void GetStiffness( realT c[6][6]) const override;

  struct ViewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    dataRepository::ViewKey bulkModulus = { "BulkModulus" };
    dataRepository::ViewKey referenceDensity = { "referenceDensity" };
    dataRepository::ViewKey referencePressure = { "referencePressure" };
    dataRepository::ViewKey fluidViscosity = { "fluidViscosity" };

    dataRepository::ViewKey fluidPressure = { "fluidPressure" };
    dataRepository::ViewKey fluidDensity = { "fluidDensity" };
  } viewKeys;




  dataRepository::view_rtype<real64>       bulkModulus()       { return GetParameterData()->getData<real64>(viewKeys.bulkModulus); }
  dataRepository::view_rtype_const<real64> bulkModulus() const { return GetParameterData()->getData<real64>(viewKeys.bulkModulus); }

  dataRepository::view_rtype<real64>       density()       { return GetParameterData()->getData<real64>(viewKeys.referenceDensity); }
  dataRepository::view_rtype_const<real64> density() const { return GetParameterData()->getData<real64>(viewKeys.referenceDensity); }


private:

  /// scalar bulk modulus parameter
  real64 m_bulkModulus;

  /// reference pressure parameter for EOS relation
  real64 m_referencePressure;

  /// reference density parameter for EOS relation
  real64 m_referenceDensity;

  /// fluid viscosity parameter
  real64 m_fluidViscosity;

  /// fluid pressure state variable
  array<real64> m_fluidPressure;

  /// fluid density state variable
  array<real64> m_fluidDensity;
};



}

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_HYPOELASTICLINEAR_HPP_ */
