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

/*
 * HypoElasticLinear.hpp
 *
 *  Created on: Aug 4, 2016
 *      Author: rrsettgast
 */

#ifndef LINEARELASTICISOTROPIC_HPP_
#define LINEARELASTICISOTROPIC_HPP_
#include "ConstitutiveBase.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const linearElasticIsotropic = "LinearElasticIsotropic";
}
}

namespace constitutive
{

class LinearElasticIsotropic : public ConstitutiveBase
{
public:
  LinearElasticIsotropic( std::string const & name, ManagedGroup * const parent );

  virtual ~LinearElasticIsotropic() override;

  virtual std::unique_ptr<ConstitutiveBase>
  DeliverClone( string const & name, ManagedGroup * const parent ) const override;

  virtual void AllocateMaterialData( dataRepository::ManagedGroup * const parent,
                                     localIndex const numConstitutivePointsPerParentIndex ) override;

  static std::string CatalogName() { return dataRepository::keys::linearElasticIsotropic; }
  virtual string GetCatalogName() override { return CatalogName(); }


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


  virtual void FillDocumentationNode() override;

  virtual void ReadXML_PostProcess() override;

  void GetStiffness( realT c[6][6]) const override;

  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr auto bulkModulus0String  = "BulkModulus0";
    static constexpr auto shearModulus0String = "ShearModulus0";
    static constexpr auto bulkModulusString  = "BulkModulus";
    static constexpr auto shearModulusString = "ShearModulus";

    static constexpr auto deviatorStressString = "DeviatorStress";
    static constexpr auto meanStressString = "MeanStress";


    dataRepository::ViewKey youngsModulus = { "YoungsModulus" };
    dataRepository::ViewKey bulkModulus = { bulkModulusString };
    dataRepository::ViewKey shearModulus = { shearModulusString };
    dataRepository::ViewKey poissonRatio = { "PoissonRatio" };
    dataRepository::ViewKey density = { "Density" };

    dataRepository::ViewKey deviatorStress = { deviatorStressString };
    dataRepository::ViewKey meanStress = { meanStressString };
  } m_linearElasticIsotropicViewKeys;

  struct groupKeyStruct : public ConstitutiveBase::groupKeyStruct
  {
  } m_linearElasticIsotropicGroupKeys;

  virtual viewKeyStruct       & viewKeys()       override { return m_linearElasticIsotropicViewKeys; }
  virtual viewKeyStruct const & viewKeys() const override { return m_linearElasticIsotropicViewKeys; }

  virtual groupKeyStruct       & groupKeys()      override { return m_linearElasticIsotropicGroupKeys; }
  virtual groupKeyStruct const & groupKeys() const override{ return m_linearElasticIsotropicGroupKeys; }



  dataRepository::view_rtype<real64>       youngsModulus()       { return GetParameterData()->getData<real64>(viewKeys().youngsModulus); }
  dataRepository::view_rtype_const<real64> youngsModulus() const { return GetParameterData()->getData<real64>(viewKeys().youngsModulus); }

  dataRepository::view_rtype<real64>       bulkModulus()       { return GetParameterData()->getData<real64>(viewKeys().bulkModulus); }
  dataRepository::view_rtype_const<real64> bulkModulus() const { return GetParameterData()->getData<real64>(viewKeys().bulkModulus); }

  dataRepository::view_rtype<real64>       shearModulus()       { return GetParameterData()->getData<real64>(viewKeys().shearModulus); }
  dataRepository::view_rtype_const<real64> shearModulus() const { return GetParameterData()->getData<real64>(viewKeys().shearModulus); }

  dataRepository::view_rtype<real64>       poissonRatio()       { return GetParameterData()->getData<real64>(viewKeys().poissonRatio); }
  dataRepository::view_rtype_const<real64> poissonRatio() const { return GetParameterData()->getData<real64>(viewKeys().poissonRatio); }

  dataRepository::view_rtype<real64>       density()       { return GetParameterData()->getData<real64>(viewKeys().density); }
  dataRepository::view_rtype_const<real64> density() const { return GetParameterData()->getData<real64>(viewKeys().density); }

  dataRepository::view_rtype<r2Sym_array>       deviatorStress()       { return GetStateData()->getData<r2Sym_array>(viewKeys().deviatorStress); }
  dataRepository::view_rtype_const<r2Sym_array> deviatorStress() const { return GetStateData()->getData<r2Sym_array>(viewKeys().deviatorStress); }

  dataRepository::view_rtype<real64_array>       meanStress()       { return GetStateData()->getData<real64_array>(viewKeys().meanStress); }
  dataRepository::view_rtype_const<real64_array> meanStress() const { return GetStateData()->getData<real64_array>(viewKeys().meanStress); }

  struct dataPointers
  {
    real64 * m_bulkModulus = nullptr;
    real64 * m_shearModulus = nullptr;
    R2SymTensor * m_deviatorStress = nullptr;
    real64 * m_meanStress = nullptr;
  } m_dataPointers;

private:
  real64 m_bulkModulus0;
  real64 m_shearModulus0;
  real64_array m_bulkModulus;
  real64_array m_shearModulus;
  real64_array m_meanStress;
  r2Sym_array m_deviatorStress;

};



}

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_HYPOELASTICLINEAR_HPP_ */
