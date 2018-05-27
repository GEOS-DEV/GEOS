/*
 * HypoElasticLinear.hpp
 *
 *  Created on: Aug 4, 2016
 *      Author: rrsettgast
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
string const lineaerEOS = "LinearEOS";
}
}

namespace constitutive
{

class LinearEOS : public ConstitutiveBase
{
public:
  LinearEOS( std::string const & name, ManagedGroup * const parent );

  virtual ~LinearEOS();

  static std::string CatalogName() { return dataRepository::keys::lineaerEOS; }


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

  struct ViewKeyStruct : public ConstitutiveBase::ViewKeyStruct
  {
    dataRepository::ViewKey bulkModulus = { "BulkModulus" };
    dataRepository::ViewKey referenceDensity = { "referenceDensity" };

    dataRepository::ViewKey meanStress = { "Pressure" };
    dataRepository::ViewKey density = { "Density" };
  } viewKeys;




  dataRepository::view_rtype<real64>       bulkModulus()       { return GetParameterData()->getData<real64>(viewKeys.bulkModulus); }
  dataRepository::view_rtype_const<real64> bulkModulus() const { return GetParameterData()->getData<real64>(viewKeys.bulkModulus); }

  dataRepository::view_rtype<real64>       density()       { return GetParameterData()->getData<real64>(viewKeys.referenceDensity); }
  dataRepository::view_rtype_const<real64> density() const { return GetParameterData()->getData<real64>(viewKeys.referenceDensity); }


  struct dataPointers
  {
    real64 * m_bulkModulus = nullptr;
  } m_dataPointers;

private:


};



}

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_HYPOELASTICLINEAR_HPP_ */
