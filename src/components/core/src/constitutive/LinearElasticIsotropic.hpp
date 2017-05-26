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
string const density("density");
string const bulkModulus("bulkModulus");
string const youngsModulus("youngsModulus");
string const poissonRatio("poissonsRatio");
string const shearModulus("shearModulus");
string const deviatorStress = "DeviatorStress";
string const meanStress = "MeanStress";
}
}

namespace constitutive
{

class LinearElasticIsotropic : public ConstitutiveBase
{
public:
  LinearElasticIsotropic( std::string const & name, ManagedGroup * const parent );

  virtual ~LinearElasticIsotropic();

  static std::string CatalogName() { return dataRepository::keys::linearElasticIsotropic; }

  virtual void StateUpdate( dataRepository::ManagedGroup const * const input,
                            dataRepository::ManagedGroup const * const parameters,
                            dataRepository::ManagedGroup * const stateVariables,
                            integer const systemAssembleFlag ) const override;

  R2SymTensor  StateUpdatePoint( R2SymTensor const & D,
                                         R2Tensor const & Rot,
                                         int32 const i,
                                         integer const systemAssembleFlag ) override;

  virtual void FillDocumentationNode( dataRepository::ManagedGroup * const group ) override;

  virtual void ReadXML_PostProcess() override;

private:
  real64 const & m_YoungsModulus;
  real64 const & m_BulkModulus;
  real64 const & m_ShearModulus;
  real64 const & m_PoissonRatio;
  real64 const & m_Density;

  r2Sym_array & m_devStress;
  real64_array & m_meanStress;

};





}

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_HYPOELASTICLINEAR_HPP_ */
