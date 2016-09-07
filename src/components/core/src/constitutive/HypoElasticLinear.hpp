/*
 * HypoElasticLinear.hpp
 *
 *  Created on: Aug 4, 2016
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_HYPOELASTICLINEAR_HPP_
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_HYPOELASTICLINEAR_HPP_
#include "ConstitutiveBase.hpp"

namespace geosx
{
namespace constitutive
{

class HypoElasticLinear : public ConstitutiveBase
{
public:
  HypoElasticLinear();
  virtual ~HypoElasticLinear();

  static std::string CatalogName() { return "HypoElasticLinear"; }


  virtual void Update( dataRepository::ManagedGroup * const parameters,
                       dataRepository::ManagedGroup * const stateVariables ) override;
};
}

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_HYPOELASTICLINEAR_HPP_ */
