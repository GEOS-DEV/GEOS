/*
 * ConstitutiveBase.hpp
 *
 *  Created on: Jul 28, 2016
 *      Author: rrsettgast
 */

#ifndef COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEBASE_HPP_
#define COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEBASE_HPP_

#include "dataRepository/WrapperCollection.hpp""

namespace geosx
{

namespace constitutive
{

class ConstitutiveBase
{
public:
  ConstitutiveBase();
  virtual ~ConstitutiveBase();

  virtual void Update( dataRepository::WrapperCollection * const parameters,
                       dataRepository::WrapperCollection * const stateVariables ) = 0;

};


}
} /* namespace geosx */

#endif /* COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEBASE_HPP_ */
