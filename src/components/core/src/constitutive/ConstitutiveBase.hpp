/*
 * ConstitutiveBase.hpp
 *
 *  Created on: Jul 28, 2016
 *      Author: rrsettgast
 */

#ifndef COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEBASE_HPP_
#define COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEBASE_HPP_

namespace geosx
{

namespace dataRepository
{
class WrapperCollection;
}

namespace constitutive
{

class ConstitutiveBase
{
public:
  ConstitutiveBase();
  virtual ~ConstitutiveBase();

  virtual void Update( WrapperCollection * const parameters,
                       WrapperCollection * const stateVariables ) = 0;

private:
  dataRepository::WrapperCollection * m_wrapperCollection = nullptr;

};


}
} /* namespace geosx */

#endif /* COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEBASE_HPP_ */
