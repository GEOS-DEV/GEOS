/*
 * ConstitutiveBase.hpp
 *
 *  Created on: Jul 28, 2016
 *      Author: rrsettgast
 */

#ifndef COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEBASE_HPP_
#define COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEBASE_HPP_

#include "dataRepository/WrapperCollection.hpp"
#include "ObjectCatalog.hpp"

namespace geosx
{

namespace constitutive
{

class ConstitutiveBase : public dataRepository::WrapperCollection
{
public:
  ConstitutiveBase( std::string const & name,
                    WrapperCollection * const parent );
  virtual ~ConstitutiveBase();

  virtual void Registration( dataRepository::WrapperCollection * const );

  virtual void Update( dataRepository::WrapperCollection * const parameters,
                       dataRepository::WrapperCollection * const stateVariables ) = 0;

  using CatalogInterface = cxx_utilities::CatalogInterface< ConstitutiveBase, std::string const &, WrapperCollection * const >;
  static CatalogInterface::CatalogType& GetCatalog();

};


}
} /* namespace geosx */

#endif /* COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEBASE_HPP_ */
