/*
 * ConstitutiveBase.hpp
 *
 *  Created on: Jul 28, 2016
 *      Author: rrsettgast
 */

#ifndef COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEBASE_HPP_
#define COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEBASE_HPP_

#include "../dataRepository/ManagedGroup.hpp"
#include "ObjectCatalog.hpp"

namespace geosx
{

namespace constitutive
{

class ConstitutiveBase : public dataRepository::ManagedGroup
{
public:
  ConstitutiveBase( std::string const & name,
                    ManagedGroup * const parent );
  virtual ~ConstitutiveBase();

  virtual void Registration( dataRepository::ManagedGroup * const );

  virtual void Update( dataRepository::ManagedGroup * const parameters,
                       dataRepository::ManagedGroup * const stateVariables ) = 0;

  using CatalogInterface = cxx_utilities::CatalogInterface< ConstitutiveBase, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();

};


}
} /* namespace geosx */

#endif /* COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEBASE_HPP_ */
