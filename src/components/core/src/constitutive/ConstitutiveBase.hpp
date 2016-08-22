/*
 * ConstitutiveBase.hpp
 *
 *  Created on: Jul 28, 2016
 *      Author: rrsettgast
 */

#ifndef COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEBASE_HPP_
#define COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEBASE_HPP_

#include "../dataRepository/SynchronizedGroup.hpp"
#include "ObjectCatalog.hpp"

namespace geosx
{

namespace constitutive
{

class ConstitutiveBase : public dataRepository::SynchronizedGroup
{
public:
  ConstitutiveBase( std::string const & name,
                    SynchronizedGroup * const parent );
  virtual ~ConstitutiveBase();

  virtual void Registration( dataRepository::SynchronizedGroup * const );

  virtual void Update( dataRepository::SynchronizedGroup * const parameters,
                       dataRepository::SynchronizedGroup * const stateVariables ) = 0;

  using CatalogInterface = cxx_utilities::CatalogInterface< ConstitutiveBase, std::string const &, SynchronizedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();

};


}
} /* namespace geosx */

#endif /* COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEBASE_HPP_ */
