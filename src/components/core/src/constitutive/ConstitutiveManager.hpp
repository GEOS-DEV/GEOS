/*
 * ConstitutiveManager.hpp
 *
 *  Created on: Aug 4, 2016
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEMANAGER_HPP_
#include "dataRepository/WrapperCollection.hpp"
#include "ConstitutiveBase.hpp"

namespace geosx
{
namespace constitutive
{


class ConstitutiveManager : public dataRepository::WrapperCollection
{
public:
  ConstitutiveManager() = delete;

  ConstitutiveManager( std::string const & name,
                       WrapperCollection * const parent );

  void ReadXMLInput();

  ~ConstitutiveManager();
};

}
} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEMANAGER_HPP_ */
