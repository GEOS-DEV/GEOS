/*
 * ConstitutiveManager.hpp
 *
 *  Created on: Aug 4, 2016
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEMANAGER_HPP_
#include "../dataRepository/ManagedGroup.hpp"
#include "ConstitutiveBase.hpp"

namespace geosx
{
namespace constitutive
{


class ConstitutiveManager : public dataRepository::ManagedGroup
{
public:
  ConstitutiveManager() = delete;

  ConstitutiveManager( std::string const & name,
                       ManagedGroup * const parent );

  void ReadXMLInput();

  ~ConstitutiveManager();
};

}
} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEMANAGER_HPP_ */
