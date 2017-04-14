/*
 * ConstitutiveBase.hpp
 *
 *  Created on: Jul 28, 2016
 *      Author: rrsettgast
 */

#ifndef CONSTITUTIVEBASE_HPP_
#define CONSTITUTIVEBASE_HPP_

#include "common/DataTypes.hpp"
#include "ObjectCatalog.hpp"
#include "../../../cxx-utilities/src/src/DocumentationNode.hpp"
#include "dataRepository/ManagedGroup.hpp"

namespace geosx
{
namespace constitutive
{


class ConstitutiveBase : public dataRepository::ManagedGroup
{
public:


  ConstitutiveBase( std::string const & name,
                    ManagedGroup * const parent  );

  virtual ~ConstitutiveBase();

  virtual void StateUpdate( dataRepository::ManagedGroup const * const input,
                            dataRepository::ManagedGroup const * const parameters,
                            dataRepository::ManagedGroup * const stateVariables,
                            integer const systemAssembleFlag ) const = 0;


  virtual void FillDocumentationNode( dataRepository::ManagedGroup * const group ) = 0;


  using CatalogInterface = cxx_utilities::CatalogInterface< ConstitutiveBase, std::string const &, ManagedGroup * const >;
  static typename CatalogInterface::CatalogType& GetCatalog();

};




}
}
#endif
