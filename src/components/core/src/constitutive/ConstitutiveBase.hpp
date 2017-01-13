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
#include "../dataRepository/ManagedGroup.hpp"

namespace geosx
{
namespace constitutive
{

class ConstitutiveBase
{
public:
  ConstitutiveBase( std::string const & name );
  virtual ~ConstitutiveBase();

//  template< typename LEAFCLASS >
  virtual void StateUpdate( dataRepository::ManagedGroup const * const input,
                            dataRepository::ManagedGroup const * const parameters,
                            dataRepository::ManagedGroup * const stateVariables,
                            integer const systemAssembleFlag ) const = 0;

//  inline void AddToGlobalSystem(  );


//  template< typename LEAFCLASS >

//  virtual void FillDocumentationNode( dataRepository::ManagedGroup * const group );


  using CatalogInterface = cxx_utilities::CatalogInterface< ConstitutiveBase, std::string const & >;
  static typename CatalogInterface::CatalogType& GetCatalog();

};




}
}
#endif
