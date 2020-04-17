/*
 * Dummy.cpp
 *
 *  Created on: Apr 14, 2020
 *      Author: settgast
 */

#include "Dummy.hpp"

namespace geosx
{
namespace constitutive
{

Dummy::Dummy( string const & name,
              Group * const parent ):
  ConstitutiveBase( name, parent )
{}

Dummy::~Dummy()
{}

void Dummy::DeliverClone( string const & name,
                          Group * const parent,
                          std::unique_ptr< ConstitutiveBase > & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique< Dummy >( name, parent );
  }
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, Dummy, std::string const &, dataRepository::Group * const )

} // constitutive
} /* namespace geosx */
