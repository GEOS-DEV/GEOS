/*
 * solidSelector.hpp
 *
 *  Created on: Mar 22, 2019
 *      Author: settgast1
 */

#ifndef SRC_CORECOMPONENTS_CONSTITUTIVE_SOLID_SOLIDSELECTOR_HPP_
#define SRC_CORECOMPONENTS_CONSTITUTIVE_SOLID_SOLIDSELECTOR_HPP_

#include "LinearElasticIsotropic.hpp"

namespace geosx
{
namespace constitutive
{


template< typename LAMBDA >
bool constitutiveUpdatePassThru( constitutive::ConstitutiveBase * const constitutiveRelation,
                                 LAMBDA && lambda )
{
  bool rval = true;
  if( constitutiveRelation->GetCatalogName()==LinearElasticIsotropic::CatalogName() )
  {
    lambda( static_cast<LinearElasticIsotropic & >( *constitutiveRelation) );
  }
  else
  {
    rval = false;
  }

  return rval;
}

}
}

#endif /* SRC_CORECOMPONENTS_CONSTITUTIVE_SOLID_SOLIDSELECTOR_HPP_ */
