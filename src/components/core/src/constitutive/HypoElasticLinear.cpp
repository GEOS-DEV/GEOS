/*
 * HypoElasticLinear.cpp
 *
 *  Created on: Aug 4, 2016
 *      Author: rrsettgast
 */

#include "HypoElasticLinear.hpp"

namespace geosx
{
namespace constitutive
{


HypoElasticLinear::HypoElasticLinear()
{
  // TODO Auto-generated constructor stub

}

HypoElasticLinear::~HypoElasticLinear()
{
  // TODO Auto-generated destructor stub
}


void HypoElasticLinear::Update( dataRepository::WrapperCollection * const parameters,
                                dataRepository::WrapperCollection * const stateVariables )
{

}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, HypoElasticLinear, std::string const &, WrapperCollection * const )
}
} /* namespace geosx */
