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


void HypoElasticLinear::Update( dataRepository::ManagedGroup * const parameters,
                                dataRepository::ManagedGroup * const stateVariables )
{

}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, HypoElasticLinear, std::string const &, ManagedGroup * const )
}
} /* namespace geosx */
