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


void HypoElasticLinear::Update( dataRepository::SynchronizedGroup * const parameters,
                                dataRepository::SynchronizedGroup * const stateVariables )
{

}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, HypoElasticLinear, std::string const &, SynchronizedGroup * const )
}
} /* namespace geosx */
