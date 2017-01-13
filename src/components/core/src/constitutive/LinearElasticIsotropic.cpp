/*
 * HypoElasticLinear.cpp
 *
 *  Created on: Aug 4, 2016
 *      Author: rrsettgast
 */

#include "LinearElasticIsotropic.hpp"

namespace geosx
{
using namespace dataRepository;
namespace constitutive
{


LinearElasticIsotropic::LinearElasticIsotropic( std::string const & name ):
    ConstitutiveBase(name)
{
  // TODO Auto-generated constructor stub

}

LinearElasticIsotropic::~LinearElasticIsotropic()
{
  // TODO Auto-generated destructor stub
}


void LinearElasticIsotropic::StateUpdate( dataRepository::ManagedGroup const * const input,
                                                  dataRepository::ManagedGroup const * const parameters,
                                                  dataRepository::ManagedGroup * const stateVariables,
                                                  integer const systemAssembleFlag ) const
{

  index_t numberOfMaterialPoints = stateVariables->size();
  ViewWrapper<real64_array>::rtype_const K = parameters->getData<real64_array>(std::string("BulkModulus"));
  ViewWrapper<real64_array>::rtype_const G = parameters->getData<real64_array>(std::string("ShearModulus"));

  ViewWrapper<real64_array>::rtype meanStress = stateVariables->getData<real64_array>(std::string("MeanStress"));
  ViewWrapper<real64_array>::rtype S11 = stateVariables->getData<real64_array>(std::string("S11"));
  ViewWrapper<real64_array>::rtype S22 = stateVariables->getData<real64_array>(std::string("S22"));
  ViewWrapper<real64_array>::rtype S33 = stateVariables->getData<real64_array>(std::string("S33"));
  ViewWrapper<real64_array>::rtype S23 = stateVariables->getData<real64_array>(std::string("S23"));
  ViewWrapper<real64_array>::rtype S13 = stateVariables->getData<real64_array>(std::string("S13"));
  ViewWrapper<real64_array>::rtype S12 = stateVariables->getData<real64_array>(std::string("S12"));

  ViewWrapper<real64_array>::rtype_const D11 = input->getData<real64_array>(std::string("D11"));
  ViewWrapper<real64_array>::rtype_const D22 = input->getData<real64_array>(std::string("D22"));
  ViewWrapper<real64_array>::rtype_const D33 = input->getData<real64_array>(std::string("D33"));
  ViewWrapper<real64_array>::rtype_const D23 = input->getData<real64_array>(std::string("D23"));
  ViewWrapper<real64_array>::rtype_const D13 = input->getData<real64_array>(std::string("D13"));
  ViewWrapper<real64_array>::rtype_const D12 = input->getData<real64_array>(std::string("D12"));

  for( index_t i=0 ; i<numberOfMaterialPoints ; ++i )
  {
    real volumeStrain = ( D11[i] + D22[i] + D33[i] );
    meanStress[i] += volumeStrain * K[i];

    S11[i] += ( D11[i] - volumeStrain/3.0 ) * 2.0 * G[i];
    S22[i] += ( D22[i] - volumeStrain/3.0 ) * 2.0 * G[i];
    S33[i] += ( D33[i] - volumeStrain/3.0 ) * 2.0 * G[i];
    S23[i] += ( D23[i] ) * 2.0 * G[i];
    S13[i] += ( D13[i] ) * 2.0 * G[i];
    S12[i] += ( D12[i] ) * 2.0 * G[i];

  }

  if ( systemAssembleFlag == 1 )
  {
    ViewWrapper<real64_array>::rtype K11 = stateVariables->getData<real64_array>(std::string("K11"));
    ViewWrapper<real64_array>::rtype K22 = stateVariables->getData<real64_array>(std::string("K22"));
    ViewWrapper<real64_array>::rtype K33 = stateVariables->getData<real64_array>(std::string("K33"));
    ViewWrapper<real64_array>::rtype K23 = stateVariables->getData<real64_array>(std::string("K23"));
    ViewWrapper<real64_array>::rtype K13 = stateVariables->getData<real64_array>(std::string("K13"));
    ViewWrapper<real64_array>::rtype K12 = stateVariables->getData<real64_array>(std::string("K12"));
    ViewWrapper<real64_array>::rtype K44 = stateVariables->getData<real64_array>(std::string("K44"));
    ViewWrapper<real64_array>::rtype K55 = stateVariables->getData<real64_array>(std::string("K55"));
    ViewWrapper<real64_array>::rtype K66 = stateVariables->getData<real64_array>(std::string("K66"));

    for( index_t i=0 ; i<numberOfMaterialPoints ; ++i )
    {
      real Stiffness[6][6] = {{ K[i]+4.0/3.0*G[i], K[i]-2.0/3.0*G[i], K[i]-2.0/3.0*G[i], 0,               0,        0 },
                              { K[i]-2.0/3.0*G[i], K[i]+4.0/3.0*G[i], K[i]-2.0/3.0*G[i], 0,               0,        0 },
                              { K[i]-2.0/3.0*G[i], K[i]-2.0/3.0*G[i], K[i]+4.0/3.0*G[i], 0,               0,        0 },
                              {                 0,                 0,                 0, 2.0*G[i],        0,        0 },
                              {                 0,                 0,                 0,        0, 2.0*G[i],        0 },
                              {                 0,                 0,                 0,        0 ,       0, 2.0*G[i] } };

      K11[i] = Stiffness[0][0];
      K22[i] = Stiffness[1][1];
      K33[i] = Stiffness[2][2];
      K44[i] = Stiffness[3][3];
      K55[i] = Stiffness[4][4];
      K66[i] = Stiffness[5][5];
      K23[i] = Stiffness[1][2];
      K13[i] = Stiffness[0][2];
      K12[i] = Stiffness[0][1];
    }
  }


}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, LinearElasticIsotropic, std::string const & )
}
} /* namespace geosx */
