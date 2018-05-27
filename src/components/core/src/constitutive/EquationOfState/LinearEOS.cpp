/*
 * HypoElasticLinear.cpp
 *
 *  Created on: Aug 4, 2016
 *      Author: rrsettgast
 */

#include "LinearEOS.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace cxx_utilities;
namespace constitutive
{


static inline void UpdateStatePoint( R2SymTensor const & D,
                                     R2Tensor const & Rot,
                                     localIndex const i,
                                     void * dataPtrs,
                                     integer const systemAssembleFlag )
{


}


LinearEOS::LinearEOS( std::string const & name, ManagedGroup * const parent ):
  ConstitutiveBase(name, parent )
{
  // TODO Auto-generated constructor stub

}

LinearEOS::~LinearEOS()
{
  // TODO Auto-generated destructor stub
}

void LinearEOS::FillDocumentationNode()
{

  DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName(this->CatalogName());
  docNode->setSchemaType("Node");
  docNode->setShortDescription("Linear Elastic Isotropic Constitutive Relation");

  ManagedGroup * parameterData = this->GetGroup( groupKeys.ParameterData.Key() );
  DocumentationNode * const parameterDocNode = parameterData->getDocumentationNode();
  parameterDocNode->setSchemaType("Node");
  parameterDocNode->setShortDescription("Parameters for Linear Elastic Isotropic Constitutive Relation");



  parameterDocNode->AllocateChildNode( viewKeys.bulkModulus.Key(),
                                       viewKeys.bulkModulus.Key(),
                                       -1,
                                       "real64",
                                       "real64",
                                       "Elastic Bulk Modulus",
                                       "Elastic Bulk Modulus",
                                       "-1",
                                       "",
                                       1,
                                       1,
                                       0 );


  parameterDocNode->AllocateChildNode( viewKeys.referenceDensity.Key(),
                                       viewKeys.referenceDensity.Key(),
                                       -1,
                                       "real64",
                                       "real64",
                                       "referenceDensity",
                                       "referenceDensity",
                                       "-1",
                                       "",
                                       1,
                                       1,
                                       0 );


  ManagedGroup * stateData     = this->GetGroup( groupKeys.StateData.Key() );
  DocumentationNode * const stateDocNode = stateData->getDocumentationNode();
  stateDocNode->setSchemaType("Node");
  stateDocNode->setShortDescription("State for Linear Elastic Isotropic Constitutive Relation");

  stateDocNode->AllocateChildNode( viewKeys.meanStress.Key(),
                                   viewKeys.meanStress.Key(),
                                   -1,
                                   "real64_array",
                                   "real64_array",
                                   "Stress",
                                   "Stress",
                                   "0",
                                   "",
                                   1,
                                   0,
                                   0 );
}

void LinearEOS::ReadXML_PostProcess()
{
  ManagedGroup * parameterData = this->GetParameterData();
  real64 & K  = *( parameterData->getData<real64>(viewKeys.bulkModulus) );

  if( K >= 0.0 )
  {
    ++numConstantsSpecified;
  }
  if( K <= 0.0 )
  {
    string const message = "An invalid value of K ("+std::to_string(K)+") is specified";
    GEOS_ERROR(message);
  }
}

void LinearEOS::SetParamStatePointers( void *& data )
{

  this->m_dataPointers.m_bulkModulus = this->bulkModulus();

  data = reinterpret_cast<void*>(&m_dataPointers);
}

ConstitutiveBase::UpdateFunctionPointer
LinearEOS::GetStateUpdateFunctionPointer()
{
  return UpdateStatePoint;
}


void LinearEOS::StateUpdate( dataRepository::ManagedGroup const * const input,
                                          dataRepository::ManagedGroup const * const parameters,
                                          dataRepository::ManagedGroup * const stateVariables,
                                          integer const systemAssembleFlag ) const
{

//  localIndex numberOfMaterialPoints = stateVariables->size();
//  ViewWrapper<real64_array>::rtype_const K = parameters->getData<real64_array>(std::string("BulkModulus"));
//
//  ViewWrapper<real64_array>::rtype mean_stress = stateVariables->getData<real64_array>(std::string("MeanStress"));
//
//  for( localIndex i=0 ; i<numberOfMaterialPoints ; ++i )
//  {
//    real volumeStrain = ( D11[i] + D22[i] + D33[i] );
//    mean_stress[i] += volumeStrain * K[i];
//
//
//  }
//
//  if ( systemAssembleFlag == 1 )
//  {
//    ViewWrapper<real64_array>::rtype K11 = stateVariables->getData<real64_array>(std::string("K11"));
//    ViewWrapper<real64_array>::rtype K22 = stateVariables->getData<real64_array>(std::string("K22"));
//    ViewWrapper<real64_array>::rtype K33 = stateVariables->getData<real64_array>(std::string("K33"));
//    ViewWrapper<real64_array>::rtype K23 = stateVariables->getData<real64_array>(std::string("K23"));
//    ViewWrapper<real64_array>::rtype K13 = stateVariables->getData<real64_array>(std::string("K13"));
//    ViewWrapper<real64_array>::rtype K12 = stateVariables->getData<real64_array>(std::string("K12"));
//    ViewWrapper<real64_array>::rtype K44 = stateVariables->getData<real64_array>(std::string("K44"));
//    ViewWrapper<real64_array>::rtype K55 = stateVariables->getData<real64_array>(std::string("K55"));
//    ViewWrapper<real64_array>::rtype K66 = stateVariables->getData<real64_array>(std::string("K66"));
//
//    for( localIndex i=0 ; i<numberOfMaterialPoints ; ++i )
//    {
//      real Stiffness[6][6] = {
//        { K[i]+4.0/3.0*G[i], K[i]-2.0/3.0*G[i], K[i]-2.0/3.0*G[i], 0,               0,        0 },
//        { K[i]-2.0/3.0*G[i], K[i]+4.0/3.0*G[i], K[i]-2.0/3.0*G[i], 0,               0,        0 },
//        { K[i]-2.0/3.0*G[i], K[i]-2.0/3.0*G[i], K[i]+4.0/3.0*G[i], 0,               0,        0 },
//        {                 0,                 0,                 0, 2.0*G[i],        0,        0 },
//        {                 0,                 0,                 0,        0, 2.0*G[i],        0 },
//        {                 0,                 0,                 0,        0,       0, 2.0*G[i] }
//      };
//
//      K11[i] = Stiffness[0][0];
//      K22[i] = Stiffness[1][1];
//      K33[i] = Stiffness[2][2];
//      K44[i] = Stiffness[3][3];
//      K55[i] = Stiffness[4][4];
//      K66[i] = Stiffness[5][5];
//      K23[i] = Stiffness[1][2];
//      K13[i] = Stiffness[0][2];
//      K12[i] = Stiffness[0][1];
//    }
//  }


}

R2SymTensor LinearEOS::StateUpdatePoint( R2SymTensor const & D,
                                         R2Tensor const & Rot,
                                         localIndex const i,
                                         integer const systemAssembleFlag )
{
  real64 volumeStrain = D.Trace();
  meanStress()[i] += volumeStrain * bulkModulus()[0];

  R2SymTensor temp = D;
  temp.PlusIdentity(-volumeStrain / 3.0);
  temp *= 2.0 * shearModulus()[0];
  deviatorStress()[i] += temp;


  temp.QijAjkQlk(deviatorStress()[i],Rot);
  deviatorStress()[i] = temp;

  temp.PlusIdentity(meanStress()[i]);
  return temp;
}

void LinearEOS::GetStiffness( realT c[6][6]) const
{
  real64 G = *shearModulus();
  real64 Lame = *bulkModulus() - 2.0/3.0 * G;
  c[0][0] = Lame + 2 * G;
  c[0][1] = Lame;
  c[0][2] = Lame;
  c[0][3] = 0.0;
  c[0][4] = 0.0;
  c[0][5] = 0.0;

  c[1][0] = Lame;
  c[1][1] = Lame + 2 * G;
  c[1][2] = Lame;
  c[1][3] = 0.0;
  c[1][4] = 0.0;
  c[1][5] = 0.0;

  c[2][0] = Lame;
  c[2][1] = Lame;
  c[2][2] = Lame + 2 * G;
  c[2][3] = 0.0;
  c[2][4] = 0.0;
  c[2][5] = 0.0;

  c[3][0] = 0.0;
  c[3][1] = 0.0;
  c[3][2] = 0.0;
  c[3][3] = G;
  c[3][4] = 0.0;
  c[3][5] = 0.0;

  c[4][0] = 0.0;
  c[4][1] = 0.0;
  c[4][2] = 0.0;
  c[4][3] = 0.0;
  c[4][4] = G;
  c[4][5] = 0.0;

  c[5][0] = 0.0;
  c[5][1] = 0.0;
  c[5][2] = 0.0;
  c[5][3] = 0.0;
  c[5][4] = 0.0;
  c[5][5] = G;
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, LinearEOS, std::string const &, ManagedGroup * const )
}
} /* namespace geosx */
