/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  m_parameterData.RegisterViewWrapper( viewKeys.bulkModulus.Key(), &m_bulkModulus, 0 );
  m_parameterData.RegisterViewWrapper( viewKeys.referenceDensity.Key(), &m_referenceDensity, 0 );
  m_parameterData.RegisterViewWrapper( viewKeys.referencePressure.Key(), &m_referencePressure, 0 );
  m_parameterData.RegisterViewWrapper( viewKeys.fluidViscosity.Key(), &m_fluidViscosity, 0 );
  m_stateData.RegisterViewWrapper( viewKeys.fluidDensity.Key(), &m_fluidDensity, 0 );
  m_stateData.RegisterViewWrapper( viewKeys.fluidPressure.Key(), &m_fluidPressure, 0 );
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

  ManagedGroup * parameterData = this->GetGroup( groupKeys().ParameterData.Key() );
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

  parameterDocNode->AllocateChildNode( viewKeys.referencePressure.Key(),
                                       viewKeys.referencePressure.Key(),
                                       -1,
                                       "real64",
                                       "real64",
                                       "referencePressure",
                                       "referencePressure",
                                       "-1",
                                       "",
                                       1,
                                       1,
                                       0 );

  parameterDocNode->AllocateChildNode( viewKeys.fluidViscosity.Key(),
                                       viewKeys.fluidViscosity.Key(),
                                       -1,
                                       "real64",
                                       "real64",
                                       "fluidViscosity",
                                       "fluidViscosity",
                                       "0.001",
                                       "",
                                       1,
                                       1,
                                       0 );


  ManagedGroup * stateData     = this->GetGroup( groupKeys().StateData.Key() );
  DocumentationNode * const stateDocNode = stateData->getDocumentationNode();
  stateDocNode->setSchemaType("Node");
  stateDocNode->setShortDescription("State for Linear Elastic Isotropic Constitutive Relation");

  stateDocNode->AllocateChildNode( viewKeys.fluidPressure.Key(),
                                   viewKeys.fluidPressure.Key(),
                                   -1,
                                   "real64_array",
                                   "real64_array",
                                   "pressure",
                                   "pressure",
                                   "0",
                                   "",
                                   1,
                                   0,
                                   0 );

  stateDocNode->AllocateChildNode( viewKeys.fluidDensity.Key(),
                                   viewKeys.fluidDensity.Key(),
                                   -1,
                                   "real64_array",
                                   "real64_array",
                                   "density",
                                   "density",
                                   "0",
                                   "",
                                   1,
                                   0,
                                   0 );

}

void LinearEOS::ReadXML_PostProcess()
{
  if( m_bulkModulus <= 0.0 )
  {
    string const message = "An invalid value of bulk modulus ("+std::to_string(m_bulkModulus)+") is specified";
    GEOS_ERROR(message);
  }

  if( m_referenceDensity <= 0.0 )
   {
     string const message = "An invalid value of reference density ("+std::to_string(m_referenceDensity)+") is specified";
     GEOS_ERROR(message);
   }
}


void LinearEOS::SetParamStatePointers( void *& data )
{
//
//  this->m_dataPointers.m_bulkModulus = this->bulkModulus();
//
//  data = reinterpret_cast<void*>(&m_dataPointers);
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
//  real64 volumeStrain = D.Trace();
//  meanStress()[i] += volumeStrain * bulkModulus()[0];
//
//  R2SymTensor temp = D;
//  temp.PlusIdentity(-volumeStrain / 3.0);
//  temp *= 2.0 * shearModulus()[0];
//  deviatorStress()[i] += temp;
//
//
//  temp.QijAjkQlk(deviatorStress()[i],Rot);
//  deviatorStress()[i] = temp;
//
//  temp.PlusIdentity(meanStress()[i]);
//  return temp;
}

void LinearEOS::EquationOfStatePressureUpdate( real64 const & dRho,
                                               localIndex const i,
                                               real64 & dP,
                                               real64 & dPdRho )
{
  dP = m_bulkModulus * log( ( m_fluidDensity[i] + dRho ) / m_fluidDensity[i] );
//  m_fluidDensity[i] += dRho;
//  m_fluidPressure[i] = m_bulkModulus * log( m_fluidDensity[i] / m_referenceDensity );
  dPdRho = m_bulkModulus / m_fluidDensity[i];

}

void LinearEOS::EquationOfStateDensityUpdate( real64 const & dP,
                                              localIndex const i,
                                              real64 & dRho,
                                              real64 & dRho_dP )
{

  //dRho = m_fluidDensity[i] * ( exp(dP / m_bulkModulus) - 1) ;
  // use taylor expansion to avoid loss of precision
  dRho = m_fluidDensity[i] * ( dP / m_bulkModulus) ;//* ( 1 + 0.5 * dP / m_bulkModulus );

//  m_fluidPressure[i] += dP;
//  m_fluidDensity[i] = m_referenceDensity * exp( ( m_fluidPressure[i] - m_referencePressure ) / m_bulkModulus );
  dRho_dP = m_fluidDensity[i] / m_bulkModulus;
}

void LinearEOS::GetStiffness( realT c[6][6]) const
{
  c[0][0] = m_bulkModulus;
  c[0][1] = 0.0;
  c[0][2] = 0.0;
  c[0][3] = 0.0;
  c[0][4] = 0.0;
  c[0][5] = 0.0;

  c[1][0] = 0.0;
  c[1][1] = m_bulkModulus;
  c[1][2] = 0.0;
  c[1][3] = 0.0;
  c[1][4] = 0.0;
  c[1][5] = 0.0;

  c[2][0] = 0.0;
  c[2][1] = 0.0;
  c[2][2] = m_bulkModulus;
  c[2][3] = 0.0;
  c[2][4] = 0.0;
  c[2][5] = 0.0;

  c[3][0] = 0.0;
  c[3][1] = 0.0;
  c[3][2] = 0.0;
  c[3][3] = 0.0;
  c[3][4] = 0.0;
  c[3][5] = 0.0;

  c[4][0] = 0.0;
  c[4][1] = 0.0;
  c[4][2] = 0.0;
  c[4][3] = 0.0;
  c[4][4] = 0.0;
  c[4][5] = 0.0;

  c[5][0] = 0.0;
  c[5][1] = 0.0;
  c[5][2] = 0.0;
  c[5][3] = 0.0;
  c[5][4] = 0.0;
  c[5][5] = 0.0;
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, LinearEOS, std::string const &, ManagedGroup * const )
}
} /* namespace geosx */
