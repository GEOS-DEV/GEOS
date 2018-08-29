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

/*
 * HypoElasticLinear.cpp
 *
 *  Created on: Aug 4, 2016
 *      Author: rrsettgast
 */

#include "LinearElasticIsotropic.hpp"

#ifdef USE_ATK
#include "slic/slic.hpp"
#endif

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

  LinearElasticIsotropic::dataPointers * castedDataPtrs = reinterpret_cast<LinearElasticIsotropic::dataPointers *>(dataPtrs);
  real64 volumeStrain = D.Trace();
  castedDataPtrs->m_meanStress[i] += volumeStrain * castedDataPtrs->m_bulkModulus[0];

  R2SymTensor temp = D;
  temp.PlusIdentity( -volumeStrain / 3.0 );
  temp *= 2.0 * castedDataPtrs->m_shearModulus[0];
  castedDataPtrs->m_deviatorStress[i] += temp;


  temp.QijAjkQlk( castedDataPtrs->m_deviatorStress[i], Rot );
  castedDataPtrs->m_deviatorStress[i] = temp;

  temp.PlusIdentity( castedDataPtrs->m_meanStress[i] );
}


LinearElasticIsotropic::LinearElasticIsotropic( std::string const & name, ManagedGroup * const parent ):
  ConstitutiveBase( name, parent )
{
  this->RegisterViewWrapper( viewKeyStruct::bulkModulus0String, &m_bulkModulus0, 0 );
  this->RegisterViewWrapper( viewKeyStruct::bulkModulusString, &m_bulkModulus, 0 );
  this->RegisterViewWrapper( viewKeyStruct::density0String, &m_density0, 0 );
  this->RegisterViewWrapper( viewKeyStruct::densityString, &m_density, 0 );
  this->RegisterViewWrapper( viewKeyStruct::shearModulus0String, &m_shearModulus0, 0 );
  this->RegisterViewWrapper( viewKeyStruct::shearModulusString, &m_shearModulus, 0 );


  this->RegisterViewWrapper( viewKeyStruct::deviatorStressString, &m_deviatorStress, 0 );
  this->RegisterViewWrapper( viewKeyStruct::meanStressString, &m_meanStress, 0 );

}

LinearElasticIsotropic::~LinearElasticIsotropic()
{}


std::unique_ptr<ConstitutiveBase>
LinearElasticIsotropic::DeliverClone( string const & name,
                                      ManagedGroup * const parent ) const
{
  std::unique_ptr<LinearElasticIsotropic>
  newConstitutiveRelation = std::make_unique<LinearElasticIsotropic>( name, parent );

  newConstitutiveRelation->m_bulkModulus0 = m_bulkModulus0;
  newConstitutiveRelation->m_bulkModulus = m_bulkModulus;
  newConstitutiveRelation->m_density0 = m_density0;
  newConstitutiveRelation->m_density = m_density;
  newConstitutiveRelation->m_shearModulus0 = m_shearModulus0;
  newConstitutiveRelation->m_shearModulus = m_shearModulus;

  newConstitutiveRelation->m_meanStress = m_meanStress;
  newConstitutiveRelation->m_deviatorStress = m_deviatorStress;

  std::unique_ptr<ConstitutiveBase> rval = std::move( newConstitutiveRelation );

  return rval;
}

void LinearElasticIsotropic::AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                                       localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );
  m_bulkModulus.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_deviatorStress.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_density.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_meanStress.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_shearModulus.resize( parent->size(), numConstitutivePointsPerParentIndex );

  m_bulkModulus = m_bulkModulus0;
  m_density = m_density0;
  m_shearModulus = m_shearModulus0;

}

void LinearElasticIsotropic::FillDocumentationNode()
{

  DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName( this->CatalogName());
  docNode->setSchemaType( "Node" );
  docNode->setShortDescription( "Linear Elastic Isotropic Constitutive Relation" );

  docNode->AllocateChildNode( viewKeys().youngsModulus.Key(),
                              viewKeys().youngsModulus.Key(),
                              -1,
                              "real64",
                              "real64",
                              "Young's Elastic Modulus",
                              "Young's Elastic Modulus",
                              "-1",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeyStruct::bulkModulus0String,
                              viewKeyStruct::bulkModulus0String,
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

  docNode->AllocateChildNode( viewKeyStruct::shearModulus0String,
                              viewKeyStruct::shearModulus0String,
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

  docNode->AllocateChildNode( viewKeys().poissonRatio.Key(),
                              viewKeys().poissonRatio.Key(),
                              -1,
                              "real64",
                              "real64",
                              "Elastic Poisson's Ratio",
                              "Elastic Poisson's Ratio",
                              "-1",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeyStruct::density0String,
                              viewKeyStruct::density0String,
                              -1,
                              "real64",
                              "real64",
                              "density",
                              "density",
                              "Required",
                              "",
                              1,
                              1,
                              0 );



}

void LinearElasticIsotropic::ReadXML_PostProcess()
{
  real64 & nu = *( getData<real64>( viewKeys().poissonRatio ) );
  real64 & E  = *( getData<real64>( viewKeys().youngsModulus ) );
  real64 & K  = m_bulkModulus0;
  real64 & G  = m_shearModulus0;

  int numConstantsSpecified = 0;
  if( nu >= 0.0 )
  {
    ++numConstantsSpecified;
  }
  if( E >= 0.0 )
  {
    ++numConstantsSpecified;
  }
  if( K >= 0.0 )
  {
    ++numConstantsSpecified;
  }
  if( G >= 0.0 )
  {
    ++numConstantsSpecified;
  }

  if( numConstantsSpecified == 2 )
  {
    if( nu >= 0.0 && E >= 0.0 )
    {
      K = E / (3 * ( 1 - 2*nu ) );
      G = E / (2 * ( 1 + nu ) );
    }
    else if( !( K >= 0.0 && G >= 0.0 ) )
    {
      string const message = "A specific pair of elastic constants is required. Either (K,G) or (E,nu)";
      GEOS_ERROR( message );
    }
    else
    {
      E = 9 * K * G / ( 3 * K + G );
      nu = ( 3 * K - 2 * G ) / ( 2 * ( 3 * K + G ) );
    }
  }
  else
  {
    string const message = std::to_string( numConstantsSpecified ) + " Elastic Constants Specified. Must specify 2 constants!";
    GEOS_ERROR( message );
  }
}

void LinearElasticIsotropic::SetParamStatePointers( void *& data )
{

  this->m_dataPointers.m_bulkModulus = &m_bulkModulus0;
  this->m_dataPointers.m_shearModulus = &m_shearModulus0;
  this->m_dataPointers.m_meanStress = m_meanStress.data();
  this->m_dataPointers.m_deviatorStress = m_deviatorStress.data();

  data = reinterpret_cast<void*>(&m_dataPointers);
}

ConstitutiveBase::UpdateFunctionPointer
LinearElasticIsotropic::GetStateUpdateFunctionPointer()
{
  return UpdateStatePoint;
}


void LinearElasticIsotropic::StateUpdate( dataRepository::ManagedGroup const * const input,
                                          dataRepository::ManagedGroup const * const parameters,
                                          dataRepository::ManagedGroup * const stateVariables,
                                          integer const systemAssembleFlag ) const
{

  localIndex numberOfMaterialPoints = stateVariables->size();
  ViewWrapper<real64_array>::rtype_const K = parameters->getData<real64_array>( std::string( "BulkModulus" ));
  ViewWrapper<real64_array>::rtype_const G = parameters->getData<real64_array>( std::string( "ShearModulus" ));

  ViewWrapper<real64_array>::rtype mean_stress = stateVariables->getData<real64_array>( std::string( "MeanStress" ));
  ViewWrapper<real64_array>::rtype S11 = stateVariables->getData<real64_array>( std::string( "S11" ));
  ViewWrapper<real64_array>::rtype S22 = stateVariables->getData<real64_array>( std::string( "S22" ));
  ViewWrapper<real64_array>::rtype S33 = stateVariables->getData<real64_array>( std::string( "S33" ));
  ViewWrapper<real64_array>::rtype S23 = stateVariables->getData<real64_array>( std::string( "S23" ));
  ViewWrapper<real64_array>::rtype S13 = stateVariables->getData<real64_array>( std::string( "S13" ));
  ViewWrapper<real64_array>::rtype S12 = stateVariables->getData<real64_array>( std::string( "S12" ));

  ViewWrapper<real64_array>::rtype_const D11 = input->getData<real64_array>( std::string( "D11" ));
  ViewWrapper<real64_array>::rtype_const D22 = input->getData<real64_array>( std::string( "D22" ));
  ViewWrapper<real64_array>::rtype_const D33 = input->getData<real64_array>( std::string( "D33" ));
  ViewWrapper<real64_array>::rtype_const D23 = input->getData<real64_array>( std::string( "D23" ));
  ViewWrapper<real64_array>::rtype_const D13 = input->getData<real64_array>( std::string( "D13" ));
  ViewWrapper<real64_array>::rtype_const D12 = input->getData<real64_array>( std::string( "D12" ));

  for( localIndex i=0 ; i<numberOfMaterialPoints ; ++i )
  {
    real64 volumeStrain = ( D11[i] + D22[i] + D33[i] );
    mean_stress[i] += volumeStrain * K[i];

    S11[i] += ( D11[i] - volumeStrain/3.0 ) * 2.0 * G[i];
    S22[i] += ( D22[i] - volumeStrain/3.0 ) * 2.0 * G[i];
    S33[i] += ( D33[i] - volumeStrain/3.0 ) * 2.0 * G[i];
    S23[i] += ( D23[i] ) * 2.0 * G[i];
    S13[i] += ( D13[i] ) * 2.0 * G[i];
    S12[i] += ( D12[i] ) * 2.0 * G[i];

  }

  if( systemAssembleFlag == 1 )
  {
    ViewWrapper<real64_array>::rtype K11 = stateVariables->getData<real64_array>( std::string( "K11" ));
    ViewWrapper<real64_array>::rtype K22 = stateVariables->getData<real64_array>( std::string( "K22" ));
    ViewWrapper<real64_array>::rtype K33 = stateVariables->getData<real64_array>( std::string( "K33" ));
    ViewWrapper<real64_array>::rtype K23 = stateVariables->getData<real64_array>( std::string( "K23" ));
    ViewWrapper<real64_array>::rtype K13 = stateVariables->getData<real64_array>( std::string( "K13" ));
    ViewWrapper<real64_array>::rtype K12 = stateVariables->getData<real64_array>( std::string( "K12" ));
    ViewWrapper<real64_array>::rtype K44 = stateVariables->getData<real64_array>( std::string( "K44" ));
    ViewWrapper<real64_array>::rtype K55 = stateVariables->getData<real64_array>( std::string( "K55" ));
    ViewWrapper<real64_array>::rtype K66 = stateVariables->getData<real64_array>( std::string( "K66" ));

    for( localIndex i=0 ; i<numberOfMaterialPoints ; ++i )
    {
      real64 Stiffness[6][6] = {
        { K[i]+4.0/3.0*G[i], K[i]-2.0/3.0*G[i], K[i]-2.0/3.0*G[i], 0, 0, 0 },
        { K[i]-2.0/3.0*G[i], K[i]+4.0/3.0*G[i], K[i]-2.0/3.0*G[i], 0, 0, 0 },
        { K[i]-2.0/3.0*G[i], K[i]-2.0/3.0*G[i], K[i]+4.0/3.0*G[i], 0, 0, 0 },
        {                 0, 0, 0, 2.0*G[i], 0, 0 },
        {                 0, 0, 0, 0, 2.0*G[i], 0 },
        {                 0, 0, 0, 0, 0, 2.0*G[i] }
      };

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

R2SymTensor LinearElasticIsotropic::StateUpdatePoint( R2SymTensor const & D,
                                                      R2Tensor const & Rot,
                                                      localIndex const i,
                                                      localIndex const q,
                                                      integer const systemAssembleFlag )
{
  real64 volumeStrain = D.Trace();
  m_meanStress[i][q] += volumeStrain * m_bulkModulus0;

  R2SymTensor temp = D;
  temp.PlusIdentity( -volumeStrain / 3.0 );
  temp *= 2.0 * m_shearModulus0;
  m_deviatorStress[i][q] += temp;


  temp.QijAjkQlk( m_deviatorStress[i][q], Rot );
  m_deviatorStress[i][q] = temp;

  temp.PlusIdentity( m_meanStress[i][q] );
  return temp;
}

void LinearElasticIsotropic::GetStiffness( realT c[6][6] ) const
{
  real64 G = m_shearModulus0;
  real64 Lame = m_bulkModulus0 - 2.0/3.0 * G;
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


REGISTER_CATALOG_ENTRY( ConstitutiveBase, LinearElasticIsotropic, std::string const &, ManagedGroup * const )
}
} /* namespace geosx */
