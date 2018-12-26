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

namespace geosx
{
using namespace dataRepository;
using namespace cxx_utilities;
namespace constitutive
{


static inline void UpdateStatePoint( R2SymTensor const & D,
                                     R2Tensor const & Rot,
                                     localIndex const i,
                                     localIndex const q,
                                     void * dataPtrs,
                                     integer const systemAssembleFlag )
{

  LinearElasticIsotropic::dataPointers * restrict const
  castedDataPtrs = reinterpret_cast<LinearElasticIsotropic::dataPointers *>(dataPtrs);
  real64 volumeStrain = D.Trace();
  (*castedDataPtrs->m_meanStress)[i][q] += volumeStrain * castedDataPtrs->m_bulkModulus0[0];
//
//  R2SymTensor temp = D;
//  temp.PlusIdentity( -volumeStrain / 3.0 );
//  temp *= 2.0 * castedDataPtrs->m_shearModulus0[0];
//  (*castedDataPtrs->m_deviatorStress)[i][q] += temp;
//
//
//  temp.QijAjkQlk( (*castedDataPtrs->m_deviatorStress)[i][q], Rot );
//  (*castedDataPtrs->m_deviatorStress)[i][q] = temp;
//
//  temp.PlusIdentity( (*castedDataPtrs->m_meanStress)[i][q] );
}


LinearElasticIsotropic::LinearElasticIsotropic( std::string const & name, ManagedGroup * const parent ):
  ConstitutiveBase( name, parent ),
  m_bulkModulus0(),
  m_shearModulus0(),
  m_density0(),
  m_density(),
  m_bulkModulus(),
  m_shearModulus(),
  m_meanStress(),
  m_deviatorStress(),
  m_compressibility(),
  m_referencePressure(),
  m_biotCoefficient(),
  m_poreVolumeMultiplier(),
  m_dPVMult_dPressure(),
  m_poreVolumeRelation( ExponentApproximationType::Linear )
{
  RegisterViewWrapper( viewKeyStruct::density0String, &m_density0, 0 )->
      setInputFlag(InputFlags::REQUIRED)->
      setDescription("Reference Material Density");

  RegisterViewWrapper( viewKeyStruct::bulkModulus0String, &m_bulkModulus0, 0 )->
      setDefaultValue(-1)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("Elastic Bulk Modulus Parameter");

  RegisterViewWrapper( viewKeyStruct::shearModulus0String, &m_shearModulus0, 0 )->
      setDefaultValue(-1)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("Elastic Shear Modulus Parameter");


  RegisterViewWrapper<real64>( viewKeys().poissonRatio.Key() )->
      setDefaultValue(-1)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("Poisson's ratio");

  RegisterViewWrapper<real64>( viewKeys().biotCoefficient.Key() )->
      setDefaultValue(0)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("Biot's coefficient");


  RegisterViewWrapper( viewKeys().compressibility.Key(), &m_compressibility, 0 )->
      setDefaultValue(-1)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("Fluid Compressibilty");

  RegisterViewWrapper( viewKeys().referencePressure.Key(), &m_referencePressure, 0 )->
      setDefaultValue(0)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("ReferencePressure");

  RegisterViewWrapper( viewKeys().biotCoefficient.Key(), &m_biotCoefficient, 0 )->
      setDefaultValue(-1)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("Young's Elastic Modulus");




  RegisterViewWrapper( viewKeyStruct::deviatorStressString, &m_deviatorStress, 0 )->
      setDescription("Stress Deviator stress");

  RegisterViewWrapper( viewKeyStruct::meanStressString, &m_meanStress, 0 )->
      setDefaultValue(-1)->
      setDescription("Young's Elastic Modulus");


  RegisterViewWrapper( viewKeyStruct::poreVolumeMultiplierString, &m_poreVolumeMultiplier, 0 )->
      setDefaultValue(-1)->
      setDescription("");

  RegisterViewWrapper( viewKeyStruct::dPVMult_dPresString, &m_dPVMult_dPressure, 0 )->
      setDefaultValue(-1)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("");

  RegisterViewWrapper( viewKeyStruct::bulkModulusString, &m_bulkModulus, 0 )->
      setDefaultValue(-1)->
      setDescription("Elastic Bulk Modulus Field");

  RegisterViewWrapper( viewKeyStruct::densityString, &m_density, 0 )->
      setDefaultValue(-1)->
      setDescription("Material Density");

  RegisterViewWrapper( viewKeyStruct::shearModulusString, &m_shearModulus, 0 )->
      setDefaultValue(-1)->
      setDescription("Elastic Shear Modulus");

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

  newConstitutiveRelation->m_compressibility   = this->m_compressibility;
  newConstitutiveRelation->m_referencePressure  = this->m_referencePressure;
  newConstitutiveRelation->m_biotCoefficient   = this->m_biotCoefficient;

  newConstitutiveRelation->m_poreVolumeRelation  = this->m_poreVolumeRelation;

  return std::move(newConstitutiveRelation);
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
  m_poreVolumeMultiplier.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dPVMult_dPressure.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_poreVolumeMultiplier = 1.0;

  m_bulkModulus = m_bulkModulus0;
  m_density = m_density0;
  m_shearModulus = m_shearModulus0;

}

void LinearElasticIsotropic::ReadXML_PostProcess()
{
  real64 & nu = getReference<real64>( viewKeys().poissonRatio );
  real64 & E  = getReference<real64>( viewKeys().youngsModulus );
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
    m_compressibility = 1 / K;
  }
  else if (m_compressibility <= 0)
  {
    string const message = std::to_string( numConstantsSpecified ) + " Elastic Constants Specified. Must specify 2 constants!";
    GEOS_ERROR( message );
  }
}

void LinearElasticIsotropic::FinalInitializationPreSubGroups(ManagedGroup * const parent)
{
  m_poreVolumeRelation.SetCoefficients( m_referencePressure, 1.0, m_compressibility );
}

void LinearElasticIsotropic::SetParamStatePointers( void *& data )
{

  this->m_dataPointers.m_bulkModulus0 = &m_bulkModulus0;
  this->m_dataPointers.m_shearModulus0 = &m_shearModulus0;
  this->m_dataPointers.m_bulkModulus = &m_bulkModulus;
  this->m_dataPointers.m_shearModulus = &m_shearModulus;
  this->m_dataPointers.m_meanStress = &m_meanStress;
  this->m_dataPointers.m_deviatorStress = &m_deviatorStress;

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
  array_view<real64> const & K = parameters->getReference<real64_array>( std::string( "BulkModulus" ));
  array_view<real64> const & G = parameters->getReference<real64_array>( std::string( "ShearModulus" ));

  array_view<real64>& mean_stress = stateVariables->getReference<real64_array>( std::string( "MeanStress" ));
  array_view<real64>& S11 = stateVariables->getReference<real64_array>( std::string( "S11" ));
  array_view<real64>& S22 = stateVariables->getReference<real64_array>( std::string( "S22" ));
  array_view<real64>& S33 = stateVariables->getReference<real64_array>( std::string( "S33" ));
  array_view<real64>& S23 = stateVariables->getReference<real64_array>( std::string( "S23" ));
  array_view<real64>& S13 = stateVariables->getReference<real64_array>( std::string( "S13" ));
  array_view<real64>& S12 = stateVariables->getReference<real64_array>( std::string( "S12" ));

  array_view<real64> const & D11 = input->getReference<real64_array>( std::string( "D11" ));
  array_view<real64> const & D22 = input->getReference<real64_array>( std::string( "D22" ));
  array_view<real64> const & D33 = input->getReference<real64_array>( std::string( "D33" ));
  array_view<real64> const & D23 = input->getReference<real64_array>( std::string( "D23" ));
  array_view<real64> const & D13 = input->getReference<real64_array>( std::string( "D13" ));
  array_view<real64> const & D12 = input->getReference<real64_array>( std::string( "D12" ));

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
    array_view<real64>& K11 = stateVariables->getReference<real64_array>( std::string( "K11" ));
    array_view<real64>& K22 = stateVariables->getReference<real64_array>( std::string( "K22" ));
    array_view<real64>& K33 = stateVariables->getReference<real64_array>( std::string( "K33" ));
    array_view<real64>& K23 = stateVariables->getReference<real64_array>( std::string( "K23" ));
    array_view<real64>& K13 = stateVariables->getReference<real64_array>( std::string( "K13" ));
    array_view<real64>& K12 = stateVariables->getReference<real64_array>( std::string( "K12" ));
    array_view<real64>& K44 = stateVariables->getReference<real64_array>( std::string( "K44" ));
    array_view<real64>& K55 = stateVariables->getReference<real64_array>( std::string( "K55" ));
    array_view<real64>& K66 = stateVariables->getReference<real64_array>( std::string( "K66" ));

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
