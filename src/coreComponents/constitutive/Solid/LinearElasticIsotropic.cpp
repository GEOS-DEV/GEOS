/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

/**
 *  @file LinearElasticIsotropic.cpp
 */

#include "LinearElasticIsotropic.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace cxx_utilities;
namespace constitutive
{



LinearElasticIsotropic::LinearElasticIsotropic( std::string const & name, ManagedGroup * const parent ):
  SolidBase( name, parent ),
  m_defaultBulkModulus(),
  m_defaultShearModulus(),
  m_bulkModulus(),
  m_shearModulus()
{
  RegisterViewWrapper( viewKeyStruct::bulkModulus0String, &m_defaultBulkModulus, 0 )->
    setApplyDefaultValue(-1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Elastic Bulk Modulus Parameter");

  RegisterViewWrapper( viewKeyStruct::shearModulus0String, &m_defaultShearModulus, 0 )->
    setApplyDefaultValue(-1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Elastic Shear Modulus Parameter");


  RegisterViewWrapper<real64>( viewKeyStruct::youngsModulus0String )->
    setApplyDefaultValue(-1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Elastic Young's Modulus.");

  RegisterViewWrapper<real64>( viewKeyStruct::poissonRatioString )->
    setApplyDefaultValue(-1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Poisson's ratio");



  RegisterViewWrapper( viewKeyStruct::bulkModulusString, &m_bulkModulus, 0 )->
    setApplyDefaultValue(-1)->
    setDescription("Elastic Bulk Modulus Field");

  RegisterViewWrapper( viewKeyStruct::shearModulusString, &m_shearModulus, 0 )->
    setApplyDefaultValue(-1)->
    setDescription("Elastic Shear Modulus");

}


LinearElasticIsotropic::~LinearElasticIsotropic()
{}


void
LinearElasticIsotropic::DeliverClone( string const & name,
                                      ManagedGroup * const parent,
                                      std::unique_ptr<ConstitutiveBase> & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique<LinearElasticIsotropic>( name, parent );
  }
  SolidBase::DeliverClone( name, parent, clone );
  LinearElasticIsotropic * const newConstitutiveRelation = dynamic_cast<LinearElasticIsotropic *>(clone.get());


  newConstitutiveRelation->m_defaultBulkModulus = m_defaultBulkModulus;
  newConstitutiveRelation->m_bulkModulus = m_bulkModulus;
  newConstitutiveRelation->m_defaultDensity = m_defaultDensity;
  newConstitutiveRelation->m_density = m_density;
  newConstitutiveRelation->m_defaultShearModulus = m_defaultShearModulus;
  newConstitutiveRelation->m_shearModulus = m_shearModulus;

  newConstitutiveRelation->m_meanStress = m_meanStress;
  newConstitutiveRelation->m_deviatorStress = m_deviatorStress;
}

void LinearElasticIsotropic::AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                          localIndex const numConstitutivePointsPerParentIndex )
{
  SolidBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );
  m_bulkModulus.resize( parent->size() );
  m_shearModulus.resize( parent->size() );

  m_bulkModulus = m_defaultBulkModulus;
  m_shearModulus = m_defaultShearModulus;

}

void LinearElasticIsotropic::PostProcessInput()
{
  real64 & nu = getReference<real64>( viewKeyStruct::poissonRatioString );
  real64 & E  = getReference<real64>( viewKeyStruct::youngsModulus0String );
  real64 & K  = m_defaultBulkModulus;
  real64 & G  = m_defaultShearModulus;

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
}

void LinearElasticIsotropic::StateUpdatePoint( localIndex const k,
                                               localIndex const q,
                                               R2SymTensor const & D,
                                               R2Tensor const & Rot,
                                               integer const updateStiffnessFlag )
{
  real64 volumeStrain = D.Trace();
  m_meanStress[k][q] += volumeStrain * m_bulkModulus[k];

  R2SymTensor temp = D;
  temp.PlusIdentity( -volumeStrain / 3.0 );
  temp *= 2.0 * m_shearModulus[k];
  m_deviatorStress[k][q] += temp;


  temp.QijAjkQlk( m_deviatorStress[k][q], Rot );
  m_deviatorStress[k][q] = temp;
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, LinearElasticIsotropic, std::string const &, ManagedGroup * const )
}
} /* namespace geosx */
