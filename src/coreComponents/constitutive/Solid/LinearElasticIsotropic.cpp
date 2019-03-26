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
  m_poreVolumeRelation()
{
  RegisterViewWrapper( viewKeyStruct::density0String, &m_density0, 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Reference Material Density");

  RegisterViewWrapper( viewKeyStruct::bulkModulus0String, &m_bulkModulus0, 0 )->
    setApplyDefaultValue(-1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Elastic Bulk Modulus Parameter");

  RegisterViewWrapper( viewKeyStruct::shearModulus0String, &m_shearModulus0, 0 )->
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

  RegisterViewWrapper( viewKeyStruct::biotCoefficientString, &m_biotCoefficient, 0 )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Biot's coefficient");

  RegisterViewWrapper( viewKeyStruct::compressibilityString, &m_compressibility, 0 )->
    setApplyDefaultValue(-1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Fluid Compressibilty");

  RegisterViewWrapper( viewKeyStruct::referencePressureString, &m_referencePressure, 0 )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("ReferencePressure");

  RegisterViewWrapper( viewKeyStruct::deviatorStressString, &m_deviatorStress, 0 )->
    setPlotLevel(PlotLevel::LEVEL_0)->
    setDescription("Stress Deviator stress");

  RegisterViewWrapper( viewKeyStruct::meanStressString, &m_meanStress, 0 )->
    setApplyDefaultValue(-1)->
    setPlotLevel(PlotLevel::LEVEL_0)->
    setDescription("Mean stress");


  RegisterViewWrapper( viewKeyStruct::poreVolumeMultiplierString, &m_poreVolumeMultiplier, 0 )->
    setApplyDefaultValue(-1)->
    setDescription("");

  RegisterViewWrapper( viewKeyStruct::dPVMult_dPresString, &m_dPVMult_dPressure, 0 )->
    setApplyDefaultValue(-1)->
    setDescription("");

  RegisterViewWrapper( viewKeyStruct::bulkModulusString, &m_bulkModulus, 0 )->
    setApplyDefaultValue(-1)->
    setDescription("Elastic Bulk Modulus Field");

  RegisterViewWrapper( viewKeyStruct::densityString, &m_density, 0 )->
    setApplyDefaultValue(-1)->
    setDescription("Material Density");

  RegisterViewWrapper( viewKeyStruct::shearModulusString, &m_shearModulus, 0 )->
    setApplyDefaultValue(-1)->
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
  m_bulkModulus.resize( parent->size() );
  m_shearModulus.resize( parent->size() );
  m_deviatorStress.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_density.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_meanStress.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_poreVolumeMultiplier.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dPVMult_dPressure.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_poreVolumeMultiplier = 1.0;

  m_bulkModulus = m_bulkModulus0;
  m_density = m_density0;
  m_shearModulus = m_shearModulus0;

}

void LinearElasticIsotropic::PostProcessInput()
{
  real64 & nu = getReference<real64>( viewKeyStruct::poissonRatioString );
  real64 & E  = getReference<real64>( viewKeyStruct::youngsModulus0String );
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
  m_poreVolumeRelation.SetCoefficients( m_referencePressure, 1.0, m_compressibility );
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

void LinearElasticIsotropic::GetStiffness( localIndex const k, real64 c[6][6] ) const
{
  real64 const G = m_shearModulus[k];
  real64 const Lame = m_bulkModulus[k] - 2.0/3.0 * G;
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
