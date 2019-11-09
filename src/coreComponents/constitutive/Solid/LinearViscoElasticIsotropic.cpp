/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 *  @file LinearViscoElasticIsotropic.cpp
 */

#include "LinearViscoElasticIsotropic.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace cxx_utilities;
namespace constitutive
{



LinearViscoElasticIsotropic::LinearViscoElasticIsotropic( std::string const & name, Group * const parent ):
  SolidBase( name, parent ),
  m_compressibility(),
  m_viscosity(),
  m_defaultBulkModulus(),
  m_defaultShearModulus(),
  m_bulkModulus(),
  m_shearModulus(),
  m_postProcessed(false)
{
  registerWrapper( viewKeyStruct::defaultBulkModulusString, &m_defaultBulkModulus, 0 )->
    setApplyDefaultValue(-1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Elastic Bulk Modulus Parameter");

  registerWrapper( viewKeyStruct::defaultShearModulusString, &m_defaultShearModulus, 0 )->
    setApplyDefaultValue(-1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Elastic Shear Modulus Parameter");

  registerWrapper<real64>( viewKeyStruct::defaultYoungsModulusString )->
    setApplyDefaultValue(-1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Elastic Young's Modulus.");

  registerWrapper<real64>( viewKeyStruct::defaultPoissonRatioString )->
    setApplyDefaultValue(-1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Poisson's ratio");

  registerWrapper( viewKeyStruct::bulkModulusString, &m_bulkModulus, 0 )->
    setApplyDefaultValue(-1)->
    setDescription("Elastic Bulk Modulus Field");

  registerWrapper( viewKeyStruct::shearModulusString, &m_shearModulus, 0 )->
    setApplyDefaultValue(-1)->
    setDescription("Elastic Shear Modulus");

  registerWrapper( viewKeyStruct::compressibilityString, &m_compressibility, 0 )->
    setApplyDefaultValue(-1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Rock Compressibilty");

  registerWrapper( viewKeyStruct::referencePressureString, &m_referencePressure, 0 )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("ReferencePressure");

  registerWrapper( viewKeyStruct::viscosityString, &m_viscosity, 0 )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Material Viscosity");
}


LinearViscoElasticIsotropic::~LinearViscoElasticIsotropic()
{}


void
LinearViscoElasticIsotropic::DeliverClone( string const & name,
                                      Group * const parent,
                                      std::unique_ptr<ConstitutiveBase> & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique<LinearViscoElasticIsotropic>( name, parent );
  }
  SolidBase::DeliverClone( name, parent, clone );
  LinearViscoElasticIsotropic * const newConstitutiveRelation = dynamic_cast<LinearViscoElasticIsotropic *>(clone.get());

  newConstitutiveRelation->m_compressibility      = m_compressibility;
  newConstitutiveRelation->m_viscosity      = m_viscosity;
  newConstitutiveRelation->m_defaultBulkModulus = m_defaultBulkModulus;
  newConstitutiveRelation->m_bulkModulus = m_bulkModulus;
  newConstitutiveRelation->m_defaultDensity = m_defaultDensity;
  newConstitutiveRelation->m_density = m_density;
  newConstitutiveRelation->m_defaultShearModulus = m_defaultShearModulus;
  newConstitutiveRelation->m_shearModulus = m_shearModulus;

  newConstitutiveRelation->m_stress = m_stress;
}

void LinearViscoElasticIsotropic::AllocateConstitutiveData( dataRepository::Group * const parent,
                                          localIndex const numConstitutivePointsPerParentIndex )
{
  SolidBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );
  m_bulkModulus.resize( parent->size() );
  m_shearModulus.resize( parent->size() );

  m_bulkModulus = m_defaultBulkModulus;
  m_shearModulus = m_defaultShearModulus;

}

void LinearViscoElasticIsotropic::PostProcessInput()
{

  if( !m_postProcessed )
  {
    real64 & nu = getReference<real64> (viewKeyStruct::defaultPoissonRatioString);
    real64 & E  = getReference<real64> (viewKeyStruct::defaultYoungsModulusString);
    real64 & K  = m_defaultBulkModulus;
    real64 & G  = m_defaultShearModulus;

    string errorCheck( "( ");
    int numConstantsSpecified = 0;
    if( nu >= 0.0 )
    {
      ++numConstantsSpecified;
      errorCheck += "nu, ";
    }
    if( E >= 0.0 )
    {
      ++numConstantsSpecified;
      errorCheck += "E, ";
    }
    if( K >= 0.0 )
    {
      ++numConstantsSpecified;
      errorCheck += "K, ";
    }
    if( G >= 0.0 )
    {
      ++numConstantsSpecified;
      errorCheck += "G, ";
    }
    errorCheck += ")";

    GEOS_ERROR_IF( numConstantsSpecified != 2,
                   "A specific pair of elastic constants is required. Either (K,G) or (E,nu). "<<
                   "You have specified "<<errorCheck );

    if( nu >= 0.0 && E >= 0.0 )
    {
      K = E / (3 * ( 1 - 2*nu ) );
      G = E / (2 * ( 1 + nu ) );
    }
    else if( nu >= 0.0 && G >= 0.0 )
    {
      E = 2 * G * ( 1 + nu );
      K = E / (3 * ( 1 - 2*nu ) );
    }
    else if( nu >= 0 && K >= 0.0 )
    {
      E = 3 * K * ( 1 - 2 * nu );
      G = E / ( 2 * ( 1 + nu ) );
    }
    else if( E >= 0.0 && K >=0 )
    {
      nu = 0.5 * ( 1 - E /  ( 3 * K ) );
      G = E / ( 2 * ( 1 + nu ) );
    }
    else if( E >= 0.0 && G >= 0 )
    {
      nu = 0.5 * E / G - 1.0;
      K = E / (3 * ( 1 - 2*nu ) );
    }
    else if( K >= 0.0 && G >= 0.0)
    {
      E = 9 * K * G / ( 3 * K + G );
      nu = ( 3 * K - 2 * G ) / ( 2 * ( 3 * K + G ) );
    }
    else
    {
      GEOS_ERROR( "invalid specification for default elastic constants. "<<errorCheck<<" has been specified.");
    }
  }

  m_compressibility = 1 / m_defaultBulkModulus;
  if (m_compressibility <= 0)
  {
    string const message = "Must specify 2 elastic constants!";
    GEOS_ERROR( message );
  }

  m_postProcessed = true;
}

void LinearViscoElasticIsotropic::StateUpdatePoint( localIndex const k,
                                                   localIndex const q,
                                                   R2SymTensor const & D,
                                                   R2Tensor const & Rot,
                                                   real64 const dt,
                                                   integer const GEOSX_UNUSED_ARG( updateStiffnessFlag ) )
{
  real64 meanStresIncrement = D.Trace();

  R2SymTensor temp = D;
  temp.PlusIdentity( -meanStresIncrement / 3.0 );
  R2SymTensor deviatorStrain = temp;
  temp *= 2.0 * m_shearModulus[k];
  meanStresIncrement *= m_bulkModulus[k];
  temp.PlusIdentity( meanStresIncrement );

  m_stress[k][q] += temp;

  // store elastic stress and add viscous stress into total stress
  m_elasticStress[k][q] = m_stress[k][q];
  temp.QijAjkQlk( m_elasticStress[k][q], Rot );
  m_elasticStress[k][q] = temp;
  m_stress[k][q] += m_viscosity / dt * deviatorStrain;

  temp.QijAjkQlk( m_stress[k][q], Rot );
  m_stress[k][q] = temp;
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, LinearViscoElasticIsotropic, std::string const &, Group * const )
}
} /* namespace geosx */
