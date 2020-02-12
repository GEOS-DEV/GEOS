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
  LinearElasticIsotropic( name, parent ),
  m_viscosity()
{
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
  LinearElasticIsotropic::DeliverClone( name, parent, clone );
  LinearViscoElasticIsotropic * const newConstitutiveRelation = dynamic_cast<LinearViscoElasticIsotropic *>(clone.get());

  newConstitutiveRelation->m_viscosity      = m_viscosity;
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

  m_elasticStress[k][q] += temp;

  // store elastic stress and add viscous stress into total stress
  m_stress[k][q] = m_elasticStress[k][q];
  m_stress[k][q] += m_viscosity / dt * deviatorStrain;
  temp.QijAjkQlk( m_elasticStress[k][q], Rot );
  m_elasticStress[k][q] = temp;

  temp.QijAjkQlk( m_stress[k][q], Rot );
  m_stress[k][q] = temp;
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, LinearViscoElasticIsotropic, std::string const &, Group * const )
}
} /* namespace geosx */
