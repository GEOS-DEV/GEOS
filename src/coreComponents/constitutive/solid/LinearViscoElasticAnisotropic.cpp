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
 *  @file LinearViscoElasticAnisotropic.cpp
 */

#include "LinearViscoElasticAnisotropic.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace cxx_utilities;
namespace constitutive
{

LinearViscoElasticAnisotropic::LinearViscoElasticAnisotropic( std::string const & name, Group * const parent ):
  LinearElasticAnisotropic( name, parent ),
  m_viscosity()
{
  registerWrapper( viewKeyStruct::viscosityString, &m_viscosity, 0 )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Material Viscosity");
}


LinearViscoElasticAnisotropic::~LinearViscoElasticAnisotropic()
{}


void
LinearViscoElasticAnisotropic::DeliverClone( string const & name,
                                      Group * const parent,
                                      std::unique_ptr<ConstitutiveBase> & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique<LinearViscoElasticAnisotropic>( name, parent );
  }
  LinearElasticAnisotropic::DeliverClone( name, parent, clone );
  LinearViscoElasticAnisotropic * const newConstitutiveRelation = dynamic_cast<LinearViscoElasticAnisotropic *>(clone.get());

  newConstitutiveRelation->m_viscosity      = m_viscosity;
}

void LinearViscoElasticAnisotropic::AllocateConstitutiveData( dataRepository::Group * const parent,
                                                            localIndex const numConstitutivePointsPerParentIndex )
{
  LinearElasticAnisotropic::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );
  m_elasticStress.resize( parent->size(), numConstitutivePointsPerParentIndex, 6 );
}

void LinearViscoElasticAnisotropic::StateUpdatePoint( localIndex const k,
                                                 localIndex const q,
                                                 R2SymTensor const & D,
                                                 R2Tensor const & Rot,
                                                 real64 const dt,
                                                 integer const GEOSX_UNUSED_PARAM( updateStiffnessFlag ) )
{
  R2SymTensor T;
  real64 * const GEOSX_RESTRICT Tdata = T.Data();
  real64 const * const GEOSX_RESTRICT Ddata = D.Data();
//  real64 const (&c)[6][6] = m_stiffness[k].m_data;

  Tdata[0] += m_c00[k]*Ddata[0] + m_c01[k]*Ddata[2] + m_c02[k]*Ddata[5] + m_c03[k]*(2*Ddata[4]) + m_c04[k]*(2*Ddata[3]) + m_c05[k]*(2*Ddata[1]);
  Tdata[2] += m_c10[k]*Ddata[0] + m_c11[k]*Ddata[2] + m_c12[k]*Ddata[5] + m_c13[k]*(2*Ddata[4]) + m_c14[k]*(2*Ddata[3]) + m_c15[k]*(2*Ddata[1]);
  Tdata[5] += m_c20[k]*Ddata[0] + m_c21[k]*Ddata[2] + m_c22[k]*Ddata[5] + m_c23[k]*(2*Ddata[4]) + m_c24[k]*(2*Ddata[3]) + m_c25[k]*(2*Ddata[1]);
  Tdata[4] += m_c30[k]*Ddata[0] + m_c31[k]*Ddata[2] + m_c32[k]*Ddata[5] + m_c33[k]*(2*Ddata[4]) + m_c34[k]*(2*Ddata[3]) + m_c35[k]*(2*Ddata[1]);
  Tdata[3] += m_c40[k]*Ddata[0] + m_c41[k]*Ddata[2] + m_c42[k]*Ddata[5] + m_c43[k]*(2*Ddata[4]) + m_c44[k]*(2*Ddata[3]) + m_c45[k]*(2*Ddata[1]);
  Tdata[1] += m_c50[k]*Ddata[0] + m_c51[k]*Ddata[2] + m_c52[k]*Ddata[5] + m_c53[k]*(2*Ddata[4]) + m_c54[k]*(2*Ddata[3]) + m_c55[k]*(2*Ddata[1]);

  R2SymTensor deviatorStrain = D;
  deviatorStrain.PlusIdentity( -D.Trace() / 3.0 );

  R2SymTensor temp2 = m_elasticStress[ k ][ q ];
  temp2 += T;

  T.QijAjkQlk( temp2, Rot );

  for ( localIndex i = 0; i < 6; ++i )
  {
    m_elasticStress( k, q, i ) = T.Data()[ i ];
  }

  temp2 += m_viscosity / dt * deviatorStrain;

  T.QijAjkQlk( temp2, Rot );

  for ( localIndex i = 0; i < 6; ++i )
  {
    m_stress( k, q, i ) = T.Data()[ i ];
  }
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, LinearViscoElasticAnisotropic, std::string const &, Group * const )
}
} /* namespace geosx */
