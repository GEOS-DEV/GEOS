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
 *  @file LinearElasticAnisotropic.cpp
 */

#include "LinearElasticAnisotropic.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace cxx_utilities;
namespace constitutive
{



LinearElasticAnisotropic::LinearElasticAnisotropic( std::string const & name, ManagedGroup * const parent ):
  SolidBase( name, parent ),
  m_stiffness0{},
  m_stiffness{}
{

  RegisterViewWrapper( viewKeyStruct::c11, &(m_stiffness0.m_data[0][0]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c12, &(m_stiffness0.m_data[0][1]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c13, &(m_stiffness0.m_data[0][2]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c14, &(m_stiffness0.m_data[0][3]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c15, &(m_stiffness0.m_data[0][4]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c16, &(m_stiffness0.m_data[0][5]), 0 )->setInputFlag(InputFlags::REQUIRED);

  RegisterViewWrapper( viewKeyStruct::c21, &(m_stiffness0.m_data[1][0]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c22, &(m_stiffness0.m_data[1][1]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c23, &(m_stiffness0.m_data[1][2]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c24, &(m_stiffness0.m_data[1][3]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c25, &(m_stiffness0.m_data[1][4]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c26, &(m_stiffness0.m_data[1][5]), 0 )->setInputFlag(InputFlags::REQUIRED);

  RegisterViewWrapper( viewKeyStruct::c31, &(m_stiffness0.m_data[2][0]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c32, &(m_stiffness0.m_data[2][1]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c33, &(m_stiffness0.m_data[2][2]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c34, &(m_stiffness0.m_data[2][3]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c35, &(m_stiffness0.m_data[2][4]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c36, &(m_stiffness0.m_data[2][5]), 0 )->setInputFlag(InputFlags::REQUIRED);

  RegisterViewWrapper( viewKeyStruct::c41, &(m_stiffness0.m_data[3][0]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c42, &(m_stiffness0.m_data[3][1]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c43, &(m_stiffness0.m_data[3][2]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c44, &(m_stiffness0.m_data[3][3]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c45, &(m_stiffness0.m_data[3][4]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c46, &(m_stiffness0.m_data[3][5]), 0 )->setInputFlag(InputFlags::REQUIRED);

  RegisterViewWrapper( viewKeyStruct::c51, &(m_stiffness0.m_data[4][0]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c52, &(m_stiffness0.m_data[4][1]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c53, &(m_stiffness0.m_data[4][2]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c54, &(m_stiffness0.m_data[4][3]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c55, &(m_stiffness0.m_data[4][4]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c56, &(m_stiffness0.m_data[4][5]), 0 )->setInputFlag(InputFlags::REQUIRED);

  RegisterViewWrapper( viewKeyStruct::c61, &(m_stiffness0.m_data[5][0]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c62, &(m_stiffness0.m_data[5][1]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c63, &(m_stiffness0.m_data[5][2]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c64, &(m_stiffness0.m_data[5][3]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c65, &(m_stiffness0.m_data[5][4]), 0 )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( viewKeyStruct::c66, &(m_stiffness0.m_data[5][5]), 0 )->setInputFlag(InputFlags::REQUIRED);

  RegisterViewWrapper( viewKeyStruct::stiffness0String, &m_stiffness0, 0 )->
//    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Elastic Stiffness Tensor");

  RegisterViewWrapper( viewKeyStruct::stiffnessString, &m_stiffness, 0 )->
//    setApplyDefaultValue(-1)->
    setDescription("Elastic Modulus Field");

}


LinearElasticAnisotropic::~LinearElasticAnisotropic()
{}


void
LinearElasticAnisotropic::DeliverClone( string const & name,
                                      ManagedGroup * const parent,
                                      std::unique_ptr<ConstitutiveBase> & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique<LinearElasticAnisotropic>( name, parent );
  }
  SolidBase::DeliverClone( name, parent, clone );
  LinearElasticAnisotropic * const newConstitutiveRelation = dynamic_cast<LinearElasticAnisotropic *>(clone.get());


//  newConstitutiveRelation->m_stiffness0 = m_stiffness0;
  newConstitutiveRelation->m_stiffness = m_stiffness;
}

void LinearElasticAnisotropic::AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                          localIndex const numConstitutivePointsPerParentIndex )
{
  SolidBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );
  m_stiffness.resize( parent->size() );
  m_stiffness = m_stiffness0;
}

void LinearElasticAnisotropic::PostProcessInput()
{
  m_stiffness = m_stiffness0;
}

void LinearElasticAnisotropic::StateUpdatePoint( localIndex const k,
                                                 localIndex const q,
                                                 R2SymTensor const & D,
                                                 R2Tensor const & Rot,
                                                 integer const systemAssembleFlag )
{
  R2SymTensor T;
  real64 * const restrict Tdata = T.Data();
  real64 const * const restrict Ddata = D.Data();
  real64 (&c)[6][6] = m_stiffness[k].m_data;

  Tdata[0] += c[0][0]*Ddata[0] + c[0][1]*Ddata[2] + c[0][2]*Ddata[5] + c[0][3]*Ddata[4] + c[0][4]*Ddata[3] + c[0][5]*Ddata[1];
  Tdata[2] += c[1][0]*Ddata[0] + c[1][1]*Ddata[2] + c[1][2]*Ddata[5] + c[1][3]*Ddata[4] + c[1][4]*Ddata[3] + c[1][5]*Ddata[1];
  Tdata[5] += c[2][0]*Ddata[0] + c[2][1]*Ddata[2] + c[2][2]*Ddata[5] + c[2][3]*Ddata[4] + c[2][4]*Ddata[3] + c[2][5]*Ddata[1];
  Tdata[4] += c[3][0]*Ddata[0] + c[3][1]*Ddata[2] + c[3][2]*Ddata[5] + c[3][3]*Ddata[4] + c[3][4]*Ddata[3] + c[3][5]*Ddata[1];
  Tdata[3] += c[4][0]*Ddata[0] + c[4][1]*Ddata[2] + c[4][2]*Ddata[5] + c[4][3]*Ddata[4] + c[4][4]*Ddata[3] + c[4][5]*Ddata[1];
  Tdata[1] += c[5][0]*Ddata[0] + c[5][1]*Ddata[2] + c[5][2]*Ddata[5] + c[5][3]*Ddata[4] + c[5][4]*Ddata[3] + c[5][5]*Ddata[1];


  m_meanStress[k][q] = ( Tdata[0] + Tdata[2] + Tdata[5] ) / 3.0;
  T.PlusIdentity( -m_meanStress[k][q] );
  m_deviatorStress[k][q] += T;

  T.QijAjkQlk( m_deviatorStress[k][q], Rot );
  m_deviatorStress[k][q] = T;

}

void LinearElasticAnisotropic::GetStiffness( localIndex const k, real64 c[6][6] ) const
{
  for( int i=0 ; i<6 ; ++i )
  {
    for( int j=0 ; j<6 ; ++j )
    {
      c[i][j] = m_stiffness0.m_data[i][j];
    }
  }
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, LinearElasticAnisotropic, std::string const &, ManagedGroup * const )
}
} /* namespace geosx */
